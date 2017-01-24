#include <cmath>
#include <iostream>
#include <memory>

#include <gismo.h>

namespace gismo
{

std::unique_ptr<gsBSpline<>> loadBoundaryCurve(int n , int degree)
{
    gsKnotVector<real_t> kv(0 , 1 , n - degree - 1 , degree + 1);

    gsBSplineBasis<real_t> basis(kv);

    gsMatrix<> greville = basis.anchors();
    gsMatrix<> coefs(greville.cols() , 2);

    for (index_t col = 0 ; col != greville.cols() ; ++col) {
        real_t t = greville(0, col);

        coefs(col, 0) = sin(2. * M_PI * t);
        coefs(col, 1) = cos(2. * M_PI * t);
    }

    return std::unique_ptr<gsBSpline<real_t>>
           (new gsBSpline<>(basis , coefs));
}

std::unique_ptr<gsBSpline<>> loadBoundaryFromFile(const std::string &filename)
{
    gsFileData<> file(filename);

    if (file.has< gsBSpline<> >()) {
        gsBSpline<> *pBSpline;
        pBSpline = file.getFirst< gsBSpline<> >();
        return std::unique_ptr<gsBSpline<>>(pBSpline);
    }
    else
    {
        gsWarn << "Could not load boundary file : " << filename << "!\n\n";
    }

    return std::unique_ptr<gsBSpline<>>();
}

void saveParameterization( const gsTensorBSpline<2,real_t>& surface , std::string& filename , std::string postfix )
{
    // write gismo xml file
    gsFileData<real_t> file;
    file << surface;
    file.save( filename + postfix );
    
    // write paraview surfaces
    gsWriteParaview( surface , filename + postfix , 2000 );
    
    gsMesh<> controlNet;
    surface.controlNet( controlNet );    
    gsWriteParaview( controlNet , filename + "_net" + postfix );
}


bool isClockwise( const gsBSpline<>& curve )
{
    const gsMatrix<>& coefs = curve.coefs();
    
    double orientation = 0.;
    for( int i = 1 ; i < coefs.cols() ; ++i )
    {
        const gsMatrix<>& diff = coefs.col(i) - coefs.col(i-1);
        orientation += diff(1) * diff(2);
    }
    
    return orientation > 0;
}

void fillBoundaryData(const gsBSpline<> &boundary
                      , int n
                      , int m
                      , gsMatrix< gsMatrix<real_t, 2, 1> > &X
                      , real_t start
                      , bool bRectangle )
{
    X.resize(n, m);

    // fill in the allready know boundary points

    const int numberOfPoints = 2 * (m + n) - 2;

    gsMatrix<real_t> boundaryPoints(2 , numberOfPoints);

    gsMatrix<real_t> parameterRange = boundary.parameterRange();
    
    real_t sign = 1.;
    if( !isClockwise(boundary) )
        sign = -1;
    
    gsMatrix<real_t> times( 1 , numberOfPoints );
    if( !isClockwise(boundary) )
        times = uniformPointGrid<real_t>( parameterRange.block(0,0,1,1) , parameterRange.block(0,1,1,1) , numberOfPoints );
    else
        times = uniformPointGrid<real_t>( parameterRange.block(0,1,1,1) , parameterRange.block(0,0,1,1) , numberOfPoints );
    
    real_t a = parameterRange(0);
    real_t b = parameterRange(1);
    int iStart = static_cast<int>( start / (b-a) * static_cast<real_t>(numberOfPoints));

    if ( iStart != 0 )
    {
        // TODO: Replace
        gsInfo << "Internal: Parameter start in fillBoundary is unimplemented\n";
    }

    
    
    boundary.eval_into(times , boundaryPoints);

    for (int i = 0 ; i < n - 1 ; ++i)
        X(i, 0) = boundaryPoints.col(i);

    for (int j = 0 ; j < m - 1 ; ++j)
        X(n - 1, j) = boundaryPoints.col(n - 1 + j);

    for (int i = n - 1 ; i > 0  ; --i)
        X(i, m - 1) = boundaryPoints.col(n - 1 + m - 1 + n - 1 - i);

    for (int j = m - 1 ; j > 0 ; --j)
        X(0, j) = boundaryPoints.col(n - 1 + m - 1 + n - 1 + m - 1 - j);

    //TODO: How to choose X(0,0) ? depending on b(0) or b(end) ?
    X(0,0) = 0.5 * boundaryPoints.col(0) + 0.5 * boundaryPoints.col( numberOfPoints -1 );
}

gsTensorBSpline<2, real_t> bilinearlyBlendedCoons(const gsBSpline<> &boundary
, int n , int m , int degree , bool bRectangle )
{
    gsMatrix<gsMatrix<real_t, 2, 1>> X(n, m);

    // fill in the allready know boundary points
    fillBoundaryData(boundary, n , m , X , 0. , bRectangle );

    const int numberOfPoints = 2 * (m + n) - 4;

    // Ref: Discrete Coons patches
    // 1. Background (2)  Bilinearly blended coons patch

    gsMatrix<real_t, 2, 4> B;
    B << X(0, 0) , X(0, m - 1) , X(n - 1, 0) , X(n - 1, m - 1);

    real_t px , py;

    for (int j = 1 ; j < m - 1 ; ++j)
        for (int i = 1 ; i < n - 1 ; ++i) {
            px = static_cast<real_t>(i) / static_cast<real_t>(n);
            py = static_cast<real_t>(j) / static_cast<real_t>(m);

            gsMatrix<real_t, 4, 1> p;

            p	<< (1. - px)*(1. - py) , (1. - px)* py
                ,     px *(1. - py) ,     px *py;

            X(i, j)  = (1. - px) * X(0, j) + px * X(n - 1, j)
                       + (1. - py) * X(i, 0) + py * X(i, m - 1)
                       - B * p;
        }

    // now we have created all control points, next we construct
    // the Tensor BSpline

    gsKnotVector<real_t> kvx(0. , 1. , n - degree + 1 , degree - 1);
    gsKnotVector<real_t> kvy(0. , 1. , m - degree + 1 , degree - 1);


    gsMatrix<real_t> coefs(m * n , 2);

    int k = 0;
    for (int j = 0 ; j < m ; ++j)
        for (int i = 0 ; i < n ; ++i) {
            coefs.row(k++) = X(i, j).transpose();
        }

    gsTensorBSplineBasis<2, real_t> basis(kvx , kvy);
    // use std::move  (but I excpect the compiler to be clever enought to optimize this...)
    return gsTensorBSpline<2, real_t>(basis , coefs);
}

int fromMatrixToVector(int i , int j , int n , int)
{
    return i + j * n;
}

gsTensorBSpline<2, real_t> minimumPrinciple(const gsBSpline<> &boundary
, int n , int m , int degree , double start , bool bRectangle)
{
    gsMatrix< gsMatrix<real_t, 2, 1> > X;
    fillBoundaryData(boundary , n , m , X , start , bRectangle );

    // we want to access the points in the inside.
    // the matrix is stored column major, i.e. a11, a21, a31

    const size_t numberOfBoundaryPoints =     2 * (n + m - 2);
    const size_t numberOfInnerPoints    = (n - 2) * (m - 2);
    const size_t numberOfEquations      = (n - 1) * (m - 1);

    gsMatrix<real_t, Eigen::Dynamic, 1> b_x(numberOfEquations, 1);
    gsMatrix<real_t, Eigen::Dynamic, 1> b_y(numberOfEquations, 1);

    // assembly step, maybe we can use gsAssembly
    gsSparseMatrix< real_t > A(numberOfEquations , numberOfInnerPoints);

    std::vector< Eigen::Triplet<real_t> > tripletList;
    tripletList.reserve(4 * numberOfEquations); // only a approximation

    for (int j = 0 ; j < m - 1 ; ++j)
        for (int i = 0 ; i < n - 1 ; ++i) {
            // equation index
            int iEquation = i + j * (n - 1);
            b_x(iEquation) = 0.;
            b_y(iEquation) = 0.;

            // index refering to inner points
            int i0 = fromMatrixToVector(i - 1   , j - 1   , n - 2 , m - 2);
            int i1 = fromMatrixToVector(i - 1 + 1 , j - 1   , n - 2 , m - 2);
            int i2 = fromMatrixToVector(i - 1 + 1 , j - 1 + 1 , n - 2 , m - 2);
            int i3 = fromMatrixToVector(i - 1   , j - 1 + 1 , n - 2 , m - 2);

            // detect boundary points
            if (i > 0 && j > 0 && i < n - 1 && j < m - 1)
                tripletList.push_back(Eigen::Triplet<real_t>(iEquation , i0 , 1.));
            else {
                b_x(iEquation) -= X(i, j)(0);
                b_y(iEquation) -= X(i, j)(1);
            }

            if (i + 1 > 0 && j > 0 && i + 1 < n - 1 && j < m - 1)
                tripletList.push_back(Eigen::Triplet<real_t>(iEquation , i1 , -1.));
            else {
                b_x(iEquation) += X(i + 1, j)(0);
                b_y(iEquation) += X(i + 1, j)(1);
            }


            if (i + 1 > 0 && j + 1 > 0 && i + 1 < n - 1 && j + 1 < m - 1)
                tripletList.push_back(Eigen::Triplet<real_t>(iEquation , i2 , 1.));
            else {
                b_x(iEquation) -= X(i + 1, j + 1)(0);
                b_y(iEquation) -= X(i + 1, j + 1)(1);
            }

            if (i > 0 && j + 1 > 0 && i < n - 1 && j + 1 < m - 1)
                tripletList.push_back(Eigen::Triplet<real_t>(iEquation , i3 , -1.));
            else {
                b_x(iEquation) += X(i, j + 1)(0);
                b_y(iEquation) += X(i, j + 1)(1);
            }
        }

    A.setFromTriplets(tripletList.begin() , tripletList.end());

    // Ref: G Farin, 2. The minimum principle
    // we use QR-decomposition to solve the least square problem min || Ax-b ||

    gsEigenAdaptor<real_t>::SparseQR QR(A);
    gsMatrix<real_t> inner_x = QR.solve(b_x);
    gsMatrix<real_t> inner_y = QR.solve(b_y);



    // now we have created all control points, next we construct
    // the Tensor BSpline
    gsMatrix<real_t> coefs(m * n , 2);

    // using map would be better
    // Eigen::Map< Eigen::MatrixX<real_t> >
    // maybe we use use a higher dimensional system
    // including also the boundary values. Then no composition is nessecary anymore.

    int k = 0;
    int k_in = 0;
    for (int j = 0 ; j < m ; ++j)
        for (int i = 0 ; i < n ; ++i) {
            if (i > 0 && j > 0 && i < n - 1 && j < m - 1) {
                coefs(k, 0) = inner_x(k_in);
                coefs(k, 1) = inner_y(k_in++);
            } else
                coefs.row(k) = X(i, j).transpose();

            ++k;
        }

    gsKnotVector<real_t> kvx(0. , 1. , n - degree + 1 , degree - 1);
    gsKnotVector<real_t> kvy(0. , 1. , m - degree + 1 , degree - 1);

    gsTensorBSplineBasis<2, real_t> basis(kvx , kvy);
    return gsTensorBSpline<2, real_t>(basis , coefs);

}

gsTensorBSpline<2, real_t> unidirectionalInterpolation(gsBSpline<> boundary
        , int cOnBoundary
        , int cHorizontal
        , int degree
        , bool bRectangle )
{
    if (cHorizontal < 2) {
        gsInfo << "cHorizontal must be >= 2 \n\n";
        return gsTensorBSpline<2, real_t>();
    }

    gsMatrix<real_t> times(1 , 2 * cOnBoundary);

    for (int i = 0 ; i < cOnBoundary ; ++i)
        times(0 , i) = 0.5 * i / static_cast<double>(cOnBoundary - 1);

    for (int i = 0 ; i < cOnBoundary ; ++i)
        times(0 , cOnBoundary + i) = 1. - 0.5 * i / static_cast<double>(cOnBoundary - 1);

    gsMatrix<real_t> coefs(2 , 2 * cOnBoundary);

    boundary.eval_into(times, coefs);

    gsKnotVector<real_t> kvx(0. , 1. , cOnBoundary - degree + 1 , degree - 1);
    gsKnotVector<real_t> kvy(0. , 1. , 2 - degree + 1            , degree - 1);

    gsTensorBSplineBasis<2, real_t> basis(kvx , kvy);
    gsTensorBSpline<2, real_t> nurbs(basis , coefs.transpose());

    // refine, to get 'cHorizontal' many inner control points

    for (int i = 1 ; i < cHorizontal - 1 ; ++i)
        nurbs.insertKnot(static_cast<double>(i) / (cHorizontal - 1) , 1 , 1);

    return nurbs;
}

typedef gsMatrix<real_t, 1, 2> TVec2;

bool isInHalfspace(const TVec2 &a , TVec2 b , TVec2 p)
{
    p = p - a;
    b = b - a;
    return b(0, 0) * p(0, 1) - b(0, 1) * p(0, 0) >= 0.;
}


gsMatrix<real_t> findKernelPoint(gsMatrix<real_t> p)
{

    // idea:
    // consider the half space, spanned by two connected points of the
    // polygon.
    // Delete all points which are on the wrong side,
    // add sufficient many points to create a valid polygon which lies
    // inside the original polygon

    // 1. we need to know if the points are arranged clockwise or counter-clockwise
    // we use the schoelace formular

    //auto X = [&](int i){ return polygon(i,0);   };
    //auto Y = [&](int i){ return polygon(i,1);   };
    //auto P = [&](int i){ return polygon.row(i); };

    real_t area = 0.;

    int n = p.cols();

    for (int i = 0 ; i < p.cols() ; ++i)
        area += 0.5 * (p(i, 0) * p(i + 1 % n, 1) - p(i + 1 % n, 0) * p(i, 1));

    bool isClockwise = area > 0.;

    gsInfo << "Clockwise: " << (isClockwise ? "yes" : "no") << "\n";

    if (!isClockwise) {
        p.colwise() = p.colwise().reverse().eval();
        isClockwise = true;
    }

    int i = 0;

    while (i < n) {
        int j = i + 1 % n;

        // consider a line l through the points p_i , p_j.
        // Check if the point k = j+1 , ... i-1 % n
        // lies under or over the line l.


        int k = j + 1;
        while (k % n != i) {



            ++k;
        }

        ++i;
    }
}

gsTensorBSpline<2, real_t> naivePolarChart(gsBSpline<> boundary
        , int numberOfPolarLines
        , int degree
        , bool bRectangle )
{
    // first of, we need a periodic boundary
    // boundary.setPeriodic();
    // boundary.basis().enforceOuterKnotsPeriodic();

    // we need to find a point in the middle,
    // say the average of all controlpoints.
    // (Maybe in 'Cormen' better algorithms can be found.

    // to interpolate the boundary exactly, we need to ensure just to
    // use the boundary. A higher level of detail can be reached by
    // knot insertion afterwards.

    /*
    gsMatrix<real_t,1,2> center;
    center << 0. , 0.;

    for( size_t i = 0 ; i < boundary.coefsSize() ; ++i )
    {
        center += boundary.coef( i );
    }

    center /= static_cast<real_t>( boundary.coefsSize() );
    */

    // bilinearly blended coons:
    //

    /*
    gsMatrix<real_t,1,2> center;
    center << 0. , 0.;

    gsMatrix<real_t> times(1,8);
    for( int i = 0 ; i < 8 ; ++i )
        times(i) = i / 8.;

    gsMatrix<real_t> b(2,8);
    boundary.eval_into( times , b );

    center.transpose() = 0.5  * ( b.col(1) + b.col(3) + b.col(5) + b.col(7) )
           - 0.25 * ( b.col(0) + b.col(2) + b.col(4) + b.col(6) );

    */


    gsMatrix<real_t, 1, 2> center;
    center << 0.5 , 0.5;

    gsMatrix<real_t> coefs(boundary.coefsSize() * numberOfPolarLines , 2);

    const double dNumOfPolarLines = static_cast<double>(numberOfPolarLines);
    index_t k = 0;
    for (double j = 0 ; j < boundary.coefsSize() ; ++j) {
        for (double i = 0 ; i < numberOfPolarLines ; ++i) {
            double p = i / (dNumOfPolarLines - 1.);
            coefs.row(k++) =  p * boundary.coef(j) + (1.-p) * center;
        }
    }

    gsKnotVector<real_t> kvx(0. , 1. , numberOfPolarLines   - degree + 1 , degree - 1);
    gsKnotVector<real_t> kvy(0. , 1. , boundary.coefsSize() - degree + 1 , degree - 1);

    gsTensorBSplineBasis<2, real_t> basis(kvx , kvy);
    return gsTensorBSpline<2, real_t>(basis , coefs);
}




} //end namespace gismo
