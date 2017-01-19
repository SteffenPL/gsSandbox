#include <cmath>
#include <iostream>
#include <memory>

#include <gismo.h>

namespace gismo
{

std::unique_ptr<gsBSpline<> > loadBoundaryCurve(int n , int degree);

std::unique_ptr<gsBSpline<> > loadBoundaryFromFile(const std::string &filename);

/*
 * Save the surface in various fileformats: .xml , Paraview: Surface/Controlnet
 */
void saveParameterization( const gsTensorNurbs<2,real_t>& surface , std::string& filename , std::string postfix = "" );

bool isClockwise( const gsBSpline<>& curve );

/**
 * n : number of points in x dir
 * m : number of points in y dir
 * X : return matrix of control points, where the boundary points are set and the interior is zero
 * start : X(0,0) = boundary( start )
 * bRectangle : Make sure to map the points boundary( t = {0,0.25,0.5,0.75} ) onto corners of X.
 */
void fillBoundaryData(const gsBSpline<> &boundary
, int n
, int m
, gsMatrix< gsMatrix<real_t, 2, 1> > &X
, real_t start = 0.
, bool bRectangle = false
);

gsTensorNurbs<2, real_t> bilinearlyBlendedCoons(const gsBSpline<> &boundary
, int n , int m , int degree , bool bRectangle = false);

gsTensorNurbs<2, real_t> naivePolarChart(gsBSpline<> boundary
, int numberOfPolarLines
, int degree
, bool bRectangle = false
);


gsTensorNurbs<2, real_t> unidirectionalInterpolation(gsBSpline<> boundary
, int cOnBoundary
, int cHorizontal
, int degree
, bool bRectangle = false );

/**
 * convert 2d index to 1d (array) index
 */
int fromMatrixToVector(int i , int j , int n , int);

gsTensorNurbs<2, real_t> minimumPrinciple(const gsBSpline<> &boundary
, int n , int m , int degree , double start = 0. , bool bRectangle = false);



gsTensorNurbs<2, real_t> centricPotentialFlow(const gsBSpline<> &curve ,
                                              int numberOfPolarLines ,
                                              int numberOfPolarCircles ,
                                              real_t midPoint_x ,
                                              real_t midPoint_y);




} //end namespace gismo
