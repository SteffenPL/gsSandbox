#include "OptParameterization.h"

// how to include Eigen correctly here?
#include <Eigen/Dense>

#include <limits>

namespace gismo
{
    
template< typename T >
T determinant(const gsMatrix<T, 2, 2> &J)
{
    return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
}
    
    
// geometric measures
template< typename T >
T integrateFunctional( const gsTensorNurbs<2,T>& surface, T(*F)( const gsMatrix<T,2,2>&), bool onParameterSpace )
{
    /*
    gsVector<T> p1(2);
    gsVector<T> p2(2);
    p1 << T(0.) , T(0.);
    p2 << T(1.) , T(1.);    
    gsMatrix<T> collocationPoints = uniformPointGrid( p1 , p2  );*/

    
    const int quadratureDegree = 6;
    
    gsMatrix<T> collocationPoints( 2, pow(quadratureDegree,2) );
    
    // Source: https://pomax.github.io/bezierinfo/legendre-gauss.html
    double weights[6] = {0.3607615730481386,0.3607615730481386,0.4679139345726910,
            0.4679139345726910,0.1713244923791704,0.1713244923791704};
        
    double points[6] = {0.6612093864662645,-0.6612093864662645,-0.2386191860831969,
            0.2386191860831969,-0.9324695142031521,0.9324695142031521};
    
    
    const gsKnotVector<T>& knotsX = surface.knots(0);
    const gsKnotVector<T>& knotsY = surface.knots(1);
    
    //gsInfo << surface.coefs() << " G " <<greville<< "\n";
    gsMatrix<T> DSurface;
    T integral = T(0.);
        
    gsMatrix<T> pos(2,1);
    pos << T(0.) , T(0.);
    typename gsKnotVector<T>::uiterator nextX , nextY;
    for( typename gsKnotVector<T>::uiterator itX = knotsX.ubegin() ; itX != knotsX.uend()-1 ; ++itX )
    {
        pos(0) = *itX;
        nextX = itX + 1;
        T dx = *nextX - *itX;
        for( typename gsKnotVector<T>::uiterator itY = knotsY.ubegin() ; itY != knotsY.uend()-1 ; ++itY )
        {
            pos(1) = *itY;
            nextY = itY + 1;
            T dy = *nextY - *itY;
            
            for( int i = 0 ; i < quadratureDegree ; ++i )
                for( int j = 0; j < quadratureDegree ; ++j )
                {
                    collocationPoints(0, quadratureDegree*j + i ) = pos(0) + dx*(0.5 + 0.5*points[i]);
                    collocationPoints(1, quadratureDegree*j + i ) = pos(1) + dy*(0.5 + 0.5*points[j]);
                }    
                
            surface.deriv_into( collocationPoints , DSurface );            
            // now the jacobians are represented as columns of the matrix
            // 'DSurface' as ( df1/d1 , df1/d2 , df2 /d1 , df2/d2 )^T
                                    
            T det = T(1.);
            for( int i = 0 ; i < quadratureDegree ; ++i )
            {
                for( int j = 0; j < quadratureDegree ; ++j )
                {
                    int index = quadratureDegree*j+i;
                    gsMatrix<T,2,2> Jac;
                    Jac << DSurface(0,index) , DSurface(1,index) , DSurface(2,index) , DSurface(3,index);
                    
                    if( !onParameterSpace )
                    {
                        det = Jac.determinant();
                        integral +=  0.5*dx*dy*weights[i]*weights[j] * fabs(det) * F( Jac );
                    }
                    else
                    {
                        integral +=  0.5*dx*dy*weights[i]*weights[j] * F( Jac );
                    }
                }
            }
            
        }
    }
    
    return integral;
}


template<typename T>
OptParameterization<T>::OptParameterization( bool forcePositiveDet ):
    m_numElemX(0),
    m_numElemY(0),
    m_functional( areaOrthogonalityMeasure<T> ),
    m_bForcePositiveDet( forcePositiveDet )
{}

//** default destructor */
template<typename T>
OptParameterization<T>::~OptParameterization()
{}


template<typename T>
void OptParameterization<T>::setParameterization( const gsTensorNurbs<2,T>& parameterization )
{
    m_param = parameterization;
    
    m_numElemX = m_param.knots(0).uSize();
    m_numElemY = m_param.knots(1).uSize();
    
    const T inf = std::numeric_limits<T>::infinity(); //or should be use max instead?
    
    this->m_numDesignVars  = 2*(m_numElemX-2) * (m_numElemY-2);
    this->m_desLowerBounds = gsMatrix<T>::Constant( this->m_numDesignVars , 1 , -inf );
    this->m_desUpperBounds = gsMatrix<T>::Constant( this->m_numDesignVars , 1 ,  inf );
    
    if( m_bForcePositiveDet )
    {            
        this->m_numConstraints =  4 * (m_numElemX-1) * (m_numElemY-1);
        this->m_numConJacNonZero = 0;
        this->m_conLowerBounds = gsMatrix<T>::Constant( this->m_numConstraints , 1 ,   0. );
        this->m_conUpperBounds = gsMatrix<T>::Constant( this->m_numConstraints , 1 ,  inf );
    }
    else
    {
        this->m_numConstraints =   0;
        this->m_numConJacNonZero = 0;
        this->m_conLowerBounds = gsMatrix<T>::Constant( this->m_numConstraints , 1 , -inf );
        this->m_conUpperBounds = gsMatrix<T>::Constant( this->m_numConstraints , 1 ,  inf );
    }
        
    
    int k = 0;    
    this->m_curDesign.resize( this->m_numDesignVars , 1 );
    for( int j = 1 ; j < m_numElemY-1 ; ++j )
        for( int i = 1 ; i < m_numElemX-1 ; ++i )
        {
            this->m_curDesign( k++ ) = m_param.coef( i + m_numElemX * j , 0 );
            this->m_curDesign( k++ ) = m_param.coef( i + m_numElemX * j , 1 );
        }
}

//TODO: Use move schematics
template<typename T>
const gsTensorNurbs<2,T>& OptParameterization<T>::getParameterization()
{
    //updateParameterization( gsAsConstVector<T>( this->m_curDesign.data() , this->m_curDesign.cols()* this->m_curDesign.rows() ) );
    
    int k = 0;
    
    for( int j = 1 ; j < m_numElemY-1 ; ++j )
        for( int i = 1 ; i < m_numElemX-1 ; ++i )
        {
            m_param.coef( i + m_numElemX*j , 0 ) = this->m_curDesign( k++ );
            m_param.coef( i + m_numElemX*j , 1 ) = this->m_curDesign( k++ );
        }  
         
    //m_param.coefs().asVector() = this->m_curDesign;
    return m_param;
}

template< typename T >
void OptParameterization<T>::updateParameterization(const gsAsConstVector<T>& u ) const
{
    int k = 0;
    
    for( int j = 1 ; j < m_numElemY-1 ; ++j )
        for( int i = 1 ; i < m_numElemX-1 ; ++i )
        {
            m_param.coef( i + m_numElemX*j , 0 ) = u( k++ );
            m_param.coef( i + m_numElemX*j , 1 ) = u( k++ );
        }    
}

template<typename T>
T areaOrthogonalityMeasure( const gsMatrix<T,2,2>& J )
{
    gsMatrix<T,2,2> g = J.transpose()*J;
    return g(0,0) * g(1,1);
}

template<typename T>
T liaoFunctional( const gsMatrix<T,2,2>& J )
{
    gsMatrix<T,2,2> g = J.transpose()*J;
    return g(0,0) * g(0,0) + 
    2 * g(1,0) * g(1,0) + 
    g(1,1) * g(1,1) ;
}

template<typename T>
T winslowFunctional( const gsMatrix<T,2,2>& J )
{
    return 
    ( pow( J(0,0) , 2 ) + 
    pow( J(1,0) , 2 ) +
    pow( J(1,1) , 2 ) +
    pow( J(0,1) , 2 ) ) / fabs( ( J(0,0)*J(1,1) - J(1,0)*J(0,1) ) );
}

template<typename T>
T detMeasure( const gsMatrix<T,2,2>& J )
{
    return ( (J(0,0)*J(1,1) - J(0,1)*J(1,0))-1 ) * ((J(0,0)*J(1,1) - J(0,1)*J(1,0))-1) ;
}

template<typename T>
T contMechanics( const gsMatrix<T,2,2>& J )
{
    T c_1,c_2;
    c_1 = T(0.);
    c_2 = T(1.);
    gsMatrix<T,2,2> g = J.transpose()*J;
    T det = J.determinant();
    return c_1*( g.trace() - 1. - 2.*det ) + c_2*( det - 1. )*( det-1. ); 
}

template<typename T>
T OptParameterization<T>::evalObj( const gsAsConstVector<T> & u ) const
{
    // we interpret the current state u as a vector containing the positions of the 
    // control points. The Mapping is done by the standard aligment in memory.
    // Therefore depends on the current gismo implementations.
    
    //const_cast<gsTensorNurbs<2,T>&>(m_param).setCoefs( gsAsMatrix<T>( (double*)u.data() , m_param.coefs().rows() , m_param.coefs().cols() ) );
    
    updateParameterization(u);
    return integrateFunctional<T>( m_param , m_functional );
}

/*
/// \brief Returns the gradient of the objective function at design value
/// \a u
/// By default it uses finite differences, overriding it should provide exact gradient.
template<typename T>
void OptParameterization<T>::gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{
    const_cast<gsTensorNurbs<2,T>&>(m_param).setCoefs( gsAsMatrix<T>( (double*)u.data() , m_param.coefs().rows() , m_param.coefs().cols() ) );
    
    result(0) = 2*m_param.coefs()(0,0);
}
*/

/// \brief Returns values of the constraints at design value \a u
template<typename T>
void OptParameterization<T>::evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{
    if( m_bForcePositiveDet )
    {
        updateParameterization(u);
        
        const gsKnotVector<T>& knotsX = m_param.knots(0);
        const gsKnotVector<T>& knotsY = m_param.knots(1);
        typename gsKnotVector<T>::uiterator nextX , nextY , itX , itY;
        
        gsMatrix<T> collocationPoints( 2 , 4*(m_numElemX-1) * (m_numElemY-1) );
        T dx , dy;
        int k = 0;
        for( itX = knotsX.ubegin() ; itX != knotsX.uend()-1 ; ++itX )
        {
            nextX = itX+1;
            dx = *nextX - *itX;
            for( itY = knotsY.ubegin() ; itY != knotsY.uend()-1 ; ++itY )
            {
                nextY = itY+1;
                dy = *nextY - *itY;
                
                nextY = itY+1;
                collocationPoints(0,k) = *itX + 0.2*dx;
                collocationPoints(1,k) = *itY + 0.2*dy;
                ++k;
                collocationPoints(0,k) = *itX + 0.8*dx;
                collocationPoints(1,k) = *itY + 0.2*dy;
                ++k;
                collocationPoints(0,k) = *itX + 0.2*dx;
                collocationPoints(1,k) = *itY + 0.8*dy;
                ++k;
                collocationPoints(0,k) = *itX + 0.8*dx;
                collocationPoints(1,k) = *itY + 0.8*dy;
                ++k;
            }
        }
        
        gsMatrix<T> D(4, collocationPoints.cols() );
        m_param.deriv_into( collocationPoints , D );
        result = D.row(0).array() * D.row(3).array() - D.row(2).array() * D.row(1).array();
    }
    // can we use rows().array() here instead?
}

/// \brief Returns Jacobian of the constraints at design value \a u.
/// Format of \a result is sparse, complying to \a m_conJacRows
/// and \a m_conJacCols
template<typename T>
void OptParameterization<T>::jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{}

/// \brief Returns Hessian Lagrangian of the constraints at design value
/// \a u
/*template<typename T>
void OptParameterization<T>::hessLagr_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{GISMO_NO_IMPLEMENTATION }
*/

} // end namespace gismo
