#include "OptParametrization.h"

#include <limits>

#include <gsTensor/gsGridIterator.h>

namespace gismo
{
    
template< int d , typename T >
T determinant(const gsMatrix<T, d, d> &J)
{
    return J.determinant();
}    

template< int d, typename T>
T areaOrthogonalityMeasure( const gsMatrix<T,d,d>& J )
{
    gsMatrix<T,d,d> g = J.transpose()*J;
    return g.diagonal().prod();
}

template< int d, typename T>
T liaoFunctional( const gsMatrix<T,d,d>& J )
{
    gsMatrix<T,d,d> g = J.transpose()*J;
    return g.squaredNorm();
}

template< int d, typename T>
T winslowFunctional( const gsMatrix<T,d,d>& J )
{
    gsMatrix<T,d,d> g = J.transpose()*J;
    return g.trace() / g.determinant();
}

template< int d, typename T>
T detMeasure( const gsMatrix<T,d,d>& J )
{
    return pow( J.determinant() - 1.  , 2 );
}

template<typename T>
T contMechanics( const gsMatrix<T,2,2>& J )
{
    // reference: Handbook of grid generation, 
    // ed. Joe F. Thompson, Bharat K. Soni, Nigel P. Weatherill
    T c_1,c_2;
    c_1 = T(0.);
    c_2 = T(1.);
    gsMatrix<T,2,2> g = J.transpose()*J;
    T det = J.determinant();
    return c_1*( g.trace() - 1. - 2.*det ) + c_2*( det - 1. )*( det-1. ); 
}
    

template< int d, typename T>
OptParametrization<d,T>::OptParametrization():
    m_functional( liaoFunctional<d,T> ),
    m_bPolar(false)
{
}

//** default destructor */
template< int d, typename T>
OptParametrization<d,T>::~OptParametrization()
{}


template< int d, typename T>
void OptParametrization<d,T>::setParametrization( const gsTensorBSpline<d,T>& parametrization )
{
    GISMO_ASSERT( d == parametrization.geoDim() , "OptParametrization: geoDim must match parDim" );
    
    //m_bPolar = polar;
    m_param = parametrization;
    
    for( int k = 0; k < d ; ++k )
        m_controlNetDim(k) = m_param.knots(k).size() - m_param.basis().degree(k) - 1;
    
    const T inf = std::numeric_limits<T>::infinity(); 
    
    // each control point has d components
    this->m_numDesignVars  = d;
        
    
    // next we multiply the number of inner knots in each dimension
    for( int k = 0 ; k < d ; ++k )
        this->m_numDesignVars *= m_controlNetDim(k)-2;

    this->m_desLowerBounds = gsMatrix<T>::Constant( this->m_numDesignVars , 1 , -inf );
    this->m_desUpperBounds = gsMatrix<T>::Constant( this->m_numDesignVars , 1 ,  inf );        
        
    // build up basis functions for the determinant
    
    m_det = gsTensorBSpline<d,T>( m_param.basis() , gsMatrix<T>::Zero( m_param.coefs().rows() , 1 ) );
    
    for( index_t k = 0; k < d ; ++k )
        m_det.basis().degreeIncrease( m_det.basis().degree(k)-1 , k );
    
    m_det.basis().anchor_into( m_detAnchorPoints );
    m_det.basis().collocationMatrix( m_detAnchorPoints , m_detCollocationMatrix );
    
    
    
    // init quadrature rule    
    m_quadRule.setNodes( m_positiveDetConstrainsPerElements * gsVector<int>::Ones(d,1) );
    

    
    //if( m_positiveDetConstrainsPerElements > 0 )
    //{            
        unsigned numberOfControlPoints      = 1;
        unsigned numberOfInnerControlPoints = 1;
        
        // find the number of total knots (including the boundary)
        for( int k = 0 ; k < d ; ++k )
        {
            numberOfControlPoints      *= m_controlNetDim(k);
            numberOfInnerControlPoints *= m_controlNetDim(k) - 2;
        }
                
        this->m_numConstraints =  pow(m_positiveDetConstrainsPerElements,d) * (numberOfControlPoints - numberOfInnerControlPoints);
        this->m_numConJacNonZero = 0;
        this->m_conLowerBounds = gsMatrix<T>::Constant( this->m_numConstraints , 1 ,   0. );
        this->m_conUpperBounds = gsMatrix<T>::Constant( this->m_numConstraints , 1 ,  inf );
        
        m_detEvaluationPoints.resize( d , this->m_numConstraints );
        

        gsTensorDomainIterator<T,d> iter( m_param.basis() );
        gsMatrix<T> nodes;
        gsVector<T> weights;
        
        int i = 0;
        for(  ; iter.good() ; iter.next() )
        {
            m_quadRule.mapTo( iter.lowerCorner() , iter.upperCorner() , nodes , weights );
            
            for( int j = 0 ; j < nodes.cols() ; ++j )
            {
                m_detEvaluationPoints.col(i++) = nodes.col(j);
            }
        }

    //}
    //else
    //{
    //   this->m_numConstraints =   0;
    //    this->m_numConJacNonZero = 0;
    //    this->m_conLowerBounds = gsMatrix<T>::Constant( this->m_numConstraints , 1 , -inf );
    //    this->m_conUpperBounds = gsMatrix<T>::Constant( this->m_numConstraints , 1 ,  inf );
    //}
        
    
    
    // we need to map all control points, onto a 1-D vector
    // and the other way around, the vector back to the control points
    
    gsGridIterator<unsigned,gsGridIteratorMode::CUBE,d> iter( gsVector<unsigned,d>::Ones() , m_controlNetDim - gsVector<unsigned,d>::Ones() );
    
    this->m_curDesign.resize( this->m_numDesignVars , 1);

    for( index_t i = 0 ; iter ; ++iter , ++i )
    {
        this->m_curDesign.block(d*i,0,d,1) = m_param.coef( fromTensorIndex( *iter , m_controlNetDim ) ).transpose();
        gsInfo  << "iter: " << *iter << "\n";
    }
}

//TODO: Use move schematics
template< int d, typename T>
const gsTensorBSpline<d,T>& OptParametrization<d,T>::getParametrization()
{
    updateParametrization( gsAsConstVector<T>( this->m_curDesign.data() , this->m_curDesign.cols()* this->m_curDesign.rows() ) );
    
    return m_param;
}

template< int d, typename T>
void OptParametrization<d,T>::updateParametrization(const gsAsConstVector<T>& u ) const
{   
    gsGridIterator<unsigned,gsGridIteratorMode::CUBE,d> iter( gsVector<unsigned,d>::Ones() , m_controlNetDim - gsVector<unsigned,d>::Ones() );
    
    for( index_t i = 0 ; iter ; ++iter , ++i )
    {
        m_param.coef( fromTensorIndex( *iter , m_controlNetDim ) ) = u.block(d*i,0 , d,1).transpose();
    }
}

template< int d , typename T>
T OptParametrization<d,T>::integrateFunctional( GeometricFunctional fnc ) const
{
    gsTensorDomainIterator<T,d> iter( m_param.basis() );
    gsMatrix<T> nodes;
    gsMatrix<T> jac;
    gsVector<T> funcValues( m_quadRule.numNodes() );
    gsVector<T> weights;
    
    T integral = T(0.);
    
    for(  ; iter.good() ; iter.next() )
    {
        m_quadRule.mapTo( iter.lowerCorner() , iter.upperCorner() , nodes , weights );
        
        for( int i = 0 ; i < nodes.cols() ; ++i )
            funcValues(i) = fnc( m_param.deriv( nodes.col(i) )->reshape(d,d) );
        
        integral += funcValues.dot( weights );
    }
    
    return integral;
}


template< int d, typename T>
T OptParametrization<d,T>::evalObj( const gsAsConstVector<T> & u ) const
{
    // we overwrite the control points of param with the values form 'u' and then compute the integral.
    updateParametrization(u);
    return integrateFunctional( m_functional );
}

/*
/// \brief Returns the gradient of the objective function at design value
/// \a u
/// By default it uses finite differences, overriding it should provide exact gradient.
template<typename T>
void OptParametrization<T>::gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{
    const_cast<gsTensorNurbs<2,T>&>(m_param).setCoefs( gsAsMatrix<T>( (double*)u.data() , m_param.coefs().rows() , m_param.coefs().cols() ) );
    
    result(0) = 2*m_param.coefs()(0,0);
}
*/

/// \brief Returns values of the constraints at design value \a u
template< int d, typename T>
void OptParametrization<d,T>::evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{
    
    /*
    if( m_positiveDetConstrainsPerElements > 0 )
    {
        updateParametrization(u);
        
        const gsMatrix<T>& points = this->m_detEvaluationPoints;
        
        for( index_t i = 0 ; i != points.cols() ; ++i )
            result(i) = m_param.deriv( points.col(i) )->reshape(d,d).determinant();
    }*/
}

/// \brief Returns Jacobian of the constraints at design value \a u.
/// Format of \a result is sparse, complying to \a m_conJacRows
/// and \a m_conJacCols
template< int d, typename T>
void OptParametrization<d,T>::jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{
    // at the moment only a full finite difference matrix is returned.
}

/// \brief Returns Hessian Lagrangian of the constraints at design value
/// \a u
/*template<typename T>
void OptParametrization<T>::hessLagr_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
{GISMO_NO_IMPLEMENTATION }
*/

} // end namespace gismo
