#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsTensorNurbs.h>

#include "gsIpopt/gsOptProblem.h"

#include <functional>

#pragma once

namespace gismo
{
    
    template<typename T>
    T determinant( const gsMatrix<T,2,2>& J );
    
    template<typename T>
    T areaOrthogonalityMeasure( const gsMatrix<T,2,2>& J );
    
    template<typename T>
    T liaoFunctional( const gsMatrix<T,2,2>& J );
    
    template<typename T>
    T winslowFunctional( const gsMatrix<T,2,2>& J );
    
    template<typename T>
    T contMechanics( const gsMatrix<T,2,2>& J );
    
    template<typename T>
    T detMeasure( const gsMatrix<T,2,2>& J );
    
    // deprecatred
    template<typename T>
    void appendFunctional( gsTensorNurbs<2,real_t>& surface , T(*F)(const gsMatrix<T,2,2>&) )
    {
        gsMatrix<real_t> greville = surface.basis().anchors();
        
        gsMatrix<real_t> jac;
        surface.deriv_into( greville , jac );
        
        surface.coefs().conservativeResize( Eigen::NoChange , surface.coefs().cols()+1 );
        
        for( index_t i = 0 ; i < greville.cols() ; ++i )
        {
            surface.coef(i,2) = F( gsAsMatrix<T>( jac.col(i).data() , 2, 2) );        
        }
    }
    
    template< typename T >
    T integrateFunctional( const gsTensorNurbs<2,T>& surface, T(*F)( const gsMatrix<T,2,2>&) , bool onParameterSpace = true);
    

template< typename T > 
class OptParameterization : public  gsOptProblem<T>
{
public:
    typedef T(*GeometricFunctional)( const gsMatrix<T,2,2>& );
    
public:
    //** default constructor */
    OptParameterization( bool forcePositiveDet = false );    
    
    //** default destructor */
    virtual ~OptParameterization();
    
public:
    
    void setParameterization( const gsTensorNurbs<2,T>& parameterization );
    const gsTensorNurbs<2,T>& getParameterization();
    
private:
    void updateParameterization(const gsAsConstVector<T>& u ) const;
public:
    
    void setFunctional( GeometricFunctional fnc ){ m_functional = fnc; }
    const GeometricFunctional& getFunctional() const { return m_functional; }
    
    void setForcePositiveDeterminate( bool value = true ) 
    { m_bForcePositiveDet = value; }
    
    /// \brief Returns the gradient value of the objective function at design
    /// value \a u
    virtual T evalObj( const gsAsConstVector<T> & u ) const;
    
    /// \brief Returns the gradient of the objective function at design value
    /// \a u
    /// By default it uses finite differences, overriding it should provide exact gradient.
    //virtual void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const;
    
    /// \brief Returns values of the constraints at design value \a u
    virtual void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const;
    
    /// \brief Returns Jacobian of the constraints at design value \a u.
    /// Format of \a result is sparse, complying to \a m_conJacRows
    /// and \a m_conJacCols
    virtual void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const;
    
    /// \brief Returns Hessian Lagrangian of the constraints at design value
    /// \a u
    virtual void hessLagr_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {GISMO_NO_IMPLEMENTATION }
    
    
private:
    mutable gsTensorNurbs<2,T>  m_param;
    unsigned int m_numElemX;
    unsigned int m_numElemY;
    
    GeometricFunctional m_functional;
    bool m_bForcePositiveDet;
};


}

// implementation 
#include "OptParameterization.hpp"
