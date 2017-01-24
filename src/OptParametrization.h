#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsAssembler/gsGaussRule.h>

#include "gsIpopt/gsOptProblem.h"

#pragma once

namespace gismo
{
    
template< int d, typename T>
T determinant( const gsMatrix<T,d,d>& J );

template< int d, typename T>
T areaOrthogonalityMeasure( const gsMatrix<T,d,d>& J );

template< int d, typename T>
T liaoFunctional( const gsMatrix<T,d,d>& J );

template< int d, typename T>
T winslowFunctional( const gsMatrix<T,d,d>& J );

template< typename T>
T contMechanics( const gsMatrix<T,2,2>& J );

template< int d, typename T>
T detMeasure( const gsMatrix<T,d,d>& J );


// 
// Warning: Currently the class works only for T being double (not float!)
// since it uses the class gsOptProblem, which has the same restriction.
//
// We need geoDim == d
template< int d , typename T  > 
class OptParametrization : public  gsOptProblem<T>
{
public:
    typedef T(*GeometricFunctional)( const gsMatrix<T,d,d>& );
    
public:
    //** default constructor */
    OptParametrization();    
    
    //** default destructor */
    virtual ~OptParametrization();
    
    
    //gsTensorBSpline<d,T> optimize( const gsTensorBSpline<d,T>& chart , polar = false );
    
public:
    
    void setParametrization( const gsTensorBSpline<d,T>& parametrization );
    
    const gsTensorBSpline<d,T>& getParametrization();
        
private:
    void updateParametrization(const gsAsConstVector<T>& u ) const;
public:
    
    void setFunctional( GeometricFunctional fnc ){ m_functional = fnc; }
    const GeometricFunctional& getFunctional() const { return m_functional; }
    
    T integrateFunctional( GeometricFunctional fnc ) const;
    
    /*void forcePositiveDeterminate( bool value = true ) 
    { 
        forcePositiveDeterminate( value ? 2 : 0);
    }
    
    void forcePositiveDeterminate( int value ) 
    { 
        m_positiveDetConstrainsPerElements = value; 
    }*/
    
    //bool isPolar() const;
    //bool setPolar( bool polar = true );
    
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
    void calculateDetTensorBSpline();
    
private:
    mutable gsTensorBSpline<d,T>  m_param;
    mutable gsTensorBSpline<d,T>  m_det;
    
    GeometricFunctional m_functional;
       
    
    gsGaussRule<T> m_quadRule;
    gsVector<unsigned,d> m_controlNetDim;
    
    //int  m_positiveDetConstrainsPerElements;
    //gsMatrix<T> m_detEvaluationPoints;
    
    gsMatrix<T> m_detCollocationMatrix;
    gsVector<T> m_detAnchorPoints;
    
    bool m_bCounterclockwise;
    bool m_bPolar;
};


}

// implementation 
#include "OptParametrization.hpp"
