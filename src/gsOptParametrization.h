#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsMatrix/gsSparseSolver.h>

#include "gsIpopt/gsOptProblem.h"

#include <Eigen/Sparse>

#pragma once

namespace gismo
{

/// \brief Computes the determinant of J.
/// \a J: dxd matrix
template< int d, typename T>
T determinant( const gsMatrix<T,d,d>& J );

/// \brief Computes the product of the metric tensors diagonal.
/// \a J: dxd matrix
template< int d, typename T>
T areaOrthogonalityMeasure( const gsMatrix<T,d,d>& J );

/// \brief Computes the Frobenius norm of the metric tensor.
/// \a J: dxd matrix
template< int d, typename T>
T liaoFunctional( const gsMatrix<T,d,d>& J );

/// \brief Computes the Winslow functional $ \frac{tr( G )}{\sqrt{det(G)}}$ for the metric tensor $G = J^T J$.
/// \a J: dxd matrix
template< int d, typename T>
T winslowFunctional( const gsMatrix<T,d,d>& J );

/// \brief Computes the cell deformation measure given in Chapter 33-9 of "The Handbook of Grid Generation" by Jeo F. Thompson, Bharat K. Soni, Nigel P. Weatherill. This functional measures is a stress term from elasticity theory.
/// \a J: dxd matrix
template< typename T>
T contMechanics( const gsMatrix<T,2,2>& J );

/** \brief For a given function $f : \mathbb{R}^{d\times d} \to \mathbb{R}$ we aim to minimize
 *  $\int f( \operatorname{D} S ) \operatorname{d} \xi$. Applications are the improvement of a
 *  given parametrization with respect to some measure.
 * 
 *  The functional will be optimized w.r.t. the inner control points of the parametrization.
 *  To ensure a locally bijective mapping, the inequality constrains $d_{ij} > 0$ for the 
 * coefficients of the B-Spline representation of $det \operatorname{D} S$ are imposed.
 * 
 *  A minimal example how to use this class:
 * 
 * gsTensorBSpline<d,T> spline = ...
 * gsOptParametrization<d,T> opt;
 * opt.setFunctional();
 * opt.setParametrization( spline );
 * spline = opt.solve();
 * 
 *  In addition the degree of quadrature can be changed.
 *  Options for IPOPT can be set in "<GISMO_DIR>/filedata/options/ipopt.opt".
 * 
 *  This class depends on the extension gsOptProblem, therefore the CMAKE option GISMO_WITH_IPOPT must be active.
 *  
 *  The interface is in an early state and could change in future.
 * 
 *  \tparam d: codimension (must be equal to the geometrical dimension)
 *  \tparam T: coefficient type (currently IPOPT uses double internally, this could lead to incompabilities!)
 */
template< int d , typename T  > 
class gsOptParametrization : public  gsOptProblem<T>
{
public:
    
    /// \brief User defined functionals have to be of this type.
    typedef T(*GeometricFunctional)( const gsMatrix<T,d,d>& );
    
public:
    //** default constructor */
    gsOptParametrization();    
    
    //** default destructor */
    virtual ~gsOptParametrization();
    
public:
    
    /// \brief Defines the initial parametrization. 
    ///  - The geometrical dimension must be equal to 'd'!
    ///  - Open boundary knots are assumed (if not the boundary curve may change)
    ///  - The map should be nearly orientation preserving.
    void setInitialParametrization( const gsTensorBSpline<d,T>& parametrization );
    
    /// \brief Returns the current state (or result) of the optimization.
    const gsTensorBSpline<d,T>& getParametrization();
        
private:
    
    /// \brief internal class. Copies the state vector into the TensorBSpline for further computations.
    void updateParametrization(const gsAsConstVector<T>& u ) const;
    
public:
    
    /// \brief 
    void setFunctional( GeometricFunctional fnc ){ m_functional = fnc; }
    const GeometricFunctional& getFunctional() const { return m_functional; }
    
    /// \brief Defines order of quadrature for the computation of the functional.
    /// The integration is done for each element.
    void setQuadratureDegreee( int degree = 3 );
    
    /// \brief Defines order of quadrature for each dimension of the computation of the functional.
    /// The integration is done for each element. 
    void setQuadratureDegreee( const gsVector<int> degree );
    
    /// \brief Returns the values of the 'fnc' integrated over the codomain of the current parametrization.
    /// The function will be evaluated the the Jacobians of the parametrization.
    //
    /// \a GeometricFunctional fnc: A function pointer of the type T (fnct*) ( const gsMatrix<T,d,d>& J )
    T integrateFunctional( GeometricFunctional fnc ) const;
    
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
    /// Currently implemented via finite differences.
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
    
    
    gsSparseMatrix<T>   m_detCollocationMatrix;
    gsEigenSparseLU<T>  m_solver;
    gsMatrix<T> m_detAnchorPoints;
    
    //bool m_bCounterclockwise;
    //bool m_bPolar;
};


}

// implementation 
#include "gsOptParametrization.hpp"
