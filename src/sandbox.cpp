#include <gismo.h>
#include <Eigen/Dense>
#include <iostream>

#include "algorithms.h"

using namespace std;
using namespace Eigen;
using namespace gismo;

int main( int , char** )
{
    gsKnotVector<real_t> kv(0 , 1 , 10 , 4 );
    
    gsBSplineBasis<real_t> basis(kv);
    
    int n = 14; 
    gsMatrix<> coefs( n , 2);
    gsMatrix<> times( 1 , n);
    
    for (index_t i = 0; i < n ; ++i) {
        double t = i/(n-1.);
        times(i) = t;
        coefs(i, 0) = sin(2. * M_PI * t);
        coefs(i, 1) = cos(2. * M_PI * t);
    }
    
    gsBSpline<> curve( kv , coefs );
    
    gsWriteParaview( curve , "atemp/atmpfile");
    
    double dt = 0.01;
    
    for( int step = 0; step < 100; ++step )
    {
        gsMatrix<>::uPtr hess = curve.deriv2( times );
        
        curve.coefs() += dt * hess->transpose();
        
        gsWriteParaview( curve , "atemp/atmpfile" + std::to_string(step) );
    }
    
    return 0;
}


/*
 t e*mplate< typename T>
 class gsParameterization
 {
 
 public:
     
     void loadBoundaryFromFile( const std::string& fileName );
     void setBoundary( const gsBSpline<> boundary );
     
     const gsBSpline<>& getBoundary() const ;
     gsBSpline<>& getBoundary();
     
     private:
         
         gsBSpline<>        m_boundary;
         gsTensorNurbs<2> m_surface;
         
         
         
         };
         
         #include "gsIpopt/gsOptProblem.h"
         
         template< typename T>
         class gsOpt : public gsOptProblem<T>
         {
         public:
             gsMatrix<T> A;
             gsVector<T> b;
             T           c;
             
             
             public:
                 gsOpt(int n , gsMatrix<T> _A , gsVector<T> _b , T _c):
                 gsOptProblem<T>(),
                 A(_A),
                 b(_b),
                 c(_c) 
                 {
                 this->m_numDesignVars = n;
                 this->m_numConstraints = 0;
                 this->m_numConJacNonZero = 0;
                 
                 const T inf = std::numeric_limits<T>::infinity(); //or should be use max instead?        
                 this->m_desLowerBounds = gsMatrix<T>::Constant( this->m_numDesignVars , 1 , 0. );
                 this->m_desUpperBounds = gsMatrix<T>::Constant( this->m_numDesignVars , 1 ,  inf );
                 this->m_conLowerBounds = gsMatrix<T>::Constant( this->m_numDesignVars , 1 , -inf );
                 this->m_conUpperBounds = gsMatrix<T>::Constant( this->m_numDesignVars , 1 ,  inf );
                 
                 this->m_curDesign = gsMatrix<T>::Constant( n , 1 , T(0) );
                 
                 }
                 
                 virtual ~gsOpt(){}
                 
                 virtual T evalObj( const gsAsConstVector<T>& u ) const
                 {
                 //return T(0.5) * (u.transpose() * A * u).value()  + (b.transpose() * u).value() + c;
                 return sin( 2.* u.value() );
                 }
                 
                 virtual void gradObj_into( const gsAsConstVector<T>& u , gsAsVector<T>& result ) const
                 {
                 //result = A * u + b;
                 
                 result(0,0) = 2.* cos( u.value() );
                 }
                 
                 virtual void evalCon_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const {}
                 
                 virtual void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const {}
                 
                 virtual void hessLagr_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const {}
                 
                 
                 };
                 */
