#include <gismo.h>
#include <iostream>

#include "gsOptParametrization.h"

using namespace std;
using namespace Eigen;
using namespace gismo;

template< int d , typename T>
T myFunc( const gsMatrix<T,d,d>& jac )
{
    return jac.determinant();
}

int main( int , char** )
{
    auto param = gsNurbsCreator<real_t>::BSplineCube(3);
    param->uniformRefine(1);
    
    gsOptParametrization<3, real_t > opt;
    
    opt.setInitialParametrization( *param );
    opt.setFunctional( myFunc<3,real_t> );    
    opt.solve();
    
    auto paramOpt = opt.getParametrization();
    
    
    gsMesh<> controlNet, controlNetOpt; 
    param->controlNet( controlNet ); 
    paramOpt.controlNet( controlNetOpt );
    
    gsWriteParaview( controlNet    , "cube_net_000" );
    gsWriteParaview( controlNetOpt , "cube_net_001" );
    
    return 0;
}
