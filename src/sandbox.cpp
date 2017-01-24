#include <gismo.h>

#include <iostream>

#include "OptParametrization.h"

using namespace std;
using namespace Eigen;
using namespace gismo;

template< typename T>
T myFunc( const gsMatrix<T,2,2>& jac )
{
    return jac.determinant();
}

int main( int , char** )
{
    //auto param = gsNurbsCreator<real_t>::BSplineLShape_p2C0();
    auto param = gsNurbsCreator<real_t>::BSplineCube(3);
    
    gsFileData<real_t> fd;
    fd << *param;
    fd.save("cube");
    
    param->uniformRefine(1);
    
    OptParametrization<3, real_t > opt;
    
    opt.setParametrization( *param );
    
    gsInfo << "Area : " << opt.integrateFunctional( winslowFunctional<3,real_t> );
    
    gsWriteParaview( *param , "Cube" , 1000 );
    
    opt.solve();
    
    auto param_new = opt.getParametrization();
    
    gsWriteParaview( param_new , "Cube_new" , 1000 );
    
    gsFileData<real_t> fd2;
    fd2 << *param;
    fd2.save("cube_new");
    
    gsMesh<> controlNet;
    param->controlNet( controlNet );    
    gsWriteParaview( controlNet , "cube_net" );
    
    return 0;
}
