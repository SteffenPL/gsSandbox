#include <cmath>
#include <iostream>
#include <memory>

#include <gismo.h>

#include <boost/filesystem.hpp>

#include "OptParameterization.h"

#include "algorithms.h"

using namespace gismo;


int main(int argc, char* argv[])
{
    std::string input("");
    std::string output("");
    std::string functional("");
    bool plot = false;
    int cIterations = 1;
    int iVariation = -1;
    bool bForcePositiveDeterminat = false;
    
    gsCmdLine cmd("Optimizes a given parameterization w.r.t. a functional.");
    cmd.addString("i", "input", "Name of the input file.", input);
    cmd.addString("o", "output", "Name of the output file.", output);
    cmd.addString("f", "func", "Name of the functional to optimize", functional);
    cmd.addInt("v", "var", "Generate a scalar field where the i.th control point is varied", iVariation);
    cmd.addSwitch("forcePositiveDet", "Force positive determinate as boundary condition.", bForcePositiveDeterminat );
    cmd.addSwitch("plot", "Plot surface using ParaView" , plot );
    cmd.addInt("c","iter","Number of iterations. (Used to generate intermediate outputs.)", cIterations );

    cmd.getValues(argc,argv);

    // load domain from input file
    gsFileData<> domainFile;
    domainFile.read( input );
    
    if( !domainFile.has< gsTensorNurbs<2> >() )
    {
        gsInfo  << "Cannot open file '"
        << input << "'.";
        return -1;
    }
    
    std::unique_ptr< gsTensorNurbs<2>> pNurbs( domainFile.getFirst< gsTensorNurbs<2> >() );
    
    if( pNurbs == NULL )
    {
        gsInfo  << "Cannot find a TensorNurbs<2,real_t> domain in file '"
        << input << "'.";
        return -1;
    }
    
    gsTensorNurbs<2>& surface = *pNurbs;
  
    
    // create folder for the output
    boost::filesystem::create_directory( boost::filesystem::path( output ) );
    std::string outputPath = output + "/" + output;
    
    // write basis to paraview file
    gsWriteParaview( surface.basis() , outputPath + "_basis" ); 
    
    
    // init optimization
    OptParameterization<real_t> opt;
    opt.setForcePositiveDeterminate( bForcePositiveDeterminat ); 
    
    if( functional == "areaOrth" )
        opt.setFunctional( areaOrthogonalityMeasure<real_t> );
    else if( functional == "liao" )
        opt.setFunctional( liaoFunctional<real_t> );
    else if( functional == "winslow" )
        opt.setFunctional( winslowFunctional<real_t> );
    else if( functional == "contMechanics" )
        opt.setFunctional( contMechanics<real_t> );
    else if( functional == "det" )
        opt.setFunctional( detMeasure<real_t> );
    else
    {
        cIterations = -1;
        gsInfo << "No functional given (-f).\n";
    }
    
    // we want to plot the effect of variation
    if( iVariation > -1 )
        {
        // original surface with small pertubations
        gsTensorNurbs<2,real_t> dsurface = surface; 


        // codim(2) , dim(3) surface with the value of the functuanal as additional component    
        gsTensorNurbs<2,real_t> variationField = dsurface; 
            
        variationField.uniformRefine(4);
        gsMatrix<real_t> greville = variationField.basis().anchors();
        variationField.coefs().conservativeResize( Eigen::NoChange , variationField.coefs().cols()+1 );
            
        for( index_t i = 0 ; i < greville.cols() ; ++i )
        {
            dsurface.coef(iVariation,0) = greville(0,i);
            dsurface.coef(iVariation,1) = greville(1,i);
            variationField.coef(i,2) = integrateFunctional( dsurface , opt.getFunctional() ); 
        }

        saveParameterization( variationField , outputPath , "_functional" );
    }
    
    // save intermediate results
    saveParameterization( surface , outputPath , "_0" );
    
    for( int i = 0 ; i < cIterations ; ++i )
    {
        
        opt.setParameterization( surface );
        opt.solve();
        surface = opt.getParameterization();
        
        // save intermediate results
        saveParameterization( surface , outputPath , "_" + std::to_string(i+1) );
    }
    
    // save without postfix
    saveParameterization( surface , outputPath );
    
    return 0;
}
