#include <cmath>
#include <iostream>
#include <memory>

#include <gismo.h>

#include <boost/filesystem.hpp>

#include "gsOptParametrization.h"

#include "algorithms.h"

using namespace gismo;


int main(int argc, char* argv[])
{
    std::string input("");
    std::string output("");
    std::string functional("");
    std::string quadratureRule("gauss");
    bool plot = false;
    int cIterations = 1;
    int iVariation = -1;
    bool bForcePositiveDeterminat = false;
    bool bPolar = false;
    
    gsCmdLine cmd("Optimizes a given parameterization w.r.t. a functional.");
    cmd.addString("i", "input", "Name of the input file.", input);
    cmd.addString("o", "output", "Name of the output file.", output);
    cmd.addString("f", "func", "Name of the functional to optimize", functional);
    cmd.addString("q", "quadrule", "Quadrature rule for computing the functional", quadratureRule);
    cmd.addInt("v", "var", "Generate a scalar field where the i.th control point is varied", iVariation);
    cmd.addSwitch("forcePositiveDet", "Force positive determinate as boundary condition.", bForcePositiveDeterminat );
    cmd.addSwitch("polar","If the parameterization is similar to polar coordinates with the north side contracted onto one point",bPolar);
    cmd.addSwitch("plot", "Plot surface using ParaView" , plot );
    cmd.addInt("c","iter","Number of iterations. (Used to generate intermediate outputs.)", cIterations );

    cmd.getValues(argc,argv);

    // load domain from input file
    gsFileData<> domainFile;
    domainFile.read( input );
    
    std::unique_ptr< gsTensorBSpline<2>> pTensorBSpline;
    
    if( domainFile.has< gsTensorBSpline<2> >() )
    {
        pTensorBSpline = std::unique_ptr< gsTensorBSpline<2>>( domainFile.getFirst< gsTensorBSpline<2> >() );
        
        if( pTensorBSpline == NULL )
        {
            gsInfo  << "Cannot find a gsTensorBSpline<2,real_t> domain in file '"
            << input << "'.";
            return -1;
        }
        
    }
    else if( domainFile.has< gsTensorNurbs<2> >() )
    {
        std::unique_ptr< gsTensorNurbs<2> >pNurbs;
        
        pNurbs = std::unique_ptr< gsTensorNurbs<2>>( domainFile.getFirst< gsTensorNurbs<2> >() );
        
        if( pNurbs == NULL )
        {
            gsInfo  << "Cannot find a TensorBSpline<2,real_t> domain in file '"
            << input << "'.";
            return -1;
        }
        
        gsInfo << "Warning: Converting Tensor Nurbs to BSpline, i.e. ignoring all weights.\n";
        
        pTensorBSpline = std::unique_ptr< gsTensorBSpline<2> >( new gsTensorBSpline<2>( pNurbs->knots(0) , pNurbs->knots(1) , give( pNurbs->coefs()) ) ); 
    }
    else
    {
        gsInfo  << "Cannot open file  '"
        << input << "' (no TensorBSpline or TensorNurbs found).";
        return -1;
    }
    
    
    gsTensorBSpline<2>& surface = *pTensorBSpline; // TODO: remove later, possible memory leak if copied in a wrong way
    
    
    // create folder for the output
    boost::filesystem::create_directory( boost::filesystem::path( output ) );
    std::string outputPath = output + "/" + output;
    
    // write basis to paraview file
    gsWriteParaview( surface.basis() , outputPath + "_basis" ); 
    
    
    // init optimization
    gsOptParametrization<2,real_t> opt;
    //opt.forcePositiveDeterminate( bForcePositiveDeterminat ); 
        
    if( functional == "areaOrth" )
        opt.setFunctional( areaOrthogonalityMeasure<2,real_t> );
    else if( functional == "liao" )
        opt.setFunctional( liaoFunctional<2,real_t> );
    else if( functional == "winslow" )
        opt.setFunctional( winslowFunctional<2,real_t> );
    else if( functional == "contMechanics" )
        opt.setFunctional( contMechanics<real_t> );
    else if( functional == "det" )
        opt.setFunctional( detMeasure<2,real_t> );
    else
    {
        cIterations = -1;
        gsInfo << "No functional given (-f).\n";
    }
    
    // we want to plot the effect of variation
    /*if( iVariation > -1 )
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
    }*/
    
    // save intermediate results
    saveParameterization( surface , outputPath , "_0" );
    
    for( int i = 0 ; i < cIterations ; ++i )
    {
        
        opt.setInitialParametrization( surface );//, bPolar );
        opt.solve();
        surface = opt.getParametrization();
        
        // save intermediate results
        saveParameterization( surface , outputPath , "_" + std::to_string(i+1) );
    }
    
    // save without postfix
    saveParameterization( surface , outputPath );
    
    return 0;
}
