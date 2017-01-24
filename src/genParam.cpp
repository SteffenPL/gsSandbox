#include <cmath>
#include <iostream>
#include <memory>

#include <gismo.h>

//#include "extensions/gsIpopt/gsOptProblem.h"

#include "gsOptParametrization.h"

#include <boost/filesystem.hpp>

#include "algorithms.h"

using namespace gismo;


int main(int argc, char* argv[])
{
    int n = 8;
    int m = 8;
    int degree = 3;
    bool plot = false;
    bool bOpt = false;
    bool bRectangle = false;
    bool plotDet = false;
    real_t start = 0.;

    std::string input("");
    std::string output("");
    std::string algorithm("min");
    bool nonPeriodicBoundary = false;

    int nDetValues = 20;

    gsCmdLine cmd("Tutorial on gsTensorBSpline class.");
    cmd.addInt   ("n", "n", "Number of basis function in one direction"  , n);
    cmd.addInt   ("m", "m", "Number of basis function in other direction", m);
    cmd.addInt   ("d", "degree", "Degree of a surface", degree);
    cmd.addString("i", "input", "Name of the input file.", input);
    cmd.addString("o", "output", "Name of the output file.", output);
    cmd.addString("a", "algorithm", "Algorithms: bil , min , ..." , algorithm );
    cmd.addSwitch("plot", "Plot surface using ParaView" , plot );
    cmd.addSwitch("opt", "Optimize the surface parameterization using IOPT" , bOpt );
    cmd.addSwitch("plotDet", "Append the determinate as a third coordinate to the surface NURBS" , plotDet );
    cmd.addReal("s" , "start" , "The surface map should map the point (0,0) onto boundary(start)" , start );
    cmd.addSwitch("rectangle" , "Consider the boundary curve to have corners at t={0,0.25,0.5,0.75}" , bRectangle );
    cmd.addSwitch("nonPeriodic" , "Do not force the boundary to be periodic, i.e. no call of .setPeriodic()" , nonPeriodicBoundary );
    cmd.addInt( "c" , "det" , "number of det values ((todo: better desc))" , nDetValues );
    cmd.getValues(argc,argv);

    
    
    // Adjust values to the minimum required
    degree = math::max(0, degree    );
    n      = math::max(n, degree + 1);
    m      = math::max(m, degree + 1);
    gsInfo << "----------------------\n\n"
              << "n: " << n << "\n\n"
              << "m: " << m << "\n\n"
              << "degree: " << degree << "\n\n"
              << "input: "  << input  << "\n\n"
              << "output: " << output << "\n\n"
              << "----------------------\n\n";

    // 0. input: Border Nurbs curve:
    // finally we shoud read in some curve from a file,
    // since tools to create appropriate files are not avaliable for free
    // we just create a simple border curve here
    int curve_n   = 10;
    int curve_degree = 2;

    std::unique_ptr<gsBSpline<real_t>> pCurve;

    //findKernelPoint( curve.coefs() );

    if( input == "" )
        pCurve = loadBoundaryCurve( curve_n , curve_degree );
    else
        pCurve = loadBoundaryFromFile( input );

    gsBSpline<>& curve = *pCurve;

    if( !nonPeriodicBoundary )
        curve.setPeriodic();

    //curve.basis().enforceOuterKnotsPeriodic();

    gsTensorBSpline<2,real_t> surface;

    if( algorithm == "bil" )
        surface = bilinearlyBlendedCoons( curve , n , m , degree , bRectangle );
    else if( algorithm == "min" )
        surface = minimumPrinciple( curve , n , m , degree , start , bRectangle );
    else if( algorithm == "naivePolar" )
        surface = naivePolarChart( curve , m , degree , bRectangle );
    else if( algorithm == "uni" )
        surface = unidirectionalInterpolation( curve , n , m , degree , bRectangle );

    // 5. saving surface, basis and control net to a file
    if (output != "")
    {
        
        // create folder for output stuff (too many files)
        boost::filesystem::create_directory( boost::filesystem::path( output ) );
        std::string prefix = output + "/" + output;
        
        std::string surfaceOutput = prefix;
        std::string curveOutput = prefix + "Curve";
        std::string detOutput = prefix + "Det";

        gsFileData<> fd;
        fd << surface;
        fd.save( surfaceOutput );

        gsFileData<> fd2;
        fd2 << curve;
        fd2.save( curveOutput );


        // write det:
        gsJacDetField<real_t> detFnct( surface );
        gsMultiPatch<> patch(surface);
        
        gsField<> detField( patch ,detFnct );
        
        gsInfo  << "Writing the surface (with determinate coordinate) to a paraview file: "
        << detOutput << ".vts"
        << "\n\n";
        
        gsWriteParaview(detField, detOutput , 5000 );


        if( plot )
        {       
            gsInfo << "Writing the surface to a paraview file: "
            << surfaceOutput << ".vts"
            << "\n\n";

            gsMesh<> controlNet;
            surface.controlNet( controlNet );
            
            gsWriteParaview( controlNet , prefix + "Mesh" );
            
            gsInfo << "Writing the boundary curve to a paraview file: "
                << curveOutput << ".vtp"
                << "\n\n";

            gsWriteParaview( curve , curveOutput , 200 );
            //system( ("paraview " + curveOutput + ".vtp").c_str() );
}
    }

    return 0;
}
