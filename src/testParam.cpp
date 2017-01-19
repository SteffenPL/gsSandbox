#include <gismo.h>

#include <memory>

#include <boost/filesystem.hpp>


using namespace gismo;

int main( int argc , char** argv )
{
    std::string domainFilePath;
    std::string outputName = "poisson";
    bool bPolar = false;
    bool bPlot = false;
    
    int refine = 0;
    int cDegreeSteps = 1;

    std::string boundaryExpr = "sin(x) + cos(y)";
    std::string rhsExpr   = "sin(x) + cos(y)";
    std::string exactExpr = "sin(x) + cos(y)";
    std::string boundaryMethod = "elimination";

    std::string singularitySide = "w";
    
    
    // 0. Parse Input
    gsCmdLine cmd("Solves Poisson's equation with an isogeometric discretization.");

    cmd.addString( "i" , "input" , "Gismo .xml file containing the parameterization of the domain as a tensor nurbs curve." , domainFilePath );
    cmd.addString( "o" , "output" , "Output filename." , outputName );
    cmd.addSwitch( "polar" , "If the domain is parameterizised using a polar coordinate system, different boundary conditions apply." , bPolar );
    cmd.addString("s","singularity", "Side where the singularity of the polar map lies. (n,s,w,e)" , singularitySide );
    cmd.addSwitch("plot","Plot the solution using paraview." , bPlot );
    cmd.addString("b","boundary","Function expression in x,y to describe the boundary values", boundaryExpr );
    cmd.addString("r","rhs","Function expression in x,y to describe the rhs", rhsExpr );
    cmd.addInt("k","refine","Number of refinement steps.",refine);
    cmd.addInt("d","degree","Number of degree incrementations",cDegreeSteps );
    cmd.addString("e","exact","Function expression in x,y to describe the exact solution",exactExpr);
    cmd.addString("t","trick","Method used to introduce boundary conditions (nitsche,elimination)",boundaryMethod);

    cmd.getValues( argc , argv );

    // create folder for output stuff (too many files)
    boost::filesystem::create_directory( boost::filesystem::path( outputName ) );
    std::string prefix = outputName + "/" + outputName;
    
    
    // open output file
    std::ofstream globalErrorFile( prefix + "_error.csv");
    
    // use first line as a header
    const char csvDelimiter = ',';
    globalErrorFile   << "step"       << csvDelimiter
                << "degree"     << csvDelimiter
                << "errorL1"    << csvDelimiter
                << "errorL2"    << csvDelimiter
                << "errorLinf"  << std::endl;
    
    // 1. set up PDE
    
    // boundary conditions
    gsFunctionExpr<> exactSol( exactExpr , 2 );
    gsFunctionExpr<> boundaryFunction( boundaryExpr, 2 );
    gsFunctionExpr<> rhsFunction( rhsExpr , 2);
    
    // use the values from the exact solution by default as boundary values
    if( exactExpr != "" )
        boundaryFunction = gsFunctionExpr<>( exactExpr , 2 );
    
    gsFunction<>& boundaryValues = boundaryFunction;
    gsFunction<>& rhsValues = rhsFunction;
    
    // load domain from input file
    gsFileData<> domainFile;
    domainFile.read( domainFilePath );
    
    if( !domainFile.has< gsTensorNurbs<2> >() )
    {
        gsInfo  << "Cannot open file '"
        << domainFilePath << "'.";
        return -1;
    }
    
    std::unique_ptr< gsTensorNurbs<2>> pNurbs( domainFile.getFirst< gsTensorNurbs<2> >() );
    
    if( pNurbs == NULL )
    {
        gsInfo  << "Cannot find a TensorNurbs<2,real_t> domain in file '"
        << domainFilePath << "'.";
        return -1;
    }
    
    gsTensorNurbs<2>& nurbs = *pNurbs;
    
    for( int i = 0 ; i <= refine ; ++i )
    {  
        gsMultiPatch<> patch( nurbs );

        // create dirichlet boundary conditions for all boundaries
        gsBoundaryConditions<> bcInfo;

        gsConstantFunction<> zero( 0., 2);

        if( !bPolar )
        {
            for( gsMultiPatch<>::const_biterator
                bit = patch.bBegin() ; bit != patch.bEnd() ; ++bit )
            {
                bcInfo.addCondition( *bit
                        , condition_type::dirichlet
                        , &boundaryFunction );
            }
        }
        else
        {
            boundary::side singularity;
            boundary::side dirichletSide;
            boundary::side periodicSide1;
            boundary::side periodicSide2;
            
            
            if( singularitySide=="n" || singularitySide=="north" )
            {
                singularity = boundary::side::north;
                dirichletSide = boundary::side::south;
                periodicSide1 = boundary::side::east;
                periodicSide2 = boundary::side::west;
            }
            if( singularitySide=="s" || singularitySide == "south" )
            {
                singularity = boundary::side::south;
                dirichletSide = boundary::side::north;
                periodicSide1 = boundary::side::east;
                periodicSide2 = boundary::side::west;
            }
            if( singularitySide=="w" || singularitySide == "west" )
            {
                singularity = boundary::side::west;
                dirichletSide = boundary::side::east;
                periodicSide1 = boundary::side::north;
                periodicSide2 = boundary::side::south;
            }
            if( singularitySide=="e" || singularitySide == "east" )
            {
                singularity = boundary::side::east;
                dirichletSide = boundary::side::west;
                periodicSide1 = boundary::side::north;
                periodicSide2 = boundary::side::south;
            }
            
            bcInfo.addCondition( singularity 
                                , condition_type::neumann
                                , &zero );

            bcInfo.addCondition( dirichletSide
                                , condition_type::dirichlet
                                , &boundaryValues );

            patch.gsBoxTopology::addInterface(0 ,boxSide(periodicSide1) ,
                                            0 , boxSide(periodicSide2));
        }

        // init poisson pde
        gsPoissonPde<> pde( patch , bcInfo , rhsValues );

        // create basis for the element functions
        gsMultiBasis<> basis( patch );
              
        for( int p = 0 ; p < cDegreeSteps ; ++p )
        {
            // 2. setup solver

            gsPoissonAssembler<> assembler;
            gsOptionList options = assembler.defaultOptions();
            //Use Nitsche's method for Dirichlet boundaries
            if( boundaryMethod == "elimination" )
            {
                options.setInt("DirichletStrategy", dirichlet::elimination);
                gsInfo << "Using elimination for Dirichlet boundaries.\n";
            }
            else
            {
                options.setInt("DirichletStrategy", dirichlet::nitsche);
                gsInfo<<"Using Nitsche's method for Dirichlet boundaries.\n";
            }
            

            assembler.initialize( pde , basis , options );

            // generate matrix and load vector
            gsInfo << "Assembling...\n";
            assembler.assemble();


            // init cg solver
            gsInfo << "Solving...\n";

            gsEigenSparseLU<real_t> solver( assembler.matrix() );
            gsMatrix<> solVector = solver.solve( assembler.rhs()    );

            // construct the solution as a scalar field

            gsMultiPatch<> mpsol;
            assembler.constructSolution( solVector , mpsol );

            gsField<> sol( assembler.patches() , mpsol );


            gsInfo << "Sol: " << mpsol << "\n";


            // 3. return output

            
            gsInfo << "Writing to Paraview...\n";
            gsWriteParaview( sol , prefix + "_" + std::to_string( i ) , 2000 );
        
            
            gsMesh<> controlNet;
            nurbs.controlNet( controlNet );
            
            gsWriteParaview( controlNet , prefix + "_net_" + std::to_string(i) );
            
            
            if( !exactExpr.empty() )
            {
                gsAbsError<real_t> error( sol.function(0) , nurbs , exactSol );
                
                gsField<> errorField( assembler.patches() , error );
                
                gsInfo << "Plotting absolute error in Paraview...\n";
                gsWriteParaview( errorField , prefix + "_error_" + std::to_string( i ) , 2000 );
                
                gsNormL<1> norm1(   sol , exactSol );
                gsNormL<2> norm2(   sol , exactSol );
                gsNormL<0> normInf( sol , exactSol );
                
                
                globalErrorFile   << i        << csvDelimiter 
                            << basis.degree()       << csvDelimiter
                            << norm1.compute()      << csvDelimiter
                            << norm2.compute()      << csvDelimiter
                            << normInf.compute()    << std::endl;  
            }
        
            // 4. refinement
            
            
            //basis[0].degreeDecrease(1);
            //basis[1].degreeDecrease(1);
        }      
        
        nurbs.uniformRefine(1);
        
    }
    
    
    // 3. return output
    
    //write local file
    
    //plot solution in paraview
    int result = 0;
    
    if( bPlot )
        result = system( std::string("paraview " + prefix + "_" + std::to_string( refine ) + "0.vts &").c_str() );
    
    gsInfo << "Test is done. Exiting.\n";

    return 0;
}
