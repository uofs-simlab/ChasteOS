##Change accordingly


#Save the current directory location
startingDir=`pwd`

#Set Chaste location
CHASTE_LOCATION=/home/jcervi/release_3.4

#Set the location of the testing file (this needs to be a project in Chaste)
FILE_LOCATION=/projects/NiedererAKS/test

#Set the folder location for the solution data output generated by the script (doesn't need to be changed)
export CHASTE_TEST_OUTPUT=/home/jcervi/release_3.4/NiedererAKS

#Set the name of the C++ file to be generated
fileName=NiedererAKS
fileNameCplusplus=$fileName.hpp
#Set the name of the results file displaying digits matched
fileResults=NiedererAKS.txt

#Maximum number of digits that are displayed in the results file
numDigits=10

#Number of time steps (we are outputting the solution every 0.1 ms)
numTimeSteps=21

#The timesteps we want to test (currently PDE and ODE timesteps are the same)
dt="1e-3" # "5e-3 4e-3"
#The mesh sizes we want to test (must be n*countRef-1 in size, for n a non-negative integer). Same mesh spacing
#is used for x and y spacing
mesh="101" #  "101 201" 

for timeStep in `echo $dt` #For each timestep listed in dt
do
	for meshSize in `echo $mesh` #For each of the meshes listed
	do
		#Move to the starting directory
		cd $startingDir

		#Sizes of the 2 meshes and the reference mesh and mesh values minus 1, used for extracting data
		#Note that the mesh1 must be n*countRef-1 in size, for n a non-negative integer
		#Note that mesh2 is 2*count1 -1 in size (eg mesh1 has 101 points, mesh2 has 201 points)
		meshRef=101
		countRef=$(($meshRef-1))
		mesh1=$meshSize 
		count1=$(($mesh1-1))
		mesh2=20001
		count2=$(($mesh2-1)) 

		#These values are all set in the auto-generated C++ file
		#ODE and PDE timesteps
		timeStepODE=$timeStep
		timeStepPDE=$timeStep
		#Duration in time of the simulation
		timeDuration=40.0
		#Length of the 1D problem
		xLength=1.0 
		#For some of the cell models available in Chaste, see heart/build/debug/src/odes/cellml
		cellModelfileName=TenTusscher2006Epi #LuoRudy1991
		cellModelConstructor=CellTenTusscher2006EpiFromCellML #CellLuoRudy1991FromCellML
		#ODE solver method
		#For some of the ODE methods available in Chaste, see odes/src/solver
		odeSolver=RungeKutta4IvpOdeSolver #EulerIvpOdeSolver #RungeKutta2IvpOdeSolver

		#Save all of the values being used in the results file
		echo -e "
		Results auto-generated from runReferenceSolution1D_continuous_IC.sh
		ODE timestep = $timeStepODE
		PDE timestep = $timeStepPDE
		Simulation duration = $timeDuration
		1D interval = 0 to $xLength
		Cell model = $cellModelConstructor
		ODE method = $odeSolver
		Number of mesh points = $mesh1
		Initial conditions = Voltage = (voltage+100*(1-sin(x)))

		Results:
		" >> $fileResults

		#Move to the Chaste folder
		cd $CHASTE_LOCATION

		#Generate the C++ file
		##################################
		echo -e "/*

		Copyright (C) University of Oxford, 2005-2011

		University of Oxford means the Chancellor, Masters and Scholars of the
		University of Oxford, having an administrative office at Wellington
		Square, Oxford OX1 2JD, UK.

		This file is part of Chaste.

		Chaste is free software: you can redistribute it and/or modify it
		under the terms of the GNU Lesser General Public License as published
		by the Free Software Foundation, either version 2.1 of the License, or
		(at your option) any later version.

		Chaste is distributed in the hope that it will be useful, but WITHOUT
		ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
		FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
		License for more details. The offer of Chaste under the terms of the
		License is subject to the License being interpreted in accordance with
		English Law and subject to any action against the University of Oxford
		being under the jurisdiction of the English Courts.

		You should have received a copy of the GNU Lesser General Public License
		along with Chaste. If not, see <http://www.gnu.org/licenses/>.

		*/


		/** This file has been auto-generated by runReferenceSolution1D.sh */

		#ifndef TEST1DREFERENCESOLUTIONGENERATOREXAMPLE_HPP_
		#define TEST1DREFERENCESOLUTIONGENERATOREXAMPLE_HPP_

		#include <cxxtest/TestSuite.h>
		#include <boost/archive/text_oarchive.hpp>
		#include <boost/archive/text_iarchive.hpp>
		#include <vector>
		#include \"MonodomainProblem.hpp\"
		#include \"ZeroStimulusCellFactory.hpp\"
		#include \"AbstractCardiacCellFactory.hpp\"
		#include \"$cellModelfileName.hpp\"
		#include \"PlaneStimulusCellFactory.hpp\"
		#include \"TetrahedralMesh.hpp\"
		#include \"PetscTools.hpp\"
		#include \"PetscSetupAndFinalize.hpp\"
		#include \"PropagationPropertiesCalculator.hpp\"
		#include \"$odeSolver.hpp\"
		#include \"CardiacSimulationArchiver.hpp\"

		// stimulate a block of cells (an interval in 1d, a block in a corner in 2d)
		 template<unsigned DIM>
		class BlockCellFactory : public AbstractCardiacCellFactory<1>
		{

		private:
		    boost::shared_ptr<SimpleStimulus> mpStimulus;

		public:
		    BlockCellFactory()
		      : AbstractCardiacCellFactory<1>(),
			mpStimulus(new SimpleStimulus(-50000.0,2,0))
		    {
		    }

		    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
		    {
			//double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
			double x = pNode->rGetLocation()[0];

			boost::shared_ptr<$odeSolver> p_solver(new $odeSolver()); 

			if ( (x<0.15+1e-6))
			{

			    return new $cellModelConstructor(p_solver, mpStimulus);
			}
			else
			{
			    return new $cellModelConstructor(p_solver, mpZeroStimulus);
			} 


		    }
		};

		//Class to generate the solutions
		class $fileName : public CxxTest::TestSuite
		{
		public:
		    void TestFindReferenceSolutions_example() //throw(Exception)
		    {
		      	//Length of mesh
			double xLength=$xLength;
			//Grid spacing
			double numPoints=$count1;
			double h = xLength/numPoints;
			//ODE timestep (must be less than or equal to dtPDE)
			double dtODE=$timeStepODE;
			//PDE timestep
			double dtPDE=$timeStepPDE;

			//Duration of simulation in ms (only solving first half)
			double timeDuration=$timeDuration/2;
			double numTimeSteps=$numTimeSteps -1;
			//Timestep to output
			double dtOut = (timeDuration*2)/numTimeSteps;
	
			HeartConfig::Instance()->SetSimulationDuration(2*timeDuration); //ms
			//Output with 16 digits
			HeartConfig::Instance()->SetVisualizerOutputPrecision(16);
	
			// Reference solution with h1 x spacing
			HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dtODE, dtPDE, dtOut);
			DistributedTetrahedralMesh<1,1> mesh1;
			mesh1.ConstructRegularSlabMesh(h, xLength);

			HeartConfig::Instance()->SetOutputDirectory(\"Reference_results_h1\");
			HeartConfig::Instance()->SetOutputFilenamePrefix(\"results\");
		 	HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);

			BlockCellFactory<1> cell_factory1;

			//HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75));   // unit: mS/cm
      			//HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1.0));   // unit: mS/cm

			
			//HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400.0); // X
                        //HeartConfig::Instance()->SetCapacitance(1.0); // C_m

	 	  	// Godunov Operator splitting
        		HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplittingMonodomainSolver();
			HeartConfig::Instance()->SetUseReactionDiffusionThirdOrderPalindromeOperatorSplittingMonodomainSolver();




			MonodomainProblem<1> monodomain_problem1( &cell_factory1 );
			monodomain_problem1.SetMesh(&mesh1);
			monodomain_problem1.Initialise();
			monodomain_problem1.Solve(); 


			HeartEventHandler::Headings();
		        HeartEventHandler::Report();


		    }
		}; 


		#endif /* TEST1DREFERENCESOLUTIONGENERATOREXAMPLE_HPP_ */
		" > $CHASTE_LOCATION/$FILE_LOCATION/$fileNameCplusplus

		##################################

		#Build and execute the C++ file
                scons ts=$FILE_LOCATION/$fileNameCplusplus

		cd $startingDir
		#qsub runChastePBS_gidaspow_1D_TenTusscher_Monodomain.sh

	done
done





