#!/bin/bash
#Jessica Cervi
#University of Saskatchewan
#September 2016
#Script to generate reference solutions using Chaste 3.4 for LR1 cell model



#Change here accordingly

#Save the current directory location
startingDir=`pwd`

#Set Chaste location
CHASTE_LOCATION=/Users/jessicacervi/Chaste

#Set the location of the testing file (this needs to be a project in Chaste)
FILE_LOCATION=projects/ThirdOrderOS/test

#Set the folder location for the solution data output generated by the script (doesn't need to be changed)
export CHASTE_TEST_OUTPUT=/Users/jessicacervi/Chaste/3DAKSBidomainLR1

#Set the name of the C++ file to be generated
fileName=3DAKSBidomainLR1
fileNameCplusplus=$fileName.hpp
#Set the name of the results file displaying digits matched
fileResults=Results_3DAKSBidomainLR1.txt

#Maximum number of digits that are displayed in the results file
numDigits=10

#Number of time steps (we are outputting the solution every 0.1 ms)
numTimeSteps=21

#The timesteps we want to test (currently PDE and ODE timesteps are the same)
#change back to 10-4 order
dt="1e-2"
#The mesh sizes we want to test (must be n*countRef-1 in size, for n a non-negative integer). Same mesh spacing
#is used for x and y spacing

#change back to 1001
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
timeDuration=40
#Length of the 1D problem
xLength=1.0
#For some of the cell models available in Chaste, see heart/build/debug/src/odes/cellml
cellModelfileName=LuoRudy1991
cellModelConstructor=CellLuoRudy1991FromCellML


#ODE solver method
#For some of the ODE methods available in Chaste, see odes/src/solver
odeSolver=RungeKutta3IvpOdeSolver

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
#include \"BidomainProblem.hpp\"
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

class BlockCellFactory : public AbstractCardiacCellFactory<3>
{
private:
boost::shared_ptr<ZeroStimulus> mpStimulus;

public:
BlockCellFactory()
: AbstractCardiacCellFactory<3>(),
mpStimulus(new ZeroStimulus())
{
}

AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
{
double x = pNode->rGetLocation()[0];
double y = pNode->rGetLocation()[1];
double z = pNode->rGetLocation()[2];

boost::shared_ptr<$odeSolver> p_solver(new $odeSolver());
AbstractCardiacCell* cardiac_cell =  new $cellModelConstructor(p_solver, this->mpZeroStimulus);
double voltage = cardiac_cell->GetVoltage();
cardiac_cell->SetVoltage(voltage+100*(1-sin(x*y*z)));
return cardiac_cell;
}
};



//Class to generate the solutions
class $fileName : public CxxTest::TestSuite
{
public:
void TestFindReferenceSolutions_example() throw(Exception)
{
//Length of mesh
double xLength=$xLength;
double yLength = xLength;
double zLength = xLength;
//Grid spacing
double numPoints=$count1;
double h = xLength/numPoints;
//ODE timestep (must be less than or equal to dtPDE)
double dtODE=$timeStepODE;
//PDE timestep
double dtPDE=$timeStepPDE;
HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);

//Duration of simulation in ms (only solving first half)
double timeDuration=$timeDuration;
double numTimeSteps=$numTimeSteps -1;
//Timestep to output
double dtOut = (timeDuration*2)/numTimeSteps;

HeartConfig::Instance()->SetSimulationDuration(timeDuration*2); //ms
HeartConfig::Instance()->SetOutputFilenamePrefix(\"results\");
//Output with 16 digits
HeartConfig::Instance()->SetVisualizerOutputPrecision(16);


// Reference solution with h1 x spacing
HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dtODE, dtPDE, dtOut);
DistributedTetrahedralMesh<3,3> mesh1;
mesh1.ConstructRegularSlabMesh(h, xLength, yLength, zLength);
HeartConfig::Instance()->SetOutputDirectory(\"Reference_results_h1\");
BlockCellFactory<3> cell_factory1;


// Third-Order Operator splitting
HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplitting();
HeartConfig::Instance()->SetUseReactionDiffusionThirdOrderPalindromeOperatorSplittingBidomainSolver();



BidomainProblem<3> bidomain_problem1( &cell_factory1 );
bidomain_problem1.SetMesh(&mesh1);
bidomain_problem1.Initialise();
bidomain_problem1.Solve();

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


done
done







