/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/* 
Megan E. Marsh, Saeed Torabi Ziaratgahi, and Raymond J. Spiteri
Numerical Simulation Laboratory 
University of Saskatchewan 
Aril 2015
Partial support provided by research grants from the National Science and Engineering Research Council (NSERC) of Canada and the MITACS/Mprime Canadian Network of Centres of Excellence.
*/

#ifndef TESTOPERATORSPLITTINGGODUNOVBIDOMAINSOLVER_HPP_
#define TESTOPERATORSPLITTINGGODUNOVBIDOMAINSOLVER_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "BidomainProblem.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PropagationPropertiesCalculator.hpp"

/* 
Megan E. Marsh, Raymond J. Spiteri 
Numerical Simulation Laboratory 
University of Saskatchewan 
February 2012
Partial support provided by research grants from the National Science and Engineering Research Council (NSERC) of Canada and the MITACS/Mprime Canadian Network of Centres of Excellence.
*/

//Use a continuous initial condition of V_rest - 100*(1-sin(x))
template<unsigned DIM>
class BlockCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BlockCellFactory()
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
        assert(DIM<3);
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
	    boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver());
	    AbstractCardiacCell* cardiac_cell =  new CellLuoRudy1991FromCellML(p_solver, this->mpZeroStimulus);
	    double voltage = cardiac_cell->GetVoltage();
	    cardiac_cell->SetVoltage(voltage+100*(1-sin(x)));
	    return cardiac_cell;

    }
};

class TestOperatorSplittingGodunovBidomainSolver : public CxxTest::TestSuite
{
public:

    // The operator splitting and normal methods should agree closely with very small dt and h, but this takes
    // too long to run in the continuous build (see instead TestOperatorSplittingBiomainSolverLong)
    //
    // Here we run on a fine (as opposed to v fine) mesh and with a normal dt, and check that the solutions
    // are near.
    void TestComparisonOnNormalMeshes() throw(Exception)
    {
		ReplicatableVector final_voltage_normal;
		ReplicatableVector final_voltage_reference;
		ReplicatableVector final_voltage_operator_splitting;
		ReplicatableVector final_voltage_operator_splitting_halved;
		ReplicatableVector final_voltage_operator_splitting_normal_halved;

		HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
		HeartConfig::Instance()->SetOutputFilenamePrefix("results");

		double h = 0.01;
		// Reference solution
		HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(1e-5, 1e-5, 0.1);
		TetrahedralMesh<1,1> mesh;
		mesh.ConstructRegularSlabMesh(h, 1.0);
		HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplitting_reference");
		HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-6);
		BlockCellFactory<1> cell_factory;

		BidomainProblem<1> bidomain_problem( &cell_factory );
		bidomain_problem.SetMesh(&mesh);
		bidomain_problem.Initialise();
		bidomain_problem.Solve();

        final_voltage_reference.ReplicatePetscVector(bidomain_problem.GetSolution());

        //Reset ODE and PDE timesteps
        double OdeDt=0.002;
        double PdeDt=0.02;
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(OdeDt, PdeDt, 0.1);
        // Normal
        {
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, 1.0);
            HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplitting_normal");
            BlockCellFactory<1> cell_factory;
            BidomainProblem<1> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();

            final_voltage_normal.ReplicatePetscVector(bidomain_problem.GetSolution());
        }

        // Operator splitting
        {
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, 1.0);
            HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplitting_splitting");
            BlockCellFactory<1> cell_factory;

            HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplitting();

            BidomainProblem<1> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();
            final_voltage_operator_splitting.ReplicatePetscVector(bidomain_problem.GetSolution());
        }

        // Operator splitting with halved timestep
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(OdeDt/2, PdeDt/2, 0.1);
        {

            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, 1.0);
            HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplitting_splitting_halved");
            BlockCellFactory<1> cell_factory;

            HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplittingBidomainSolver();
	    HeartConfig::Instance()->SetUseGodunovReactionDiffusionOperatorSplittingBidomainSolver();

            BidomainProblem<1> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();
            final_voltage_operator_splitting_halved.ReplicatePetscVector(bidomain_problem.GetSolution());
        }
        {
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, 1.0);
            HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplitting_splitting_normal_halved");
            BlockCellFactory<1> cell_factory;

            HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplittingBidomainSolver();
	    HeartConfig::Instance()->SetUseGodunovReactionDiffusionOperatorSplittingBidomainSolver();

            BidomainProblem<1> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();
            final_voltage_operator_splitting_normal_halved.ReplicatePetscVector(bidomain_problem.GetSolution());
        }

	//Test for first-order convergence
	unsigned int solution_size = final_voltage_reference.GetSize();

	for(unsigned int index = 0;index<solution_size;index+=2)//points)
	{
		TS_ASSERT_DELTA(2,fabs(final_voltage_normal[index]-final_voltage_reference[index])/fabs(final_voltage_operator_splitting_normal_halved[index]-final_voltage_reference[index]),0.4);
		TS_ASSERT_DELTA(2,fabs(final_voltage_operator_splitting[index]-final_voltage_reference[index])/fabs(final_voltage_operator_splitting_halved[index]-final_voltage_reference[index]),0.4);
	}

        bool some_node_depolarised = false;
        assert(final_voltage_normal.GetSize()==final_voltage_operator_splitting.GetSize());
        for(unsigned j=0; j<final_voltage_normal.GetSize(); j++)
        {
            // this tolerance means the wavefronts are not on top of each other, but not too far
            // separated (as otherwise max difference between the voltages across space would be
            // greater than 80).
            double tol=25;

            TS_ASSERT_DELTA(final_voltage_normal[j], final_voltage_operator_splitting[j], tol);

            if(final_voltage_normal[j]>-80)
            {
                // shouldn't be exactly equal, as long as away from resting potential
                TS_ASSERT_DIFFERS(final_voltage_normal[j], final_voltage_operator_splitting[j]);
            }

            if(final_voltage_normal[j]>0.0)
            {
            	some_node_depolarised = true;
            }
        }
        assert(some_node_depolarised);
    }
};

#endif /* TESTOPERATORSPLITTINGGODUNOVBIDOMAINSOLVER_HPP_ */
