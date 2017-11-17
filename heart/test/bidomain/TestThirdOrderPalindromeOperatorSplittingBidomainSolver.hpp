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
Jessica Cervi and Raymond J. Spiteri
 Numerical Simulation Laboratory
 University of Saskatchewan
 Aril 2016
 Partial support provided by research grants from the National Science and Engineering Research Council (NSERC) of Canada and the MITACS/Mprime Canadian Network of Centres of Excellence.
 */

#ifndef TESTTHIRDORDERPALINDROMEOPERATORSPLITTINGBIDOMAINSOLVER_HPP_
#define TESTTHIRDORDERPALINDROMEOPERATORSPLITTINGBIDOMAINSOLVER_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "BidomainProblem.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "LuoRudy1991.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "HeunIvpOdeSolver.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include <algorithm>
#include <fstream>
#include <iterator>
#include <iostream>


// stimulate a block of cells (an interval in 1d, a block in a corner in 2d)
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
        boost::shared_ptr<HeunIvpOdeSolver> p_solver(new HeunIvpOdeSolver());
        AbstractCardiacCell* temp =  new CellLuoRudy1991FromCellML(p_solver, this->mpZeroStimulus);
        temp->SetVoltage(-85.533+100*(1-sin(x)));
        return temp;
    }
};


class TestThirdOrderPalindromeOperatorSplittingBidomainSolver : public CxxTest::TestSuite
{
public:
    
    //Test the CN method for nearness to the semi-implicit default method and for 2nd order convergence
    void TestComparisonOnNormalMeshes() throw(Exception)
    {
        ReplicatableVector final_voltage_normal;
        ReplicatableVector final_voltage_reference;
        ReplicatableVector final_voltage_operator_splitting;
        ReplicatableVector final_voltage_operator_splitting_halved;
        
        double time_duration=0.05;
        HeartConfig::Instance()->SetSimulationDuration(time_duration); //ms
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-6);
        double h = 0.01;
        double length = 1.0;
        //        		// Code used to generate Reference solution, saved for use in this test
        //            	HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(5e-8, 5e-8, time_duration);
        //            	TetrahedralMesh<1,1> mesh;
        //             	mesh.ConstructRegularSlabMesh(h, length);
        //             	HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplittingThirdOrder_reference");
        //             	BlockCellFactory<1> cell_factory;
        //
        //             BidomainProblem<1> bidomain_problem( &cell_factory );
        //             bidomain_problem.SetMesh(&mesh);
        //             bidomain_problem.Initialise();
        //             bidomain_problem.Solve();
        //
        //             final_voltage_reference.ReplicatePetscVector(bidomain_problem.GetSolution());
        
        
        double dt = 0.002;
        double PDE_timestep_factor=5;
        //Reset ODE and PDE timesteps
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dt, dt*PDE_timestep_factor, time_duration);
        // Normal
        {
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, length);
            HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplittingThirdOrder_normal");
            BlockCellFactory<1> cell_factory;
            BidomainProblem<1> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();
            
            final_voltage_normal.ReplicatePetscVector(bidomain_problem.GetSolution());
        }
        
        // Operator splitting
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dt, dt*PDE_timestep_factor, time_duration);
        {
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, length);
            HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplittingThirdOrder_splitting");
            BlockCellFactory<1> cell_factory;
            
            HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplittingBidomainSolver();
            HeartConfig::Instance()->SetUseReactionDiffusionThirdOrderPalindromeOperatorSplittingBidomainSolver();
            BidomainProblem<1> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();
            
            final_voltage_operator_splitting.ReplicatePetscVector(bidomain_problem.GetSolution());
        }
        
        // Operator splitting with halved timestep
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(dt/2, PDE_timestep_factor*dt/2, time_duration);
        {
            
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, length);
            HeartConfig::Instance()->SetOutputDirectory("BidomainCompareWithOperatorSplittingThirdOrder_splitting_halved");
            BlockCellFactory<1> cell_factory;
            HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplittingBidomainSolver();
            HeartConfig::Instance()->SetUseReactionDiffusionThirdOrderPalindromeOperatorSplittingBidomainSolver();
            BidomainProblem<1> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();
            final_voltage_operator_splitting_halved.ReplicatePetscVector(bidomain_problem.GetSolution());
        }
        
        //Reference solution, computed with commented code (dt=5e-8)
        double voltage_reference [100] =
        {
            16.20370975117123,
            15.82490467655516,
            14.94317330350061,
            13.92501077890292,
            12.90526460501189,
            11.88740316527865,
            10.86684904184003,
            9.844328812735791,
            8.820091751593694,
            7.794235532367327,
            6.766924163521837,
            5.738315090772364,
            4.708571680257752,
            3.677865036962222,
            2.646373571315049,
            1.614283150466309,
            0.5817864852026416,
            -0.4509176741089461,
            -1.483625408271839,
            -2.516128415143577,
            -3.548214603718919,
            -4.579668455320419,
            -5.610271150172954,
            -6.639800505283805,
            -7.66803079876112,
            -8.694732555165812,
            -9.719672383524253,
            -10.74261295237278,
            -11.7633131591765,
            -12.78152854803889,
            -13.79701201173075,
            -14.80951476073303,
            -15.81878756855522,
            -16.82458222541106,
            -17.82665313698597,
            -18.82475896679595,
            -19.81866419583048,
            -20.80814047426187,
            -21.7929675939258,
            -22.77293401587059,
            -23.74783687800583,
            -24.71748155968567,
            -25.68168098653246,
            -26.64025499143472,
            -27.59303002743267,
            -28.53983955922145,
            -29.48052518618569,
            -30.41493837896273,
            -31.34294248697753,
            -32.2644146145889,
            -33.17924693851479
            -34.08734730810113,
            -34.98863903002006,
            -35.88306000786559,
            -36.77056128505918,
            -37.65110463894776,
            -38.52467568935062,
            -39.39111435616871,
            -40.24997648637871,
            -41.10227254621776,
            -41.94774398141056,
            -42.7862392663989,
            -43.61773511932009,
            -44.44221091981328,
            -45.25965401565648,
            -46.07005588313255,
            -46.87340966813223,
            -47.66970852776827,
            -48.45894444389575,
            -49.24110739357462,
            -50.01618477438735,
            -50.78416103816804,
            -51.54501740442571,
            -52.29873174363007,
            -53.04527852420806,
            -53.7846288090979,
            -54.51675031611021,
            -55.2416075194106,
            -55.95916178437702,
            -56.66937152986366,
            -57.37219241233316,
            -58.06757752701714,
            -58.7554776208314,
            -59.43584131439355,
            -60.10861535620681,
            -60.77374480873649,
            -61.43117331148028,
            -62.08084331855668,
            -62.72269631382157,
            -63.35667306549475,
            -63.98271368497076,
            -64.60075834753307,
            -65.21074701928956,
            -65.81261366815758,
            -66.4063293031725,
            -66.99177511094047,
            -67.56858037999291,
            -68.13939240993179,
            -68.70184132767298,
            -69.18087946887874,
            -69.38276750917774
        };
        int check = 0;
        int solution_size = 202;//final_voltage_reference.GetSize();
        int ref_index=0;
        for(int index = 0;index<solution_size-2;index+=2)//points
        {
            
            //If the value hasn't changed, don't test for convergence
            if(fabs(final_voltage_operator_splitting[index]-final_voltage_operator_splitting_halved[index])>5e-4)
            {
                check++;
                TS_ASSERT_DELTA(4,fabs(final_voltage_operator_splitting[index]-voltage_reference[ref_index])/fabs(voltage_reference[ref_index]-final_voltage_operator_splitting_halved[index]),0.5);
                //TS_ASSERT_DELTA(4,fabs(final_voltage_operator_splitting[index]-final_voltage_reference[index])/fabs(final_voltage_reference[index]-final_voltage_operator_splitting_halved[index]),1.0);
                
            }
            ref_index++;
        }
        //Make sure at least 5 points were tested for convergence
        TS_ASSERT(check>5);
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

#endif /* TESTTHIRDORDERPALINDROMEOPERATORSPLITTINGBIDOMAINSOLVER_HPP_ */
