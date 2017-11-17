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
 NOVEMBER 2016
 */

#include "ThirdOrderPalindromeOperatorSplittingMonodomainSolver.hpp"
#include "MonodomainAssembler.hpp"
#include "TimeStepper.hpp"
//#include "PetscMatTools.hpp"
//#include "PdeSimulationTime.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ThirdOrderPalindromeOperatorSplittingMonodomainSolver<ELEMENT_DIM, SPACE_DIM>::SetupLinearSystem(
                                                                                            Vec currentSolution,
                                                                                            bool computeMatrix)
{
    
    // first linear system to be solved
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);
    
    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix and stiffness matrix)
    /////////////////////////////////////////
    if(computeMatrix)
    {
        mpMonodomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        mpMonodomainAssembler->AssembleMatrix();
        
        
        MassMatrixAssembler<ELEMENT_DIM,SPACE_DIM> mass_matrix_assembler(this->mpMesh, HeartConfig::Instance()->GetUseMassLumping());
        mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        mass_matrix_assembler.Assemble();
        
        MonodomainStiffnessMatrixAssembler<ELEMENT_DIM,SPACE_DIM> stiffness_matrix_assembler(this->mpMesh, this->mpMonodomainTissue);
        stiffness_matrix_assembler.SetMatrixToAssemble(mAiMatrix);
        stiffness_matrix_assembler.Assemble();
        
        this->mpLinearSystem->FinaliseLhsMatrix();
        PetscMatTools::Finalise(mMassMatrix);
        PetscMatTools::Finalise(mAiMatrix);
    }
    
    
    //Assemble and solve the first linear system
    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);
    
    //////////////////////////////////////////
    //Set up b = M z1
    /////////////////////////////////////////
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    
    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
    
    // dist stripe for z1
    DistributedVector dist_vec_matrix_based_z1 = p_factory->CreateDistributedVector(mVecForConstructingRhs1);
    
    
    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm  = HeartConfig::Instance()->GetCapacitance();
    
    
    for (DistributedVector::Iterator index = dist_vec_matrix_based_z1.Begin();
         index!= dist_vec_matrix_based_z1.End();
         ++index)
    {
        double V = distributed_current_solution[index];
        // in the main solver, the nodal ionic current and stimuli is computed and used.
        // However in operator splitting, this part of the solve is diffusion only, no reaction terms
        //double F = - Am*this->mpMonodomainTissue->rGetIionicCacheReplicated()[index.Global]
        //           - this->mpMonodomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global];
        
        //double time_linear_system = PdeSimulationTime::GetPdeTimeStepInverse();
        
        //solve for delta_t/6
        //dist_vec_matrix_based_z1[index] = Am*Cm*V*time_linear_system*1/6;
        
        dist_vec_matrix_based_z1[index] = Am*Cm*V*(PdeSimulationTime::GetPdeTimeStepInverse()*(1/6));
        
    }
    
    dist_vec_matrix_based_z1.Restore();
    
    
    
    //First RHS to be solved for the SDIRK2O3 method
    MatMult(mMassMatrix, mVecForConstructingRhs1, this->mpLinearSystem->rGetRhsVector());
    
    HeartEventHandler::EndEvent(HeartEventHandler::ASSEMBLE_RHS);
    
    /////////////////////////////////////////
    // apply Neumann boundary conditions
    /////////////////////////////////////////
    
    mpNeumannSurfaceTermsAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
    mpNeumannSurfaceTermsAssembler->AssembleVector();
    
    
    if(computeMatrix)
    {
        this->mpLinearSystem->FinaliseLhsMatrix();
    }
    
    this->mpLinearSystem->FinaliseRhsVector();
    
    this->FinaliseLinearSystem(currentSolution);
    
    this->mpLinearSystem->AssembleFinalLinearSystem();
    
    ////////////////////////////////
    
    Vec solution = currentSolution;
    
    Vec intermediate_solution;
    
    if (computeMatrix)
    {
        this->mpLinearSystem->ResetKspSolver();
    }
    intermediate_solution = this->mpLinearSystem->Solve(solution);
    
    // Avoid memory leaks
    if (solution != currentSolution)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        //VecDestroy(solution);
        PetscTools::Destroy(solution);
        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
    }
    
    //Set up the second linear system
    Vec intermediateSolution = intermediate_solution;
    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution_inter = p_factory->CreateDistributedVector(intermediateSolution);
    
    
    // dist stripe for z2
    DistributedVector dist_vec_matrix_based_z2 = p_factory->CreateDistributedVector(mVecForConstructingRhs2);
    
    
    
    for (DistributedVector::Iterator index = dist_vec_matrix_based_z2.Begin();
         index!= dist_vec_matrix_based_z2.End();
         ++index)
    {
        double V_star = distributed_current_solution_inter[index];
        
        dist_vec_matrix_based_z2[index] = -1*((1.0-2*mGamma)*V_star);
    }
    
    dist_vec_matrix_based_z2.Restore();
    
    //Create a temporary vector
    Vec temp;
    VecDuplicate(mVecForConstructingRhs1, &temp);
    
    // b = M z1 + Ai z2
    //////////////////////////////////////////
    MatMult(mMassMatrix, mVecForConstructingRhs1, temp);
    //Add M z1 + Ai z2
    MatMultAdd(mAiMatrix, mVecForConstructingRhs2, temp, this->mpLinearSystem->rGetRhsVector());
    
    
    
    /////////////////////////////////////////
    // apply Neumann boundary conditions
    /////////////////////////////////////////
    // mpBidomainNeumannSurfaceTermAssembler->ResetBoundaryConditionsContainer(this->mpBoundaryConditions); // as the BCC can change
    mpNeumannSurfaceTermsAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
    mpNeumannSurfaceTermsAssembler->AssembleVector();
    
    
    //finalise
    this->mpLinearSystem->FinaliseRhsVector();
    this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);
    
    
    
    this->mpLinearSystem->FinaliseRhsVector();
    // VecDestroy(intermediate_solution);
    PetscTools::Destroy(intermediate_solution);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ThirdOrderPalindromeOperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec currentSolution)
{
    
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ThirdOrderPalindromeOperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::FollowingSolveLinearSystem(Vec currentSolution)
{
    
}


//ok so the linear system should be assembled using SDIRK203, at least for the forward time marching

//Define my Solve function here
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Vec ThirdOrderPalindromeOperatorSplittingMonodomainSolver<ELEMENT_DIM, SPACE_DIM>::SolveOS(Vec currentSolution)

{
    double time = PdeSimulationTime::GetTime();
    double dt = PdeSimulationTime::GetPdeTimeStep();
     bool SolvingCellSystem[6] = {TRUE, FALSE, TRUE,FALSE, TRUE, FALSE};
     double ThirdOrderTimeSteps[6] = {0.919661523017399857,0.268330095781759925,-0.187991618799159782,-0.187991618799159782, 0.268330095781759925, 0.919661523017399857};
    int BackwardTimeIntegrationIndex = 12;
    
    
    Vec solution=currentSolution;
    Vec next_solution;
    
    //Start solving the substeps
    for (TimeSubStepIndex =0; TimeSubStepIndex < 6; ++TimeSubStepIndex)
    {
        //If we are solving the cell system
        if (SolvingCellSystem[TimeSubStepIndex])
        {
            
            //The true here is to solve the whole cell system (including the voltage)
            this->mpMonodomainTissue->SolveCellSystems(solution, time, time+(ThirdOrderTimeSteps[TimeSubStepIndex]*dt), true);
            
            // We only need to clean up memory if we are NOT on the first PDE time step,
            // as someone else cleans up the mInitialCondition vector in higher classes.
            //                if (solution != initialSolution)
            //                {
            //                    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
            //                    PetscTools::Destroy(solution);
            //                    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
            //                }
            
            
        }
        // Otherwise we solve the solve the PDEs system
        
        //find a way to set the time here
        else
        {
            if (TimeSubStepIndex == BackwardTimeIntegrationIndex)
            {
                bool compute_matrix = TRUE;
                
                this->SetupLinearSystem(solution, compute_matrix);
                
                this->FinaliseLinearSystem(solution);
                
                if (compute_matrix)
                {
                    this->mpLinearSystem->ResetKspSolver();
                }
            }
            next_solution = this->mpLinearSystem->Solve(solution);
            
        }
    }
    time = time+dt;
    
    // Avoid memory leaks
    //        if (solution != initialSolution)
    //        {
    //            HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    //            PetscTools::Destroy(solution);
    //            HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
    //        }
    solution = next_solution;
    
    return solution;
    
    
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ThirdOrderPalindromeOperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }
    
    // call base class version...
    AbstractLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>::InitialiseForSolve(initialSolution);
    
    //..then do a bit extra
    if(HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
        this->mpLinearSystem->SetAbsoluteTolerance(HeartConfig::Instance()->GetAbsoluteTolerance());
    }
    else
    {
        NEVER_REACHED;
        // re-implement when needed
        //this->mpLinearSystem->SetRelativeTolerance(HeartConfig::Instance()->GetRelativeTolerance());
    }
    
    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());
    this->mpLinearSystem->SetMatrixIsSymmetric(true);
    this->mpLinearSystem->SetUseFixedNumberIterations(HeartConfig::Instance()->GetUseFixedNumberIterationsLinearSolver(), HeartConfig::Instance()->GetEvaluateNumItsEveryNSolves());
    
    // initialise matrix-based RHS vector and matrix, and use the linear
    // system rhs as a template
    Vec& r_template = this->mpLinearSystem->rGetRhsVector();
    VecDuplicate(r_template, &mVecForConstructingRhs1);
    VecDuplicate(r_template, &mVecForConstructingRhs2);
    PetscInt ownership_range_lo;
    PetscInt ownership_range_hi;
    VecGetOwnershipRange(r_template, &ownership_range_lo, &ownership_range_hi);
    PetscInt local_size = ownership_range_hi - ownership_range_lo;
    PetscTools::SetupMat(mMassMatrix, this->mpMesh->GetNumNodes(), this->mpMesh->GetNumNodes(),
                         this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                         local_size, local_size);
    PetscTools::SetupMat(mAiMatrix, this->mpMesh->GetNumNodes(), this->mpMesh->GetNumNodes(),
                         this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                         local_size, local_size);
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ThirdOrderPalindromeOperatorSplittingMonodomainSolver<ELEMENT_DIM, SPACE_DIM>::ThirdOrderPalindromeOperatorSplittingMonodomainSolver(
                                                                                                                 AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                                                                                                 MonodomainTissue<ELEMENT_DIM, SPACE_DIM>* pTissue,
                                                                                                                 BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions )    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>(pMesh),
mpBoundaryConditions(pBoundaryConditions), mpMonodomainTissue(pTissue)
{
    assert(pTissue);
    assert(pBoundaryConditions);
    this->mMatrixIsConstant = true;
    
    mpMonodomainAssembler = new MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpMonodomainTissue);
    mpNeumannSurfaceTermsAssembler = new NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,1>(pMesh,pBoundaryConditions);
    
    // Tell tissue there's no need to replicate ionic caches
    pTissue->SetCacheReplication(false);
    mVecForConstructingRhs1 = NULL;
    mVecForConstructingRhs2 = NULL;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ThirdOrderPalindromeOperatorSplittingMonodomainSolver<ELEMENT_DIM, SPACE_DIM>::~ThirdOrderPalindromeOperatorSplittingMonodomainSolver()
{
    delete mpMonodomainAssembler;
    delete mpNeumannSurfaceTermsAssembler;
    
    if(mVecForConstructingRhs1)
    {
        PetscTools::Destroy(mVecForConstructingRhs1);
        if(mVecForConstructingRhs2)
        {
            PetscTools::Destroy(mVecForConstructingRhs2);
        }
        PetscTools::Destroy(mMassMatrix);
        PetscTools::Destroy(mAiMatrix);
    }
}

/////////////////////////////////////////////////////
//explicit instantiation
/////////////////////////////////////////////////////

template class ThirdOrderPalindromeOperatorSplittingMonodomainSolver<1,1>;
template class ThirdOrderPalindromeOperatorSplittingMonodomainSolver<1,2>;
template class ThirdOrderPalindromeOperatorSplittingMonodomainSolver<1,3>;
template class ThirdOrderPalindromeOperatorSplittingMonodomainSolver<2,2>;
template class ThirdOrderPalindromeOperatorSplittingMonodomainSolver<3,3>;
