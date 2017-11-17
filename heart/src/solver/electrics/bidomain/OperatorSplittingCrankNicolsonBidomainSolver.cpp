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

#include "OperatorSplittingCrankNicolsonBidomainSolver.hpp"
#include "BidomainAssembler.hpp"
#include "BidomainWithBathAssembler.hpp"
#include "PetscMatTools.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OperatorSplittingCrankNicolsonBidomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem != NULL)
    {
        return;
    }

    AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(initialSolution);

    // initialise matrix-based RHS vector and matrix, and use the linear
    // system rhs as a template
    Vec& r_template = this->mpLinearSystem->rGetRhsVector();
    VecDuplicate(r_template, &mVecForConstructingRhs1);
    VecDuplicate(r_template, &mVecForConstructingRhs2);
    PetscInt ownership_range_lo;
    PetscInt ownership_range_hi;
    VecGetOwnershipRange(r_template, &ownership_range_lo, &ownership_range_hi);
    PetscInt local_size = ownership_range_hi - ownership_range_lo;
    PetscTools::SetupMat(mMassMatrix, 2*this->mpMesh->GetNumNodes(), 2*this->mpMesh->GetNumNodes(), 
                         2*this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                         local_size, local_size); 
    PetscTools::SetupMat(mAiMatrix, 2*this->mpMesh->GetNumNodes(), 2*this->mpMesh->GetNumNodes(), 
                         2*this->mpMesh->CalculateMaximumNodeConnectivityPerProcess(),
                         local_size, local_size);

}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OperatorSplittingCrankNicolsonBidomainSolver<ELEMENT_DIM,SPACE_DIM>::SetupLinearSystem(
        Vec currentSolution,
        bool computeMatrix)
{
    assert(this->mpLinearSystem->rGetLhsMatrix() != NULL);
    assert(this->mpLinearSystem->rGetRhsVector() != NULL);
    assert(currentSolution != NULL);

    /////////////////////////////////////////
    // set up LHS matrix (and mass matrix and stiffness matrix)
    /////////////////////////////////////////
    if (computeMatrix)
    {
        this->mpBidomainAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        this->mpBidomainAssembler->AssembleMatrix();

        // the BidomainMassMatrixAssembler deals with the mass matrix
        // for both bath and nonbath problems 
        assert(SPACE_DIM==ELEMENT_DIM);
        BidomainMassMatrixAssembler<SPACE_DIM> mass_matrix_assembler(this->mpMesh);
        mass_matrix_assembler.SetMatrixToAssemble(mMassMatrix);
        mass_matrix_assembler.Assemble();

        BidomainStiffnessMatrixAssembler<ELEMENT_DIM,SPACE_DIM> stiffness_matrix_assembler(this->mpMesh, this->mpBidomainTissue);
        stiffness_matrix_assembler.SetMatrixToAssemble(mAiMatrix);
        stiffness_matrix_assembler.Assemble();

        this->mpLinearSystem->SwitchWriteModeLhsMatrix();
        PetscMatTools::Finalise(mMassMatrix);
        PetscMatTools::Finalise(mAiMatrix);
    }


    HeartEventHandler::BeginEvent(HeartEventHandler::ASSEMBLE_RHS);

    //////////////////////////////////////////  
    //Set up b = M z1 + Ai z2
    /////////////////////////////////////////  
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();

    // dist stripe for the current Voltage
    DistributedVector distributed_current_solution = p_factory->CreateDistributedVector(currentSolution);
    DistributedVector::Stripe distributed_current_solution_vm(distributed_current_solution, 0);

    // dist stripe for z1
    DistributedVector dist_vec_matrix_based_z1 = p_factory->CreateDistributedVector(mVecForConstructingRhs1);
    DistributedVector::Stripe dist_vec_z1_matrix_based_vm(dist_vec_matrix_based_z1, 0);
    DistributedVector::Stripe dist_vec_z1_matrix_based_phie(dist_vec_matrix_based_z1, 1);

    double Am = HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio();
    double Cm  = HeartConfig::Instance()->GetCapacitance();
    
    // dist stripe for z2
    DistributedVector dist_vec_matrix_based_z2 = p_factory->CreateDistributedVector(mVecForConstructingRhs2);
    DistributedVector::Stripe dist_vec_z2_matrix_based_vm(dist_vec_matrix_based_z2, 0);
    DistributedVector::Stripe dist_vec_z2_matrix_based_phie(dist_vec_matrix_based_z2, 1);

    if(!(this->mBathSimulation))
    {
        for (DistributedVector::Iterator index = dist_vec_matrix_based_z1.Begin();
             index!= dist_vec_matrix_based_z1.End();
             ++index)
        {
            double V = distributed_current_solution_vm[index];
	    //In Operator Splitting, the nodal ionic current and stimuli are computed with the cell system
           // double F = - Am*this->mpBidomainTissue->rGetIionicCacheReplicated()[index.Global]
           //            - this->mpBidomainTissue->rGetIntracellularStimulusCacheReplicated()[index.Global];

            dist_vec_z1_matrix_based_vm[index] = Am*Cm*V*PdeSimulationTime::GetPdeTimeStepInverse(); //+ F;
            dist_vec_z1_matrix_based_phie[index] = 0.0;
        }
        
        for (DistributedVector::Iterator index = dist_vec_matrix_based_z2.Begin();
             index!= dist_vec_matrix_based_z2.End();
             ++index)
        {
            double V = distributed_current_solution_vm[index];
            dist_vec_z2_matrix_based_vm[index] = -0.5*V;
            dist_vec_z2_matrix_based_phie[index] = -V;
        }
    }
    else
    {
	printf("Error, Crank-Nicolson operator splitting not set up for a bath.\n");
	    exit(-1);
    }

    dist_vec_matrix_based_z1.Restore();
    dist_vec_matrix_based_z2.Restore();
    
    //Create a temporary vector
    Vec temp;
    VecDuplicate(mVecForConstructingRhs1, &temp);

    //////////////////////////////////////////
    // b = M z1 + Ai z2
    //////////////////////////////////////////
    MatMult(mMassMatrix, mVecForConstructingRhs1, temp);
    //Add M z1 + Ai z2
    MatMultAdd(mAiMatrix, mVecForConstructingRhs2, temp, this->mpLinearSystem->rGetRhsVector());

    // assembling RHS is not finished yet, as Neumann bcs are added below, but
    // the event will be begun again inside this->mpBidomainAssembler->AssembleVector(); 
    HeartEventHandler::EndEvent(HeartEventHandler::ASSEMBLE_RHS);

    /////////////////////////////////////////
    // apply Neumann boundary conditions
    /////////////////////////////////////////
    mpBidomainNeumannSurfaceTermAssembler->ResetBoundaryConditionsContainer(this->mpBoundaryConditions); // as the BCC can change
    mpBidomainNeumannSurfaceTermAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
    mpBidomainNeumannSurfaceTermAssembler->AssembleVector();


    /////////////////////////////////////////
    // apply correction term
    /////////////////////////////////////////
    if(mpBidomainCorrectionTermAssembler)
    {
        mpBidomainCorrectionTermAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false/*don't zero vector!*/);
        // don't need to set current solution
        mpBidomainCorrectionTermAssembler->AssembleVector();
    }

    this->mpLinearSystem->FinaliseRhsVector();

    this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

    if(this->mBathSimulation)
    {
        this->mpLinearSystem->FinaliseLhsMatrix();
        this->FinaliseForBath(computeMatrix,true);
    }

    if(computeMatrix)
    {
        this->mpLinearSystem->FinaliseLhsMatrix();
    }
    this->mpLinearSystem->FinaliseRhsVector();
    //this->mpLinearSystem->DisplayMatrix();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OperatorSplittingCrankNicolsonBidomainSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec currentSolution)
{
    double time = PdeSimulationTime::GetTime();
    double dt = PdeSimulationTime::GetPdeTimeStep();
    //The true here is to solve the whole cell system (including the voltage)
    this->mpBidomainTissue->SolveCellSystems(currentSolution, time, time+dt/2.0, true);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OperatorSplittingCrankNicolsonBidomainSolver<ELEMENT_DIM,SPACE_DIM>::FollowingSolveLinearSystem(Vec currentSolution)
{
    // solve cell models for second half timestep
    double time = PdeSimulationTime::GetTime();
    double dt = PdeSimulationTime::GetPdeTimeStep();
    //The true here is to solve the whole cell system (including the voltage)
    this->mpBidomainTissue->SolveCellSystems(currentSolution, time + dt/2.0, time+dt, true);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OperatorSplittingCrankNicolsonBidomainSolver<ELEMENT_DIM,SPACE_DIM>::OperatorSplittingCrankNicolsonBidomainSolver(
        bool bathSimulation,
        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
        BidomainTissue<SPACE_DIM>* pTissue,
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions)
    : AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>(bathSimulation,pMesh,pTissue,pBoundaryConditions)
{

    // Tell tissue there's no need to replicate ionic caches
    pTissue->SetCacheReplication(false);
    mVecForConstructingRhs1 = NULL;
    mVecForConstructingRhs2 = NULL;

    // create assembler
     if(bathSimulation)
     {
    	 std::cout<<"not set up for CN operator splitting."<<std::endl;
         mpBidomainAssembler = new BidomainWithBathAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpBidomainTissue);
     }
     else
     {
         mpBidomainAssembler = new BidomainAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpBidomainTissue,BidomainAssembler_helper::CRANKNICOLSON);
     }


     mpBidomainNeumannSurfaceTermAssembler = new BidomainNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM>(pMesh,pBoundaryConditions);

     if(HeartConfig::Instance()->GetUseStateVariableInterpolation())
     {
         mpBidomainCorrectionTermAssembler
             = new BidomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,this->mpBidomainTissue);
         //We are going to need those caches after all
         pTissue->SetCacheReplication(true);
     }
     else
     {
         mpBidomainCorrectionTermAssembler = NULL;
     }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OperatorSplittingCrankNicolsonBidomainSolver<ELEMENT_DIM,SPACE_DIM>::~OperatorSplittingCrankNicolsonBidomainSolver()
{
    delete mpBidomainAssembler;
    delete mpBidomainNeumannSurfaceTermAssembler;

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
    
    if(mpBidomainCorrectionTermAssembler)
    {
        delete mpBidomainCorrectionTermAssembler;
    }
}

///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class OperatorSplittingCrankNicolsonBidomainSolver<1,1>;
template class OperatorSplittingCrankNicolsonBidomainSolver<2,2>;
template class OperatorSplittingCrankNicolsonBidomainSolver<3,3>;

