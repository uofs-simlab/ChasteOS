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
 November 2016
 */

#ifndef THIRDORDERPALINDROMEOPERATORSPLITTINGMONODOMAINSOLVER_HPP
#define THIRDORDERPALINDROMEOPERATORSPLITTINGMONODOMAINSOLVER_HPP


#include "MonodomainTissue.hpp"
#include "MonodomainAssembler.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "MassMatrixAssembler.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"
#include "MonodomainStiffnessMatrixAssembler.hpp"

/**
 Fix description
 */


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//might have a problem here, FROM WHERE AM I INHERITING????
class ThirdOrderPalindromeOperatorSplittingMonodomainSolver : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>
{
private:
    
    /** Coefficient for SDIRK2O2 method */
    double mGamma;
    
    /** Boundary conditions */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* mpBoundaryConditions;
    
    /** Monodomain tissue class (collection of cells, and conductivities) */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* mpMonodomainTissue;
    
    /** The monodomain assembler, used to set up the LHS matrix */
    MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainAssembler;
    
    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
    NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,1>* mpNeumannSurfaceTermsAssembler;
    
    /** The mass matrix, used to computing the RHS vector */
    Mat mMassMatrix;
    
    /** Ai matrix, used to computing the RHS vector */
    Mat mAiMatrix;
    
    /**
     *  The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b, this vector is z1 where b = M z1 + Ai z2.
     */
    Vec mVecForConstructingRhs1;
    
    /**
     *  The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b, this vector is z2 where b = M z1 + Ai z2.
     */
    Vec mVecForConstructingRhs2;
    
    
    /**
     *  Implementation of SetupLinearSystem() which uses the assembler to compute the
     *  LHS matrix, but sets up the RHS vector using the mass-matrix (constructed
     *  using a separate assembler) multiplied by a vector
     *
     *  @param currentSolution  Solution at current time
     *  @param computeMatrix  Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);
    
    /**
     *  Called before setting up the linear system, used to solve the cell models for first half timestep (step (i) above)
     *  @param currentSolution the latest solution vector
     */
    void PrepareForSetupLinearSystem(Vec currentSolution);
    
    /**
     *  Called after solving the linear system, used to solve the cell models for second half timestep (step (iii) above)
     *  @param currentSolution the latest solution vector (ie the solution of the linear system).
     */
    void FollowingSolveLinearSystem(Vec currentSolution);
    
    
    //New part, for third-order OS
    /**
     * Called when we are ready to solve the splitted system
     */
    Vec SolveOS(Vec currentSolution);
    
    
    
    
    /**
     * Whether we are solving the cell system or not
     */
    bool SolvingCellSystem;
    
    /**
     * When using third-order OS, the substeps for the ODE and PDE system have to be the following
     */
    double ThirdOrderTimeSteps;
    
    /**
     * Index for the time-step array
     */
    int TimeSubStepIndex;
    
    /**
     * Index for the whether solving the cell system of not
     */
    int BackwardTimeIntegrationIndex;
    
    /*
     *Stepper
     */
    TimeStepper stepper();
    
    //end new part
    
    
public:
    
    /** Overloaded InitialiseForSolve() which calls base version but also
     *  initialises #mMassMatrix and #mVecForConstructingRhs.
     *
     *  @param initialSolution  initial solution
     */
    void InitialiseForSolve(Vec initialSolution);
    
    
    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    ThirdOrderPalindromeOperatorSplittingMonodomainSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                                MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                                                BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions);
    
    /**
     *  Destructor
     */
    ~ThirdOrderPalindromeOperatorSplittingMonodomainSolver();
};




#endif /* THIRDORDERPALINDROMEOPERATORSPLITTINGMONODOMAINSOLVER_HPP_ */
