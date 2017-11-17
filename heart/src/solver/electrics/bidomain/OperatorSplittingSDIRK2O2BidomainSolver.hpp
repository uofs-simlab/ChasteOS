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

#ifndef OPERATORSPLITTINGSDIRK2O2BIDOMAINSOLVER_HPP_
#define OPERATORSPLITTINGSDIRK2O2BIDOMAINSOLVER_HPP_


#include "UblasIncludes.hpp"

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "AbstractBidomainSolver.hpp"
#include "HeartConfig.hpp"
#include "BidomainAssembler.hpp"
#include "BidomainMassMatrixAssembler.hpp"
#include "BidomainStiffnessMatrixAssembler.hpp"
#include "BidomainCorrectionTermAssembler.hpp"
#include "BidomainNeumannSurfaceTermAssembler.hpp"
#include "TimeStepper.hpp"


#include <boost/numeric/ublas/vector_proxy.hpp>

/**
 *  A bidomain solver that uses Strang operator splitting of the diffusion (conductivity) term and the reaction
 *  (ionic current) term, instead of solving the full reaction-diffusion PDE. This does NOT refer to operator splitting
 *  of the two PDEs in the bidomain equations. For details see for example Sundnes et al "Computing the Electrical
 *  Activity of the Heart". This solves the PDEs using the SDIRK2O2 method (second order).
 *
 *  The algorithm is, for solving from t=T to T+dt.
 *
 *  (i)   Solve ODEs   dV/dt = Iionic, du/dt = f(u,V)  for t=T to T+dt/2     
 *        [giving updated V (internally, and in solution vector) and updated state variables]
 *  (ii)  Solve PDE    dV/dt = div (sigma_i grad V) + div (sigma_i grad phi_e), div (sigma_i grad V) + div ((sigma_i+sigma_e) grad phi_e) = 0  
 *        for t=T to dt         [using V from step i, --> updated V]
 *  (iii) Solve ODEs   dV/dt = Iionic, du/dt = f(u,V)  for t=T+dt/2 to T+dt  [using V from step ii, --> final V]
 *
 *  Notes
 *   (a)  Stages (iii) and (i) can normally be solved together in one go, except just before/after printing the voltage to file.
 *        However for simplicity of code this has not been implemented
 *   (b)  Therefore, the effective ODE timestep will be:  min(ode_dt, pde_dt/2), where ode_dt and pde_dt are those
 *        given via HeartConfig.
 *   (c)  This solver is FOR COMPARING ACCURACY, NOT PERFORMANCE. It has not been optimised and may or may not
 *        perform well in parallel.
 *   (d)  Usage: add the follwoing lines to your original code:
 *        HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplittingBidomainSolver();
 * 	  HeartConfig::Instance()->SetUseSDIRK2O2ReactionDiffusionOperatorSplittingBidomainSolver(); 
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class OperatorSplittingSDIRK2O2BidomainSolver : public AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>
{
private:

    /** Coefficient for SDIRK2O2 method */
    double mGamma;
  
    /** Mass matrix, used to computing the RHS vector (actually: mass-matrix in
     *  voltage-voltage block, zero elsewhere
     */
    Mat mMassMatrix;
    
    /** Ai matrix, used to computing the RHS vector (actually: Ai-matrix in
     *  voltage-voltage block, zero elsewhere
     */
    Mat mAiMatrix;

    /** 
     *  The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b, this vector is z1 where b = M z1 + Ai z2.
     */
    Vec mVecForConstructingRhs1;
    
    /** 
     *  The vector multiplied by the Ai matrix. Ie, if the linear system to
     *  be solved is Ax=b, this vector is z2 where b = M z1 + Ai z2.
     */
    Vec mVecForConstructingRhs2;
    
    /** The bidomain assembler, used to set up the LHS matrix */
    BidomainAssembler<ELEMENT_DIM,SPACE_DIM>* mpBidomainAssembler;

    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
     BidomainNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpBidomainNeumannSurfaceTermAssembler;

    /**
     * If using state variable interpolation, points to an assembler to use in
     * computing the correction term to apply to the RHS.
     */
    BidomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpBidomainCorrectionTermAssembler;

     /** 
     *  The linear system that will be set up and solved first as part of the
     *  PDE solve
     */
     //LinearSystem* mpLinearSystem_first_solve;

     /**
     * If using state variable interpolation, points to an assembler to use in 
     * computing the correction term to apply to the RHS.
     */
     //BidomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpBidomainCorrectionTermAssembler;


    /** Overloaded InitialiseForSolve() which calls base version but also
     *  initialises mMassMatrix and mVecForConstructingRhs
     * 
     *  @param initialSolution initial solution
     */
    void InitialiseForSolve(Vec initialSolution);

    /** 
     *  Implementation of SetupLinearSystem() which uses the assembler to compute the
     *  LHS matrix, but sets up the RHS vector using the mass-matrix (constructed 
     *  using a separate assembler) multiplied by a vector
     * 
     *  @param currentSolution Solution at current time
     *  @param computeMatrix Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);

    /**
     *  Called before setting up the linear system, used to solve the cell models for first half timestep (step (i) above)
     *  @param currentSolution the latest solution vector
     */
    //void PrepareForSetupLinearSystem(Vec currentSolution);

    /**
     *  Called after solving the linear system, used to solve the cell models for second half timestep (step (iii) above)
     *  @param currentSolution the latest solution vector (ie the solution of the linear system).
     */
    //void FollowingSolveLinearSystem(Vec currentSolution);
    
    /**
     * Called when we are ready to solve the splitted system
     */
    void SolveStrang(Vec currentSolution);
    
    /**
     * Whether we are solving the cell system or not
     */
    bool SolvingCellSystem;
    
    /**
     * When using third-order OS, the substeps for the ODE and PDE system have to be the following
     */
    double StrangTimeSteps;
    
    /**
     * Index for the time-step array
     */
    int TimeSubStepIndex;

    
    /*
     *Stepper
     */
    TimeStepper stepper();




public:
    /**
     * Constructor
     *
     * @param bathSimulation Whether the simulation involves a perfusing bath
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    OperatorSplittingSDIRK2O2BidomainSolver(bool bathSimulation,
                              AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                              BidomainTissue<SPACE_DIM>* pTissue,
                              BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions);

    ~OperatorSplittingSDIRK2O2BidomainSolver();
    
  //  void FinaliseLinearSystem_first_solve(Vec existingSolution);

};


#endif /*OPERATORSPLITTINGSDIRK2O2BIDOMAINSOLVER_HPP_*/

