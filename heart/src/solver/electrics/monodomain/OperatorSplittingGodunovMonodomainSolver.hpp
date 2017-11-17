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

#ifndef OPERATORSPLITTINGGODUNOVMONODOMAINSOLVER_HPP_
#define OPERATORSPLITTINGGODUNOVMONODOMAINSOLVER_HPP_

#include "AbstractDynamicLinearPdeSolver.hpp"
#include "MassMatrixAssembler.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"
#include "MonodomainCorrectionTermAssembler.hpp"
#include "MonodomainTissue.hpp"
#include "MonodomainAssembler.hpp"

/**
 *  A monodomain solver that uses Godunov operator splitting of the diffusion (conductivity) term and the reaction
 *  (ionic current) term, instead of solving the full reaction-diffusion PDE. This does NOT refer to operator splitting
 *  of the two PDEs in the bidomain equations. For details see for example Sundnes et al "Computing the Electrical
 *  Activity of the Heart". This solves the PDEs using the BE method (first order).
 *
 *  The algorithm is, for solving from t=T to T+dt.
 *
 *  (i)   Solve ODEs   dV/dt = Iionic, du/dt = f(u,V)  for t=T to T+dt
 *        [giving updated V (internally, and in solution vector) and updated state variables]
 *  (ii)  Solve PDE    dV/dt = div (sigma_i grad V) + div (sigma_i grad phi_e), div (sigma_i grad V) + div ((sigma_i+sigma_e) grad phi_e) = 0  
 *        for t=T to dt         [using V from step i, --> updated V]
 *
 *  Notes
 *   (a)  This solver is FOR COMPARING ACCURACY, NOT PERFORMANCE. It has not been optimised and may or may not
 *        perform well in parallel.
 *   (b)  Usage: add the follwoing lines to your original code:
 *        HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplittingMonodomainSolver();
 * 	  HeartConfig::Instance()->SetUseGodunovReactionDiffusionOperatorSplittingMonodomainSolver();   
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class OperatorSplittingGodunovMonodomainSolver
  : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,1>
{
private:

    /** Monodomain tissue class (collection of cells, and conductivities) */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* mpMonodomainTissue;

    /** Boundary conditions */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* mpBoundaryConditions;

    /** The monodomain assembler, used to set up the LHS matrix */
    MonodomainAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainAssembler;

    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
    NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,1>* mpNeumannSurfaceTermsAssembler;

    /**
     * If using state variable interpolation, points to an assembler to use in
     * computing the correction term to apply to the RHS.
     */
    MonodomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainCorrectionTermAssembler;

    /** The mass matrix, used to computing the RHS vector */
    Mat mMassMatrix;

    /** The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b (excluding surface integrals), this vector is z where b=Mz.
     */
    Vec mVecForConstructingRhs;


    /**
     *  Implementation of SetupLinearSystem() which uses the assembler to compute the
     *  LHS matrix, but sets up the RHS vector using the mass-matrix (constructed
     *  using a separate assembler) multiplied by a vector
     *
     *  @param currentSolution  Solution at current time
     *  @param computeMatrix  Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);

public:
    /**
     *  Overloaded PrepareForSetupLinearSystem() methods which
     *  gets the cell models to solve themselves
     *
     *  @param currentSolution solution at current time
     */
    void PrepareForSetupLinearSystem(Vec currentSolution);

    /**
     *  Overloaded InitialiseForSolve
     *
     *  @param initialSolution initial solution
     */
    virtual void InitialiseForSolve(Vec initialSolution);

    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    OperatorSplittingGodunovMonodomainSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                     MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                     BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions);

    /**
     *  Destructor
     */
    virtual ~OperatorSplittingGodunovMonodomainSolver();
};



#endif /*OPERATORSPLITTINGGODUNOVMONODOMAINSOLVER_HPP_*/
