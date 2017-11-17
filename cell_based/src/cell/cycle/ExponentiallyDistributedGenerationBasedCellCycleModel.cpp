/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "ExponentiallyDistributedGenerationBasedCellCycleModel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

ExponentiallyDistributedGenerationBasedCellCycleModel::ExponentiallyDistributedGenerationBasedCellCycleModel()
    : AbstractSimpleGenerationBasedCellCycleModel(),
      mRate(1.0/mTransitCellG1Duration)
{
}

ExponentiallyDistributedGenerationBasedCellCycleModel::ExponentiallyDistributedGenerationBasedCellCycleModel(const ExponentiallyDistributedGenerationBasedCellCycleModel& rModel)
   :  AbstractSimpleGenerationBasedCellCycleModel(rModel),
      mRate(rModel.mRate)
{
    /*
     * Initialize only the member variable defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* ExponentiallyDistributedGenerationBasedCellCycleModel::CreateCellCycleModel()
{
    return new ExponentiallyDistributedGenerationBasedCellCycleModel(*this);
}

void ExponentiallyDistributedGenerationBasedCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (    mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>()
            || mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() )
    {
        // Generate an exponential random number with mScale
        mG1Duration = p_gen->ExponentialRandomDeviate(mRate);
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

double ExponentiallyDistributedGenerationBasedCellCycleModel::GetRate() const
{
    return mRate;
}

void ExponentiallyDistributedGenerationBasedCellCycleModel::SetRate(double rate)
{
    mRate = rate;

    // These are now set to the average value of the G1Duration
    SetTransitCellG1Duration(1.0/rate);
    SetStemCellG1Duration(1.0/rate);
}

void ExponentiallyDistributedGenerationBasedCellCycleModel::SetStemCellG1Duration(double stemCellG1Duration)
{
    assert(stemCellG1Duration > 0.0);
    mStemCellG1Duration = stemCellG1Duration;

    mRate = 1.0/stemCellG1Duration;
}

void ExponentiallyDistributedGenerationBasedCellCycleModel::SetTransitCellG1Duration(double transitCellG1Duration)
{
    assert(transitCellG1Duration > 0.0);
    mTransitCellG1Duration = transitCellG1Duration;

    mRate = 1.0/transitCellG1Duration;
}

void ExponentiallyDistributedGenerationBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
     *rParamsFile << "\t\t\t<Rate>" << mRate << "</Rate>\n";

    // Call method on direct parent class
    AbstractSimpleGenerationBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ExponentiallyDistributedGenerationBasedCellCycleModel)
