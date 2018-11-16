/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE 
   This software is licensed under the Eclipse Public License 2.0. 
   Please see the README and LICENSE files for more information.
*/

#pragma once
#include "ISolutionStrategy.h"
#include "../Tasks/TaskSelectPrimalFixedNLPPointsFromSolutionPool.h"
#include "../Tasks/TaskAddIntegerCuts.h"
#include "../Tasks/TaskFindInteriorPoint.h"

#include "../Tasks/TaskBase.h"
#include "../Tasks/TaskSequential.h"
#include "../Tasks/TaskGoto.h"
#include "../Tasks/TaskConditional.h"

#include "../Tasks/TaskInitializeOriginalProblem.h"
#include "../Tasks/TaskInitializeIteration.h"
#include "../Tasks/TaskTerminate.h"

#include "../Tasks/TaskInitializeDualSolver.h"
#include "../Tasks/TaskCreateDualProblem.h"

#include "../Tasks/TaskExecuteSolutionLimitStrategy.h"
#include "../Tasks/TaskExecuteRelaxationStrategy.h"

#include "../Tasks/TaskPrintIterationReport.h"

#include "../Tasks/TaskSolveIteration.h"
#include "../Tasks/TaskPresolve.h"

#include "../Tasks/TaskCheckAbsoluteGap.h"
#include "../Tasks/TaskCheckIterationError.h"
#include "../Tasks/TaskCheckIterationLimit.h"
#include "../Tasks/TaskCheckObjectiveStagnation.h"
#include "../Tasks/TaskCheckConstraintTolerance.h"
#include "../Tasks/TaskCheckRelativeGap.h"
#include "../Tasks/TaskCheckTimeLimit.h"
#include "../Tasks/TaskCheckObjectiveGapNotMet.h"

#include "../Tasks/TaskInitializeLinesearch.h"
#include "../Tasks/TaskSelectHyperplanePointsLinesearch.h"
#include "../Tasks/TaskSelectHyperplanePointsIndividualLinesearch.h"
#include "../Tasks/TaskSelectHyperplanePointsSolution.h"
#include "../Tasks/TaskSelectHyperplanePointsProjection.h"
#include "../Tasks/TaskAddHyperplanes.h"

#include "../Tasks/TaskSelectPrimalCandidatesFromSolutionPool.h"
#include "../Tasks/TaskSelectPrimalCandidatesFromLinesearch.h"
#include "../Tasks/TaskSelectPrimalCandidatesFromNLP.h"

#include "../Tasks/TaskUpdateInteriorPoint.h"
#include "../Tasks/TaskUpdateNonlinearObjectiveByLinesearch.h"

#include "../Tasks/TaskSolveFixedDualProblem.h"
#include "../Tasks/TaskAddIntegerCuts.h"

#include "SHOTSettings.h"
#include "../ProcessInfo.h"

class SolutionStrategySingleTree : public ISolutionStrategy
{
  public:
    SolutionStrategySingleTree(OSInstance *osInstance);
    virtual ~SolutionStrategySingleTree();

    virtual bool solveProblem();
    virtual void initializeStrategy();

  protected:
};
