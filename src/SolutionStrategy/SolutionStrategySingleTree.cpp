/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE 
   This software is licensed under the Eclipse Public License 2.0. 
   Please see the README and LICENSE files for more information.
*/

#include "SolutionStrategySingleTree.h"

SolutionStrategySingleTree::SolutionStrategySingleTree(OSInstance *osInstance)
{
    ProcessInfo::getInstance().createTimer("ProblemInitialization", " - problem initialization");
    ProcessInfo::getInstance().createTimer("InteriorPointSearch", " - interior point search");

    ProcessInfo::getInstance().createTimer("DualProblemsRelaxed", "   - solving relaxed problems");
    ProcessInfo::getInstance().createTimer("DualStrategy", " - dual strategy");
    ProcessInfo::getInstance().createTimer("DualProblemsDiscrete", "   - solving MIP problems");
    ProcessInfo::getInstance().createTimer("DualCutGenerationRootSearch", "   - performing root search for cuts");
    ProcessInfo::getInstance().createTimer("DualObjectiveLiftRootSearch", "   - performing root search for objective lift");

    ProcessInfo::getInstance().createTimer("PrimalStrategy", " - primal strategy");
    ProcessInfo::getInstance().createTimer("PrimalBoundStrategyNLP", "   - solving NLP problems");
    ProcessInfo::getInstance().createTimer("PrimalBoundStrategyRootSearch", "   - performing root searches");

    auto solver = static_cast<ES_InteriorPointStrategy>(Settings::getInstance().getIntSetting("ESH.InteriorPoint.Solver", "Dual"));
    auto solverMIP = static_cast<ES_MIPSolver>(Settings::getInstance().getIntSetting("MIP.Solver", "Dual"));

    TaskBase *tFinalizeSolution = new TaskSequential();

    TaskBase *tInitMIPSolver = new TaskInitializeDualSolver(solverMIP, true);
    ProcessInfo::getInstance().tasks->addTask(tInitMIPSolver, "InitMIPSolver");

    auto MIPSolver = ProcessInfo::getInstance().MIPSolver;

    TaskBase *tInitOrigProblem = new TaskInitializeOriginalProblem(osInstance);
    ProcessInfo::getInstance().tasks->addTask(tInitOrigProblem, "InitOrigProb");

    if (Settings::getInstance().getIntSetting("CutStrategy", "Dual") == (int)ES_HyperplaneCutStrategy::ESH && (ProcessInfo::getInstance().originalProblem->getObjectiveFunctionType() != E_ObjectiveFunctionType::Quadratic || ProcessInfo::getInstance().originalProblem->getNumberOfNonlinearConstraints() != 0))
    {
        TaskBase *tFindIntPoint = new TaskFindInteriorPoint();
        ProcessInfo::getInstance().tasks->addTask(tFindIntPoint, "FindIntPoint");
    }

    TaskBase *tCreateDualProblem = new TaskCreateDualProblem(MIPSolver);
    ProcessInfo::getInstance().tasks->addTask(tCreateDualProblem, "CreateDualProblem");

    TaskBase *tInitializeLinesearch = new TaskInitializeLinesearch();
    ProcessInfo::getInstance().tasks->addTask(tInitializeLinesearch, "InitializeLinesearch");

    TaskBase *tInitializeIteration = new TaskInitializeIteration();
    ProcessInfo::getInstance().tasks->addTask(tInitializeIteration, "InitIter");

    TaskBase *tAddHPs = new TaskAddHyperplanes(MIPSolver);
    ProcessInfo::getInstance().tasks->addTask(tAddHPs, "AddHPs");

    TaskBase *tExecuteRelaxStrategy = new TaskExecuteRelaxationStrategy(MIPSolver);
    ProcessInfo::getInstance().tasks->addTask(tExecuteRelaxStrategy, "ExecRelaxStrategyInitial");

    if (static_cast<ES_MIPPresolveStrategy>(Settings::getInstance().getIntSetting("MIP.Presolve.Frequency", "Dual")) != ES_MIPPresolveStrategy::Never)
    {
        TaskBase *tPresolve = new TaskPresolve(MIPSolver);
        ProcessInfo::getInstance().tasks->addTask(tPresolve, "Presolve");
    }

    TaskBase *tSolveIteration = new TaskSolveIteration(MIPSolver);
    ProcessInfo::getInstance().tasks->addTask(tSolveIteration, "SolveIter");

    TaskBase *tSelectPrimSolPool = new TaskSelectPrimalCandidatesFromSolutionPool();
    ProcessInfo::getInstance().tasks->addTask(tSelectPrimSolPool, "SelectPrimSolPool");
    dynamic_cast<TaskSequential *>(tFinalizeSolution)->addTask(tSelectPrimSolPool);

    TaskBase *tPrintIterReport = new TaskPrintIterationReport();
    ProcessInfo::getInstance().tasks->addTask(tPrintIterReport, "PrintIterReport");

    TaskBase *tCheckAbsGap = new TaskCheckAbsoluteGap("FinalizeSolution");
    ProcessInfo::getInstance().tasks->addTask(tCheckAbsGap, "CheckAbsGap");

    TaskBase *tCheckRelGap = new TaskCheckRelativeGap("FinalizeSolution");
    ProcessInfo::getInstance().tasks->addTask(tCheckRelGap, "CheckRelGap");

    TaskBase *tCheckTimeLim = new TaskCheckTimeLimit("FinalizeSolution");
    ProcessInfo::getInstance().tasks->addTask(tCheckTimeLim, "CheckTimeLim");

    TaskBase *tCheckIterError = new TaskCheckIterationError("FinalizeSolution");
    ProcessInfo::getInstance().tasks->addTask(tCheckIterError, "CheckIterError");

    TaskBase *tCheckConstrTol = new TaskCheckConstraintTolerance("FinalizeSolution");
    ProcessInfo::getInstance().tasks->addTask(tCheckConstrTol, "CheckConstrTol");

    TaskBase *tCheckObjectiveGapNotMet = new TaskCheckObjectiveGapNotMet("FinalizeSolution");
    ProcessInfo::getInstance().tasks->addTask(tCheckObjectiveGapNotMet, "CheckObjGapNotMet");

    ProcessInfo::getInstance().tasks->addTask(tInitializeIteration, "InitIter");

    if (static_cast<ES_HyperplaneCutStrategy>(Settings::getInstance().getIntSetting("CutStrategy", "Dual")) == ES_HyperplaneCutStrategy::ESH)
    {
        TaskBase *tUpdateInteriorPoint = new TaskUpdateInteriorPoint();
        ProcessInfo::getInstance().tasks->addTask(tUpdateInteriorPoint, "UpdateInteriorPoint");

        if (static_cast<ES_RootsearchConstraintStrategy>(Settings::getInstance().getIntSetting(
                "ESH.Linesearch.ConstraintStrategy", "Dual")) == ES_RootsearchConstraintStrategy::AllAsMaxFunct)
        {
            TaskBase *tSelectHPPts = new TaskSelectHyperplanePointsLinesearch();
            ProcessInfo::getInstance().tasks->addTask(tSelectHPPts, "SelectHPPts");
        }
        else
        {
            TaskBase *tSelectHPPts = new TaskSelectHyperplanePointsIndividualLinesearch();
            ProcessInfo::getInstance().tasks->addTask(tSelectHPPts, "SelectHPPts");
        }
    }
    else if (static_cast<ES_HyperplaneCutStrategy>(Settings::getInstance().getIntSetting("CutStrategy", "Dual")) == ES_HyperplaneCutStrategy::ECP)
    {
        TaskBase *tSelectHPPts = new TaskSelectHyperplanePointsSolution();
        ProcessInfo::getInstance().tasks->addTask(tSelectHPPts, "SelectHPPts");
    }
    else
    {
        TaskBase *tSelectHPPts = new TaskSelectHyperplanePointsProjection();
        ProcessInfo::getInstance().tasks->addTask(tSelectHPPts, "SelectHPPts");
    }

    TaskBase *tGoto = new TaskGoto("AddHPs");
    ProcessInfo::getInstance().tasks->addTask(tGoto, "Goto");

    if (Settings::getInstance().getIntSetting("FixedInteger.CallStrategy", "Primal") && ProcessInfo::getInstance().originalProblem->getNumberOfNonlinearConstraints() > 0 && ProcessInfo::getInstance().originalProblem->getNumberOfDiscreteVariables() > 0)
    {
        TaskBase *tSelectPrimFixedNLPSolPool = new TaskSelectPrimalFixedNLPPointsFromSolutionPool();
        ProcessInfo::getInstance().tasks->addTask(tSelectPrimFixedNLPSolPool, "SelectPrimFixedNLPSolPool");
        dynamic_cast<TaskSequential *>(tFinalizeSolution)->addTask(tSelectPrimFixedNLPSolPool);

        TaskBase *tSelectPrimNLPCheck = new TaskSelectPrimalCandidatesFromNLP();
        ProcessInfo::getInstance().tasks->addTask(tSelectPrimNLPCheck, "SelectPrimNLPCheck");
        dynamic_cast<TaskSequential *>(tFinalizeSolution)->addTask(tSelectPrimNLPCheck);

        ProcessInfo::getInstance().tasks->addTask(tCheckAbsGap, "CheckAbsGap");
        ProcessInfo::getInstance().tasks->addTask(tCheckRelGap, "CheckRelGap");
    }

    ProcessInfo::getInstance().tasks->addTask(tFinalizeSolution, "FinalizeSolution");
}

SolutionStrategySingleTree::~SolutionStrategySingleTree()
{
}

bool SolutionStrategySingleTree::solveProblem()
{
    TaskBase *nextTask;

    while (ProcessInfo::getInstance().tasks->getNextTask(nextTask))
    {
        Output::getInstance().outputInfo("┌─── Started task:  " + nextTask->getType());
        nextTask->run();
        Output::getInstance().outputInfo("└─── Finished task: " + nextTask->getType());
    }

    return (true);
}

void SolutionStrategySingleTree::initializeStrategy()
{
}
