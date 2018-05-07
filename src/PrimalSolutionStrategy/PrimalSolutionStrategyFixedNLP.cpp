/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE 
   This software is licensed under the Eclipse Public License 2.0. 
   Please see the README and LICENSE files for more information.
*/

#include "PrimalSolutionStrategyFixedNLP.h"

PrimalSolutionStrategyFixedNLP::PrimalSolutionStrategyFixedNLP()
{
    originalNLPTime = Settings::getInstance().getDoubleSetting("FixedInteger.Frequency.Time", "Primal");
    originalNLPIter = Settings::getInstance().getIntSetting("FixedInteger.Frequency.Iteration", "Primal");

    switch (static_cast<ES_PrimalNLPSolver>(Settings::getInstance().getIntSetting("FixedInteger.Solver", "Primal")))
    {
    case (ES_PrimalNLPSolver::CuttingPlane):
    {
        ProcessInfo::getInstance().usedPrimalNLPSolver = ES_PrimalNLPSolver::CuttingPlane;
        NLPSolver = new NLPSolverCuttingPlaneRelaxed();
        break;
    }
    case (ES_PrimalNLPSolver::Ipopt):
    {
        ProcessInfo::getInstance().usedPrimalNLPSolver = ES_PrimalNLPSolver::Ipopt;
        NLPSolver = new NLPSolverIpoptRelaxed();
        break;
    }
#ifdef HAS_GAMS
    case (ES_PrimalNLPSolver::GAMS):
    {
        ProcessInfo::getInstance().usedPrimalNLPSolver = ES_PrimalNLPSolver::GAMS;
        NLPSolver = new NLPSolverGAMS();
        break;
    }
#endif
    default:
        Output::getInstance().Output::getInstance().outputError("Error in solver definition for primal NLP solver. Check option 'Primal.FixedInteger.Solver'.");
        throw new ErrorClass("Error in solver definition for primal NLP solver. Check option 'Primal.FixedInteger.Solver'.");

        throw std::logic_error("Unknown PrimalNLPSolver setting.");
    }

    NLPSolver->setProblem(ProcessInfo::getInstance().originalProblem->getProblemInstance());

    if (Settings::getInstance().getBoolSetting("FixedInteger.CreateInfeasibilityCut", "Primal"))
    {
        if (static_cast<ES_HyperplaneCutStrategy>(Settings::getInstance().getIntSetting("CutStrategy", "Dual")) == ES_HyperplaneCutStrategy::ESH)
        {
            if (static_cast<ES_RootsearchConstraintStrategy>(Settings::getInstance().getIntSetting(
                    "ESH.Linesearch.ConstraintStrategy", "Dual")) == ES_RootsearchConstraintStrategy::AllAsMaxFunct)
            {
                taskSelectHPPts = new TaskSelectHyperplanePointsLinesearch();
            }
            else
            {
                taskSelectHPPts = new TaskSelectHyperplanePointsIndividualLinesearch();
            }
        }
        else
        {
            taskSelectHPPts = new TaskSelectHyperplanePointsSolution();
        }
    }

    this->originalIterFrequency = Settings::getInstance().getIntSetting("FixedInteger.Frequency.Iteration", "Primal");
    this->originalTimeFrequency = Settings::getInstance().getDoubleSetting("FixedInteger.Frequency.Time", "Primal");
}

PrimalSolutionStrategyFixedNLP::~PrimalSolutionStrategyFixedNLP()
{
    delete taskSelectHPPts;
    delete NLPSolver;

    discreteVariableIndexes.clear();
    testedPoints.clear();
    fixPoint.clear();
}

bool PrimalSolutionStrategyFixedNLP::runStrategy()
{
    auto currIter = ProcessInfo::getInstance().getCurrentIteration();

    NLPSolver->initializeProblem();

    int numVars = ProcessInfo::getInstance().originalProblem->getNumberOfVariables();

    auto discreteVariableIndexes = ProcessInfo::getInstance().originalProblem->getDiscreteVariableIndices();
    auto realVariableIndexes = ProcessInfo::getInstance().originalProblem->getRealVariableIndices();

    bool isSolved;

    vector<PrimalFixedNLPCandidate> testPts;

    // Fix variables
    auto varTypes = ProcessInfo::getInstance().originalProblem->getVariableTypes();

    if (ProcessInfo::getInstance().primalFixedNLPCandidates.size() == 0)
    {
        ProcessInfo::getInstance().solutionStatistics.numberOfIterationsWithoutNLPCallMIP++;
        return (false);
    }

    if (testedPoints.size() > 0)
    {
        for (int j = 0; j < ProcessInfo::getInstance().primalFixedNLPCandidates.size(); j++)
        {
            for (int i = 0; i < testedPoints.size(); i++)
            {
                if (UtilityFunctions::isDifferentRoundedSelectedElements(
                        ProcessInfo::getInstance().primalFixedNLPCandidates.at(j).point, testedPoints.at(i),
                        discreteVariableIndexes))
                {
                    testPts.push_back(ProcessInfo::getInstance().primalFixedNLPCandidates.at(j));
                    testedPoints.push_back(ProcessInfo::getInstance().primalFixedNLPCandidates.at(j).point);
                    break;
                }
            }
        }
    }
    else
    {
        testPts.push_back(ProcessInfo::getInstance().primalFixedNLPCandidates.at(0));
        testedPoints.push_back(ProcessInfo::getInstance().primalFixedNLPCandidates.at(0).point);
    }

    if (testPts.size() == 0)
    {
        ProcessInfo::getInstance().solutionStatistics.numberOfIterationsWithoutNLPCallMIP++;
        return (false);
    }

    auto lbs = NLPSolver->getVariableLowerBounds();
    auto ubs = NLPSolver->getVariableUpperBounds();

    for (int j = 0; j < testPts.size(); j++)
    {
        auto oldPrimalBound = ProcessInfo::getInstance().getPrimalBound();
        double timeStart = ProcessInfo::getInstance().getElapsedTime("Total");
        std::vector<double> fixedVariableValues(discreteVariableIndexes.size());

        int sizeOfVariableVector = NLPSolver->NLPProblem->getNumberOfVariables();

        std::vector<int> startingPointIndexes(sizeOfVariableVector);
        std::vector<double> startingPointValues(sizeOfVariableVector);

        // Sets the fixed values for discrete variables
        for (int k = 0; k < discreteVariableIndexes.size(); k++)
        {
            int currVarIndex = discreteVariableIndexes.at(k);

            auto tmpSolPt = UtilityFunctions::round(testPts.at(j).point.at(currVarIndex));

            fixedVariableValues.at(k) = tmpSolPt;

            // Sets the starting point to the fixed value
            if (Settings::getInstance().getBoolSetting("FixedInteger.Warmstart", "Primal"))
            {
                startingPointIndexes.at(currVarIndex) = currVarIndex;
                startingPointValues.at(currVarIndex) = tmpSolPt;
            }
        }

        if (Settings::getInstance().getBoolSetting("FixedInteger.Warmstart", "Primal"))
        {
            for (int k = 0; k < realVariableIndexes.size(); k++)
            {
                int currVarIndex = realVariableIndexes.at(k);

                if (NLPSolver->isObjectiveFunctionNonlinear() && currVarIndex == NLPSolver->getObjectiveFunctionVariableIndex())
                {
                    continue;
                }

                auto tmpSolPt = testPts.at(j).point.at(currVarIndex);

                startingPointIndexes.at(currVarIndex) = currVarIndex;
                startingPointValues.at(currVarIndex) = tmpSolPt;
            }
        }

        NLPSolver->setStartingPoint(startingPointIndexes, startingPointValues);

        NLPSolver->fixVariables(discreteVariableIndexes, fixedVariableValues);

        if (Settings::getInstance().getBoolSetting("Debug.Enable", "Output"))
        {
            std::string filename = Settings::getInstance().getStringSetting("Debug.Path", "Output") + "/primalnlp" + to_string(ProcessInfo::getInstance().getCurrentIteration()->iterationNumber) + "_" + to_string(j);
            NLPSolver->saveProblemToFile(filename + ".txt");
            NLPSolver->saveOptionsToFile(filename + ".osrl");
        }

        auto solvestatus = NLPSolver->solveProblem();

        NLPSolver->unfixVariables();
        ProcessInfo::getInstance().solutionStatistics.numberOfProblemsFixedNLP++;

        double timeEnd = ProcessInfo::getInstance().getElapsedTime("Total");

        std::string sourceDesc;
        switch (testPts.at(j).sourceType)
        {
        case E_PrimalNLPSource::FirstSolution:
            sourceDesc = "SOLPT ";
            break;
        case E_PrimalNLPSource::FeasibleSolution:
            sourceDesc = "FEASPT";
            break;
        case E_PrimalNLPSource::InfeasibleSolution:
            sourceDesc = "UNFEAS";
            break;
        case E_PrimalNLPSource::SmallestDeviationSolution:
            sourceDesc = "SMADEV";
            break;
        case E_PrimalNLPSource::FirstSolutionNewDualBound:
            sourceDesc = "NEWDB";
            break;
        default:
            break;
        }

        //std::string solExpr = ((boost::format("(%.3f s)") % (timeEnd - timeStart)).str());

        if (solvestatus == E_NLPSolutionStatus::Feasible || solvestatus == E_NLPSolutionStatus::Optimal)
        {
            double tmpObj = NLPSolver->getObjectiveValue();
            auto variableSolution = NLPSolver->getSolution();

            if (ProcessInfo::getInstance().originalProblem->isObjectiveFunctionNonlinear())
            {
                variableSolution.push_back(tmpObj);
            }

            auto mostDevConstr = ProcessInfo::getInstance().originalProblem->getMostDeviatingConstraint(
                variableSolution);

            if (Settings::getInstance().getBoolSetting("FixedInteger.Frequency.Dynamic", "Primal"))
            {
                int iters = max(ceil(Settings::getInstance().getIntSetting("FixedInteger.Frequency.Iteration", "Primal") * 0.98),
                                originalNLPIter);

                if (iters > max(0.1 * this->originalIterFrequency, 1.0))
                    Settings::getInstance().updateSetting("FixedInteger.Frequency.Iteration", "Primal", iters);

                double interval = max(
                    0.9 * Settings::getInstance().getDoubleSetting("FixedInteger.Frequency.Time", "Primal"),
                    originalNLPTime);

                if (interval > 0.1 * this->originalTimeFrequency)
                    Settings::getInstance().updateSetting("FixedInteger.Frequency.Time", "Primal", interval);

                ProcessInfo::getInstance().addPrimalSolutionCandidate(variableSolution, E_PrimalSolutionSource::NLPFixedIntegers, currIter->iterationNumber);
            }

            Output::getInstance().outputIterationDetail(ProcessInfo::getInstance().solutionStatistics.numberOfProblemsFixedNLP,
                                                        ("NLP" + sourceDesc),
                                                        ProcessInfo::getInstance().getElapsedTime("Total"),
                                                        currIter->numHyperplanesAdded,
                                                        currIter->totNumHyperplanes,
                                                        ProcessInfo::getInstance().getDualBound(),
                                                        ProcessInfo::getInstance().getPrimalBound(),
                                                        ProcessInfo::getInstance().getAbsoluteObjectiveGap(),
                                                        ProcessInfo::getInstance().getRelativeObjectiveGap(),
                                                        tmpObj,
                                                        mostDevConstr.idx,
                                                        mostDevConstr.value,
                                                        E_IterationLineType::PrimalNLP);
        }
        else
        {
            double tmpObj = NLPSolver->getObjectiveValue();

            // Utilize the solution point for adding a cutting plane / supporting hyperplane
            std::vector<SolutionPoint> solutionPoints(1);

            auto variableSolution = NLPSolver->getSolution();

            if (ProcessInfo::getInstance().originalProblem->isObjectiveFunctionNonlinear())
            {
                variableSolution.push_back(tmpObj);
            }

            auto mostDevConstr = ProcessInfo::getInstance().originalProblem->getMostDeviatingConstraint(variableSolution);

            if (Settings::getInstance().getBoolSetting("FixedInteger.CreateInfeasibilityCut", "Primal"))
            {
                SolutionPoint tmpSolPt;
                tmpSolPt.point = variableSolution;
                tmpSolPt.objectiveValue = ProcessInfo::getInstance().originalProblem->calculateOriginalObjectiveValue(
                    variableSolution);
                tmpSolPt.iterFound = ProcessInfo::getInstance().getCurrentIteration()->iterationNumber;
                tmpSolPt.maxDeviation = mostDevConstr;

                solutionPoints.at(0) = tmpSolPt;

                if (static_cast<ES_HyperplaneCutStrategy>(Settings::getInstance().getIntSetting(
                        "CutStrategy", "Dual")) == ES_HyperplaneCutStrategy::ESH)
                {
                    if (static_cast<ES_RootsearchConstraintStrategy>(Settings::getInstance().getIntSetting(
                            "ESH.Linesearch.ConstraintStrategy", "Dual")) == ES_RootsearchConstraintStrategy::AllAsMaxFunct)
                    {
                        static_cast<TaskSelectHyperplanePointsLinesearch *>(taskSelectHPPts)->run(solutionPoints);
                    }
                    else
                    {
                        static_cast<TaskSelectHyperplanePointsIndividualLinesearch *>(taskSelectHPPts)->run(solutionPoints);
                    }
                }
                else
                {
                    static_cast<TaskSelectHyperplanePointsSolution *>(taskSelectHPPts)->run(solutionPoints);
                }
            }

            if (Settings::getInstance().getBoolSetting("FixedInteger.Frequency.Dynamic", "Primal"))
            {
                int iters = ceil(Settings::getInstance().getIntSetting("FixedInteger.Frequency.Iteration", "Primal") * 1.02);

                if (iters < 10 * this->originalIterFrequency)
                    Settings::getInstance().updateSetting("FixedInteger.Frequency.Iteration", "Primal", iters);

                double interval = 1.1 * Settings::getInstance().getDoubleSetting("FixedInteger.Frequency.Time", "Primal");

                if (interval < 10 * this->originalTimeFrequency)
                    Settings::getInstance().updateSetting("FixedInteger.Frequency.Time", "Primal", interval);

                Output::getInstance().outputInfo(
                    "     Duration:  " + to_string(timeEnd - timeStart) + " s. New interval: " + to_string(interval) + " s or " + to_string(iters) + " iters.");
            }

            Output::getInstance().outputIterationDetail(ProcessInfo::getInstance().solutionStatistics.numberOfProblemsFixedNLP,
                                                        ("NLP" + sourceDesc),
                                                        ProcessInfo::getInstance().getElapsedTime("Total"),
                                                        currIter->numHyperplanesAdded,
                                                        currIter->totNumHyperplanes,
                                                        ProcessInfo::getInstance().getDualBound(),
                                                        ProcessInfo::getInstance().getPrimalBound(),
                                                        ProcessInfo::getInstance().getAbsoluteObjectiveGap(),
                                                        ProcessInfo::getInstance().getRelativeObjectiveGap(),
                                                        NAN,
                                                        mostDevConstr.idx,
                                                        mostDevConstr.value,
                                                        E_IterationLineType::PrimalNLP);

            if (Settings::getInstance().getBoolSetting("HyperplaneCuts.UseIntegerCuts", "Dual") && ProcessInfo::getInstance().originalProblem->getNumberOfIntegerVariables() == 0)
            {
                //Add integer cut.

                auto binVars = ProcessInfo::getInstance().originalProblem->getBinaryVariableIndices();

                if (binVars.size() > 0)
                {
                    std::vector<int> elements;

                    for (int i = 0; i < binVars.size(); i++)
                    {
                        if (testPts.at(j).point.at(binVars.at(i)) > 0.99)
                        {
                            elements.push_back(binVars.at(i));
                        }
                    }
                    ProcessInfo::getInstance().integerCutWaitingList.push_back(elements);
                }
            }
        }

        ProcessInfo::getInstance().solutionStatistics.numberOfIterationsWithoutNLPCallMIP = 0;
        ProcessInfo::getInstance().solutionStatistics.timeLastFixedNLPCall = ProcessInfo::getInstance().getElapsedTime("Total");
    }

    ProcessInfo::getInstance().primalFixedNLPCandidates.clear();

    return (true);
}
