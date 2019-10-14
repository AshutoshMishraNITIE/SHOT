/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE 
   This software is licensed under the Eclipse Public License 2.0. 
   Please see the README and LICENSE files for more information.
*/

#include "TaskSelectHyperplanePointsProjection.h"

TaskSelectHyperplanePointsProjection::TaskSelectHyperplanePointsProjection()
{
}

TaskSelectHyperplanePointsProjection::~TaskSelectHyperplanePointsProjection()
{
}

void TaskSelectHyperplanePointsProjection::run()
{
    this->run(ProcessInfo::getInstance().getPreviousIteration()->solutionPoints);
}

void TaskSelectHyperplanePointsProjection::run(vector<SolutionPoint> solPoints)
{
    ProcessInfo::getInstance().startTimer("DualCutGenerationRootSearch");

    int addedHyperplanes = 0;

    auto currIter = ProcessInfo::getInstance().getCurrentIteration(); // The unsolved new iteration

    auto originalProblem = ProcessInfo::getInstance().originalProblem;

    auto constrSelFactor = Settings::getInstance().getDoubleSetting("ECP.ConstraintSelectionFactor", "Dual");

    int maxProjections = Settings::getInstance().getIntSetting("PECP.MaxIterations", "Dual");

    double constrTolProjection = Settings::getInstance().getDoubleSetting("PECP.ConstraintTolerance", "Dual");

    for (int i = 0; i < solPoints.size(); i++)
    {
        if (addedHyperplanes >= Settings::getInstance().getIntSetting("HyperplaneCuts.MaxPerIteration", "Dual"))
        {
            ProcessInfo::getInstance().stopTimer("DualCutGenerationRootSearch");
            return;
        }

        auto projectedPoint = std::vector<double>(solPoints.at(i).point);
        auto mostDevConstraint = originalProblem->getMostDeviatingConstraint(projectedPoint);
        auto gradient = originalProblem->calculateConstraintFunctionGradient(mostDevConstraint.idx, projectedPoint);

        if (mostDevConstraint.idx == -1 || !currIter->isMIP() || mostDevConstraint.value <= constrTolProjection)
        {
            Hyperplane hyperplane;
            hyperplane.sourceConstraintIndex = mostDevConstraint.idx;
            hyperplane.generatedPoint = projectedPoint;

            if (solPoints.at(i).isRelaxedPoint)
            {
                hyperplane.source = E_HyperplaneSource::MIPCallbackRelaxed;
            }
            else if (i == 0 && currIter->isMIP())
            {
                hyperplane.source = E_HyperplaneSource::MIPOptimalSolutionPoint;
            }
            else if (currIter->isMIP())
            {
                hyperplane.source = E_HyperplaneSource::MIPSolutionPoolSolutionPoint;
            }
            else
            {
                hyperplane.source = E_HyperplaneSource::LPRelaxedSolutionPoint;
            }

            ProcessInfo::getInstance().hyperplaneWaitingList.push_back(hyperplane);

            addedHyperplanes++;
        }
        else
        {
            int p = 0;

            for (p = 0; p < maxProjections; p++)
            {
                std::vector<double> newProjectedPoint(projectedPoint);

                double scalarProduct = 0;

                for (int k = 0; k < gradient->number; k++)
                {
                    scalarProduct += gradient->values[k] * gradient->values[k];
                }

                double factor = mostDevConstraint.value / scalarProduct;

                bool isValid = true;

                auto discreteVariableIndexes = originalProblem->getDiscreteVariableIndices();

                for (int k = 0; k < gradient->number; k++)
                {
                    if (std::find(discreteVariableIndexes.begin(), discreteVariableIndexes.end(), gradient->indexes[k]) != discreteVariableIndexes.end())
                        continue;

                    newProjectedPoint[gradient->indexes[k]] = newProjectedPoint[gradient->indexes[k]] - factor * gradient->values[k];

                    auto variableLB = originalProblem->getVariableLowerBound(gradient->indexes[k]);
                    auto variableUB = originalProblem->getVariableUpperBound(gradient->indexes[k]);
                    auto currPointComponent = newProjectedPoint[gradient->indexes[k]];

                    if (currPointComponent < variableLB || currPointComponent > variableUB)
                    {
                        Output::getInstance().outputAlways("         Projected point is outside of variable bounds: " + std::to_string(variableLB) + " < " + std::to_string(currPointComponent) + " < " + std::to_string(variableUB));
                        isValid = false;
                    }
                }

                if (!isValid)
                {
                    ProcessInfo::getInstance().solutionStatistics.totalNumberOfFailedProjections++;
                    break;
                }

                auto newMostDevConstraint = originalProblem->getMostDeviatingConstraint(newProjectedPoint);

                if (newMostDevConstraint.value <= constrTolProjection)
                    break;

                auto newGradient = originalProblem->calculateConstraintFunctionGradient(newMostDevConstraint.idx, newProjectedPoint);

                double checkValue = newMostDevConstraint.value;

                for (int k = 0; k < gradient->number; k++)
                {
                    checkValue += newGradient->values[k] * (solPoints.at(i).point[gradient->indexes[k]] - newProjectedPoint[gradient->indexes[k]]);
                }

                if (checkValue <= 0)
                {
                    ProcessInfo::getInstance().solutionStatistics.totalNumberOfFailedProjections++;
                    Output::getInstance().outputAlways("         Projected point does not remove solution point: " + std::to_string(checkValue));
                    break;
                }
                else
                {
                    Output::getInstance().outputWarning("         Projected point does remove solution point: " + std::to_string(checkValue));
                }

                projectedPoint = newProjectedPoint;
                mostDevConstraint = newMostDevConstraint;
                gradient = newGradient;

                Output::getInstance().outputInfo("         New projected point with max function value: " + std::to_string(mostDevConstraint.value));

                ProcessInfo::getInstance().solutionStatistics.totalNumberOfSuccessfulProjections++;
            }

            if (p > 0)
                Output::getInstance().outputAlways("         Number of projections performed: " + std::to_string(p));

            mostDevConstraint = originalProblem->getMostDeviatingConstraint(projectedPoint);

            Hyperplane hyperplane;
            hyperplane.sourceConstraintIndex = mostDevConstraint.idx;
            hyperplane.generatedPoint = projectedPoint;

            if (solPoints.at(i).isRelaxedPoint)
            {
                hyperplane.source = E_HyperplaneSource::MIPCallbackRelaxed;
            }
            else if (i == 0 && currIter->isMIP())
            {
                hyperplane.source = E_HyperplaneSource::MIPOptimalSolutionPoint;
            }
            else if (currIter->isMIP())
            {
                hyperplane.source = E_HyperplaneSource::MIPSolutionPoolSolutionPoint;
            }
            else
            {
                hyperplane.source = E_HyperplaneSource::LPRelaxedSolutionPoint;
            }

            ProcessInfo::getInstance().hyperplaneWaitingList.push_back(hyperplane);

            addedHyperplanes++;
        }
    }

    ProcessInfo::getInstance().stopTimer("DualCutGenerationRootSearch");
}

std::string TaskSelectHyperplanePointsProjection::getType()
{
    std::string type = typeid(this).name();
    return (type);
}
