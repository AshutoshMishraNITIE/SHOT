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

        if (mostDevConstraint.value < 0)
        {
        }
        else
        {
            int p = 0;

            for (p = 0; p < 3; p++)
            {
                std::vector<double> newProjectedPoint(projectedPoint);

                double scalarProduct = 0;

                for (int k = 0; k < gradient->number; k++)
                {
                    scalarProduct += gradient->values[k] * gradient->values[k];
                }

                double factor = mostDevConstraint.value / scalarProduct;

                bool isValid = true;

                for (int k = 0; k < gradient->number; k++)
                {
                    newProjectedPoint[gradient->indexes[k]] = newProjectedPoint[gradient->indexes[k]] - factor * gradient->values[k];

                    auto variableLB = originalProblem->getVariableLowerBound(gradient->indexes[k]);
                    auto variableUB = originalProblem->getVariableUpperBound(gradient->indexes[k]);
                    auto currPointComponent = newProjectedPoint[gradient->indexes[k]];

                    if (currPointComponent < variableLB || currPointComponent > variableUB)
                    {
                        Output::getInstance().outputWarning("Projected point is outside of variable bounds: " + std::to_string(variableLB) + " < " + std::to_string(currPointComponent) + " < " + std::to_string(variableUB));
                        isValid = false;
                        break;
                    }
                }

                if (!isValid)
                    break;

                auto newMostDevConstraint = originalProblem->getMostDeviatingConstraint(newProjectedPoint);

                if (newMostDevConstraint.value <= Settings::getInstance().getDoubleSetting("ConstraintTolerance", "Termination"))
                    break;

                auto newGradient = originalProblem->calculateConstraintFunctionGradient(newMostDevConstraint.idx, newProjectedPoint);

                double checkValue = newMostDevConstraint.value;

                for (int k = 0; k < gradient->number; k++)
                {
                    checkValue += newGradient->values[k] * (solPoints.at(i).point[gradient->indexes[k]] - newProjectedPoint[gradient->indexes[k]]);
                }

                if (checkValue <= 0)
                {
                    Output::getInstance().outputWarning("Projected point does not remove solution point: " + std::to_string(checkValue));
                    break;
                }
                else
                {
                    Output::getInstance().outputDebug("Projected point does remove solution point: " + std::to_string(checkValue));
                }

                projectedPoint = newProjectedPoint;
                mostDevConstraint = newMostDevConstraint;
                gradient = newGradient;

                Output::getInstance().outputDebug("New projected point with max function value: " + std::to_string(mostDevConstraint.value));
            }

            Output::getInstance().outputWarning("Number of projections performed: " + std::to_string(p));

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
