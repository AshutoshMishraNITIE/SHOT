﻿/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE
   This software is licensed under the Eclipse Public License 2.0.
   Please see the README and LICENSE files for more information.
*/

#include "NLPSolverCuttingPlaneMinimax.h"

#include "../Output.h"
#include "../Report.h"
#include "../Settings.h"
#include "../Timing.h"
#include "../Utilities.h"
#include "../DualSolver.h"
#include "../MIPSolver/IMIPSolver.h"

#include <functional>
#include <map>

#include "boost/math/tools/minima.hpp"

namespace SHOT
{

class MinimizationFunction
{
private:
    VectorDouble firstPt;
    VectorDouble secondPt;
    ProblemPtr NLPProblem;

public:
    MinimizationFunction(VectorDouble ptA, VectorDouble ptB, ProblemPtr prob)
    {
        firstPt = ptA;
        secondPt = ptB;
        NLPProblem = prob;
    }

    double operator()(const double x)
    {
        int length = firstPt.size();
        VectorDouble ptNew(length);

        for(int i = 0; i < length; i++)
        {
            ptNew.at(i) = x * firstPt.at(i) + (1 - x) * secondPt.at(i);
        }

        auto maxDev = NLPProblem->getMaxNumericConstraintValue(ptNew, NLPProblem->nonlinearConstraints);

        return (maxDev.normalizedValue);
    }
};

NLPSolverCuttingPlaneMinimax::NLPSolverCuttingPlaneMinimax(EnvironmentPtr envPtr, ProblemPtr problem)
    : INLPSolver(envPtr), sourceProblem(problem)
{
    auto solver = static_cast<ES_MIPSolver>(env->settings->getSetting<int>("MIP.Solver", "Dual"));

#ifdef HAS_CPLEX
    if(solver == ES_MIPSolver::Cplex)
    {
        LPSolver = std::make_unique<MIPSolverCplex>(env);
        env->output->outputDebug("Cplex selected as MIP solver for minimax solver.");
    }
#endif

#ifdef HAS_GUROBI
    if(solver == ES_MIPSolver::Gurobi)
    {
        LPSolver = std::make_unique<MIPSolverGurobi>(env);
        env->output->outputDebug("Gurobi selected as MIP solver for minimax solver.");
    }
#endif

#ifdef HAS_CBC
    if(solver == ES_MIPSolver::Cbc)
    {
        LPSolver = std::make_unique<MIPSolverCbc>(env);
        env->output->outputDebug("Cbc selected as MIP solver for minimax solver.");
    }
#endif

    env->output->outputDebug("Creating LP problem for minimax solver");
    createProblem(LPSolver.get(), sourceProblem);
    env->output->outputDebug("LP problem for minimax solver created");

    LPSolver->activateDiscreteVariables(false);
    LPSolver->initializeSolverSettings();
}

NLPSolverCuttingPlaneMinimax::~NLPSolverCuttingPlaneMinimax() = default;

void NLPSolverCuttingPlaneMinimax::saveProblemToFile([[maybe_unused]] std::string fileName) {}

E_NLPSolutionStatus NLPSolverCuttingPlaneMinimax::solveProblemInstance()
{
    int numVar = sourceProblem->properties.numberOfVariables;

    // Sets the maximal number of iterations
    int maxIter = env->settings->getSetting<int>("ESH.InteriorPoint.CuttingPlane.IterationLimit", "Dual");
    double termObjTolAbs
        = env->settings->getSetting<double>("ESH.InteriorPoint.CuttingPlane.TerminationToleranceAbs", "Dual");
    double termObjTolRel
        = env->settings->getSetting<double>("ESH.InteriorPoint.CuttingPlane.TerminationToleranceRel", "Dual");
    double constrSelFactor
        = env->settings->getSetting<double>("ESH.InteriorPoint.CuttingPlane.ConstraintSelectionFactor", "Dual");
    boost::uintmax_t maxIterSubsolver
        = env->settings->getSetting<int>("ESH.InteriorPoint.CuttingPlane.IterationLimitSubsolver", "Dual");
    int bitPrecision = env->settings->getSetting<int>("ESH.InteriorPoint.CuttingPlane.BitPrecision", "Dual");

    // currSol is the current LP solution, and prevSol the previous one
    VectorDouble currSol, prevSol;

    double lambda; // Variable in the linesearch minimization
    double mu; // Objective value for the linesearch minimization value

    // Corresponds to the difference between the LP solution objective value and
    // the objective found in the linesearch minimization procedure
    double maxObjDiffAbs = SHOT_DBL_MAX;
    double maxObjDiffRel = SHOT_DBL_MAX;

    double LPObjVar = 0.0;

    E_NLPSolutionStatus statusCode = E_NLPSolutionStatus::Error;

    int numHyperAdded = 0;
    int numHyperTot = 0;
    for(int i = 0; i <= maxIter; i++)
    {
        boost::uintmax_t maxIterSubsolverTmp = maxIterSubsolver;

        // Saves the LP problem to file if in debug mode
        if(env->settings->getSetting<bool>("Debug.Enable", "Output"))
        {
            std::stringstream ss;
            ss << env->settings->getSetting<std::string>("Debug.Path", "Output");
            ss << "/lpminimax";
            ss << i;
            ss << ".lp";
            LPSolver->writeProblemToFile(ss.str());
        }

        // Solves the problem and obtains the solution
        auto solStatus = LPSolver->solveProblem();
        env->solutionStatistics.numberOfProblemsMinimaxLP++;

        if(solStatus == E_ProblemSolutionStatus::Infeasible)
        {
            statusCode = E_NLPSolutionStatus::Infeasible;
            break;
        }
        else if(solStatus == E_ProblemSolutionStatus::Error)
        {
            statusCode = E_NLPSolutionStatus::Error;
            break;
        }
        else if(solStatus == E_ProblemSolutionStatus::Unbounded)
        {
            statusCode = E_NLPSolutionStatus::Unbounded;
            break;
        }
        else if(solStatus == E_ProblemSolutionStatus::TimeLimit)
        {
            statusCode = E_NLPSolutionStatus::TimeLimit;
        }
        else if(solStatus == E_ProblemSolutionStatus::IterationLimit)
        {
            statusCode = E_NLPSolutionStatus::IterationLimit;
        }
        else
        {
            statusCode = E_NLPSolutionStatus::Optimal;
        }

        auto LPVarSol = LPSolver->getVariableSolution(0);
        LPObjVar = LPSolver->getObjectiveValue();

        // Saves the LP solution to file if in debug mode
        if(env->settings->getSetting<bool>("Debug.Enable", "Output"))
        {
            std::stringstream ss;
            ss << env->settings->getSetting<std::string>("Debug.Path", "Output");
            ss << "/lpminimaxsolpt";
            ss << i;
            ss << ".txt";
            Utilities::saveVariablePointVectorToFile(LPVarSol, variableNames, ss.str());
        }

        if(std::isnan(LPObjVar))
        {
            statusCode = E_NLPSolutionStatus::Error;
            break;
        }

        numHyperAdded = 0;

        if(i == 0) // No linesearch minimization in first iteration, just add cutting plane in LP solution point
        {
            currSol = LPVarSol;
            lambda = -1; // For reporting purposes only
            mu = LPObjVar;
            env->report->outputIterationDetailHeaderMinimax();
        }
        else
        {
            MinimizationFunction funct(LPVarSol, prevSol, sourceProblem);

            // Solves the minization problem wrt lambda in [0, 1]
            auto minimizationResult
                = boost::math::tools::brent_find_minima(funct, 0.0, 1.0, bitPrecision, maxIterSubsolverTmp);

            lambda = minimizationResult.first; // The value for the line search parameter
            mu = minimizationResult.second; // The objective value

            // Calculates the corresponding solution point
            for(size_t j = 0; j < LPVarSol.size(); j++)
            {
                currSol.at(j) = lambda * LPVarSol.at(j) + (1 - lambda) * prevSol.at(j);
            }

            // The difference between linesearch and LP objective values
            maxObjDiffAbs = std::abs(mu - LPObjVar);
            maxObjDiffRel = maxObjDiffAbs / ((1e-10) + std::abs(LPObjVar));

            // Saves the LP solution to file if in debug mode
            if(env->settings->getSetting<bool>("Debug.Enable", "Output"))
            {
                std::stringstream ss;
                ss << env->settings->getSetting<std::string>("Debug.Path", "Output");
                ss << "/lpminimaxlinesearchsolpt";
                ss << i;
                ss << ".txt";
                Utilities::saveVariablePointVectorToFile(currSol, variableNames, ss.str());
            }
        }

        env->report->outputIterationDetailMinimax((i + 1), "LP", env->timing->getElapsedTime("Total"), numHyperAdded,
            numHyperTot, LPObjVar, mu, maxObjDiffAbs, maxObjDiffRel);

        if(mu < 0 && (maxObjDiffAbs < termObjTolAbs || maxObjDiffRel < termObjTolRel))
        {
            statusCode = E_NLPSolutionStatus::Optimal;
            break;
        }

        // Gets the most deviated constraints
        auto constraintValues = sourceProblem->getFractionOfDeviatingNonlinearConstraints(
            currSol, SHOT_DBL_MIN, constrSelFactor, LPObjVar);

        for(auto& NCV : constraintValues)
        {
            // Contains the coefficient and variable index for the terms in the generated cut
            std::map<int, double> elements;

            double constant = NCV.normalizedValue;
            auto gradient = NCV.constraint->calculateGradient(currSol, true);

            for(auto& G : gradient)
            {
                int variableIndex = G.first->index;
                double coefficient = G.second;

                auto element = elements.emplace(variableIndex, coefficient);

                if(!element.second)
                {
                    // Element already exists for the variable
                    element.first->second += coefficient;
                }

                constant = constant - coefficient * currSol.at(variableIndex);
            }

            // Adding the objective term
            elements.emplace(numVar, -1.0);

            // Adds the linear constraint
            if(LPSolver->addLinearConstraint(elements, constant,
                   "minimax_" + std::to_string(NCV.constraint->index) + "_" + std::to_string(numHyperTot))
                > 0)
            {
                numHyperTot++;
                numHyperAdded++;

                if(mu >= 0 && env->settings->getSetting<bool>("ESH.InteriorPoint.CuttingPlane.Reuse", "Dual")
                    && NCV.constraint->properties.convexity == E_Convexity::Convex)
                {
                    auto tmpPoint = currSol;
                    tmpPoint.pop_back();

                    Hyperplane hyperplane;
                    hyperplane.sourceConstraint = NCV.constraint;
                    hyperplane.sourceConstraintIndex = NCV.constraint->index;
                    hyperplane.generatedPoint = tmpPoint;
                    hyperplane.source = E_HyperplaneSource::InteriorPointSearch;

                    env->dualSolver->hyperplaneWaitingList.push_back(hyperplane);
                }
            }
        }

        prevSol = currSol;

        if(i == maxIter - 1)
        {
            statusCode = E_NLPSolutionStatus::IterationLimit;
            break;
        }

        if(numHyperAdded == 0)
        {
            statusCode = (objectiveValue > 0) ? E_NLPSolutionStatus::Infeasible : E_NLPSolutionStatus::Feasible;
            break;
        }
    }

    if(currSol.size() > 0)
    {
        // Removes the objective variable
        currSol.pop_back();

        solution = currSol;
        objectiveValue = LPObjVar;
    }

    return (statusCode);
}

double NLPSolverCuttingPlaneMinimax::getSolution(int i) { return (solution.at(i)); }

VectorDouble NLPSolverCuttingPlaneMinimax::getSolution() { return (solution); }

double NLPSolverCuttingPlaneMinimax::getObjectiveValue() { return (objectiveValue); }

// Creates the minimax problem in the LP solver
bool NLPSolverCuttingPlaneMinimax::createProblem(IMIPSolver* destination, ProblemPtr sourceProblem)
{
    // Now creating the variables
    bool variablesInitialized = true;

    for(auto& V : sourceProblem->allVariables)
    {
        variablesInitialized = variablesInitialized
            && destination->addVariable(V->name, V->properties.type, V->lowerBound, V->upperBound);

        if(env->settings->getSetting<bool>("Debug.Enable", "Output"))
        {
            variableNames.push_back(V->name);
        }
    }

    // Auxiliary objective variable for minimax problem
    double objUpperBound = env->settings->getSetting<double>("ESH.InteriorPoint.MinimaxObjectiveUpperBound", "Dual");

    variablesInitialized = variablesInitialized
        && destination->addVariable("shot_mmobjvar", E_VariableType::Real, -1e+10 + 1, objUpperBound);

    if(env->settings->getSetting<bool>("Debug.Enable", "Output"))
    {
        variableNames.push_back("shot_mmobjvar");
    }

    if(!variablesInitialized)
        return false;

    // Now creating the objective function

    bool objectiveInitialized = true;

    objectiveInitialized = objectiveInitialized && destination->initializeObjective();

    objectiveInitialized = objectiveInitialized
        && destination->addLinearTermToObjective(1.0, sourceProblem->properties.numberOfVariables);

    objectiveInitialized = objectiveInitialized && destination->finalizeObjective(true);

    if(!objectiveInitialized)
        return false;

    // Now creating the constraints

    bool constraintsInitialized = true;

    for(auto& C : env->problem->linearConstraints)
    {
        constraintsInitialized = constraintsInitialized && destination->initializeConstraint();

        if(C->properties.hasLinearTerms)
        {
            for(auto& T : C->linearTerms)
            {
                constraintsInitialized = constraintsInitialized
                    && destination->addLinearTermToConstraint(T->coefficient, T->variable->index);
            }
        }

        constraintsInitialized
            = constraintsInitialized && destination->finalizeConstraint(C->name, C->valueLHS, C->valueRHS);
    }

    if(!constraintsInitialized)
        return false;

    bool problemFinalized = destination->finalizeProblem();

    return (problemFinalized);
}

void NLPSolverCuttingPlaneMinimax::fixVariables(VectorInteger variableIndexes, VectorDouble variableValues)
{
    LPSolver->fixVariables(variableIndexes, variableValues);
}

void NLPSolverCuttingPlaneMinimax::unfixVariables() { LPSolver->unfixVariables(); }

void NLPSolverCuttingPlaneMinimax::setStartingPoint(
    [[maybe_unused]] VectorInteger variableIndexes, [[maybe_unused]] VectorDouble variableValues)
{
}

bool NLPSolverCuttingPlaneMinimax::isObjectiveFunctionNonlinear() { return (false); }

int NLPSolverCuttingPlaneMinimax::getObjectiveFunctionVariableIndex() { return (SHOT_INT_MAX); }

VectorDouble NLPSolverCuttingPlaneMinimax::getVariableLowerBounds()
{
    return (sourceProblem->getVariableLowerBounds());
}

VectorDouble NLPSolverCuttingPlaneMinimax::getVariableUpperBounds()
{
    return (sourceProblem->getVariableUpperBounds());
}

void NLPSolverCuttingPlaneMinimax::updateVariableLowerBound(int variableIndex, double bound)
{
    LPSolver->updateVariableLowerBound(variableIndex, bound);
}

void NLPSolverCuttingPlaneMinimax::updateVariableUpperBound(int variableIndex, double bound)
{
    LPSolver->updateVariableUpperBound(variableIndex, bound);
}

void NLPSolverCuttingPlaneMinimax::clearStartingPoint() {}

void NLPSolverCuttingPlaneMinimax::saveOptionsToFile([[maybe_unused]] std::string fileName) {}
} // namespace SHOT