/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE
   This software is licensed under the Eclipse Public License 2.0.
   Please see the README and LICENSE files for more information.
*/

#pragma once
#include "TaskBase.h"

#include <map>
#include <tuple>

#include "../Model/AuxiliaryVariables.h"
#include "../Model/Constraints.h"
#include "../Model/NonlinearExpressions.h"
#include "../Model/Problem.h"
#include "../Model/Terms.h"
#include "../Model/Variables.h"

namespace SHOT
{
struct Reformulation
{
    LinearConstraints linearConstraints;
    QuadraticConstraints quadraticConstraints;
    NonlinearConstraint nonlinearConstraint;

    AuxiliaryVariables reformulationVariables;
};

class TaskReformulateProblem : public TaskBase
{
public:
    TaskReformulateProblem(EnvironmentPtr envPtr);
    ~TaskReformulateProblem() override;

    void run() override;
    std::string getType() override;

private:
    bool useQuadraticConstraints = false;
    bool useQuadraticObjective = true;
    bool quadraticObjectiveRegardedAsNonlinear = false;
    bool partitionQuadraticTermsInObjective = false;
    bool partitionQuadraticTermsInConstraint = false;

    void reformulateObjectiveFunction();
    NumericConstraints reformulateConstraint(NumericConstraintPtr constraint);

    template <class T> void copyLinearTermsToConstraint(LinearTerms terms, T destination, bool reversedSigns = false);

    template <class T>
    void copyQuadraticTermsToConstraint(QuadraticTerms terms, T destination, bool reversedSigns = false);

    template <class T>
    void copyMonomialTermsToConstraint(MonomialTerms terms, T destination, bool reversedSigns = false);

    template <class T>
    void copySignomialTermsToConstraint(SignomialTerms terms, T destination, bool reversedSigns = false);

    template <class T>
    void copyLinearTermsToObjectiveFunction(LinearTerms terms, T destination, bool reversedSigns = false);

    template <class T>
    void copyQuadraticTermsToObjectiveFunction(QuadraticTerms terms, T destination, bool reversedSigns = false);

    template <class T>
    void copyMonomialTermsToObjectiveFunction(MonomialTerms terms, T destination, bool reversedSigns = false);

    template <class T>
    void copySignomialTermsToObjectiveFunction(SignomialTerms terms, T destination, bool reversedSigns = false);

    LinearTerms partitionNonlinearSum(const std::shared_ptr<ExpressionSum> source, bool reversedSigns);
    LinearTerms partitionMonomialTerms(const MonomialTerms sourceTerms, bool reversedSigns);
    LinearTerms partitionSignomialTerms(const SignomialTerms sourceTerms, bool reversedSigns);

    LinearTerms partitionNonlinearBinaryProduct(const std::shared_ptr<ExpressionSum> source, bool reversedSigns);

    std::tuple<LinearTerms, QuadraticTerms> reformulateAndPartitionQuadraticSum(
        const QuadraticTerms& quadraticTerms, bool reversedSigns, bool partitionNonBinaryTerms);
    std::tuple<LinearTerms, MonomialTerms> reformulateMonomialSum(
        const MonomialTerms& monomialTerms, bool reversedSigns);

    AuxiliaryVariablePtr getBilinearAuxiliaryVariable(VariablePtr firstVariable, VariablePtr secondVariable);

    void addBilinearMcCormickEnvelope(
        AuxiliaryVariablePtr auxVariable, VariablePtr firstVariable, VariablePtr secondVariable);

    int auxVariableCounter = 0;
    int auxConstraintCounter = 0;

    std::map<VariablePtr, Variables> integerAuxiliaryBinaryVariables;

    std::map<std::tuple<VariablePtr, VariablePtr>, AuxiliaryVariablePtr> bilinearAuxVariables;

    ProblemPtr reformulatedProblem;
};
} // namespace SHOT