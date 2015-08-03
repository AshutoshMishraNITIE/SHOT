﻿#pragma once
#include "IMILPSolver.h"
#include "MILPSolverBase.h"

class MILPSolverOsiCbc: public IMILPSolver, MILPSolverBase
{
	public:
		MILPSolverOsiCbc();
		virtual ~MILPSolverOsiCbc();

		virtual bool createLinearProblem(OptProblem *origProblem);
		virtual bool addLinearConstraint(std::vector<IndexValuePair> elements, int numNonZero, double constant);
		virtual std::vector<double> getVariableSolution();
		virtual void activateDiscreteVariables(bool activate);
		virtual bool getDiscreteVariableStatus();
		virtual E_ProblemSolutionStatus solveProblem();
		virtual E_ProblemSolutionStatus getSolutionStatus();
		virtual double getLastObjectiveValue();
		virtual double getBestObjectiveValue();

		virtual int increaseSolutionLimit(int increment);
		virtual void setSolutionLimit(int limit);
		virtual int getSolutionLimit();

		virtual void writeProblemToFile(std::string filename);

		virtual std::vector<SolutionPoint> getAllVariableSolutions();

		virtual void setTimeLimit(double seconds);

		virtual void setCutOff(double cutOff);

		virtual void addMIPStart(std::vector<double> point);
		virtual void deleteMIPStarts();
		virtual void changeConstraintToLazy(std::vector<int> constrIdxs);

		virtual void fixVariable(int varIndex, double value);
		virtual void updateVariableBound(int varIndex, double lowerBound, double upperBound);
		virtual pair<double, double> getCurrentVariableBounds(int varIndex);
		virtual void createHyperplane(int constrIdx, std::vector<double> point);
};
