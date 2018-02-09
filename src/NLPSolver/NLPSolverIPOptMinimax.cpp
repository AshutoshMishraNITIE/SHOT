#include "NLPSolverIPOptMinimax.h"

NLPSolverIPOptMinimax::NLPSolverIPOptMinimax()
{

	osolwriter = new OSoLWriter();

	NLPProblem = new OptProblemNLPMinimax();

	setInitialSettings();
}

NLPSolverIPOptMinimax::~NLPSolverIPOptMinimax()
{
	delete NLPSolver;
	delete osolwriter;
	delete NLPProblem;
}

std::vector<double> NLPSolverIPOptMinimax::getSolution()
{
	int numVar = NLPProblem->getNumberOfVariables();
	std::vector<double> tmpPoint(numVar);

	for (int i = 0; i < numVar; i++)
	{
		tmpPoint.at(i) = NLPSolverIPOptBase::getSolution(i);
	}

	if (ProcessInfo::getInstance().originalProblem->getObjectiveFunctionType() == E_ObjectiveFunctionType::Quadratic)
	{
		tmpPoint.pop_back();
	}

	// Removes the minimax objective variable from solution
	tmpPoint.pop_back();

	return (tmpPoint);
}

bool NLPSolverIPOptMinimax::createProblemInstance(OSInstance * origInstance)
{

	ProcessInfo::getInstance().outputInfo("     Creating Ipopt minimax problem.");

	dynamic_cast<OptProblemNLPMinimax*>(NLPProblem)->reformulate(origInstance);

	ProcessInfo::getInstance().outputInfo("     Ipopt minimax problem created.");

	return (true);
}

void NLPSolverIPOptMinimax::setSolverSpecificInitialSettings()
{
	auto constrTol = Settings::getInstance().getDoubleSetting("Ipopt.ConstraintViolationTolerance", "Subsolver");
	osOption->setAnotherSolverOption("constr_viol_tol", UtilityFunctions::toStringFormat(constrTol, "%.10f"), "ipopt",
			"", "double", "");

	osOption->setAnotherSolverOption("tol",
			UtilityFunctions::toStringFormat(
					Settings::getInstance().getDoubleSetting("Ipopt.RelativeConvergenceTolerance", "Subsolver"), "%.10f"),
			"ipopt", "", "double", "");

	osOption->setAnotherSolverOption("max_iter",
			to_string(Settings::getInstance().getIntSetting("Ipopt.MaxIterations", "Subsolver")), "ipopt", "",
			"integer", "");

	auto timeLimit = Settings::getInstance().getDoubleSetting("FixedInteger.TimeLimit", "Primal");
	osOption->setAnotherSolverOption("max_cpu_time", UtilityFunctions::toStringFormat(timeLimit, "%.10f"), "ipopt", "",
			"number", "");
}
