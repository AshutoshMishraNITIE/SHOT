#pragma once
#include "ILinesearchMethod.h"
#include "SHOTSettings.h"

//#include "..\OptProblems\OptProblemReformulated.h"
class LinesearchMethodBisection: public ILinesearchMethod
{
	public:
		LinesearchMethodBisection();
		virtual ~LinesearchMethodBisection();

		virtual std::pair<std::vector<double>, std::vector<double>> findZero(std::vector<double> ptA,
				std::vector<double> ptB, int Nmax, double lambdaTol, double constrTol);

		virtual std::pair<std::vector<double>, std::vector<double>> findZero(std::vector<double> ptA,
				std::vector<double> ptB, int Nmax, double lambdaTol, double constrTol, std::vector<int> constrIdxs);
	private:
		ProcessInfo* processInfo;
		SHOTSettings::Settings *settings;
};
