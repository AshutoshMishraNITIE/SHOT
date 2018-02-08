/*
 * TaskAddHyperplanes.cpp
 *
 *  Created on: Mar 28, 2015
 *      Author: alundell
 */

#include "TaskAddHyperplanes.h"

TaskAddHyperplanes::TaskAddHyperplanes(IMILPSolver *MILPSolver)
{

	itersWithoutAddedHPs = 0;

	this->MILPSolver = MILPSolver;
}

TaskAddHyperplanes::~TaskAddHyperplanes()
{
	// TODO Auto-generated destructor stub
}

void TaskAddHyperplanes::run()
{
	this->MILPSolver = MILPSolver;
	auto currIter = ProcessInfo::getInstance().getCurrentIteration(); // The unsolved new iteration

	if (!currIter->isMILP() || !Settings::getInstance().getBoolSetting("HyperplaneCuts.Delay", "Dual")
			|| !currIter->MILPSolutionLimitUpdated || itersWithoutAddedHPs > 5)
	{
		int addedHyperplanes = 0;

		for (int k = ProcessInfo::getInstance().hyperplaneWaitingList.size(); k > 0; k--)
		{
			if (addedHyperplanes >= Settings::getInstance().getIntSetting("HyperplaneCuts.MaxPerIteration", "Dual")) break;

			auto tmpItem = ProcessInfo::getInstance().hyperplaneWaitingList.at(k - 1);

			if (tmpItem.source == E_HyperplaneSource::PrimalSolutionSearchInteriorObjective)
			{
				MILPSolver->createInteriorHyperplane(tmpItem);
			}
			else
			{
				MILPSolver->createHyperplane(tmpItem);

				ProcessInfo::getInstance().addedHyperplanes.push_back(tmpItem);

				addedHyperplanes++;
			}

		}

		ProcessInfo::getInstance().hyperplaneWaitingList.clear();
		itersWithoutAddedHPs = 0;
	}
	else
	{
		itersWithoutAddedHPs++;
	}
}

std::string TaskAddHyperplanes::getType()
{
	std::string type = typeid(this).name();
	return (type);

}

