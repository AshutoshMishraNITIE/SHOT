/*
 * TaskCheckIterationLimit.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: alundell
 */

#include "TaskCheckIterationLimit.h"

TaskCheckIterationLimit::TaskCheckIterationLimit(std::string taskIDTrue)
{

	taskIDIfTrue = taskIDTrue;
}

TaskCheckIterationLimit::~TaskCheckIterationLimit()
{
	// TODO Auto-generated destructor stub
}

void TaskCheckIterationLimit::run()
{
	auto currIter = ProcessInfo::getInstance().getCurrentIteration();

	if (currIter->iterationNumber
			>= Settings::getInstance().getIntSetting("Relaxation.IterationLimit", "Dual")
					+ Settings::getInstance().getIntSetting("IterationLimit", "Termination"))
	{
		ProcessInfo::getInstance().terminationReason = E_TerminationReason::IterationLimit;
		ProcessInfo::getInstance().tasks->setNextTask(taskIDIfTrue);
	}
}

std::string TaskCheckIterationLimit::getType()
{
	std::string type = typeid(this).name();
	return (type);

}
