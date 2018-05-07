/**
   The Supporting Hyperplane Optimization Toolkit (SHOT).

   @author Andreas Lundell, Åbo Akademi University

   @section LICENSE 
   This software is licensed under the Eclipse Public License 2.0. 
   Please see the README and LICENSE files for more information.
*/

#include "TaskCalculateSolutionChangeNorm.h"

TaskCalculateSolutionChangeNorm::TaskCalculateSolutionChangeNorm()
{
}

TaskCalculateSolutionChangeNorm::~TaskCalculateSolutionChangeNorm()
{
}

void TaskCalculateSolutionChangeNorm::run()
{
	auto currIter = ProcessInfo::getInstance().getCurrentIteration();

	currIter->boundaryDistance = OSDBL_MAX;

	if (ProcessInfo::getInstance().iterations.size() < 3)
	{
		return;
	}

	if (ProcessInfo::getInstance().getCurrentIteration()->hyperplanePoints.size() == 0 || ProcessInfo::getInstance().getCurrentIteration()->isMIP())
	{
		return;
	}

	auto currIterSol = ProcessInfo::getInstance().getCurrentIteration()->hyperplanePoints.at(0);

	for (int i = ProcessInfo::getInstance().iterations.size() - 2; i >= 1; i--)
	{
		if (ProcessInfo::getInstance().iterations.size() > 0 && !ProcessInfo::getInstance().iterations.at(i).isMIP())
		{
			auto prevIterSol = ProcessInfo::getInstance().iterations.at(i).hyperplanePoints.at(0);

			double distance = 0;

			for (int j = 0; j < currIterSol.size(); j++)
			{
				distance = distance + (currIterSol.at(j) - prevIterSol.at(j)) * (currIterSol.at(j) - prevIterSol.at(j));
			}

			distance = sqrt(distance + 0.001);

			if (OSIsnan(distance)) // Checks for INF, do not remove!
			{
				currIter->boundaryDistance = OSDBL_MAX;
			}
			else
			{
				currIter->boundaryDistance = distance;
			}

			return;
		}
	}
}

std::string TaskCalculateSolutionChangeNorm::getType()
{
	std::string type = typeid(this).name();
	return (type);
}
