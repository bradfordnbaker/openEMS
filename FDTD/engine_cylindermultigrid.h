/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ENGINE_CYLINDERMULTIGRID_H
#define ENGINE_CYLINDERMULTIGRID_H

#include "engine_sse.h"

class Operator_CylinderMultiGrid;

//! The SSE engine of a cylindrical multi-grid level with its sub-grid
/*!
  The sub-grid engine is stepped by the threads of this engine (Engine_ThreadPool):
  every half step this grid, then the sub-grid, then the coupling of the two
  (SyncSubGridVoltages(), SyncSubGridCurrents()); after IterateTS() this grid is
  filled inside the split radius from the sub-grid (InterpolateFromSubGrid()).
  */
class Engine_CylinderMultiGrid : public Engine_sse
{
public:
	static Engine_CylinderMultiGrid* New(const Operator_CylinderMultiGrid* op, unsigned int numThreads = 0);
	virtual ~Engine_CylinderMultiGrid();

	virtual Engine* GetSubGridEngine() const {return m_InnerEngine;}
	//! the voltages of the sub-grid on the split radius from this grid
	virtual void SyncSubGridVoltages();
	//! the currents of this grid next to the split radius from the sub-grid
	virtual void SyncSubGridCurrents();
	virtual void InterpolateFromSubGrid();

	virtual void InterpolVoltChild2Base(unsigned int rPos);
	virtual void InterpolCurrChild2Base(unsigned int rPos);

protected:
	Engine_CylinderMultiGrid(const Operator_CylinderMultiGrid* op);
	const Operator_CylinderMultiGrid* Op_CMG;

	Engine_sse* m_InnerEngine;
};

#endif // ENGINE_CYLINDERMULTIGRID_H
