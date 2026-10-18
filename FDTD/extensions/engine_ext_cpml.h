/*
*	Copyright (C) 2025 openEMS Contributors
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

#ifndef ENGINE_EXT_CPML_H
#define ENGINE_EXT_CPML_H

#include <vector>
#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "engine_extension_dispatcher.h"

class Operator_Ext_CPML;

class Engine_Ext_CPML : public Engine_Extension
{
public:
	Engine_Ext_CPML(Operator_Ext_CPML* op_ext);
	virtual ~Engine_Ext_CPML();

	virtual void SetNumberOfThreads(int nrThread);

	virtual void DoPreVoltageUpdates() {Engine_Ext_CPML::DoPreVoltageUpdates(0);};
	virtual void DoPreVoltageUpdates(int threadID);
	virtual void DoPostVoltageUpdates() {Engine_Ext_CPML::DoPostVoltageUpdates(0);};
	virtual void DoPostVoltageUpdates(int threadID);

	virtual void DoPreCurrentUpdates() {Engine_Ext_CPML::DoPreCurrentUpdates(0);};
	virtual void DoPreCurrentUpdates(int threadID);
	virtual void DoPostCurrentUpdates() {Engine_Ext_CPML::DoPostCurrentUpdates(0);};
	virtual void DoPostCurrentUpdates(int threadID);

protected:
	template <typename EngType>
	void DoPreVoltageUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void DoPostVoltageUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void DoPreCurrentUpdatesImpl(EngType* eng, int threadID);

	template <typename EngType>
	void DoPostCurrentUpdatesImpl(EngType* eng, int threadID);

	Operator_Ext_CPML* m_Op_CPML;

	std::vector<unsigned int> m_start;
	std::vector<unsigned int> m_numX;

	// Per cell data for the two field components tangential to the slab
	// normal, packed for the update loop. ca and ck include the curl sign and
	// vi (iv), so psi is stored pre-scaled: field += psi + ck*diff
	struct Cell
	{
		FDTD_FLOAT b;
		FDTD_FLOAT ca[2];
		FDTD_FLOAT ck[2];
		FDTD_FLOAT psi[2];
	};
	std::vector<Cell> m_volt;
	std::vector<Cell> m_curr;

	// tangential components and the field component differenced along the normal
	int m_comp[2];
	int m_diff[2];

	void InitCells();
	inline size_t CellIndex(const unsigned int loc_pos[3]) const;
};

#endif // ENGINE_EXT_CPML_H
