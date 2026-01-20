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

	vector<unsigned int> m_start;
	vector<unsigned int> m_numX;

	// CPML auxiliary psi fields for E-field (voltage) updates
	// psi_volt_x stores psi fields for derivatives in y-direction
	// psi_volt_y stores psi fields for derivatives in z-direction
	// First index: field component (0=Ex, 1=Ey, 2=Ez)
	ArrayLib::ArrayNIJK<FDTD_FLOAT> psi_volt_x;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> psi_volt_y;

	// CPML auxiliary psi fields for H-field (current) updates
	// psi_curr_x stores psi fields for derivatives in y-direction
	// psi_curr_y stores psi fields for derivatives in z-direction
	// First index: field component (0=Hx, 1=Hy, 2=Hz)
	ArrayLib::ArrayNIJK<FDTD_FLOAT> psi_curr_x;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> psi_curr_y;
};

#endif // ENGINE_EXT_CPML_H
