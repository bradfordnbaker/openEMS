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

#include "engine_ext_cpml.h"
#include "operator_ext_cpml.h"
#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"
#include "tools/useful.h"

Engine_Ext_CPML::Engine_Ext_CPML(Operator_Ext_CPML* op_ext) : Engine_Extension(op_ext)
{
	m_Op_CPML = op_ext;

	//this ABC extension should have high priority
	m_Priority = ENG_EXT_PRIO_CPML;

	psi_volt_x.Init("psi_volt_x", m_Op_CPML->m_numLines);
	psi_volt_y.Init("psi_volt_y", m_Op_CPML->m_numLines);
	psi_curr_x.Init("psi_curr_x", m_Op_CPML->m_numLines);
	psi_curr_y.Init("psi_curr_y", m_Op_CPML->m_numLines);

	SetNumberOfThreads(1);
}

Engine_Ext_CPML::~Engine_Ext_CPML()
{
}

void Engine_Ext_CPML::SetNumberOfThreads(int nrThread)
{
	Engine_Extension::SetNumberOfThreads(nrThread);

	m_numX = AssignJobs2Threads(m_Op_CPML->m_numLines[0],m_NrThreads,false);
	m_start.resize(m_NrThreads,0);
	m_start.at(0)=0;
	for (size_t n=1; n<m_numX.size(); ++n)
		m_start.at(n) = m_start.at(n-1) + m_numX.at(n-1);
}

template <typename EngType>
void Engine_Ext_CPML::DoPreVoltageUpdatesImpl(EngType* eng, int threadID)
{
	// CPML doesn't need pre-voltage updates
	// Psi fields are updated in post-voltage using the newly computed fields
	UNUSED(eng);
	UNUSED(threadID);
}

void Engine_Ext_CPML::DoPreVoltageUpdates(int threadID)
{
	// Not used for CPML
	UNUSED(threadID);
}

template <typename EngType>
void Engine_Ext_CPML::DoPostVoltageUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng==NULL)
		return;

	if (threadID>=m_NrThreads)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		loc_pos[0]=lineX+m_start.at(threadID);
		pos[0] = loc_pos[0] + m_Op_CPML->m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<m_Op_CPML->m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_Op_CPML->m_StartPos[1];
			for (loc_pos[2]=0; loc_pos[2]<m_Op_CPML->m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_Op_CPML->m_StartPos[2];

				// Update psi fields and apply CPML correction for Ex component
				if (pos[1] < m_Op_CPML->m_Op->GetNumberOfLines(1,true)-1)
				{
					// dHz/dy for Ex
					FDTD_FLOAT dHz_dy = (eng->EngType::GetCurr(2,pos) -
					                     eng->EngType::GetCurr(2,pos[0],pos[1]+1,pos[2]));

					// Update psi_Ex_y
					psi_volt_x(0, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_E(0, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_volt_x(0, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_E(0, loc_pos[0], loc_pos[1], loc_pos[2]) * dHz_dy;
				}

				if (pos[2] < m_Op_CPML->m_Op->GetNumberOfLines(2,true)-1)
				{
					// dHy/dz for Ex
					FDTD_FLOAT dHy_dz = (eng->EngType::GetCurr(1,pos) -
					                     eng->EngType::GetCurr(1,pos[0],pos[1],pos[2]+1));

					// Update psi_Ex_z
					psi_volt_y(0, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_E(0, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_volt_y(0, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_E(0, loc_pos[0], loc_pos[1], loc_pos[2]) * dHy_dz;
				}

				// Apply CPML correction to Ex
				FDTD_FLOAT Ex_new = eng->EngType::GetVolt(0,pos) +
					psi_volt_x(0, loc_pos[0], loc_pos[1], loc_pos[2]) +
					psi_volt_y(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetVolt(0,pos, Ex_new);

				// Update psi fields and apply CPML correction for Ey component
				if (pos[0] < m_Op_CPML->m_Op->GetNumberOfLines(0,true)-1)
				{
					// dHx/dz for Ey (note: this is stored in psi_volt_y since z-derivative)
					FDTD_FLOAT dHx_dz = (eng->EngType::GetCurr(0,pos) -
					                     eng->EngType::GetCurr(0,pos[0],pos[1],pos[2]+1));

					// Update psi_Ey_z
					psi_volt_y(1, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_E(1, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_volt_y(1, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_E(1, loc_pos[0], loc_pos[1], loc_pos[2]) * dHx_dz;
				}

				if (pos[2] < m_Op_CPML->m_Op->GetNumberOfLines(2,true)-1)
				{
					// dHz/dx for Ey
					FDTD_FLOAT dHz_dx = (eng->EngType::GetCurr(2,pos) -
					                     eng->EngType::GetCurr(2,pos[0]+1,pos[1],pos[2]));

					// Update psi_Ey_x
					psi_volt_x(1, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_E(1, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_volt_x(1, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_E(1, loc_pos[0], loc_pos[1], loc_pos[2]) * dHz_dx;
				}

				// Apply CPML correction to Ey
				FDTD_FLOAT Ey_new = eng->EngType::GetVolt(1,pos) +
					psi_volt_x(1, loc_pos[0], loc_pos[1], loc_pos[2]) +
					psi_volt_y(1, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetVolt(1,pos, Ey_new);

				// Update psi fields and apply CPML correction for Ez component
				if (pos[0] < m_Op_CPML->m_Op->GetNumberOfLines(0,true)-1)
				{
					// dHy/dx for Ez
					FDTD_FLOAT dHy_dx = (eng->EngType::GetCurr(1,pos) -
					                     eng->EngType::GetCurr(1,pos[0]+1,pos[1],pos[2]));

					// Update psi_Ez_x
					psi_volt_x(2, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_E(2, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_volt_x(2, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_E(2, loc_pos[0], loc_pos[1], loc_pos[2]) * dHy_dx;
				}

				if (pos[1] < m_Op_CPML->m_Op->GetNumberOfLines(1,true)-1)
				{
					// dHx/dy for Ez
					FDTD_FLOAT dHx_dy = (eng->EngType::GetCurr(0,pos) -
					                     eng->EngType::GetCurr(0,pos[0],pos[1]+1,pos[2]));

					// Update psi_Ez_y
					psi_volt_y(2, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_E(2, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_volt_y(2, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_E(2, loc_pos[0], loc_pos[1], loc_pos[2]) * dHx_dy;
				}

				// Apply CPML correction to Ez
				FDTD_FLOAT Ez_new = eng->EngType::GetVolt(2,pos) +
					psi_volt_x(2, loc_pos[0], loc_pos[1], loc_pos[2]) +
					psi_volt_y(2, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetVolt(2,pos, Ez_new);
			}
		}
	}
}

void Engine_Ext_CPML::DoPostVoltageUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostVoltageUpdatesImpl, threadID);
}

template <typename EngType>
void Engine_Ext_CPML::DoPreCurrentUpdatesImpl(EngType* eng, int threadID)
{
	// CPML doesn't need pre-current updates
	// Psi fields are updated in post-current using the newly computed fields
	UNUSED(eng);
	UNUSED(threadID);
}

void Engine_Ext_CPML::DoPreCurrentUpdates(int threadID)
{
	// Not used for CPML
	UNUSED(threadID);
}

template <typename EngType>
void Engine_Ext_CPML::DoPostCurrentUpdatesImpl(EngType* eng, int threadID)
{
	if (m_Eng==NULL)
		return;

	if (threadID>=m_NrThreads)
		return;

	unsigned int pos[3];
	unsigned int loc_pos[3];

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		loc_pos[0]=lineX+m_start.at(threadID);
		pos[0] = loc_pos[0] + m_Op_CPML->m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<m_Op_CPML->m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_Op_CPML->m_StartPos[1];
			for (loc_pos[2]=0; loc_pos[2]<m_Op_CPML->m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_Op_CPML->m_StartPos[2];

				// Update psi fields and apply CPML correction for Hx component
				if ((pos[1] > 0) && (pos[1] < m_Op_CPML->m_Op->GetNumberOfLines(1,true)-1))
				{
					// dEz/dy for Hx
					FDTD_FLOAT dEz_dy = (eng->EngType::GetVolt(2,pos[0],pos[1]-1,pos[2]) -
					                     eng->EngType::GetVolt(2,pos));

					// Update psi_Hx_y
					psi_curr_x(0, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_H(0, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_curr_x(0, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_H(0, loc_pos[0], loc_pos[1], loc_pos[2]) * dEz_dy;
				}

				if ((pos[2] > 0) && (pos[2] < m_Op_CPML->m_Op->GetNumberOfLines(2,true)-1))
				{
					// dEy/dz for Hx
					FDTD_FLOAT dEy_dz = (eng->EngType::GetVolt(1,pos[0],pos[1],pos[2]-1) -
					                     eng->EngType::GetVolt(1,pos));

					// Update psi_Hx_z
					psi_curr_y(0, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_H(0, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_curr_y(0, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_H(0, loc_pos[0], loc_pos[1], loc_pos[2]) * dEy_dz;
				}

				// Apply CPML correction to Hx
				FDTD_FLOAT Hx_new = eng->EngType::GetCurr(0,pos) +
					psi_curr_x(0, loc_pos[0], loc_pos[1], loc_pos[2]) +
					psi_curr_y(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetCurr(0,pos, Hx_new);

				// Update psi fields and apply CPML correction for Hy component
				if ((pos[0] > 0) && (pos[0] < m_Op_CPML->m_Op->GetNumberOfLines(0,true)-1))
				{
					// dEx/dz for Hy
					FDTD_FLOAT dEx_dz = (eng->EngType::GetVolt(0,pos[0],pos[1],pos[2]-1) -
					                     eng->EngType::GetVolt(0,pos));

					// Update psi_Hy_z
					psi_curr_y(1, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_H(1, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_curr_y(1, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_H(1, loc_pos[0], loc_pos[1], loc_pos[2]) * dEx_dz;
				}

				if ((pos[2] > 0) && (pos[2] < m_Op_CPML->m_Op->GetNumberOfLines(2,true)-1))
				{
					// dEz/dx for Hy
					FDTD_FLOAT dEz_dx = (eng->EngType::GetVolt(2,pos[0]-1,pos[1],pos[2]) -
					                     eng->EngType::GetVolt(2,pos));

					// Update psi_Hy_x
					psi_curr_x(1, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_H(1, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_curr_x(1, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_H(1, loc_pos[0], loc_pos[1], loc_pos[2]) * dEz_dx;
				}

				// Apply CPML correction to Hy
				FDTD_FLOAT Hy_new = eng->EngType::GetCurr(1,pos) +
					psi_curr_x(1, loc_pos[0], loc_pos[1], loc_pos[2]) +
					psi_curr_y(1, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetCurr(1,pos, Hy_new);

				// Update psi fields and apply CPML correction for Hz component
				if ((pos[0] > 0) && (pos[0] < m_Op_CPML->m_Op->GetNumberOfLines(0,true)-1))
				{
					// dEy/dx for Hz
					FDTD_FLOAT dEy_dx = (eng->EngType::GetVolt(1,pos[0]-1,pos[1],pos[2]) -
					                     eng->EngType::GetVolt(1,pos));

					// Update psi_Hz_x
					psi_curr_x(2, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_H(2, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_curr_x(2, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_H(2, loc_pos[0], loc_pos[1], loc_pos[2]) * dEy_dx;
				}

				if ((pos[1] > 0) && (pos[1] < m_Op_CPML->m_Op->GetNumberOfLines(1,true)-1))
				{
					// dEx/dy for Hz
					FDTD_FLOAT dEx_dy = (eng->EngType::GetVolt(0,pos[0],pos[1]-1,pos[2]) -
					                     eng->EngType::GetVolt(0,pos));

					// Update psi_Hz_y
					psi_curr_y(2, loc_pos[0], loc_pos[1], loc_pos[2]) =
						m_Op_CPML->b_H(2, loc_pos[0], loc_pos[1], loc_pos[2]) *
						psi_curr_y(2, loc_pos[0], loc_pos[1], loc_pos[2]) +
						m_Op_CPML->a_H(2, loc_pos[0], loc_pos[1], loc_pos[2]) * dEx_dy;
				}

				// Apply CPML correction to Hz
				FDTD_FLOAT Hz_new = eng->EngType::GetCurr(2,pos) +
					psi_curr_x(2, loc_pos[0], loc_pos[1], loc_pos[2]) +
					psi_curr_y(2, loc_pos[0], loc_pos[1], loc_pos[2]);
				eng->EngType::SetCurr(2,pos, Hz_new);
			}
		}
	}
}

void Engine_Ext_CPML::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}
