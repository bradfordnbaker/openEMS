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

	psi_volt.Init("psi_volt", m_Op_CPML->m_numLines);
	psi_curr.Init("psi_curr", m_Op_CPML->m_numLines);

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

	const int d = m_Op_CPML->m_dir;
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

				unsigned int prev[3] = {pos[0], pos[1], pos[2]};
				if (pos[d]>0)
					--prev[d];

				const FDTD_FLOAT b = m_Op_CPML->c_E(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				const FDTD_FLOAT a = m_Op_CPML->c_E(1, loc_pos[0], loc_pos[1], loc_pos[2]);
				const FDTD_FLOAT k = m_Op_CPML->c_E(2, loc_pos[0], loc_pos[1], loc_pos[2]);

				// E_c has +d(H_c+2)/d(c+1) - d(H_c+1)/d(c+2), same as Engine::UpdateVoltages
				for (int c=0; c<3; ++c)
				{
					if (c==d)
						continue;
					const int e = 3 - c - d;
					const FDTD_FLOAT sign = (d==(c+1)%3) ? 1 : -1;
					FDTD_FLOAT dI = sign * (eng->EngType::GetCurr(e,pos) - eng->EngType::GetCurr(e,prev));

					FDTD_FLOAT& psi = psi_volt(c, loc_pos[0], loc_pos[1], loc_pos[2]);
					psi = b * psi + a * dI;

					eng->EngType::SetVolt(c, pos, eng->EngType::GetVolt(c,pos) +
						m_Op_CPML->m_Op->GetVI(c,pos) * (psi + k * dI));
				}
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

	const int d = m_Op_CPML->m_dir;
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = m_Op_CPML->m_Op->GetNumberOfLines(n,true);

	unsigned int pos[3];
	unsigned int loc_pos[3];

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		loc_pos[0]=lineX+m_start.at(threadID);
		pos[0] = loc_pos[0] + m_Op_CPML->m_StartPos[0];
		if (pos[0]>=numLines[0]-1)
			continue;
		for (loc_pos[1]=0; loc_pos[1]<m_Op_CPML->m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_Op_CPML->m_StartPos[1];
			if (pos[1]>=numLines[1]-1)
				continue;
			for (loc_pos[2]=0; loc_pos[2]<m_Op_CPML->m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_Op_CPML->m_StartPos[2];
				// Engine::UpdateCurrents does not touch the last line
				if (pos[2]>=numLines[2]-1)
					continue;

				unsigned int next[3] = {pos[0], pos[1], pos[2]};
				++next[d];

				const FDTD_FLOAT b = m_Op_CPML->c_H(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				const FDTD_FLOAT a = m_Op_CPML->c_H(1, loc_pos[0], loc_pos[1], loc_pos[2]);
				const FDTD_FLOAT k = m_Op_CPML->c_H(2, loc_pos[0], loc_pos[1], loc_pos[2]);

				// H_c has -d(E_c+2)/d(c+1) + d(E_c+1)/d(c+2), same as Engine::UpdateCurrents
				for (int c=0; c<3; ++c)
				{
					if (c==d)
						continue;
					const int e = 3 - c - d;
					const FDTD_FLOAT sign = (d==(c+1)%3) ? -1 : 1;
					FDTD_FLOAT dV = sign * (eng->EngType::GetVolt(e,next) - eng->EngType::GetVolt(e,pos));

					FDTD_FLOAT& psi = psi_curr(c, loc_pos[0], loc_pos[1], loc_pos[2]);
					psi = b * psi + a * dV;

					eng->EngType::SetCurr(c, pos, eng->EngType::GetCurr(c,pos) +
						m_Op_CPML->m_Op->GetIV(c,pos) * (psi + k * dV));
				}
			}
		}
	}
}

void Engine_Ext_CPML::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}
