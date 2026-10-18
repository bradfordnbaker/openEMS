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

	InitCells();

	SetNumberOfThreads(1);
}

inline size_t Engine_Ext_CPML::CellIndex(const unsigned int loc_pos[3]) const
{
	const unsigned int* n = m_Op_CPML->m_numLines;
	return ((size_t)loc_pos[0]*n[1] + loc_pos[1])*n[2] + loc_pos[2];
}

void Engine_Ext_CPML::InitCells()
{
	// Done here and not in the operator extension: other extensions (lumped
	// elements, sheets, hard sources) may still change vi/iv during their build.
	const int d = m_Op_CPML->m_dir;
	m_comp[0] = (d+1)%3;
	m_diff[0] = (d+2)%3;
	m_comp[1] = (d+2)%3;
	m_diff[1] = (d+1)%3;
	// sign of the normal derivative term in Engine::UpdateVoltages, UpdateCurrents has the opposite
	const FDTD_FLOAT sign_volt[2] = {-1, 1};

	const unsigned int* n = m_Op_CPML->m_numLines;
	m_volt.assign((size_t)n[0]*n[1]*n[2], Cell());
	m_curr.assign((size_t)n[0]*n[1]*n[2], Cell());

	Operator* op = m_Op_CPML->m_Op;
	unsigned int pos[3];
	unsigned int loc_pos[3];
	for (loc_pos[0]=0; loc_pos[0]<n[0]; ++loc_pos[0])
	{
		pos[0] = loc_pos[0] + m_Op_CPML->m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<n[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_Op_CPML->m_StartPos[1];
			for (loc_pos[2]=0; loc_pos[2]<n[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_Op_CPML->m_StartPos[2];
				Cell& v = m_volt[CellIndex(loc_pos)];
				Cell& i = m_curr[CellIndex(loc_pos)];
				v.b = m_Op_CPML->c_E(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				i.b = m_Op_CPML->c_H(0, loc_pos[0], loc_pos[1], loc_pos[2]);
				for (int t=0; t<2; ++t)
				{
					FDTD_FLOAT vi = sign_volt[t] * op->GetVI(m_comp[t], pos);
					FDTD_FLOAT iv = -sign_volt[t] * op->GetIV(m_comp[t], pos);
					v.ca[t] = vi * m_Op_CPML->c_E(1, loc_pos[0], loc_pos[1], loc_pos[2]);
					v.ck[t] = vi * m_Op_CPML->c_E(2, loc_pos[0], loc_pos[1], loc_pos[2]);
					i.ca[t] = iv * m_Op_CPML->c_H(1, loc_pos[0], loc_pos[1], loc_pos[2]);
					i.ck[t] = iv * m_Op_CPML->c_H(2, loc_pos[0], loc_pos[1], loc_pos[2]);
					v.psi[t] = 0;
					i.psi[t] = 0;
				}
			}
		}
	}
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

				Cell& cell = m_volt[CellIndex(loc_pos)];
				for (int t=0; t<2; ++t)
				{
					const int c = m_comp[t];
					FDTD_FLOAT diff = eng->EngType::GetCurr(m_diff[t],pos) - eng->EngType::GetCurr(m_diff[t],prev);
					cell.psi[t] = cell.b * cell.psi[t] + cell.ca[t] * diff;
					eng->EngType::SetVolt(c, pos, eng->EngType::GetVolt(c,pos) + cell.psi[t] + cell.ck[t] * diff);
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

				Cell& cell = m_curr[CellIndex(loc_pos)];
				for (int t=0; t<2; ++t)
				{
					const int c = m_comp[t];
					FDTD_FLOAT diff = eng->EngType::GetVolt(m_diff[t],next) - eng->EngType::GetVolt(m_diff[t],pos);
					cell.psi[t] = cell.b * cell.psi[t] + cell.ca[t] * diff;
					eng->EngType::SetCurr(c, pos, eng->EngType::GetCurr(c,pos) + cell.psi[t] + cell.ck[t] * diff);
				}
			}
		}
	}
}

void Engine_Ext_CPML::DoPostCurrentUpdates(int threadID)
{
	ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}
