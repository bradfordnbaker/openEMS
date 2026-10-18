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

	m_vec = false;
	m_numVectors = 0;

	SetNumberOfThreads(1);
}

inline size_t Engine_Ext_CPML::CellIndex(const unsigned int loc_pos[3]) const
{
	const unsigned int* n = m_Op_CPML->m_numLines;
	return ((size_t)loc_pos[0]*n[1] + loc_pos[1])*n[2] + loc_pos[2];
}

void Engine_Ext_CPML::SetEngine(Engine* eng)
{
	Engine_Extension::SetEngine(eng);
	InitCells();
}

void Engine_Ext_CPML::CalcCell(const unsigned int loc_pos[3], const unsigned int pos[3], bool isH,
                               FDTD_FLOAT& b, FDTD_FLOAT ca[2], FDTD_FLOAT ck[2]) const
{
	// sign of the normal derivative term in Engine::UpdateVoltages, UpdateCurrents has the opposite
	const FDTD_FLOAT sign_volt[2] = {-1, 1};
	const ArrayLib::ArrayNIJK<FDTD_FLOAT>& c = isH ? m_Op_CPML->c_H : m_Op_CPML->c_E;
	const Operator* op = m_Op_CPML->m_Op;

	b = c(0, loc_pos[0], loc_pos[1], loc_pos[2]);
	for (int t=0; t<2; ++t)
	{
		FDTD_FLOAT s;
		if (isH)
			s = -sign_volt[t] * op->GetIV(m_comp[t], pos[0], pos[1], pos[2]);
		else
			s = sign_volt[t] * op->GetVI(m_comp[t], pos[0], pos[1], pos[2]);
		ca[t] = s * c(1, loc_pos[0], loc_pos[1], loc_pos[2]);
		ck[t] = s * c(2, loc_pos[0], loc_pos[1], loc_pos[2]);
	}
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

	const unsigned int* n = m_Op_CPML->m_numLines;
	const unsigned int* start = m_Op_CPML->m_StartPos;
	const unsigned int Nz = m_Op_CPML->m_Op->GetNumberOfLines(2,true);

	Engine_sse* eng_sse = NULL;
	if (m_Eng && (m_Eng->GetType()==Engine::SSE))
		eng_sse = (Engine_sse*)m_Eng;
	m_vec = eng_sse && (d<2) && (start[2]==0) && (n[2]==Nz);
	// as in Engine_sse::Init(), which allocates its arrays only after the extensions are set up
	const unsigned int nv = (Nz+3)/4;

	m_volt.clear();
	m_curr.clear();
	m_volt4.clear();
	m_curr4.clear();
	m_zv.clear();
	m_lane.clear();

	unsigned int pos[3];
	unsigned int loc_pos[3];
	FDTD_FLOAT b, ca[2], ck[2];

	if (m_vec)
	{
		m_numVectors = nv;
		// padding lanes and the unused last H line keep zero coefficients
		m_volt4.assign((size_t)n[0]*n[1]*m_numVectors, Cell4());
		m_curr4.assign((size_t)n[0]*n[1]*m_numVectors, Cell4());
		for (loc_pos[0]=0; loc_pos[0]<n[0]; ++loc_pos[0])
		{
			pos[0] = loc_pos[0] + start[0];
			for (loc_pos[1]=0; loc_pos[1]<n[1]; ++loc_pos[1])
			{
				pos[1] = loc_pos[1] + start[1];
				for (unsigned int zv=0; zv<m_numVectors; ++zv)
				{
					size_t idx = ((size_t)loc_pos[0]*n[1] + loc_pos[1])*m_numVectors + zv;
					for (unsigned int lane=0; lane<4; ++lane)
					{
						pos[2] = loc_pos[2] = zv + lane*m_numVectors;
						if (pos[2]>=Nz)
							continue;
						CalcCell(loc_pos, pos, false, b, ca, ck);
						m_volt4[idx].b.f[lane] = b;
						for (int t=0; t<2; ++t)
						{
							m_volt4[idx].ca[t].f[lane] = ca[t];
							m_volt4[idx].ck[t].f[lane] = ck[t];
						}
						if (pos[2]>=Nz-1)
							continue;
						CalcCell(loc_pos, pos, true, b, ca, ck);
						m_curr4[idx].b.f[lane] = b;
						for (int t=0; t<2; ++t)
						{
							m_curr4[idx].ca[t].f[lane] = ca[t];
							m_curr4[idx].ck[t].f[lane] = ck[t];
						}
					}
				}
			}
		}
		return;
	}

	if (eng_sse)
	{
		m_zv.resize(Nz);
		m_lane.resize(Nz);
		for (unsigned int z=0; z<Nz; ++z)
		{
			m_zv[z] = z%nv;
			m_lane[z] = z/nv;
		}
	}

	m_volt.assign((size_t)n[0]*n[1]*n[2], Cell());
	m_curr.assign((size_t)n[0]*n[1]*n[2], Cell());
	for (loc_pos[0]=0; loc_pos[0]<n[0]; ++loc_pos[0])
	{
		pos[0] = loc_pos[0] + start[0];
		for (loc_pos[1]=0; loc_pos[1]<n[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + start[1];
			for (loc_pos[2]=0; loc_pos[2]<n[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + start[2];
				Cell& v = m_volt[CellIndex(loc_pos)];
				Cell& i = m_curr[CellIndex(loc_pos)];
				CalcCell(loc_pos, pos, false, v.b, v.ca, v.ck);
				CalcCell(loc_pos, pos, true, i.b, i.ca, i.ck);
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

void Engine_Ext_CPML::DoPostVoltageUpdatesSSE(Engine_sse* eng, int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	ArrayLib::ArrayENG<f4vector>& volt = *eng->f4_volt_ptr;
	ArrayLib::ArrayENG<f4vector>& curr = *eng->f4_curr_ptr;
	const int d = m_Op_CPML->m_dir;
	const unsigned int* n = m_Op_CPML->m_numLines;
	const unsigned int* start = m_Op_CPML->m_StartPos;
	unsigned int pos[3];
	unsigned int prev[3];

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		unsigned int loc_x = lineX + m_start.at(threadID);
		pos[0] = loc_x + start[0];
		for (unsigned int loc_y=0; loc_y<n[1]; ++loc_y)
		{
			pos[1] = loc_y + start[1];
			prev[0] = pos[0];
			prev[1] = pos[1];
			if ((d<2) && (pos[d]>0))
				--prev[d];

			if (m_vec)
			{
				Cell4* cell = &m_volt4[((size_t)loc_x*n[1] + loc_y)*m_numVectors];
				for (unsigned int zv=0; zv<m_numVectors; ++zv, ++cell)
				{
					for (int t=0; t<2; ++t)
					{
						v4sf diff = curr(m_diff[t], pos[0], pos[1], zv).v - curr(m_diff[t], prev[0], prev[1], zv).v;
						cell->psi[t].v = cell->b.v * cell->psi[t].v + cell->ca[t].v * diff;
						volt(m_comp[t], pos[0], pos[1], zv).v += cell->psi[t].v + cell->ck[t].v * diff;
					}
				}
				continue;
			}

			// z slab, differences along z
			Cell* cell = &m_volt[((size_t)loc_x*n[1] + loc_y)*n[2]];
			for (unsigned int loc_z=0; loc_z<n[2]; ++loc_z, ++cell)
			{
				unsigned int z = loc_z + start[2];
				unsigned int zp = z>0 ? z-1 : z;
				unsigned int zv = m_zv[z], l = m_lane[z];
				unsigned int zvp = m_zv[zp], lp = m_lane[zp];
				for (int t=0; t<2; ++t)
				{
					FDTD_FLOAT diff = curr(m_diff[t], pos[0], pos[1], zv).f[l] - curr(m_diff[t], pos[0], pos[1], zvp).f[lp];
					cell->psi[t] = cell->b * cell->psi[t] + cell->ca[t] * diff;
					volt(m_comp[t], pos[0], pos[1], zv).f[l] += cell->psi[t] + cell->ck[t] * diff;
				}
			}
		}
	}
}

void Engine_Ext_CPML::DoPostVoltageUpdates(int threadID)
{
	if (m_Eng && (m_Eng->GetType()==Engine::SSE))
		DoPostVoltageUpdatesSSE((Engine_sse*)m_Eng, threadID);
	else
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

void Engine_Ext_CPML::DoPostCurrentUpdatesSSE(Engine_sse* eng, int threadID)
{
	if (threadID>=m_NrThreads)
		return;

	ArrayLib::ArrayENG<f4vector>& volt = *eng->f4_volt_ptr;
	ArrayLib::ArrayENG<f4vector>& curr = *eng->f4_curr_ptr;
	const int d = m_Op_CPML->m_dir;
	const unsigned int* n = m_Op_CPML->m_numLines;
	const unsigned int* start = m_Op_CPML->m_StartPos;
	unsigned int numLines[3];
	for (int i=0; i<3; ++i)
		numLines[i] = m_Op_CPML->m_Op->GetNumberOfLines(i,true);
	unsigned int pos[3];
	unsigned int next[3];

	for (unsigned int lineX=0; lineX<m_numX.at(threadID); ++lineX)
	{
		unsigned int loc_x = lineX + m_start.at(threadID);
		pos[0] = loc_x + start[0];
		// Engine::UpdateCurrents does not touch the last line
		if (pos[0]>=numLines[0]-1)
			continue;
		for (unsigned int loc_y=0; loc_y<n[1]; ++loc_y)
		{
			pos[1] = loc_y + start[1];
			if (pos[1]>=numLines[1]-1)
				continue;
			next[0] = pos[0];
			next[1] = pos[1];
			if (d<2)
				++next[d];

			if (m_vec)
			{
				// the last z line has zero coefficients
				Cell4* cell = &m_curr4[((size_t)loc_x*n[1] + loc_y)*m_numVectors];
				for (unsigned int zv=0; zv<m_numVectors; ++zv, ++cell)
				{
					for (int t=0; t<2; ++t)
					{
						v4sf diff = volt(m_diff[t], next[0], next[1], zv).v - volt(m_diff[t], pos[0], pos[1], zv).v;
						cell->psi[t].v = cell->b.v * cell->psi[t].v + cell->ca[t].v * diff;
						curr(m_comp[t], pos[0], pos[1], zv).v += cell->psi[t].v + cell->ck[t].v * diff;
					}
				}
				continue;
			}

			// z slab, differences along z
			Cell* cell = &m_curr[((size_t)loc_x*n[1] + loc_y)*n[2]];
			for (unsigned int loc_z=0; loc_z<n[2]; ++loc_z, ++cell)
			{
				unsigned int z = loc_z + start[2];
				if (z>=numLines[2]-1)
					continue;
				unsigned int zv = m_zv[z], l = m_lane[z];
				unsigned int zvn = m_zv[z+1], ln = m_lane[z+1];
				for (int t=0; t<2; ++t)
				{
					FDTD_FLOAT diff = volt(m_diff[t], pos[0], pos[1], zvn).f[ln] - volt(m_diff[t], pos[0], pos[1], zv).f[l];
					cell->psi[t] = cell->b * cell->psi[t] + cell->ca[t] * diff;
					curr(m_comp[t], pos[0], pos[1], zv).f[l] += cell->psi[t] + cell->ck[t] * diff;
				}
			}
		}
	}
}

void Engine_Ext_CPML::DoPostCurrentUpdates(int threadID)
{
	if (m_Eng && (m_Eng->GetType()==Engine::SSE))
		DoPostCurrentUpdatesSSE((Engine_sse*)m_Eng, threadID);
	else
		ENG_DISPATCH_ARGS(DoPostCurrentUpdatesImpl, threadID);
}
