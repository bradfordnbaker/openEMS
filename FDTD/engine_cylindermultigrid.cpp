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

#include "engine_cylindermultigrid.h"
#include "operator_cylindermultigrid.h"

using std::cout;
using std::endl;

Engine_CylinderMultiGrid* Engine_CylinderMultiGrid::New(const Operator_CylinderMultiGrid* op, unsigned int numThreads)
{
	cout << "Create FDTD engine (cylindrical multi grid mesh using sse compression + multithreading)" << endl;
	Engine_CylinderMultiGrid* e = new Engine_CylinderMultiGrid(op);
	e->SetNumThreads( numThreads );
	e->Init();
	return e;
}

Engine_CylinderMultiGrid::Engine_CylinderMultiGrid(const Operator_CylinderMultiGrid* op) : Engine_sse(op)
{
	Op_CMG = op;
	m_InnerEngine = dynamic_cast<Engine_sse*>(op->GetInnerOperator()->CreateEngine());
	if (m_InnerEngine==NULL)
		throw std::runtime_error("Engine_CylinderMultiGrid: the sub-grid engine is not an SSE engine");
}

Engine_CylinderMultiGrid::~Engine_CylinderMultiGrid()
{
	// the threads of this engine step the sub-grid too
	Reset();
	delete m_InnerEngine;
	m_InnerEngine = NULL;
}

void Engine_CylinderMultiGrid::SyncSubGridVoltages()
{
	Engine_sse* m_InnerEng = m_InnerEngine;

	//interpolate voltages from base engine to child engine...
	unsigned int pos[3];
	pos[0] = Op_CMG->GetSplitPos()-1;
	unsigned int pos1_half = 0;
	f4vector v_null;
	v_null.f[0] = 0;
	v_null.f[1] = 0;
	v_null.f[2] = 0;
	v_null.f[3] = 0;
	for (pos[1]=0; pos[1]<numLines[1]-1; pos[1]+=2)
	{
		pos1_half = pos[1]/2;
		for (pos[2]=0; pos[2]<numVectors; ++pos[2])
		{
			ArrayLib::ArrayNIJK<f4vector>& mgEng_f4_volt = *f4_volt_ptr;
			ArrayLib::ArrayNIJK<f4vector>& innerEng_f4_volt = *m_InnerEng->f4_volt_ptr;

			//r - direczion
			innerEng_f4_volt(0, pos[0], pos1_half, pos[2]).v = v_null.v;

			//z - direction
			innerEng_f4_volt(2, pos[0], pos1_half, pos[2]).v = mgEng_f4_volt(2, pos[0], pos[1], pos[2]).v;

			//alpha - direction
			innerEng_f4_volt(1, pos[0], pos1_half, pos[2]).v  = mgEng_f4_volt(1, pos[0], pos[1]  , pos[2]).v;
			innerEng_f4_volt(1, pos[0], pos1_half, pos[2]).v += mgEng_f4_volt(1, pos[0], pos[1]+1, pos[2]).v;
		}
	}
}

void Engine_CylinderMultiGrid::SyncSubGridCurrents()
{
	InterpolCurrChild2Base(Op_CMG->GetSplitPos()-2);
}

void Engine_CylinderMultiGrid::InterpolateFromSubGrid()
{
	//interpolate child data to base mesh...
	for (unsigned int n=0; n<Op_CMG->m_Split_Pos-1; ++n)
		InterpolVoltChild2Base(n);
	for (unsigned int n=0; n<Op_CMG->m_Split_Pos-2; ++n)
		InterpolCurrChild2Base(n);
}

void Engine_CylinderMultiGrid::InterpolVoltChild2Base(unsigned int rPos)
{
	//interpolate voltages from child engine to the base engine...
	ArrayLib::ArrayNIJK<f4vector>& f4_volt = *f4_volt_ptr;
	ArrayLib::ArrayNIJK<f4vector>& innerEngine_f4_volt = *m_InnerEngine->f4_volt_ptr;

	unsigned int pos[3];
	pos[0] = rPos;
	for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
	{
		for (pos[2]=0; pos[2]<numVectors; ++pos[2])
		{
			//r - direction
			f4_volt(0, pos[0], pos[1], pos[2]).v  = Op_CMG->f4_interpol_v_2p[0][pos[1]].v * innerEngine_f4_volt(0, pos[0], Op_CMG->m_interpol_pos_v_2p[0][pos[1]], pos[2]).v
													+ Op_CMG->f4_interpol_v_2pp[0][pos[1]].v * innerEngine_f4_volt(0, pos[0], Op_CMG->m_interpol_pos_v_2pp[0][pos[1]], pos[2]).v;

			//z - direction
			f4_volt(2, pos[0], pos[1], pos[2]).v  = Op_CMG->f4_interpol_v_2p[0][pos[1]].v * innerEngine_f4_volt(2, pos[0], Op_CMG->m_interpol_pos_v_2p[0][pos[1]], pos[2]).v
													+ Op_CMG->f4_interpol_v_2pp[0][pos[1]].v * innerEngine_f4_volt(2, pos[0], Op_CMG->m_interpol_pos_v_2pp[0][pos[1]], pos[2]).v;

			//alpha - direction
			f4_volt(1, pos[0], pos[1], pos[2]).v  = Op_CMG->f4_interpol_v_2p[1][pos[1]].v * innerEngine_f4_volt(1, pos[0], Op_CMG->m_interpol_pos_v_2p[1][pos[1]], pos[2]).v
													+ Op_CMG->f4_interpol_v_2pp[1][pos[1]].v * innerEngine_f4_volt(1, pos[0], Op_CMG->m_interpol_pos_v_2pp[1][pos[1]], pos[2]).v;
		}
	}
}

void Engine_CylinderMultiGrid::InterpolCurrChild2Base(unsigned int rPos)
{
	//interpolate voltages from child engine to the base engine...
	ArrayLib::ArrayNIJK<f4vector>& f4_curr = *f4_curr_ptr;
	ArrayLib::ArrayNIJK<f4vector>& innerEngine_f4_curr = *m_InnerEngine->f4_curr_ptr;

	unsigned int pos[3];
	pos[0] = rPos;
	for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
	{
		for (pos[2]=0; pos[2]<numVectors; ++pos[2])
		{
			//r - direction
			f4_curr(0, pos[0], pos[1], pos[2]).v  = Op_CMG->f4_interpol_i_2p[0][pos[1]].v * innerEngine_f4_curr(0, pos[0], Op_CMG->m_interpol_pos_i_2p[0][pos[1]], pos[2]).v
													+ Op_CMG->f4_interpol_i_2pp[0][pos[1]].v * innerEngine_f4_curr(0, pos[0], Op_CMG->m_interpol_pos_i_2pp[0][pos[1]], pos[2]).v;

			//z - direction
			f4_curr(2, pos[0], pos[1], pos[2]).v  = Op_CMG->f4_interpol_i_2p[0][pos[1]].v * innerEngine_f4_curr(2, pos[0], Op_CMG->m_interpol_pos_i_2p[0][pos[1]], pos[2]).v
													+ Op_CMG->f4_interpol_i_2pp[0][pos[1]].v * innerEngine_f4_curr(2, pos[0], Op_CMG->m_interpol_pos_i_2pp[0][pos[1]], pos[2]).v;

			//alpha - direction
			f4_curr(1, pos[0], pos[1], pos[2]).v  = Op_CMG->f4_interpol_i_2p[1][pos[1]].v * innerEngine_f4_curr(1, pos[0], Op_CMG->m_interpol_pos_i_2p[1][pos[1]], pos[2]).v
													+ Op_CMG->f4_interpol_i_2pp[1][pos[1]].v * innerEngine_f4_curr(1, pos[0], Op_CMG->m_interpol_pos_i_2pp[1][pos[1]], pos[2]).v;
		}
	}
}
