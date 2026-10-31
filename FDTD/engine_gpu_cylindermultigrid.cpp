/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#include <stdexcept>

#include "engine_gpu_cylindermultigrid.h"
#include "operator_cylindermultigrid.h"
#include "engine_subgrid.h"

using std::cout;
using std::cerr;
using std::endl;

Engine_GPU_CylinderMultiGrid* Engine_GPU_CylinderMultiGrid::New(const Operator_CylinderMultiGrid* op, GPU_Backend* backend, GPU_Backend* parent)
{
	cout << "Create FDTD engine (GPU, cylindrical multi-grid level " << op->GetMultiGridLevel() << ")" << endl;
	Engine_GPU_CylinderMultiGrid* e = new Engine_GPU_CylinderMultiGrid(op, backend, parent);
	try
	{
		e->Init();
	}
	catch (...)
	{
		delete e;
		throw;
	}
	return e;
}

Engine_GPU_CylinderMultiGrid::Engine_GPU_CylinderMultiGrid(const Operator_CylinderMultiGrid* op, GPU_Backend* backend, GPU_Backend* parent)
	: Engine_GPU(op, backend, parent)
{
	Op_CMG = op;
	m_InnerEngine = NULL;
	m_Link = NULL;
}

Engine_GPU_CylinderMultiGrid::~Engine_GPU_CylinderMultiGrid()
{
	Reset();
}

void Engine_GPU_CylinderMultiGrid::Init()
{
	// this grid: backend, extensions, device or host fallback
	Engine_GPU::Init();

	Op_CMG->GetInterpolation(m_Interpol);

	// the sub-grid engine, created by its operator on a backend sharing our work stream
	Operator_Cylinder* inner_op = Op_CMG->GetInnerOperator();
	inner_op->SetGPUBackend(NULL, m_Backend);
	m_InnerEngine = dynamic_cast<Engine_GPU*>(inner_op->CreateEngine());
	if (m_InnerEngine==NULL)
		throw std::runtime_error("Engine_GPU_CylinderMultiGrid::Init: the sub-grid engine is not a GPU engine");

	// all levels on the device with a device coupling, or all in the host fallback
	if (!m_FieldsOnHost && !m_InnerEngine->FieldsOnHost())
		m_Link = m_Backend->CreateMultiGridLink(m_InnerEngine->GetBackend(), m_Interpol);
	if (m_Link==NULL)
	{
		if (!m_FieldsOnHost || !m_InnerEngine->FieldsOnHost())
			cerr << "Engine_GPU_CylinderMultiGrid::Init: Warning: multi-grid level " << Op_CMG->GetMultiGridLevel()
			     << " and its sub-grid run in the host fallback" << endl;
		ForceHostFallback();
	}
}

void Engine_GPU_CylinderMultiGrid::Reset()
{
	// the coupling refers to the backends of both grids
	delete m_Link;
	m_Link = NULL;
	delete m_InnerEngine;
	m_InnerEngine = NULL;
	Engine_GPU::Reset();
}

void Engine_GPU_CylinderMultiGrid::ForceHostFallback()
{
	delete m_Link;
	m_Link = NULL;
	Engine_GPU::ForceHostFallback();
	if (m_InnerEngine)
		m_InnerEngine->ForceHostFallback();
}

void Engine_GPU_CylinderMultiGrid::VoltageHalfStep()
{
	Engine_GPU::VoltageHalfStep();
	m_InnerEngine->VoltageHalfStep();
	if (m_Link)
		m_Link->SyncVoltages();
	else
		HostSyncVoltages();
}

void Engine_GPU_CylinderMultiGrid::CurrentHalfStep()
{
	Engine_GPU::CurrentHalfStep();
	m_InnerEngine->CurrentHalfStep();
	if (m_Link)
		m_Link->SyncCurrents();
	else if (m_Interpol.split_pos>=2)
		HostInterpolate(*curr_ptr, m_InnerEngine->HostCurrents(), true, m_Interpol.split_pos-2);
}

void Engine_GPU_CylinderMultiGrid::NextTimestep()
{
	Engine_GPU::NextTimestep();
	m_InnerEngine->NextTimestep();
}

void Engine_GPU_CylinderMultiGrid::FinishBatch()
{
	// the sub-grid first interpolates from its own sub-grids (see Engine_ThreadPool::IterateTS())
	m_InnerEngine->FinishBatch();
	if (m_Link)
		m_Link->InterpolateToBase();
	else
	{
		for (unsigned int r=0; r+1<m_Interpol.split_pos; ++r)
			HostInterpolate(*volt_ptr, m_InnerEngine->HostVoltages(), false, r);
		for (unsigned int r=0; r+2<m_Interpol.split_pos; ++r)
			HostInterpolate(*curr_ptr, m_InnerEngine->HostCurrents(), true, r);
	}
	Engine_GPU::FinishBatch();
}

// the coupling on the host mirrors, see SubGrid_Coupling_Basic
void Engine_GPU_CylinderMultiGrid::HostSyncVoltages()
{
	SubGrid_Coupling_Basic::SyncVoltages(*volt_ptr, m_InnerEngine->HostVoltages(), m_Interpol, numLines);
}

void Engine_GPU_CylinderMultiGrid::HostInterpolate(ArrayLib::ArrayNIJK<FDTD_FLOAT>& base, const ArrayLib::ArrayNIJK<FDTD_FLOAT>& sub, bool currents, unsigned int r)
{
	SubGrid_Coupling_Basic::Interpolate(base, sub, m_Interpol, currents, r, numLines);
}
