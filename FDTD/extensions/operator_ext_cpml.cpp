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

#include <algorithm>
#include "operator_ext_cpml.h"
#include "FDTD/operator_cylindermultigrid.h"
#include "engine_ext_cpml.h"

using namespace std;

Operator_Ext_CPML::Operator_Ext_CPML(Operator* op) : Operator_Extension(op)
{
	// Default CPML parameters, see Create_CPML()
	m_kappa_max = 1.0;
	m_alpha_max = 0.05;
	m_sigma_max = -1.0;  // Auto-calculate
	m_grading_order = 4.0;
	m_reflection_coeff = 0;  // unset, sigma_max = sigma_opt

	m_dir = 0;
	m_upper = false;

	for (int n=0; n<6; ++n)
	{
		m_BC[n]=0;
		m_Size[n]=0;
	}
	for (int n=0; n<3; ++n)
	{
		m_StartPos[n]=0;
		m_numLines[n]=0;
	}
}

Operator_Ext_CPML::~Operator_Ext_CPML()
{
	DeleteOp();
}

void Operator_Ext_CPML::SetBoundaryCondition(const int* BCs, const unsigned int size[6])
{
	for (int n=0; n<6; ++n)
	{
		m_BC[n]=BCs[n];
		m_Size[n]=size[n];
	}
}

void Operator_Ext_CPML::SetRange(const unsigned int start[3], const unsigned int stop[3])
{
	for (int n=0; n<3; ++n)
	{
		m_StartPos[n]=start[n];
		m_numLines[n]=stop[n]-start[n]+1;
	}
}

void Operator_Ext_CPML::SetCPMLParameters(double kappa_max, double alpha_max, double sigma_max,
                                           double m, double reflection_coeff)
{
	m_kappa_max = kappa_max;
	m_alpha_max = alpha_max;
	m_sigma_max = sigma_max;
	m_grading_order = m;
	m_reflection_coeff = reflection_coeff;
}

void Operator_Ext_CPML::DeleteOp()
{
	// Arrays will be automatically deleted by ArrayNIJK destructor
}

Engine_Extension* Operator_Ext_CPML::CreateEngineExtention()
{
	Engine_Ext_CPML* eng_ext = new Engine_Ext_CPML(this);
	return eng_ext;
}

void Operator_Ext_CPML::ShowStat(ostream &ostr) const
{
	Operator_Extension::ShowStat(ostr);
	ostr << " CPML parameters:" << endl;
	ostr << "  kappa_max:   " << m_kappa_max << endl;
	ostr << "  alpha_max:   " << m_alpha_max << endl;
	ostr << "  sigma_max:   " << m_sigma_max << " (auto-calc if < 0)" << endl;
	ostr << "  grading order: " << m_grading_order << endl;
	ostr << "  reflection coeff: " << m_reflection_coeff << endl;
	ostr << "  Size (x,y,z): [" << m_numLines[0] << "," << m_numLines[1] << "," << m_numLines[2] << "]" << endl;
}

void Operator_Ext_CPML::CalcCoefficients(const unsigned int pos[3], bool isH, double& b, double& a, double& k) const
{
	b = 1.0;
	a = 0.0;
	k = 0.0;

	const int n = m_dir;
	const unsigned int N = m_Op->GetNumberOfLines(n,true);
	const unsigned int s = m_upper ? N - 1 - m_Size[2*n+1] : m_Size[2*n];

	// position along the normal, H is half a cell towards +dir
	double x = m_Op->GetDiscLine(n, pos[n]);
	if (isH)
		x = 0.5 * (x + m_Op->GetDiscLine(n, std::min(pos[n]+1, N-1)));

	double width = m_upper ? m_Op->GetDiscLine(n, N-1) - m_Op->GetDiscLine(n, s)
	                       : m_Op->GetDiscLine(n, s) - m_Op->GetDiscLine(n, 0);
	double depth = m_upper ? x - m_Op->GetDiscLine(n, s)
	                       : m_Op->GetDiscLine(n, s) - x;

	// to meters, alpha lines are angles
	double scale = m_Op->GetGridDelta();
	if ((n==1) && dynamic_cast<Operator_Cylinder*>(m_Op))
		scale = m_Op->GetDiscLine(0, pos[0]) * m_Op->GetGridDelta();
	width *= scale;
	depth *= scale;

	if (width <= 0)
		return;
	if (depth < 0) depth = 0;
	if (depth > width) depth = width;

	const double m = m_grading_order;
	double sigma_max = m_sigma_max;
	if ((sigma_max < 0) && (m_reflection_coeff > 0) && (m_reflection_coeff < 1))
		sigma_max = -(m + 1.0) * log(m_reflection_coeff) / (2.0 * Z0 * width);
	else if (sigma_max < 0)
	{
		// the usual optimum 0.8*(m+1)/(Z0*delta), delta the mean cell size,
		// same as R = exp(-1.6*numCells)
		const unsigned int numCells = m_Size[2*n + m_upper];
		if (numCells == 0)
			return;
		sigma_max = 0.8 * (m + 1.0) / (Z0 * width / numCells);
	}

	double kappa = 1.0 + (m_kappa_max - 1.0) * pow(depth / width, m);
	// alpha graded linearly (Gedney), it acts only where sigma > 0
	double alpha = m_alpha_max * (width - depth) / width;
	double sigma = sigma_max * pow(depth / width, m);

	// Crank-Nicolson (bilinear) form of eps0*dpsi/dt + beta*psi = -sigma/kappa^2 * dH,
	// psi+ = b*psi- + a*(dH+ + dH-). Second order in time, unlike the recursive
	// convolution of Roden & Gedney, which loses up to 20 dB at the Courant limit.
	const double dT = m_Op->GetTimestep();
	const double beta = sigma / kappa + alpha;
	b = (2.0*EPS0 - beta*dT) / (2.0*EPS0 + beta*dT);
	a = -sigma * dT / (kappa * kappa * (2.0*EPS0 + beta*dT));
	k = 1.0 / kappa - 1.0;
}

bool Operator_Ext_CPML::BuildExtension()
{
	if (m_Op == NULL)
		return false;

	c_E.Init("c_E", m_numLines);
	c_H.Init("c_H", m_numLines);

	unsigned int pos[3];
	unsigned int loc_pos[3];
	for (loc_pos[0]=0; loc_pos[0]<m_numLines[0]; ++loc_pos[0])
	{
		pos[0] = loc_pos[0] + m_StartPos[0];
		for (loc_pos[1]=0; loc_pos[1]<m_numLines[1]; ++loc_pos[1])
		{
			pos[1] = loc_pos[1] + m_StartPos[1];
			for (loc_pos[2]=0; loc_pos[2]<m_numLines[2]; ++loc_pos[2])
			{
				pos[2] = loc_pos[2] + m_StartPos[2];
				for (int isH=0; isH<2; ++isH)
				{
					double b, a, k;
					CalcCoefficients(pos, isH, b, a, k);
					ArrayLib::ArrayNIJK<FDTD_FLOAT>& c = isH ? c_H : c_E;
					c(0, loc_pos[0], loc_pos[1], loc_pos[2]) = b;
					c(1, loc_pos[0], loc_pos[1], loc_pos[2]) = a;
					c(2, loc_pos[0], loc_pos[1], loc_pos[2]) = k;
				}
			}
		}
	}

	return true;
}

bool Operator_Ext_CPML::Create_CPML(Operator* op, const int ui_BC[6], const unsigned int ui_size[6],
                                     double kappa_max, double alpha_max, double sigma_max,
                                     double m, double reflection_coeff)
{
	int BC[6]={ui_BC[0],ui_BC[1],ui_BC[2],ui_BC[3],ui_BC[4],ui_BC[5]};
	unsigned int size[6]={ui_size[0],ui_size[1],ui_size[2],ui_size[3],ui_size[4],ui_size[5]};

	//check if mesh is large enough to support the pml
	for (int n=0; n<3; ++n)
		if ( (size[2*n]*(BC[2*n]==3)+size[2*n+1]*(BC[2*n+1]==3)) >= op->GetNumberOfLines(n,true) )
		{
			cerr << "Operator_Ext_CPML::Create_CPML: Warning: Not enough lines in direction: " << n << ", resetting to PEC" << endl;
			BC[2*n]=0;
			size[2*n]=0;
			BC[2*n+1]=0;
			size[2*n+1]=0;
		}

	//check cylindrical coord compatiblility
	Operator_Cylinder* op_cyl = dynamic_cast<Operator_Cylinder*>(op);
	if (op_cyl)
	{
		if ((BC[0]==3) && (op_cyl->GetClosedAlpha() || op_cyl->GetR0Included()))
		{
			BC[0]=0;
			size[0]=0;
			cerr << "Operator_Ext_CPML::Create_CPML: Warning: A cpml in r-min direction is not possible, resetting to PEC..." << endl;
		}
		if ( (BC[2]==3) && (op_cyl->GetClosedAlpha()) )
		{
			BC[2]=0;
			size[2]=0;
			cerr << "Operator_Ext_CPML::Create_CPML: Warning: A cpml in alpha-min direction is not possible, resetting to PEC..." << endl;
		}
	}

	//check cylindrical multi-grid compatiblility
	Operator_CylinderMultiGrid* op_cyl_mg = dynamic_cast<Operator_CylinderMultiGrid*>(op);
	if (op_cyl_mg)
	{
		if ((BC[0]==3) && (op_cyl_mg->GetR0Included()))
		{
			BC[0]=0;
			size[0]=0;
			cerr << "Operator_Ext_CPML::Create_CPML: Warning: A cpml in r-min direction is not possible, resetting to PEC..." << endl;
		}
		if ( (BC[2]==3) && (op_cyl_mg->GetClosedAlpha()) )
		{
			BC[2]=0;
			size[2]=0;
			cerr << "Operator_Ext_CPML::Create_CPML: Warning: A cpml in alpha-min direction is not possible, resetting to PEC..." << endl;
		}
	}

	// One slab per boundary. A slab only corrects the derivatives along its
	// normal, so slabs span the full mesh in the other directions and overlap
	// in the corners, where the corrections of each direction add up.
	for (int n=0; n<6; ++n)
	{
		if (BC[n]!=3)
			continue;

		const int dir = n/2;
		const bool upper = n%2;
		unsigned int start[3] = {0, 0, 0};
		unsigned int stop[3];
		for (int d=0; d<3; ++d)
			stop[d] = op->GetNumberOfLines(d,true)-1;
		if (upper)
			start[dir] = stop[dir] - size[n];
		else
			stop[dir] = size[n];

		Operator_Ext_CPML* op_ext_cpml = new Operator_Ext_CPML(op);
		op_ext_cpml->SetCPMLParameters(kappa_max, alpha_max, sigma_max, m, reflection_coeff);
		op_ext_cpml->SetBoundaryCondition(BC, size);
		op_ext_cpml->SetDirection(dir, upper);
		op_ext_cpml->SetRange(start, stop);
		op->AddExtension(op_ext_cpml);
	}

	return true;
}
