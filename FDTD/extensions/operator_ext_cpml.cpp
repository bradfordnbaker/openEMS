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

#include "operator_ext_cpml.h"
#include "FDTD/operator_cylindermultigrid.h"
#include "engine_ext_cpml.h"

using namespace std;

Operator_Ext_CPML::Operator_Ext_CPML(Operator* op) : Operator_Extension(op)
{
	// Default CPML parameters (Roden & Gedney 2000)
	m_kappa_max = 15.0;
	m_alpha_max = 0.05;
	m_sigma_max = -1.0;  // Auto-calculate
	m_grading_order = 4.0;
	m_reflection_coeff = 1e-16;

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

void Operator_Ext_CPML::CalcCPMLCoefficients(int ny, unsigned int pos[3], double Zm,
                                              double& b_E_val, double& a_E_val,
                                              double& b_H_val, double& a_H_val)
{
	// Initialize to no absorption
	b_E_val = 1.0;
	a_E_val = 0.0;
	b_H_val = 1.0;
	a_H_val = 0.0;

	double width = 0;
	double depth_v = 0;  // Depth for voltage (E-field)
	double depth_i = 0;  // Depth for current (H-field)
	int pml_dir = -1;    // Which PML region

	// Check each direction to find which PML we're in
	for (int n=0; n<3; ++n)
	{
		// Lower boundary PML
		if ((pos[n] <= m_Size[2*n]) && (m_BC[2*n]==3))
		{
			width = (m_Op->GetDiscLine(n, m_Size[2*n]) - m_Op->GetDiscLine(n, 0))
			        * m_Op->GetGridDelta();
			depth_v = width - (m_Op->GetDiscLine(n, pos[n]) - m_Op->GetDiscLine(n, 0))
			          * m_Op->GetGridDelta();

			// Cylindrical correction for alpha direction
			Operator_Cylinder* op_cyl = dynamic_cast<Operator_Cylinder*>(m_Op);
			if ((op_cyl) && (n==1))
			{
				width *= op_cyl->GetDiscDelta(0, pos[0]);
				depth_v *= op_cyl->GetDiscDelta(0, pos[0]);
			}

			// Adjust for staggered grid E-field position
			if (n == ny)
				depth_v -= m_Op->GetEdgeLength(n, pos) / 2.0;

			// H-field position
			depth_i = width - (m_Op->GetDiscLine(n, pos[n]) - m_Op->GetDiscLine(n, 0))
			          * m_Op->GetGridDelta();
			if ((op_cyl) && (n==1))
				depth_i *= op_cyl->GetDiscDelta(0, pos[0]);
			if (n != ny)
				depth_i -= m_Op->GetEdgeLength(n, pos) / 2.0;

			pml_dir = 2*n;
			break;
		}
		// Upper boundary PML
		else if ((pos[n] >= m_Op->GetNumberOfLines(n,true) - 1 - m_Size[2*n+1])
		         && (m_BC[2*n+1]==3))
		{
			width = (m_Op->GetDiscLine(n, m_Op->GetNumberOfLines(n,true)-1) -
			         m_Op->GetDiscLine(n, m_Op->GetNumberOfLines(n,true) - 1 - m_Size[2*n+1]))
			        * m_Op->GetGridDelta();
			depth_v = (m_Op->GetDiscLine(n, pos[n]) -
			           m_Op->GetDiscLine(n, m_Op->GetNumberOfLines(n,true) - 1 - m_Size[2*n+1]))
			          * m_Op->GetGridDelta();

			// Cylindrical correction
			Operator_Cylinder* op_cyl = dynamic_cast<Operator_Cylinder*>(m_Op);
			if ((op_cyl) && (n==1))
			{
				width *= op_cyl->GetDiscDelta(0, pos[0]);
				depth_v *= op_cyl->GetDiscDelta(0, pos[0]);
			}

			// Adjust for staggered grid E-field position
			if (n == ny)
				depth_v -= m_Op->GetEdgeLength(n, pos) / 2.0;

			// H-field position
			depth_i = (m_Op->GetDiscLine(n, pos[n]) -
			           m_Op->GetDiscLine(n, m_Op->GetNumberOfLines(n,true) - 1 - m_Size[2*n+1]))
			          * m_Op->GetGridDelta();
			if ((op_cyl) && (n==1))
				depth_i *= op_cyl->GetDiscDelta(0, pos[0]);
			if (n != ny)
				depth_i -= m_Op->GetEdgeLength(n, pos) / 2.0;

			pml_dir = 2*n+1;
			break;
		}
	}

	// Not in PML region
	if (pml_dir < 0)
		return;

	// Make sure depth is non-negative
	if (depth_v < 0) depth_v = 0;
	if (depth_i < 0) depth_i = 0;
	if (width <= 0) return;

	double dT = m_Op->GetTimestep();

	// Calculate graded parameters for E-field position
	double kappa_v = 1.0 + (m_kappa_max - 1.0) * pow(depth_v / width, m_grading_order);
	double alpha_v = m_alpha_max * pow((width - depth_v) / width, m_grading_order);
	double sigma_v = m_sigma_max * pow(depth_v / width, m_grading_order);

	// Calculate b_E and a_E coefficients
	double damping_v = (sigma_v / kappa_v + alpha_v);

	if (fabs(damping_v) < 1e-10)
	{
		// Taylor expansion for numerical stability
		b_E_val = 1.0 - damping_v * dT / __EPS0__;
		a_E_val = -dT / __EPS0__;
	}
	else
	{
		b_E_val = exp(-damping_v * dT / __EPS0__);
		if (fabs(alpha_v) < 1e-10)
		{
			// Simplified formula when alpha = 0
			a_E_val = (1.0 - b_E_val) / kappa_v;
		}
		else
		{
			a_E_val = sigma_v / (sigma_v * kappa_v + kappa_v * kappa_v * alpha_v) * (b_E_val - 1.0);
		}
	}

	// Calculate graded parameters for H-field position
	double kappa_i = 1.0 + (m_kappa_max - 1.0) * pow(depth_i / width, m_grading_order);
	double alpha_i = m_alpha_max * pow((width - depth_i) / width, m_grading_order);
	double sigma_i = m_sigma_max * pow(depth_i / width, m_grading_order);

	// Calculate b_H and a_H coefficients
	double damping_i = (sigma_i / kappa_i + alpha_i);

	if (fabs(damping_i) < 1e-10)
	{
		// Taylor expansion for numerical stability
		b_H_val = 1.0 - damping_i * dT / __EPS0__;
		a_H_val = -dT / __EPS0__;
	}
	else
	{
		b_H_val = exp(-damping_i * dT / __EPS0__);
		if (fabs(alpha_i) < 1e-10)
		{
			// Simplified formula when alpha = 0
			a_H_val = (1.0 - b_H_val) / kappa_i;
		}
		else
		{
			a_H_val = sigma_i / (sigma_i * kappa_i + kappa_i * kappa_i * alpha_i) * (b_H_val - 1.0);
		}
	}
}

bool Operator_Ext_CPML::BuildExtension()
{
	if (m_Op == NULL)
		return false;

	// Auto-calculate sigma_max if not specified
	if (m_sigma_max < 0)
	{
		// Calculate average PML width
		double L_avg = 0;
		int count = 0;
		for (int n=0; n<3; ++n)
		{
			if (m_BC[2*n] == 3)
			{
				L_avg += (m_Op->GetDiscLine(n, m_Size[2*n]) - m_Op->GetDiscLine(n, 0))
				         * m_Op->GetGridDelta();
				count++;
			}
			if (m_BC[2*n+1] == 3)
			{
				L_avg += (m_Op->GetDiscLine(n, m_Op->GetNumberOfLines(n,true)-1) -
				          m_Op->GetDiscLine(n, m_Op->GetNumberOfLines(n,true) - 1 - m_Size[2*n+1]))
				         * m_Op->GetGridDelta();
				count++;
			}
		}
		if (count > 0)
			L_avg /= count;
		else
			L_avg = 1.0;  // Fallback

		// Standard CPML formula: sigma_max = -(m+1) * ln(R) / (2 * Z0 * L)
		m_sigma_max = -(m_grading_order + 1.0) * log(m_reflection_coeff) / (2.0 * __Z0__ * L_avg);
	}

	// Initialize coefficient arrays
	b_E.Init("b_E", m_numLines);
	a_E.Init("a_E", m_numLines);
	b_H.Init("b_H", m_numLines);
	a_H.Init("a_H", m_numLines);

	// Loop through all cells in PML region
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

				// Calculate coefficients for each field direction
				for (int n=0; n<3; ++n)
				{
					double b_E_val, a_E_val, b_H_val, a_H_val;
					CalcCPMLCoefficients(n, pos, __Z0__, b_E_val, a_E_val, b_H_val, a_H_val);

					b_E(n, loc_pos[0], loc_pos[1], loc_pos[2]) = b_E_val;
					a_E(n, loc_pos[0], loc_pos[1], loc_pos[2]) = a_E_val;
					b_H(n, loc_pos[0], loc_pos[1], loc_pos[2]) = b_H_val;
					a_H(n, loc_pos[0], loc_pos[1], loc_pos[2]) = a_H_val;
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

	// Create CPML extensions for each boundary
	Operator_Ext_CPML* op_ext_cpml = NULL;
	unsigned int start[3];
	unsigned int stop[3];

	// X-direction PML (r-direction in cylindrical)
	if (BC[0]==3)
	{
		op_ext_cpml = new Operator_Ext_CPML(op);
		op_ext_cpml->SetCPMLParameters(kappa_max, alpha_max, sigma_max, m, reflection_coeff);
		start[0] = 0;
		start[1] = 0;
		start[2] = 0;
		stop[0] = size[0];
		stop[1] = op->GetNumberOfLines(1,true)-1;
		stop[2] = op->GetNumberOfLines(2,true)-1;
		op_ext_cpml->SetBoundaryCondition(BC, size);
		op_ext_cpml->SetRange(start, stop);
		op->AddExtension(op_ext_cpml);
	}

	if (BC[1]==3)
	{
		op_ext_cpml = new Operator_Ext_CPML(op);
		op_ext_cpml->SetCPMLParameters(kappa_max, alpha_max, sigma_max, m, reflection_coeff);
		start[0] = op->GetNumberOfLines(0,true) - 1 - size[1];
		start[1] = 0;
		start[2] = 0;
		stop[0] = op->GetNumberOfLines(0,true) - 1;
		stop[1] = op->GetNumberOfLines(1,true)-1;
		stop[2] = op->GetNumberOfLines(2,true)-1;
		op_ext_cpml->SetBoundaryCondition(BC, size);
		op_ext_cpml->SetRange(start, stop);
		op->AddExtension(op_ext_cpml);
	}

	// Y-direction PML (alpha-direction in cylindrical) - skip x-direction corners
	start[0] = (size[0]+1) * (BC[0]==3);
	stop[0] = op->GetNumberOfLines(0,true) - 1 - (size[1]+1) * (BC[1]==3);

	if (BC[2]==3)
	{
		op_ext_cpml = new Operator_Ext_CPML(op);
		op_ext_cpml->SetCPMLParameters(kappa_max, alpha_max, sigma_max, m, reflection_coeff);
		start[1] = 0;
		start[2] = 0;
		stop[1] = size[2];
		stop[2] = op->GetNumberOfLines(2,true)-1;
		op_ext_cpml->SetBoundaryCondition(BC, size);
		op_ext_cpml->SetRange(start, stop);
		op->AddExtension(op_ext_cpml);
	}

	if (BC[3]==3)
	{
		op_ext_cpml = new Operator_Ext_CPML(op);
		op_ext_cpml->SetCPMLParameters(kappa_max, alpha_max, sigma_max, m, reflection_coeff);
		start[1] = op->GetNumberOfLines(1,true) - 1 - size[3];
		start[2] = 0;
		stop[1] = op->GetNumberOfLines(1,true) - 1;
		stop[2] = op->GetNumberOfLines(2,true)-1;
		op_ext_cpml->SetBoundaryCondition(BC, size);
		op_ext_cpml->SetRange(start, stop);
		op->AddExtension(op_ext_cpml);
	}

	// Z-direction PML - skip x and y-direction corners
	start[1] = (size[2]+1) * (BC[2]==3);
	stop[1] = op->GetNumberOfLines(1,true) - 1 - (size[3]+1) * (BC[3]==3);

	if (BC[4]==3)
	{
		op_ext_cpml = new Operator_Ext_CPML(op);
		op_ext_cpml->SetCPMLParameters(kappa_max, alpha_max, sigma_max, m, reflection_coeff);
		start[2] = 0;
		stop[2] = size[4];
		op_ext_cpml->SetBoundaryCondition(BC, size);
		op_ext_cpml->SetRange(start, stop);
		op->AddExtension(op_ext_cpml);
	}

	if (BC[5]==3)
	{
		op_ext_cpml = new Operator_Ext_CPML(op);
		op_ext_cpml->SetCPMLParameters(kappa_max, alpha_max, sigma_max, m, reflection_coeff);
		start[2] = op->GetNumberOfLines(2,true) - 1 - size[5];
		stop[2] = op->GetNumberOfLines(2,true) - 1;
		op_ext_cpml->SetBoundaryCondition(BC, size);
		op_ext_cpml->SetRange(start, stop);
		op->AddExtension(op_ext_cpml);
	}

	return true;
}
