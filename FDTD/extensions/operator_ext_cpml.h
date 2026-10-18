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

#ifndef OPERATOR_EXT_CPML_H
#define OPERATOR_EXT_CPML_H

#include "FDTD/operator.h"
#include <ostream>
#include <string>
#include "operator_extension.h"

#include "tools/arraylib/array_nijk.h"

//! Operator extension implementing Convolutional Perfectly Matched Layer (CPML)
/*!
  CPML uses recursive convolution with auxiliary psi fields based on the formulation
  by Roden & Gedney (2000). This implementation is compatible with the EC-FDTD method.
  The priority for this extension should be high to ensure proper field absorption.
  */
class Operator_Ext_CPML : public Operator_Extension
{
	friend class Engine_Ext_CPML;
public:
	virtual ~Operator_Ext_CPML();

	//! Returns always true, Create_CPML method will take care of creating a valid pml for the cylindrical fdtd
	virtual bool IsCylinderCoordsSave(bool closedAlpha, bool R0_included) const { UNUSED(closedAlpha); UNUSED(R0_included); return true;}

	//! Returns always true if base grid, Create_CPML will create proper child pml extensions.
	virtual bool IsCylindricalMultiGridSave(bool child) const {if (child) return false; return true;}

	virtual bool IsMPISave() const {return true;}

	void SetBoundaryCondition(const int* BCs, const unsigned int size[6]);

	void SetRange(const unsigned int start[3], const unsigned int stop[3]);

	//! Set the slab normal (0,1,2 = x,y,z) and whether it is the upper boundary
	void SetDirection(int dir, bool upper) {m_dir=dir; m_upper=upper;}

	//! Set the CPML parameters
	/*!
		Define the CPML absorption parameters.
		Parameters:
			kappa_max: Maximum kappa value (typically 15.0)
			alpha_max: Maximum alpha value (typically 0.05)
			sigma_max: Maximum sigma value (-1 for auto-calculate)
			m: Grading order (typically 4.0)
			reflection_coeff: Target reflection coefficient (typically 1e-16)
	*/
	virtual void SetCPMLParameters(double kappa_max, double alpha_max, double sigma_max,
	                               double m, double reflection_coeff);

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual std::string GetExtensionName() const {return std::string("Convolutional PML Extension");}

	virtual void ShowStat(std::ostream &ostr) const;

	//! Create the CPML
	static bool Create_CPML(Operator* op, const int ui_BC[6], const unsigned int ui_size[6],
	                        double kappa_max = 15.0, double alpha_max = 0.05, double sigma_max = -1.0,
	                        double m = 4.0, double reflection_coeff = 1e-16);

protected:
	Operator_Ext_CPML(Operator* op);
	int m_BC[6];
	unsigned int m_Size[6];

	unsigned int m_StartPos[3];
	unsigned int m_numLines[3];

	// CPML parameters
	double m_kappa_max;
	double m_alpha_max;
	double m_sigma_max;
	double m_grading_order;
	double m_reflection_coeff;

	//! Calculate b, a and kappa-1 correction (1/kappa-1) at a field position
	/*!
		The coefficients only depend on the position along the slab normal.
		E-fields are evaluated on the line, H-fields half a cell towards +dir.
	*/
	void CalcCoefficients(const unsigned int pos[3], bool isH, double& b, double& a, double& k) const;

	void DeleteOp();

	// CPML slab normal direction and side
	int m_dir;
	bool m_upper;

	// Coefficients for E-field (voltage) and H-field (current) updates
	// First index: 0 = b = exp(-(sigma/kappa + alpha) * dt / eps0)
	//              1 = a = sigma/(sigma*kappa + kappa^2*alpha) * (b - 1)
	//              2 = 1/kappa - 1
	ArrayLib::ArrayNIJK<FDTD_FLOAT> c_E;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> c_H;
};

#endif // OPERATOR_EXT_CPML_H
