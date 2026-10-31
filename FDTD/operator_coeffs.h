/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef OPERATOR_COEFFS_H
#define OPERATOR_COEFFS_H

#include <ostream>
#include <vector>

#include "tools/global.h"
#include "tools/constants.h"
#include "tools/array_ops.h"
#include "tools/aligned_allocator.h"
#include "tools/arraylib/array_nijk.h"
#include "tools/arraylib/array_ijk.h"
#include "tools/arraylib/array_e.h"

//! The update coefficients of the operator: vv, vi (voltage update) and ii, iv (current update) per node and direction
/*!
  The operator builds them through the Set functions, whatever the layout; an
  engine reads the layout it was written for directly, see Operator::GetCoeffs().
  */
class Operator_Coeffs
{
public:
	virtual ~Operator_Coeffs() {}

	//! Allocate the coefficients for \a numLines mesh lines, dropping any old ones
	virtual void Init(const unsigned int numLines[3]) = 0;
	//! Free the coefficients; the Get and Set functions are invalid after this
	virtual void Free() = 0;
	//! Called once the operator and its extensions are built
	virtual void Finalize() {}
	virtual void ShowStat(std::ostream& ostr) const {(void)ostr;}

	virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const = 0;
	virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const = 0;
	virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const = 0;
	virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const = 0;

	virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) = 0;
	virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) = 0;
	virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) = 0;
	virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) = 0;
};

//! One float per coefficient, node and direction (Engine, and the host side of Engine_GPU)
class Coeffs_Basic : public Operator_Coeffs
{
public:
	Coeffs_Basic();
	virtual ~Coeffs_Basic();

	virtual void Init(const unsigned int numLines[3]);
	virtual void Free();

	inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*vv)(n, x, y, z);}
	inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*vi)(n, x, y, z);}
	inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*ii)(n, x, y, z);}
	inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*iv)(n, x, y, z);}

	inline virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*vv)(n, x, y, z) = value;}
	inline virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*vi)(n, x, y, z) = value;}
	inline virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*ii)(n, x, y, z) = value;}
	inline virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*iv)(n, x, y, z) = value;}

	ArrayLib::ArrayNIJK<FDTD_FLOAT>* vv; //!< calc new voltage from old voltage
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* vi; //!< calc new voltage from old current
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* iv; //!< calc new current from old current
	ArrayLib::ArrayNIJK<FDTD_FLOAT>* ii; //!< calc new current from old voltage
};

//! Four z-lines interleaved per f4vector, the layout Coeffs_SSE_Compressed is built in
class Coeffs_SSE : public Operator_Coeffs
{
public:
	Coeffs_SSE();
	virtual ~Coeffs_SSE();

	virtual void Init(const unsigned int numLines[3]);
	virtual void Free();

	inline virtual FDTD_FLOAT GetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*f4_vv)(n, x, y, z%numVectors).f[z/numVectors];}
	inline virtual FDTD_FLOAT GetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*f4_vi)(n, x, y, z%numVectors).f[z/numVectors];}
	inline virtual FDTD_FLOAT GetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*f4_ii)(n, x, y, z%numVectors).f[z/numVectors];}
	inline virtual FDTD_FLOAT GetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const {return (*f4_iv)(n, x, y, z%numVectors).f[z/numVectors];}

	inline virtual void SetVV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*f4_vv)(n, x, y, z%numVectors).f[z/numVectors] = value;}
	inline virtual void SetVI(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*f4_vi)(n, x, y, z%numVectors).f[z/numVectors] = value;}
	inline virtual void SetII(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*f4_ii)(n, x, y, z%numVectors).f[z/numVectors] = value;}
	inline virtual void SetIV(unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value) {(*f4_iv)(n, x, y, z%numVectors).f[z/numVectors] = value;}

	unsigned int numVectors;
	ArrayLib::ArrayENG<f4vector>* f4_vv; //!< calc new voltage from old voltage
	ArrayLib::ArrayENG<f4vector>* f4_vi; //!< calc new voltage from old current
	ArrayLib::ArrayENG<f4vector>* f4_iv; //!< calc new current from old current
	ArrayLib::ArrayENG<f4vector>* f4_ii; //!< calc new current from old voltage
};

//! Coeffs_SSE, of which Finalize() keeps only the distinct coefficient sets plus an index per node (Engine_sse)
class Coeffs_SSE_Compressed : public Coeffs_SSE
{
public:
	Coeffs_SSE_Compressed();
	virtual ~Coeffs_SSE_Compressed();

	virtual void Init(const unsigned int numLines[3]);
	virtual void Free();
	//! Compress the coefficients and free the full arrays
	virtual void Finalize();
	virtual void ShowStat(std::ostream& ostr) const;

	bool IsCompressed() const {return use_compression;}

	inline virtual FDTD_FLOAT GetVV( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { if (use_compression) return f4_vv_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors]; else return Coeffs_SSE::GetVV(n,x,y,z);}
	inline virtual FDTD_FLOAT GetVI( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { if (use_compression) return f4_vi_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors]; else return Coeffs_SSE::GetVI(n,x,y,z);}
	inline virtual FDTD_FLOAT GetII( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { if (use_compression) return f4_ii_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors]; else return Coeffs_SSE::GetII(n,x,y,z);}
	inline virtual FDTD_FLOAT GetIV( unsigned int n, unsigned int x, unsigned int y, unsigned int z ) const { if (use_compression) return f4_iv_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors]; else return Coeffs_SSE::GetIV(n,x,y,z);}

	inline virtual void SetVV( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value ) { if (use_compression) f4_vv_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors] = value; else Coeffs_SSE::SetVV(n,x,y,z,value);}
	inline virtual void SetVI( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value ) { if (use_compression) f4_vi_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors] = value; else Coeffs_SSE::SetVI(n,x,y,z,value);}
	inline virtual void SetII( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value ) { if (use_compression) f4_ii_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors] = value; else Coeffs_SSE::SetII(n,x,y,z,value);}
	inline virtual void SetIV( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value ) { if (use_compression) f4_iv_Compressed[n][index(x, y, z%numVectors)].f[z/numVectors] = value; else Coeffs_SSE::SetIV(n,x,y,z,value);}

	bool use_compression;
	ArrayLib::ArrayIJK<unsigned int> index;
	std::vector<f4vector,aligned_allocator<f4vector> > f4_vv_Compressed[3];
	std::vector<f4vector,aligned_allocator<f4vector> > f4_vi_Compressed[3];
	std::vector<f4vector,aligned_allocator<f4vector> > f4_iv_Compressed[3];
	std::vector<f4vector,aligned_allocator<f4vector> > f4_ii_Compressed[3];

protected:
	unsigned int m_numLines[3];
	void ClearCompressed();
};

#endif // OPERATOR_COEFFS_H
