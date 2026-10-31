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

#include <cmath>
#include <cstring>
#include <iostream>
#include <map>

#include "operator_coeffs.h"
#include "tools/global.h"

using std::cout;
using std::endl;

// ----------------------------------------------------------------------------

Coeffs_Basic::Coeffs_Basic()
{
	vv = NULL;
	vi = NULL;
	iv = NULL;
	ii = NULL;
}

Coeffs_Basic::~Coeffs_Basic()
{
	Free();
}

void Coeffs_Basic::Init(const unsigned int numLines[3])
{
	Free();
	unsigned int lines[3] = {numLines[0], numLines[1], numLines[2]};
	vv = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("vv", lines);
	vi = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("vi", lines);
	iv = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("iv", lines);
	ii = new ArrayLib::ArrayNIJK<FDTD_FLOAT>("ii", lines);
}

void Coeffs_Basic::Free()
{
	delete vv;
	delete vi;
	delete iv;
	delete ii;
	vv = NULL;
	vi = NULL;
	iv = NULL;
	ii = NULL;
}

// ----------------------------------------------------------------------------

Coeffs_SSE::Coeffs_SSE()
{
	numVectors = 0;
	f4_vv = NULL;
	f4_vi = NULL;
	f4_iv = NULL;
	f4_ii = NULL;
}

Coeffs_SSE::~Coeffs_SSE()
{
	Coeffs_SSE::Free();
}

void Coeffs_SSE::Init(const unsigned int numLines[3])
{
	Coeffs_SSE::Free();

	numVectors =  ceil((double)numLines[2]/4.0);

	f4_vv = new ArrayLib::ArrayENG<f4vector>(
		"f4_vv", {numLines[0], numLines[1], numVectors}
	);
	f4_vi = new ArrayLib::ArrayENG<f4vector>(
		"f4_vi", {numLines[0], numLines[1], numVectors}
	);
	f4_iv = new ArrayLib::ArrayENG<f4vector>(
		"f4_iv", {numLines[0], numLines[1], numVectors}
	);
	f4_ii = new ArrayLib::ArrayENG<f4vector>(
		"f4_ii", {numLines[0], numLines[1], numVectors}
	);
}

void Coeffs_SSE::Free()
{
	delete f4_vv;
	delete f4_vi;
	delete f4_iv;
	delete f4_ii;
	f4_vv = NULL;
	f4_vi = NULL;
	f4_iv = NULL;
	f4_ii = NULL;
}

// ----------------------------------------------------------------------------

//! the 24 f4vector coefficients of one node, the key for finding identical ones (Coeffs_SSE_Compressed::Finalize)
class SSE_coeff
{
public:
	SSE_coeff( f4vector vv[3], f4vector vi[3], f4vector iv[3], f4vector ii[3] );
	bool operator==( const SSE_coeff& ) const;
	bool operator!=( const SSE_coeff& ) const;
	bool operator<( const SSE_coeff& ) const;
	void print( std::ostream& stream ) const;
protected:
	f4vector m_vv[3];
	f4vector m_vi[3];
	f4vector m_iv[3];
	f4vector m_ii[3];
};

Coeffs_SSE_Compressed::Coeffs_SSE_Compressed() : Coeffs_SSE()
{
	use_compression = false;
	m_numLines[0] = m_numLines[1] = m_numLines[2] = 0;
}

Coeffs_SSE_Compressed::~Coeffs_SSE_Compressed()
{
	ClearCompressed();
}

void Coeffs_SSE_Compressed::ClearCompressed()
{
	use_compression = false;
	for (int n=0; n<3; n++)
	{
		f4_vv_Compressed[n].clear();
		f4_vi_Compressed[n].clear();
		f4_iv_Compressed[n].clear();
		f4_ii_Compressed[n].clear();
	}
}

void Coeffs_SSE_Compressed::Init(const unsigned int numLines[3])
{
	ClearCompressed();
	Coeffs_SSE::Init(numLines);
	for (int n=0; n<3; ++n)
		m_numLines[n] = numLines[n];
	unsigned int lines[3] = {numLines[0], numLines[1], numLines[2]};
	index.Init("Op_index", lines);
}

void Coeffs_SSE_Compressed::Free()
{
	ClearCompressed();
	Coeffs_SSE::Free();
}

void Coeffs_SSE_Compressed::ShowStat(std::ostream& ostr) const
{
	ostr << "SSE compression enabled\t: " << (use_compression?"yes":"no") << endl;
	ostr << "Unique SSE operators\t: " << f4_vv_Compressed->size() << endl;
	ostr << "-----------------------------------" << endl;
}

void Coeffs_SSE_Compressed::Finalize()
{
	use_compression = false;

	if (g_settings.GetVerboseLevel()>0)
		cout << "Compressing the FDTD operator... this may take a while..." << endl;

	ArrayLib::ArrayNIJK<f4vector>& f4_vv_full = *f4_vv;
	ArrayLib::ArrayNIJK<f4vector>& f4_vi_full = *f4_vi;
	ArrayLib::ArrayNIJK<f4vector>& f4_iv_full = *f4_iv;
	ArrayLib::ArrayNIJK<f4vector>& f4_ii_full = *f4_ii;

	std::map<SSE_coeff,unsigned int> lookUpMap;

	unsigned int pos[3];
	for (pos[0]=0; pos[0]<m_numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<numVectors; ++pos[2])
			{
				f4vector vv[3] = { f4_vv_full(0, pos[0], pos[1], pos[2]), f4_vv_full(1, pos[0], pos[1], pos[2]), f4_vv_full(2, pos[0], pos[1], pos[2]) };
				f4vector vi[3] = { f4_vi_full(0, pos[0], pos[1], pos[2]), f4_vi_full(1, pos[0], pos[1], pos[2]), f4_vi_full(2, pos[0], pos[1], pos[2]) };
				f4vector iv[3] = { f4_iv_full(0, pos[0], pos[1], pos[2]), f4_iv_full(1, pos[0], pos[1], pos[2]), f4_iv_full(2, pos[0], pos[1], pos[2]) };
				f4vector ii[3] = { f4_ii_full(0, pos[0], pos[1], pos[2]), f4_ii_full(1, pos[0], pos[1], pos[2]), f4_ii_full(2, pos[0], pos[1], pos[2]) };
				SSE_coeff c( vv, vi, iv, ii );

				std::map<SSE_coeff,unsigned int>::iterator it;
				it = lookUpMap.find(c);
				if (it == lookUpMap.end())
				{
					// not found -> insert
					unsigned int idx = f4_vv_Compressed[0].size();
					for (int n=0; n<3; n++)
					{
						f4_vv_Compressed[n].push_back( vv[n] );
						f4_vi_Compressed[n].push_back( vi[n] );
						f4_iv_Compressed[n].push_back( iv[n] );
						f4_ii_Compressed[n].push_back( ii[n] );
					}
					lookUpMap[c] = idx;
					index(pos[0], pos[1], pos[2]) = idx;
				}
				else
				{
					// this operator is already in the list
					index(pos[0], pos[1], pos[2]) = (*it).second;
				}
			}
		}
	}

	Coeffs_SSE::Free();
	use_compression = true;
}

// ----------------------------------------------------------------------------

SSE_coeff::SSE_coeff( f4vector vv[3], f4vector vi[3], f4vector iv[3], f4vector ii[3] )
{
	for (int n=0; n<3; n++)
	{
		m_vv[n] = vv[n];
		m_vi[n] = vi[n];
		m_iv[n] = iv[n];
		m_ii[n] = ii[n];
	}
}

bool SSE_coeff::operator==( const SSE_coeff& other ) const
{
	for (int n=0; n<3; n++)
	{
		if (memcmp( &(m_vv[n]), &(other.m_vv[n]), sizeof(f4vector) ) != 0) return false;
		if (memcmp( &(m_vi[n]), &(other.m_vi[n]), sizeof(f4vector) ) != 0) return false;
		if (memcmp( &(m_iv[n]), &(other.m_iv[n]), sizeof(f4vector) ) != 0) return false;
		if (memcmp( &(m_ii[n]), &(other.m_ii[n]), sizeof(f4vector) ) != 0) return false;
	}
	return true;
}
bool SSE_coeff::operator!=( const SSE_coeff& other ) const
{
	return !(*this == other);
}
bool SSE_coeff::operator<( const SSE_coeff& other ) const
{
	for (int n=0; n<3; n++)
	{
		for (int c=0; c<4; c++)
		{
			if (m_vv[n].f[c] > other.m_vv[n].f[c]) return false;
			if (m_vv[n].f[c] < other.m_vv[n].f[c]) return true;
			if (m_vi[n].f[c] > other.m_vi[n].f[c]) return false;
			if (m_vi[n].f[c] < other.m_vi[n].f[c]) return true;
			if (m_iv[n].f[c] > other.m_iv[n].f[c]) return false;
			if (m_iv[n].f[c] < other.m_iv[n].f[c]) return true;
			if (m_ii[n].f[c] > other.m_ii[n].f[c]) return false;
			if (m_ii[n].f[c] < other.m_ii[n].f[c]) return true;
		}
	}
	return false;
}

void SSE_coeff::print( std::ostream& stream ) const
{
	stream << "SSE_coeff: (" << endl;
	for (int n=0; n<3; n++)
	{
		stream << "n=" << n << ":" << endl;
		stream << "vv=";
		for (int c=0; c<4; c++)
			stream << m_vv[n].f[c] << " ";
		stream << endl << "vi=";
		for (int c=0; c<4; c++)
			stream << m_vi[n].f[c] << " ";
		stream << endl << "iv=";
		for (int c=0; c<4; c++)
			stream << m_iv[n].f[c] << " ";
		stream << endl << "ii=";
		for (int c=0; c<4; c++)
			stream << m_ii[n].f[c] << " ";
		stream << endl;
	}
	stream << ")" << endl;
}
