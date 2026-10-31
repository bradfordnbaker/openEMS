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

#ifndef OPERATOR_SSE_COMPRESSED_H
#define OPERATOR_SSE_COMPRESSED_H

#include "operator_sse.h"

//! The operator for Engine_SSE_Compressed: its coefficients are a Coeffs_SSE_Compressed
class Operator_SSE_Compressed : public Operator_sse
{
public:
	static Operator_SSE_Compressed* New();
	virtual ~Operator_SSE_Compressed();

	virtual Engine* CreateEngine();

protected:
	Operator_SSE_Compressed();
};

#endif // OPERATOR_SSE_Compressed_H
