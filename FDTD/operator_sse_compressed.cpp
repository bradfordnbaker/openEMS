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


#include "operator_sse_compressed.h"
#include "engine_sse_compressed.h"
#include "engine_sse.h"

using std::cout;
using std::endl;

Operator_SSE_Compressed* Operator_SSE_Compressed::New()
{
	cout << "Create FDTD operator (compressed SSE)" << endl;
	Operator_SSE_Compressed* op = new Operator_SSE_Compressed();
	op->Init();
	return op;
}

Operator_SSE_Compressed::Operator_SSE_Compressed() : Operator_sse()
{
	SetCoeffs(new Coeffs_SSE_Compressed());
}

Operator_SSE_Compressed::~Operator_SSE_Compressed()
{
}

Engine* Operator_SSE_Compressed::CreateEngine()
{
	if (!GetCoeffsAs<Coeffs_SSE_Compressed>().IsCompressed())
		//! create a default sse-engine
		m_Engine = Engine_sse::New(this);
	else
		m_Engine = Engine_SSE_Compressed::New(this);
	return m_Engine;
}
