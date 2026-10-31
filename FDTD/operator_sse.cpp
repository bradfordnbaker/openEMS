/*
*	Copyright (C) 2010 Sebastian Held (Sebastian.Held@gmx.de)
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


#include "engine_sse.h"
#include "operator_sse.h"

using std::cout;
using std::endl;

Operator_sse* Operator_sse::New()
{
	cout << "Create FDTD operator (SSE)" << endl;
	Operator_sse* op = new Operator_sse();
	op->Init();
	return op;
}

Operator_sse::Operator_sse() : Operator()
{
	SetCoeffs(new Coeffs_SSE());
}

Operator_sse::~Operator_sse()
{
}

Engine* Operator_sse::CreateEngine()
{
	//! create a special sse-engine
	m_Engine = Engine_sse::New(this);
	return m_Engine;
}
