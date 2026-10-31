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


#include "operator_multithread.h"
#include "engine_multithread.h"

using std::cout;
using std::endl;

Operator_Multithread* Operator_Multithread::New(unsigned int numThreads)
{
	cout << "Create FDTD operator (compressed SSE + multi-threading)" << endl;
	Operator_Multithread* op = new Operator_Multithread();
	op->SetNumThreads(numThreads);
	op->Init();
	return op;
}

Operator_Multithread::~Operator_Multithread()
{
}

Engine* Operator_Multithread::CreateEngine()
{
	m_Engine = Engine_Multithread::New(this, m_numThreads);
	return m_Engine;
}

Operator_Multithread::Operator_Multithread() : Operator_SSE_Compressed()
{
}
