/*
*	Copyright (C) 2010-2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_THREADPOOL_H
#define ENGINE_THREADPOOL_H

#include <boost/thread.hpp>

class Engine;

//! The worker threads of a CPU engine
/*!
  Each thread updates a range of x lines (Operator::CalcStartStopLines()) with
  Engine::UpdateVoltages() and Engine::UpdateCurrents() and runs the
  thread-aware hooks of the engine extensions, with a barrier after the main
  updates and after each extension.
  */
class Engine_ThreadPool
{
	friend class Engine_ThreadPool_Worker;
public:
	//! Threads for \a eng; \a numThreads 0 starts with one and adds threads while the speed increases, see NextInterval()
	Engine_ThreadPool(Engine* eng, unsigned int numThreads);
	~Engine_ThreadPool();

	//! Run \a iterTS timesteps and return when they are done
	void IterateTS(unsigned int iterTS);

	//! The speed of the last interval, to find the best number of threads (with 0 threads requested)
	void NextInterval(float curr_speed);

	unsigned int GetNumThreads() const {return m_numThreads;}

protected:
	void ChangeNumThreads(unsigned int numThreads);
	void StopThreads();

	void DoPreVoltageUpdates(int threadID);
	void DoPostVoltageUpdates(int threadID);
	void Apply2Voltages(int threadID);
	void DoPreCurrentUpdates(int threadID);
	void DoPostCurrentUpdates(int threadID);
	void Apply2Current(int threadID);

	Engine* m_Eng;
	boost::thread_group* m_thread_group;
	boost::barrier* m_startBarrier;
	boost::barrier* m_stopBarrier;
	boost::barrier* m_IterateBarrier;
	volatile unsigned int m_iterTS;
	volatile bool m_stopThreads;
	unsigned int m_numThreads;       //!< number of worker threads
	unsigned int m_max_numThreads;   //!< max. number of worker threads
	bool m_opt_speed;
	float m_last_speed;
};

#endif // ENGINE_THREADPOOL_H
