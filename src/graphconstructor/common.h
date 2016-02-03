#ifndef _COMMON_H_
#define _COMMON_H_

#include <atomic>
#include <vector>
#include <memory>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <tbb/concurrent_queue.h>

#include "ngramhashing/cyclichash.h"

namespace Sibelia
{
	using std::max;
	using std::min;

	extern const size_t INVALID_VERTEX;
	extern const uint32_t MAX_COUNTER;

	typedef CyclicHash<uint64_t> HashFunction;
	typedef std::unique_ptr<HashFunction> HashFunctionPtr;

	struct Task
	{
		bool isFinal;
		uint32_t piece;
		uint64_t start;
		uint64_t seqId;
		std::string str;
#ifdef _DEBUG
		static const size_t TASK_SIZE = 32;
#else
		static const size_t TASK_SIZE = 1 << 19;
#endif					
		static const size_t GAME_OVER = SIZE_MAX;
		Task() {}
		Task(uint64_t seqId, uint64_t start, uint32_t piece, bool isFinal, std::string && str) :
			seqId(seqId), start(start), piece(piece), isFinal(isFinal), str(std::move(str)) {}
	};

	typedef tbb::concurrent_bounded_queue<Task> TaskQueue;
	typedef std::unique_ptr<TaskQueue> TaskQueuePtr;
}

#endif