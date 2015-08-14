#include <set>
#include <deque>
#include <ctime>
#include <memory>
#include <bitset>
#include <numeric>
#include <cassert>
#include <fstream>
#include <iostream>

#include "vertexenumerator.h"

namespace Sibelia
{
	namespace
	{
		template<size_t CAPACITY>
		std::unique_ptr<VertexEnumerator> CreateEnumeratorImpl(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			size_t aggregationThreads,
			const std::string & tmpFileName,
			const std::string & outFileName)
		{
			size_t neededCapacity = CalculateNeededCapacity(vertexLength);
			if (CAPACITY == neededCapacity)
			{
				return std::unique_ptr<VertexEnumerator>(new VertexEnumeratorImpl<CAPACITY>(fileName,
					vertexLength,
					filterSize,
					hashFunctions,
					rounds,
					threads,
					aggregationThreads,
					tmpFileName,
					outFileName));
			}
			
			return CreateEnumeratorImpl<CAPACITY + 1>(fileName,
				vertexLength,
				filterSize,
				hashFunctions,
				rounds,
				threads,
				aggregationThreads,
				tmpFileName,
				outFileName);
		}

		template<>
		std::unique_ptr<VertexEnumerator> CreateEnumeratorImpl<MAX_CAPACITY>(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			size_t aggregationThreads,
			const std::string & tmpFileName,
			const std::string & outFileName)
		{
			return 0;
		}
	}

	std::unique_ptr<VertexEnumerator> CreateEnumerator(const std::vector<std::string> & fileName,
		size_t vertexLength,
		size_t filterSize,
		size_t hashFunctions,
		size_t rounds,
		size_t threads,
		size_t aggregationThreads,
		const std::string & tmpFileName,
		const std::string & outFileName)
	{
		return CreateEnumeratorImpl<1>(fileName,
			vertexLength,
			filterSize,
			hashFunctions,
			rounds,
			threads,
			aggregationThreads,
			tmpFileName,
			outFileName);
	}
}