#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include "common.h"
#include "compressedstring.h"
#include "ngramhashing/cyclichash.h"

namespace Sibelia
{
	typedef CyclicHash<uint64_t> HashFunction;
	typedef std::unique_ptr<HashFunction> HashFunctionPtr;

	template<size_t CAPACITY>
	class BifurcationStorage
	{
	public:
		typedef CompressedString<CAPACITY> DnaString;
		BifurcationStorage(){}

		uint64_t GetDistinctVerticesCount() const
		{
			return bifurcationKey_.size();
		}

		uint64_t GetTotalVerticesCount() const
		{
			return bifurcationKey_.size() * 2;
		}

		void Init(std::istream & bifurcationTempRead, uint64_t verticesCount, uint64_t vertexLength, size_t threads)
		{
			uint64_t bitsPower = 0;
			vertexLength_ = vertexLength;
			while (verticesCount * 8 >= (uint64_t(1) << bitsPower))
			{
				++bitsPower;
			}

			DnaString buf;
			for (size_t i = 0; i < verticesCount; i++)
			{
				buf.ReadFromFile(bifurcationTempRead);
				if (!bifurcationTempRead)
				{
					throw StreamFastaParser::Exception("Can't read from a temporary file");
				}

				bifurcationKey_.push_back(buf);
			}

			tbb::task_scheduler_init init(threads);
			tbb::parallel_sort(bifurcationKey_.begin(), bifurcationKey_.end(), DnaString::Less);
		}

		uint64_t GetId(std::string::const_iterator pos) const
		{
			DnaString bitBuf;
			bool posFound = false;
			bool negFound = false;
			uint64_t ret = INVALID_VERTEX;
			bitBuf.CopyFromString(pos, vertexLength_);
			auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), bitBuf, DnaString::Less);
			if (it != bifurcationKey_.end() && *it == bitBuf)
			{
				posFound = true;
				ret = it - bifurcationKey_.begin();
			}

			if (!posFound)
			{
				bitBuf.Clear();
				bitBuf.CopyFromReverseString(pos, vertexLength_);
				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), bitBuf, DnaString::Less);
				if (it != bifurcationKey_.end() && *it == bitBuf)
				{
					negFound = true;
					ret = it - bifurcationKey_.begin();
				}
			}

			if (negFound && !posFound && !DnaChar::IsSelfReverseCompliment(pos, vertexLength_))
			{
				ret += bifurcationKey_.size();
			}

			return ret;
		}

	private:
		size_t vertexLength_;
		std::vector<DnaString> bifurcationKey_;

	};
}

#endif