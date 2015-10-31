#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include "compressedstring.h"
#include "ngramhashing/cyclichash.h"

namespace Sibelia
{
	typedef CyclicHash<uint64_t> HashFunction;
	typedef std::unique_ptr<HashFunction> HashFunctionPtr;
	
	const static size_t INVALID_VERTEX = -1;

	template<size_t CAPACITY>
		class BifurcationStorage
		{
		public:
			typedef CompressedString<CAPACITY> DnaString;
			BifurcationStorage(){}

			uint64_t GetVerticesCount() const
			{
				return bifurcationKey_.size() + selfRevCompBifurcationKey_.size();
			}

			void Init(std::istream & bifurcationTempRead, uint64_t verticesCount, uint64_t vertexLength, size_t threads)
			{
				uint64_t bitsPower = 0;
				while (verticesCount * 8 >= (uint64_t(1) << bitsPower))
				{
					++bitsPower;
				}
				
				size_t hashFunctionNumber = 3;
				bitsPower = std::max(bitsPower, size_t(28));
				bifurcationFilter_.assign(uint64_t(1) << bitsPower, false);
				hashFunction_.resize(hashFunctionNumber);
				for (HashFunctionPtr & ptr : hashFunction_)
				{
					ptr.reset(new HashFunction(vertexLength, bitsPower));
				}

				DnaString buf;
				std::string stringBuf(vertexLength, ' ');
				for (size_t i = 0; i < verticesCount; i++)
				{
					buf.ReadFromFile(bifurcationTempRead);
					if (!bifurcationTempRead)
					{
						throw StreamFastaParser::Exception("Can't read from a temporary file");
					}

					buf.ToString(stringBuf, vertexLength);
					if (DnaChar::IsSelfReverseCompliment(stringBuf.begin(), vertexLength))
					{
						selfRevCompBifurcationKey_.push_back(buf);
					}
					else
					{
						bifurcationKey_.push_back(buf);
					}

					for (HashFunctionPtr & ptr : hashFunction_)
					{
						uint64_t hf = ptr->hash(stringBuf);
						bifurcationFilter_[hf] = true;
					}
				}

				tbb::task_scheduler_init init(threads);
				tbb::parallel_sort(bifurcationKey_.begin(), bifurcationKey_.end(), DnaString::Less);
				tbb::parallel_sort(selfRevCompBifurcationKey_.begin(), selfRevCompBifurcationKey_.end(), DnaString::Less);
			}

			const std::vector<HashFunctionPtr>& GetHashFunctions() const
			{
				return hashFunction_;
			}

		private:
			std::vector<bool> bifurcationFilter_;
			std::vector<DnaString> bifurcationKey_;
			std::vector<DnaString> selfRevCompBifurcationKey_;
			std::vector<HashFunctionPtr> hashFunction_;
		};
}

#endif