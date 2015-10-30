#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

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
			void Init(size_t vertexNumber)
			{

			}

			static uint64_t GetBits(uint64_t verticesCount) const
			{
				uint64_t bitsPower = 0;
				while (verticesCount * 8 >= (uint64_t(1) << bitsPower))
				{
					++bitsPower;
				}

				return std::max(bitsPower, size_t(28));
			}

			void ReadBifurcations(std::istream & in, uint64_t verticesCount, const std::vector<HashFunctionPtr> & hashFunction)
			{
				bifurcationKey_.reserve(verticesCount);
				std::string stringBuf(vertexLength, ' ');
				for (size_t i = 0; i < verticesCount; i++)
				{
					buf.ReadFromFile(bifurcationTempRead);
					if (!bifurcationTempRead)
					{
						throw StreamFastaParser::Exception("Can't read from a temporary file");
					}

					bifurcationKey_.push_back(buf);
					buf.ToString(stringBuf, vertexLength);
					for (HashFunctionPtr & ptr : hashFunction)
					{
						uint64_t hf = ptr->hash(stringBuf);
						bifurcationFilter[hf] = true;
					}
				}
			}

			uint64_t GetBifurcationId(std::string::const_iterator it, const std::vector<HashFunctionPtr> & posHashFunction, const std::vector<HashFunctionPtr> & negHashFunction, ) const
			{
				return 0;
			}

		private:
			std::vector<bool> bifurcationFilter_;
			std::vector<DnaString> bifurcationKey_;
			std::vector<DnaString> selfRevCompBifurcationKey_;
		};
}

#endif