#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include <unordered_set>

#include "common.h"
#include "compressedstring.h"

namespace TwoPaCo
{

	template<size_t CAPACITY>
		class BifurcationStorage
		{
		public:
			typedef CompressedString<CAPACITY> DnaString;
			BifurcationStorage(){}

			uint64_t GetDistinctVerticesCount() const
			{
				return bifurcation_.size();
			}

			void Init(std::istream & bifurcationTempRead, uint64_t verticesCount, uint64_t vertexLength, size_t threads)
			{
				DnaString buf;
				vertexLength_ = vertexLength;
				for (size_t i = 0; i < verticesCount; i++)
				{
					buf.ReadFromFile(bifurcationTempRead);
					if (!bifurcationTempRead)
					{
						throw StreamFastaParser::Exception("Can't read from a temporary file");
					}
					
					bifurcation_.insert(buf);
				}

			}

			uint64_t GetId(std::string::const_iterator pos, bool & positiveStrand) const
			{
				DnaString bitBuf;
				uint64_t ret = INVALID_VERTEX;
				bitBuf.Clear();
				bitBuf.CopyFromString(pos, vertexLength_);
				auto it = bifurcation_.find(bitBuf);
				if (it != bifurcation_.end())
				{
					positiveStrand = true;
					return reinterpret_cast<uint64_t>(&(*it));
				}

				bitBuf.Clear();
				bitBuf.CopyFromReverseString(pos, vertexLength_);
				it = bifurcation_.find(bitBuf);
				if (it != bifurcation_.end())
				{
					positiveStrand = false;
					return reinterpret_cast<uint64_t>(&(*it));
				}

				return INVALID_VERTEX;
			}
		

		private:

			class DnaStringHash
			{
			public:
				uint64_t operator()(const DnaString & value) const
				{
					return value.Hash();
				}
			};

			size_t vertexLength_;
			std::unordered_set<DnaString, DnaStringHash> bifurcation_;
		};
}

#endif