#ifndef _ASSEMBLY_EDGE_CONSTRUCTOR_H_
#define _ASSEMBLY_EDGE_CONSTRUCTOR_H_

#include "vertexenumerator.h"

namespace TwoPaCo
{
	class AssemblyEdgeConstructor
	{
	public:
		AssemblyEdgeConstructor(const std::vector<std::string> & inputFileName, const std::string & marksFileName, const VertexEnumerator & vertexEnumerator) : 				
				vertexEnumerator_(vertexEnumerator)
		{
			int64_t vertexLength = vertexEnumerator_.GetHashSeed().VertexLength();
			int64_t edgeLength = vertexLength + 1;

			std::vector<std::vector<bool> > junctionMark;
			size_t chrNumber = 0;
			ChrReader chrReader(inputFileName);
			std::unique_ptr<ConcurrentBitVector> bloomFilter = vertexEnumerator_.ReloadBloomFilter();
			for (std::string chr; chrReader.NextChr(chr); chrNumber++)
			{				
				//Init hash function
				junctionMark.push_back(std::vector<bool>(chr.size(), false));
				VertexRollingHash hash(vertexEnumerator.GetHashSeed(), chr.begin(), vertexEnumerator.GetHashSeed().HashFunctionsNumber());
				for (int64_t i = 0; i <= int64_t(chr.size()) - edgeLength; i++)
				{
					std::string vertex = chr.substr(i, vertexLength);
					//Check if the Bloom filter contains an edge
					assert(IsOutgoingEdgeInBloomFilter(*bloomFilter, hash, chr[i + edgeLength - 1]));
					hash.Update(chr[i], chr[i + vertexLength]);
					//Check that the hash values were updated correctly
					assert(hash.Assert(chr.begin() + i + 1));
				}
			}
		}
		
	private:
		const VertexEnumerator & vertexEnumerator_;
		DISALLOW_COPY_AND_ASSIGN(AssemblyEdgeConstructor);
	};


}

#endif