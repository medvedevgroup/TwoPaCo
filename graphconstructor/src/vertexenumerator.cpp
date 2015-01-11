#include <deque>
#include <memory>
#include <iostream>

#include "lib/SpookyV2.h"
#include "vertexenumerator.h"

namespace Sibelia
{	
	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t k, size_t filterSize)
	{
		size_t q = 3;
		std::vector<bool> bitVector(filterSize, false);
		for (const std::string & nowFileName: fileName)
		{
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
			{
				char ch;
				DnaString kmer;
				std::vector<SpookyHash> hash(q);
				for (SpookyHash & h : hash)
				{
					h.Init(rand(), rand());
				}
				
				for (size_t j = 0; j < k && parser.GetChar(ch); j++)
				{
					kmer.AppendBack(ch);
				}

				if (kmer.GetSize() >= k)
				{
					while (true)
					{
						std::cout << kmer.ToString();
						for (SpookyHash & h : hash)
						{
							uint64_t body = kmer.GetBody();
							uint64_t hvalue = h.Hash64(&body, sizeof(body), 0);
							bitVector[hvalue % bitVector.size()] = true;
							
						}

						if (parser.GetChar(ch))
						{
							kmer.AppendBack(ch);
							kmer.PopFront();
						}
						else
						{
							break;
						}
					}
				}
			}	
		}
	}
}