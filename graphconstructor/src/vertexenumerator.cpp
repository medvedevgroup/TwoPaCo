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
		std::vector<std::pair<uint64_t, SpookyHash> > hash(q);
		for (std::pair<uint64_t, SpookyHash> & h : hash)
		{
			h.first = rand();
		}

		std::vector<bool> bitVector(filterSize, false);
		for (const std::string & nowFileName: fileName)
		{
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
			{
				char ch;
				DnaString kmer;
				for (size_t j = 0; j < k && parser.GetChar(ch); j++)
				{
					kmer.AppendBack(ch);
				}

				if (kmer.GetSize() >= k)
				{
					while (true)
					{
						for (std::pair<uint64_t, SpookyHash> & h : hash)
						{
							uint64_t body = kmer.GetBody();
							uint64_t hvalue = h.second.Hash64(&body, sizeof(body), h.first);
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

		std::vector<uint64_t> bifurcation;
		for (const std::string & nowFileName : fileName)
		{
			bool start = true;
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); start = true)
			{
				char ch;
				DnaString kmer;
				for (size_t j = 0; j < k && parser.GetChar(ch); j++)
				{
					kmer.AppendBack(ch);
				}

				if (kmer.GetSize() >= k)
				{
					while (true)
					{
						if (start)
						{
							bifurcation.push_back(kmer.GetBody());
						}
						else
						{
							for (char ch : DnaString::LITERAL)
							{

							}
						}

						if (parser.GetChar(ch))
						{
							kmer.AppendBack(ch);
							kmer.PopFront();
							start = false;
						}
						else
						{
							bifurcation.push_back(kmer.GetBody());
							break;
						}						
					}					
				}
			}
		}
	}
}