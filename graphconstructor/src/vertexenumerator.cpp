#include <deque>
#include <memory>
#include <iostream>
#include <algorithm>

#include "lib/SpookyV2.h"
#include "vertexenumerator.h"

namespace Sibelia
{	
	namespace
	{
		void PutInBloomFilter(std::vector<bool> & bitVector, const std::vector<std::pair<uint64_t, SpookyHash> > & hash, const DnaString & item)
		{
			for (const std::pair<uint64_t, SpookyHash> & h : hash)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = h.second.Hash64(&body, sizeof(body), h.first);
				bitVector[hvalue % bitVector.size()] = true;
			}
		}

		bool IsInBloomFilter(const std::vector<bool> & bitVector, const std::vector<std::pair<uint64_t, SpookyHash> > & hash, const DnaString & item)
		{
			for (const std::pair<uint64_t, SpookyHash> & h : hash)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = h.second.Hash64(&body, sizeof(body), h.first);
				if (!bitVector[hvalue % bitVector.size()])
				{
					return false;
				}
			}

			return true;
		}
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize)
	{
		size_t q = 3;
		std::vector<std::pair<uint64_t, SpookyHash> > hash(q);
		for (std::pair<uint64_t, SpookyHash> & h : hash)
		{
			h.first = rand();
		}

		size_t edgeLength = vertexLength;
		std::vector<bool> bitVector(filterSize, false);
		for (const std::string & nowFileName: fileName)
		{
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
			{
				char ch;
				DnaString edge;
				for (size_t j = 0; j < edgeLength && parser.GetChar(ch); j++)
				{
					edge.AppendBack(ch);
				}

				if (edge.GetSize() == edgeLength)
				{
					while (true)
					{
						PutInBloomFilter(bitVector, hash, edge);
						if (parser.GetChar(ch))
						{
							edge.AppendBack(ch);
							edge.PopFront();
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
				DnaString vertex;
				for (size_t j = 0; j < vertexLength && parser.GetChar(ch); j++)
				{
					vertex.AppendBack(ch);
				}

				if (vertex.GetSize() >= vertexLength)
				{
					while (true)
					{
						if (start)
						{
							bifurcation.push_back(vertex.GetBody());
						}
						else
						{
							size_t inCount = 0;
							size_t outCount = 0;
							for (char ch : DnaString::LITERAL)
							{
								DnaString inEdge = vertex;
								DnaString outEdge = vertex;
								inEdge.AppendFront(ch);
								outEdge.AppendBack(ch);
								inCount += IsInBloomFilter(bitVector, hash, inEdge) ? 1 : 0;
								outCount += IsInBloomFilter(bitVector, hash, outEdge) ? 1 : 0;
							}

							if (inCount > 1 || outCount > 1)
							{
								bifurcation.push_back(vertex.GetBody());
							}
						}

						if (parser.GetChar(ch))
						{
							vertex.AppendBack(ch);
							vertex.PopFront();
							start = false;
						}
						else
						{
							bifurcation.push_back(vertex.GetBody());
							break;
						}						
					}					
				}
			}
		}

		std::sort(bifurcation.begin(), bifurcation.end());
		bifurcation.erase(std::unique(bifurcation.begin(), bifurcation.end()), bifurcation.end());
	}
}