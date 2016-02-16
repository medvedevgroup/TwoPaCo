#include <set>
#include <map>
#include <cassert>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "test.h"
#include "vertexenumerator.h"

namespace TwoPaCo
{
	bool IsDefinite(int ch)
	{
		return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T';
	}

	void VertexEnumeratorTest(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, size_t threads, std::ostream & log)
	{
		size_t edgeLength = vertexLength + 1;
		std::unique_ptr<TwoPaCo::VertexEnumerator> vid = CreateEnumerator(fileName, vertexLength, filterSize, 4, 1, threads, "tmp", "de_bruijn.bin");
		
		int unknownCount = CHAR_MAX; 
		typedef std::vector<int> DnaString;
		std::set<DnaString> edges;
		std::vector<DnaString> genome;
		for (const std::string & nowFileName : fileName)
		{
			bool start = true;			
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); start = true)
			{
				char ch;
				DnaString nowGenome;
				nowGenome.push_back(unknownCount++);
				while(parser.GetChar(ch))
				{
					if (IsDefinite(ch))
					{
						nowGenome.push_back(ch);
					}
					else
					{
						nowGenome.push_back(unknownCount++);
					}
				}

				nowGenome.push_back(unknownCount++);
				genome.push_back(nowGenome);
				DnaString nowGenomeReverse;
				for (DnaString::const_reverse_iterator it = nowGenome.rbegin(); it != nowGenome.rend(); ++it)
				{
					if (IsDefinite(*it))
					{
						nowGenomeReverse.push_back(DnaChar::ReverseChar(*it));
					}
					else
					{
						nowGenomeReverse.push_back(unknownCount++);
					}
				}

				genome.push_back(nowGenomeReverse);
			}
		}

		std::map<DnaString, std::set<int> > inEdge;
		std::map<DnaString, std::set<int> > outEdge;
		for (const DnaString & g : genome)
		{
			if (g.size() >= vertexLength)
			{
				for (size_t i = 0; i <= g.size() - vertexLength; i++)
				{
					DnaString vertex(g.begin() + i, g.begin() + i + vertexLength);
					if (std::count_if(vertex.begin(), vertex.end(), IsDefinite) == vertexLength)
					{
						if (i + vertexLength < g.size())
						{
							outEdge[vertex].insert(g[i + vertexLength]);
						}

						if (i > 0)
						{
							inEdge[vertex].insert(g[i - 1]);
						}
					}
				}
			}
		}

		std::set<std::string> bif;
		std::map<DnaString, std::set<int> > * edge[] = { &inEdge, &outEdge };
		for (std::map<DnaString, std::set<int> > * e : edge)
		{
			for (auto it = e->begin(); it != e->end(); ++it)
			{
				if (it->second.size() > 1)
				{					
					std::string cand(it->first.begin(), it->first.end());
					bif.insert(cand);
					bif.insert(DnaChar::ReverseCompliment(cand));
					assert(vid->GetId(cand) != INVALID_VERTEX || vid->GetId(DnaChar::ReverseCompliment(cand)) != INVALID_VERTEX);
				}
			}
		}
		
		std::cout << "TP = " << bif.size() << std::endl;
	}

	bool Runtests()
	{
	//	DnaStringTest(100000, std::cerr);
		std::stringstream ss;
		std::vector<std::string> fileName;
		
		fileName.clear();
		fileName.push_back("test.fasta");
		VertexEnumeratorTest(fileName, 4, 16, 1, ss);

		fileName.clear();
		fileName.push_back("ntest.fasta");
		for (size_t k = 14; k <= 14; k++)
		{
			VertexEnumeratorTest(fileName, k, 20, 1, ss);
		}

		fileName.clear();
		fileName.push_back("tiny.fasta");
		VertexEnumeratorTest(fileName, 25, 24, 4, ss);
		
		fileName.push_back("ntest.fasta");
		VertexEnumeratorTest(fileName, 10, 32, 4, ss);

		return true;
	}
}