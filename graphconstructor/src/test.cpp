#include <set>
#include <cassert>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "test.h"
#include "vertexenumerator.h"


namespace Sibelia
{
	void VertexEnumeratorTest(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, std::ostream & log)
	{
		std::set<std::string> edges;
		size_t edgeLength = vertexLength + 1;
		std::unique_ptr<Sibelia::VertexEnumerator> vid = CreateEnumerator(fileName, vertexLength, filterSize, 4, 4, 4, 1, "graphconstructor.tmp", "de_bruijn.bin");

		for (const std::string & nowFileName : fileName)
		{
			bool start = true;			
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); start = true)
			{
				char ch;
				std::string edge;
				for (size_t j = 0; j < edgeLength && parser.GetChar(ch); j++)
				{
					edge.push_back(ch);
				}

				if (edge.size() >= edgeLength)
				{
					while (true)
					{
						edges.insert(edge);
						edges.insert(DnaChar::ReverseCompliment(edge));
						if (parser.GetChar(ch))
						{
							edge.push_back(ch);
							edge.erase(edge.begin());
						}
						else
						{
							break;
						}
					}
				}
			}
		}

		std::set<std::string> bif;
		for (const std::string & nowFileName : fileName)
		{
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
			{
				char ch;
				std::string vertex;
				for (size_t j = 0; j < vertexLength && parser.GetChar(ch); j++)
				{
					vertex.push_back(ch);
				}

				if (vertex.size() >= vertexLength)
				{
					bif.insert(vertex);
					bif.insert(DnaChar::ReverseCompliment(vertex));
					while (true)
					{
						std::string candVertex[] = { vertex, DnaChar::ReverseCompliment(vertex) };
						for (const std::string cand : candVertex)
						{
							size_t inCount = 0;
							size_t outCount = 0;
							for (char ch : DnaChar::LITERAL)
							{
								std::string inEdge = ch + cand;
								std::string outEdge = cand + ch;
								inCount += edges.count(inEdge);
								outCount += edges.count(outEdge);
							}

							if (inCount != 1 || outCount != 1)
							{
								assert(vid->GetId(cand) != VertexEnumerator::INVALID_VERTEX || vid->GetId(DnaChar::ReverseCompliment(cand)) != VertexEnumerator::INVALID_VERTEX);
								bif.insert(cand);
							}
						}						

						if (parser.GetChar(ch))
						{
							vertex.push_back(ch);
							vertex.erase(vertex.begin());
						}
						else
						{
							bif.insert(vertex);
							bif.insert(DnaChar::ReverseCompliment(vertex));
							break;
						}
					}
				}
			}
		}

		std::cout << "TP = " << bif.size() << std::endl;
		std::vector<std::string> vidSet;
		vid->Dump(vidSet);
		std::cout << "Diff:" << std::endl;
		for (const std::string & v : vidSet)
		{
			if (bif.count(v) == 0 && bif.count(DnaChar::ReverseCompliment(v)) == 0)
			{
				std::cout << v << std::endl;
			}
		}
	}

	bool Runtests()
	{
	//	DnaStringTest(100000, std::cerr);
		std::stringstream ss;
		std::vector<std::string> fileName;
		fileName.push_back("teste.fasta");
//		VertexEnumeratorTest(fileName, 3, 10, ss);

		fileName.clear();
		fileName.push_back("test.fasta");
		VertexEnumeratorTest(fileName, 4, 16, ss);
	
		fileName.clear();
		fileName.push_back("g00.fasta");
	//	fileName.push_back("g1.fasta");
	//	fileName.push_back("g2.fasta");
	//	fileName.push_back("g3.fasta");
		for (size_t k = 5; k <= 28; k++)
		{
			VertexEnumeratorTest(fileName, k, 20, ss);
		}

		fileName.clear();
		fileName.push_back("tiny.fasta");
		VertexEnumeratorTest(fileName, 25, 24, ss);
		
		return true;
	}
}