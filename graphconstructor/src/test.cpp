#include <set>
#include <cassert>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "test.h"
#include "vertexenumerator.h"


namespace Sibelia
{
	//void DnaStringTest(size_t n, std::ostream & log)
	//{
	//	size_t size = 300;
	//	DnaString str0(size, size);
	//	for (size_t i = 0; i < n; i++)
	//	{
	//		size_t idx = rand() % str0.GetSize();			
	//		char newChar = DnaString::LITERAL[rand() % DnaString::LITERAL.size()];
	//		log << "Setting char, str(" << idx << ") = " << newChar << std::endl;
	//		str0.SetChar(idx, newChar);
	//		log << "Got char " << str0.GetChar(idx) << std::endl;
	//		assert(str0.GetChar(idx) == newChar);
	//	}
	//	
	//	std::string str2;
	//	DnaString str1(300, uint64_t(0));
	//	
	//	for (size_t i = 0; i < n; i++)
	//	{
	//		log << i << '\t';
 //			if (str1.GetSize() < str1.MaxSize() && rand() % 2)
	//		{
	//			char ch = DnaString::LITERAL[rand() % DnaString::LITERAL.size()];
	//			if (rand() % 2)
	//			{
	//				log << "Appending back";
	//				str1.AppendBack(ch);
	//				str2.push_back(ch);
	//			}
	//			else
	//			{
	//				log << "Appending front";
	//				str1.AppendFront(ch);
	//				str2.insert(str2.begin(), ch);
	//			}
	//		}
	//		else if (str1.GetSize() > 0)
	//		{
	//			if (rand() % 2)
	//			{
	//				log << "Popping back";
	//				str1.PopBack();
	//				str2.erase(str2.end() - 1);
	//			}
	//			else
	//			{
	//				log << "Popping front";
	//				str1.PopFront();
	//				str2.erase(str2.begin());
	//			}
	//		}

	//		log << "\tSize=" << str1.GetSize() << std::endl;				
	//		std::string str1b = str1.ToString();
	//		log << str1b << std::endl << str2 << std::endl;
	//		assert(str1b == str2);

	//		std::string str2r = DnaString::RevComp(str2);
	//		DnaString str1r = str1.RevComp();
	//		std::string str1br = str1r.ToString();
	//		assert(str1br == str2r);
	//		/*
	//		size_t prefix = rand() % (str1.GetSize() + 1);
	//		DnaString str1pr = str1.Prefix(prefix);
	//		std::string str2pr(str2.begin(), str2.begin() + prefix);
	//		log << "Checking prefix, size = " << prefix << std::endl;
	//		
	//		log << str1pr.ToString() << std::endl;
	//		log << str2pr << std::endl;
	//		assert(str1pr.ToString() == str2pr);*/
	//	}
	//}

	void VertexEnumeratorTest(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, std::ostream & log)
	{
		std::set<std::string> edges;
		size_t edgeLength = vertexLength + 1;
		std::unique_ptr<Sibelia::VertexEnumerator> vid = CreateEnumerator(fileName, vertexLength, filterSize, 4, 4, 4, 1, "graphconstructor.tmp");

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
						edges.insert(RevComp(edge));
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
					bif.insert(RevComp(vertex));
					while (true)
					{
						std::string candVertex[] = { vertex, RevComp(vertex) };
						for (const std::string cand : candVertex)
						{
							size_t inCount = 0;
							size_t outCount = 0;
							for (char ch : LITERAL)
							{
								std::string inEdge = ch + cand;
								std::string outEdge = cand + ch;
								inCount += edges.count(inEdge);
								outCount += edges.count(outEdge);
							}

							if (inCount != 1 || outCount != 1)
							{
								assert(vid->GetId(cand) != VertexEnumerator::INVALID_VERTEX || vid->GetId(RevComp(cand)) != VertexEnumerator::INVALID_VERTEX);
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
							bif.insert(RevComp(vertex));
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
			if (bif.count(v) == 0 && bif.count(RevComp(v)) == 0)
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
//		VertexEnumeratorTest(fileName, 4, 16, ss);
	
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