#include <set>
#include <cassert>
#include <sstream>
#include <iostream>

#include "test.h"

namespace Sibelia
{
	void DnaStringTest(size_t n, std::ostream & log)
	{
		DnaString str1;
		std::string str2;
		for (size_t i = 0; i < n; i++)
		{
			log << i << '\t';
			if (str1.GetSize() < 32 && rand() % 2)
			{
				char ch = DnaString::LITERAL[rand() % DnaString::LITERAL.size()];
				if (rand() % 2)
				{
					log << "Appending back";
					str1.AppendBack(ch);
					str2.push_back(ch);
				}
				else
				{
					log << "Appending front";
					str1.AppendFront(ch);
					str2.insert(str2.begin(), ch);
				}
			}
			else if (str1.GetSize() > 0)
			{
				if (rand() % 2)
				{
					log << "Popping back";
					str1.PopBack();
					str2.pop_back();
				}
				else
				{
					log << "Popping front";
					str1.PopFront();
					str2.erase(str2.begin());
				}
			}

			log << "\tSize=" << str1.GetSize() << std::endl;				
			std::string str1b = str1.ToString();
			log << str1b << std::endl << str2 << std::endl;
			assert(str1b == str2);
		}
	}

	void VertexEnumeratorTest(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, std::ostream & log)
	{
		std::set<std::string> vertices;
		VertexEnumerator vid(fileName, vertexLength, filterSize);
		for (const std::string & nowFileName : fileName)
		{
			bool start = true;
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); start = true)
			{
				char ch;
				std::string vertex;
				for (size_t j = 0; j < vertexLength && parser.GetChar(ch); j++)
				{
					vertex.push_back(ch);
				}

				if (vertex.size() >= vertexLength)
				{
					while (true)
					{
						vertices.insert(vertex);
						if (parser.GetChar(ch))
						{
							vertex.push_back(ch);
							vertex.erase(vertex.begin());
						}
						else
						{
							break;
						}
					}
				}
			}
		}

		for (const std::string & vertex : vertices)
		{
			size_t inCount = 0;
			size_t outCount = 0;
			for (char ch : DnaString::LITERAL)
			{
				std::string inVertex = ch + std::string(vertex.begin(), vertex.end() - 1);
				std::string outVertex = std::string(vertex.begin() + 1, vertex.end()) + ch;				
				inCount += vertices.count(inVertex);
				outCount += vertices.count(outVertex);
			}

			if (inCount != 1 || outCount != 1)
			{
				DnaString check;
				for (size_t i = 0; i < vertex.size(); i++)
				{
					check.AppendBack(vertex[i]);
				}

				assert(vid.GetId(check) != VertexEnumerator::INVALID_VERTEX);
			}
		}

	}

	bool Runtests()
	{
		std::stringstream ss;
		std::vector<std::string> fileName;
		fileName.push_back("test.fasta");
		DnaStringTest(10000, ss);
		VertexEnumeratorTest(fileName, 21, (1 << 28) + 1, ss);
		return true;
	}
}