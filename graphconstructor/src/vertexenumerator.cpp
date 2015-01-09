#include <deque>
#include <memory>

#include "vertexenumerator.h"
#include "lib/SpookyV2.h"

namespace Sibelia
{
	namespace
	{
		std::string ordinary = "ACGT";

		char MakeUp(char ch)
		{
			ch = toupper(ch);
			if (std::find(ordinary.begin(), ordinary.end(), ch) != ordinary.end())
			{
				return ch;
			}

			return ordinary[rand() % ordinary.size()];
		}
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t k, size_t filterSize)
	{/*
		size_t q = 3;
		std::vector<bool> bitVector(filterSize, false);
		for (const std::string & nowFileName: fileName)
		{
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
			{
				char ch;
				std::deque<char> roller;

				typedef CyclicHash HashFunction;
				std::vector<std::unique_ptr<HashFunction> > hashPtr(q);
				for (std::unique_ptr<HashFunction> & ptr : hashPtr)
				{
					ptr.reset(new HashFunction(k, 31));
				}
				
				for (size_t j = 0; j < k && parser.GetChar(ch); j++)
				{
					for (std::unique_ptr<HashFunction> & ptr : hashPtr)
					{
						ptr->eat(ch);
					}
					
					roller.push_back(ch);
				}

				if (roller.size() >= k)
				{
					while (true)
					{
						for (std::unique_ptr<HashFunction> & ptr : hashPtr)
						{
							bitVector[ptr->hashvalue % bitVector.size()] = true;
						}

						if (parser.GetChar(ch))
						{
							ch = MakeUp(ch);
							roller.push_back(ch);
							for (std::unique_ptr<HashFunction> & ptr : hashPtr)
							{
								ptr->update(roller.front(), roller.back());
							}

							roller.pop_front();
						}
						else
						{
							break;
						}
					}
				}
			}	
		}*/
	}
}