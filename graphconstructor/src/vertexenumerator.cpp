#include <deque>
#include <memory>

#include "vertexenumerator.h"
#include "ngramhashing/rabinkarphash.h"

namespace Sibelia
{
	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t k)
	{
		k = 4;
		size_t q = 3;
		for (const std::string & nowFileName: fileName)
		{
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
			{
				char ch;
				std::deque<char> roller;
				typedef KarpRabinHash HashFunction;
				std::vector<std::unique_ptr<HashFunction> > hashPtr(q);
				for (std::unique_ptr<HashFunction> & ptr : hashPtr)
				{
					ptr.reset(new HashFunction(k, 8));
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
						std::cout << std::string(roller.begin(), roller.end());
						for (std::unique_ptr<HashFunction> & ptr : hashPtr)
						{
							std::cout << ' ' << ptr->hashvalue;
						}

						std::cout << std::endl;
						if (parser.GetChar(ch))
						{
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
		}
	}
}