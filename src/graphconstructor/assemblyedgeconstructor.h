#ifndef _ASSEMBLY_EDGE_CONSTRUCTOR_H_
#define _ASSEMBLY_EDGE_CONSTRUCTOR_H_

#include "vertexenumerator.h"

namespace TwoPaCo
{
	class AssemblyEdgeConstructor
	{
	public:
		AssemblyEdgeConstructor() {}

		virtual ~AssemblyEdgeConstructor()
		{

		}

		DISALLOW_COPY_AND_ASSIGN(AssemblyEdgeConstructor);
	};

	template<size_t CAPACITY>
	class AssemblyEdgeConstructorImpl : public AssemblyEdgeConstructor
	{
	public:
		typedef CompressedString<CAPACITY> DnaString;

		AssemblyEdgeConstructorImpl(const std::vector<std::string> & inputFileName, const std::string & marksFileName, const VertexEnumerator & vertexEnumerator) :
				vertexEnumerator_(vertexEnumerator)
		{
			DnaString buffer;
			int64_t vertexLength = vertexEnumerator_.GetHashSeed().VertexLength();
			int64_t edgeLength = vertexLength + 1;

			size_t chrNumber = 0;
			ChrReader chrReader(inputFileName);
			std::unique_ptr<ConcurrentBitVector> bloomFilter = vertexEnumerator_.ReloadBloomFilter();

			std::ofstream fout("edges.fasta");

			for (std::string chr; chrReader.NextChr(chr); chrNumber++)
			{			
				std::string vertex;	
				std::string currentOutVertex; //Current junction vertex								
				//Init hash function				
				VertexRollingHash hash(vertexEnumerator.GetHashSeed(), chr.begin(), vertexEnumerator.GetHashSeed().HashFunctionsNumber());
				for (int64_t i = 0; i <= int64_t(chr.size()) - edgeLength + 1; i++)
				{
					buffer.Clear();
					buffer.CopyFromString(chr.begin() + i, vertexLength);
					//std::cout << buffer.ToString(vertexLength) << std::endl;
					vertex = chr.substr(i, vertexLength);

					//Check if the Bloom filter contains an edge
					assert(IsOutgoingEdgeInBloomFilter(hash, *bloomFilter, chr[i + edgeLength - 1]));
					if (i > 0)
					{
						assert(IsIngoingEdgeInBloomFilter(hash, *bloomFilter, chr[i - 1]));
					}
					
					//Check the if the vertex is a junction candidate
					std::string inEdges;
					std::string outEdges;
					if (vertexEnumerator_.GetEdges(vertex.begin(), hash, inEdges, outEdges))
					{
						//Found a junction candidate, check that the mark in the vector is set
						assert(vertexEnumerator_.GetBit(chrNumber, i) == true);
					}
					
					
					if (inEdges.length() + outEdges.length() != 0) // Real junction 
					{
						if (currentOutVertex.length() == 0)
						{
							// If start of read
							currentOutVertex = vertex; 
							fout << "<edge" << std::endl << vertex;
						}
						else
						{
							fout << vertex[vertexLength - 1] << std::endl;
							if (i !=  int64_t(chr.size()) - edgeLength + 1)
								fout << "<edge" << std::endl<<vertex;
							currentOutVertex = vertex;

						}

					}
					else
						fout << chr[i + vertexLength - 1];
					
					if (i !=  int64_t(chr.size()) - edgeLength + 1)
						//Update hash				
						hash.Update(chr[i], chr[i + vertexLength]); 
				}
				if (currentOutVertex != vertex) // If end of read is not junction
				{
					bool isCurrentVertexNotJunction = true; 
					while (isCurrentVertexNotJunction)
					{
						for (char ch : DnaChar::LITERAL)
							if (IsOutgoingEdgeInBloomFilter(hash, *bloomFilter, ch))
							{
								hash.Update(vertex[0], ch);
								vertex = vertex.substr(1, vertexLength);
								vertex += ch;
								fout << ch;
								break;
							}

						std::string inEdges;
						std::string outEdges;
						vertexEnumerator_.GetEdges(vertex.begin(), hash, inEdges, outEdges);
						if (inEdges.length() + outEdges.length() != 0)
						{
							isCurrentVertexNotJunction = false;
						
						}
					}
				}
				fout<<std::endl;
				currentOutVertex = "";
			}
		}
		
	private:
		const VertexEnumerator & vertexEnumerator_;
		DISALLOW_COPY_AND_ASSIGN(AssemblyEdgeConstructorImpl);
	};

	std::unique_ptr<AssemblyEdgeConstructor> CreateAssemblyEdgeConstructor(const std::vector<std::string> & inputFileName, const std::string & marksFileName, const VertexEnumerator & vertexEnumerator);
}

#endif
