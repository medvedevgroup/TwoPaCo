#include <set>
#include <map>
#include <random>
#include <cassert>
#include <sstream>
#include <iterator>
#include <stdexcept>
#include <algorithm>

#include "test.h"
#include "vertexenumerator.h"

namespace TwoPaCo
{
	bool IsDefinite(int ch)
	{
		return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T';
	}

	void GenerateSequence(std::random_device & rd, size_t length, std::string & out)
	{
		out.resize(length);
		std::mt19937 e2(rd());
		std::string alphabet("ACGT");
		std::uniform_int_distribution<> gen(0, alphabet.size() - 1);
		for (size_t i = 0; i < out.size(); i++)
		{
			out[i] = alphabet[gen(e2)];
		}
	}

	void MutateSequence(std::random_device & rd, const std::string & chr, double changeRate, double mutationRate, std::string & out)
	{
		out.clear();
		std::mt19937 e2(rd());
		std::string alphabet("ACGT");
		std::uniform_real_distribution<> gen(0, 1);
		for (auto ch : chr)
		{
			if (gen(e2) <= changeRate)
			{
				if (gen(e2) <= mutationRate)
				{
					out.push_back(alphabet[rand() % alphabet.size()]);
				}
				else
				{
					if (gen(e2) <= 0.5)
					{
						out.push_back(ch);
						out.push_back(alphabet[rand() % alphabet.size()]);
					}
				}
			}
			else
			{
				out.push_back(ch);
			}
		}
	}

	namespace
	{
		void FindJunctionsNaively(const std::vector<std::string> & chr,
			size_t vertexLength,
			std::set<std::string> & junction,
			std::vector<std::vector<bool> > & marks,
			std::map<std::string, std::set<char> > & inEdge,
			std::map<std::string, std::set<char> > & outEdge)
		{
			inEdge.clear();
			outEdge.clear();
			typedef std::string DnaString;
			std::set<DnaString> edges;
			std::vector<DnaString> genome;
			for (const std::string & nowChr : chr)
			{
				DnaString nowGenome;
				for (auto ch : nowChr)
				{
					nowGenome.push_back(ch);
				}

				genome.push_back(nowGenome);
				DnaString nowGenomeReverse;
				for (DnaString::const_reverse_iterator it = nowGenome.rbegin(); it != nowGenome.rend(); ++it)
				{
					nowGenomeReverse.push_back(DnaChar::ReverseChar(*it));
				}

				genome.push_back(nowGenomeReverse);
			}
			
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

			std::map<DnaString, std::set<char> > * edge[] = { &inEdge, &outEdge };
			for (std::map<DnaString, std::set<char> > * e : edge)
			{
				for (auto it = e->begin(); it != e->end(); ++it)
				{
					if (it->second.size() > 1)
					{
						std::string cand(it->first.begin(), it->first.end());
						junction.insert(cand);
						junction.insert(DnaChar::ReverseCompliment(cand));					
					}
				}
			}

			for (size_t i = 0; i < chr.size(); i++)
			{
				for (size_t pos = 0; pos < chr[i].size(); pos++)
				{
					if (pos == 0 || pos == chr[i].size() - vertexLength || junction.count(chr[i].substr(pos, vertexLength)) > 0)
					{
						marks[i][pos] = true;
					}
				}
			}
		}
	}

	bool RunTests(size_t tests, size_t filterBits, size_t length, size_t chrNumber, Range vertexSize, Range hashFunctions, Range rounds, Range threads, double changeRate, double indelRate, const std::string & temporaryDir)
	{
		const std::string temporaryFasta = temporaryDir + "/test.fa";
		const std::string temporaryEdge = temporaryDir + "/out.bin";
		std::vector<std::string> fileName;
		fileName.push_back(temporaryFasta);
		std::random_device rd;		
		for (size_t t = 0; t < tests; t++)
		{
			std::vector<std::string> chr(chrNumber);
			GenerateSequence(rd, length, chr[0]);
			for (size_t i = 1; i < chrNumber; i++)
			{
				MutateSequence(rd, chr[0], changeRate, indelRate, chr[i]);
			}

			std::ofstream test(temporaryFasta.c_str());
			if (!test)
			{
				throw std::runtime_error("Can't create a temporary file for testing");
			}

			for (size_t j = 0; j < chrNumber; ++j)
			{
				test << ">" << j << std::endl;
				test << chr[j] << std::endl;
			}
			
			std::map<std::string, std::set<char> > inEdge;
			std::map<std::string, std::set<char> > outEdge;
			for (size_t k = vertexSize.first; k < vertexSize.second; k += 2)
			{
				std::set<std::string> junctions;				
				std::vector<std::vector<bool > > naiveMarks(chrNumber);				
				for (size_t i = 0; i < chrNumber; i++)
				{
					naiveMarks[i].assign(chr[i].size(), false);
				}

				std::vector<std::vector<bool > > fastMarks(chrNumber);
				FindJunctionsNaively(chr, k, junctions, naiveMarks, inEdge, outEdge);
				for (size_t hf = hashFunctions.first; hf < hashFunctions.second; ++hf)
				{
					for (size_t r = rounds.first; r < rounds.second; ++r)
					{
						for (size_t thr = threads.first; thr < threads.second; ++thr)
						{
							std::stringstream null;
							std::unique_ptr<TwoPaCo::VertexEnumerator> vid = CreateEnumerator(fileName, k, filterBits, hf, r, thr, "tmp", temporaryEdge, null);
							std::string candInEdges;
							std::string candOutEdges;
							for (auto & vertex : junctions)
							{
								TwoPaCo::VertexRollingHash hash(vid->GetHashSeed(), vertex.begin(), 1);
								bool result = vid->GetEdges(vertex.begin(), hash, candInEdges, candOutEdges);
								if (!result)
								{
									std::cerr << "Test # " << t << " FAILED" << std::endl;
									return false;
								}
								else
								{
									for (auto ch : inEdge[vertex])
									{
										if (candInEdges.find(ch) == candInEdges.npos)
										{
											bool result = vid->GetEdges(vertex.begin(), hash, candInEdges, candOutEdges);
											std::cerr << "Test # " << t << " FAILED" << std::endl;
											return false;
										}
									}

									for (auto ch : outEdge[vertex])
									{
										if (candOutEdges.find(ch) == candOutEdges.npos)
										{
											std::cerr << "Test # " << t << " FAILED" << std::endl;
											return false;
										}
									}
								}
							}
						}												
					}
				}
			}

			std::remove(temporaryFasta.c_str());
			std::remove(temporaryEdge.c_str());
			std::cerr << "Test # " << t << " PASSED" << std::endl;
		}

		return true;
	}

}