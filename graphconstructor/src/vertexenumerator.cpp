#include <set>
#include <deque>
#include <ctime>
#include <memory>
#include <bitset>
#include <numeric>
#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_set>

#include <boost/ref.hpp>
#include <boost/lockfree/spsc_queue.hpp>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>

#include "lib/SpookyV2.h"
#include "ngramhashing/cyclichash.h"

#include "vertexenumerator.h"

namespace Sibelia
{
	const size_t VertexEnumerator::INVALID_VERTEX = -1;

	namespace
	{
		char Id(char ch)
		{
			return ch;
		}

		uint64_t hashMask;
		typedef CyclicHash<uint64_t> HashFunction;
		typedef std::unique_ptr<HashFunction> HashFunctionPtr;

		template<class F>
		bool IsInBloomFilter(const ConcurrentBitVector & filter, std::vector<HashFunctionPtr> & hf, F f, char farg, uint64_t hash0)
		{
			for (size_t i = 0; i < hf.size(); i++)
			{
				uint64_t hvalue = i == 0 ? hash0 : ((*hf[i]).*f)(farg);
				if (!filter.Get(hvalue))
				{
					return false;
				}
			}

			return true;
		}

		struct Task
		{
			bool isFinal;
			uint64_t start;
			std::string str;
#ifdef _DEBUG
			static const size_t TASK_SIZE = 36;
#else
			static const size_t TASK_SIZE = 1 << 18;
#endif			
			static const size_t GAME_OVER = SIZE_MAX;
			Task() {}
			Task(uint64_t start, bool isFinal, std::string && str) : start(start), isFinal(isFinal), str(std::move(str)) {}
		};

		const size_t QUEUE_CAPACITY = 16;
		typedef boost::lockfree::spsc_queue<Task> TaskQueue;
		typedef std::unique_ptr<TaskQueue> TaskQueuePtr;

		std::string TempFile(size_t record)
		{
			std::stringstream ss;
			ss << record << ".bin";
			return ss.str();
		}

		uint64_t NormHash(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash)
		{
			return std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
		}

		bool Within(uint64_t hvalue, uint64_t low, uint64_t high)
		{
			return hvalue >= low && hvalue <= high;
		}

		void WriteCanonicalRecord(uint64_t posHash0,
			uint64_t negHash0,
			std::string::const_iterator pos,
			size_t vertexLength,
			size_t capacity,
			char posExtend,
			char posPrev,
			std::vector<VertexEnumerator::CompressedString> & out)
		{
			size_t idx = 0;
			size_t element = 0;
			out.push_back(VertexEnumerator::CompressedString());
			std::string::const_reverse_iterator rit(pos + vertexLength);
			if (posHash0 < negHash0)
			{
				char buf[] = { posExtend, posPrev };
				out.back().StrCpy(pos, element, idx, vertexLength, Id);
				out.back().StrCpy(buf, element, idx, 2, Id);
			}
			else
			{
				char buf[] = { DnaString::Reverse(posPrev), DnaString::Reverse(posExtend) };
				out.back().StrCpy(rit, element, idx, vertexLength, DnaString::Reverse);
				out.back().StrCpy(buf, element, idx, 2, Id);
			}
		}

		DnaString Generate(std::string::const_iterator it, size_t size)
		{
			DnaString str(size, size);
			for (size_t i = 0; i < size; i++)
			{
				str.SetChar(i, *(it + i));
			}

			return str;
		}

		std::string RevComp(const std::string & str)
		{
			std::string ret;
			for (auto it = str.rbegin(); it != str.rend(); ++it)
			{
				ret.push_back(DnaString::Reverse(*it));
			}

			return ret;
		}

		void InitializeHashFunctions(const std::vector<HashFunctionPtr> & seed,
			std::vector<HashFunctionPtr> & posEdgeHash,
			std::vector<HashFunctionPtr> & negEdgeHash,
			const std::string & fragment,
			size_t length,
			size_t offset = 0)
		{
			for (size_t i = 0; i < posEdgeHash.size(); i++)
			{
				posEdgeHash[i] = HashFunctionPtr(new HashFunction(*seed[i]));
				negEdgeHash[i] = HashFunctionPtr(new HashFunction(*seed[i]));
				for (auto it = fragment.begin() + offset; it != fragment.begin() + length + offset; ++it)
				{
					posEdgeHash[i]->eat(*it);
				}

				assert(posEdgeHash[i]->hashvalue == posEdgeHash[i]->hash(fragment.substr(offset, length)));
				for (std::string::const_reverse_iterator it(fragment.begin() + length + offset); it != fragment.rend() - offset; ++it)
				{
					char c = DnaString::Reverse(*it);
					negEdgeHash[i]->eat(DnaString::Reverse(*it));
				}

				assert(negEdgeHash[i]->hashvalue == negEdgeHash[i]->hash(RevComp(fragment.substr(offset, length))));
			}
		}		

		void CandidateCheckingWorker(std::pair<uint64_t, uint64_t> bound,
			const std::vector<HashFunctionPtr> & hashFunction,
			const ConcurrentBitVector & bitVector,
			size_t vertexLength,
			size_t & totalRecords,
			TaskQueue & taskQueue,			
			std::ofstream & outFile,
			boost::mutex & outMutex,
			std::unique_ptr<StreamFastaParser::Exception> & error)
		{
			uint64_t low = bound.first;
			uint64_t high = bound.second;
			std::vector<VertexEnumerator::CompressedString> output;
			size_t capacity = DnaString::CalculateCapacity(vertexLength + 2);
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{
					output.clear();
					if (task.start == Task::GAME_OVER)
					{
						break;
					}

					if (task.str.size() < vertexLength)
					{
						continue;
					}

					if (task.start == 0)
					{
						uint64_t posHash0 = hashFunction[0]->hash(task.str.substr(vertexLength));
						uint64_t negHash0 = hashFunction[0]->hash(DnaString::RevComp(task.str.substr(vertexLength)));
						if (Within(std::min(negHash0, posHash0), low, high))
						{
							WriteCanonicalRecord(posHash0, negHash0, task.str.begin(), vertexLength, capacity, 'A', 'A', output);
							WriteCanonicalRecord(posHash0, negHash0, task.str.begin(), vertexLength, capacity, 'A', 'C', output);
						}
					}

					if (task.isFinal)
					{
						uint64_t posHash0 = hashFunction[0]->hash(task.str.substr(task.str.size() - vertexLength, vertexLength));
						uint64_t negHash0 = hashFunction[0]->hash(DnaString::RevComp(task.str.substr(task.str.size() - vertexLength, vertexLength)));
						if (Within(std::min(negHash0, posHash0), low, high))
						{
							WriteCanonicalRecord(posHash0, negHash0, task.str.end() - vertexLength, vertexLength, capacity, 'A', 'A', output);
							WriteCanonicalRecord(posHash0, negHash0, task.str.end() - vertexLength, vertexLength, capacity, 'A', 'C', output);
						}
					}

					if (task.str.size() >= vertexLength + 2)
					{
						std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
						std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
						InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
						for (size_t pos = 1;; ++pos)
						{
							char posPrev = task.str[pos - 1];
							char posExtend = task.str[pos + vertexLength];
							if (Within(std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue), low, high))
							{
								size_t inCount = 0;
								size_t outCount = 0;
								for (int i = 0; i < DnaString::LITERAL.size() && inCount < 2 && outCount < 2; i++)
								{
									char nextCh = DnaString::LITERAL[i];
									char revNextCh = DnaString::Reverse(nextCh);
									if (nextCh == posPrev)
									{
										++inCount;
									}
									else
									{
										uint64_t posHash0 = posVertexHash[0]->hash_prepend(nextCh);
										uint64_t negHash0 = negVertexHash[0]->hash_extend(revNextCh);
										assert(posHash0 == posVertexHash[0]->hash(std::string(1, nextCh) + task.str.substr(pos, vertexLength)));
										assert(negHash0 == negVertexHash[0]->hash(RevComp(std::string(1, nextCh) + task.str.substr(pos, vertexLength))));
										if ((posHash0 <= negHash0 && IsInBloomFilter(bitVector, posVertexHash, &HashFunction::hash_prepend, nextCh, posHash0)) ||
											(posHash0 > negHash0 && IsInBloomFilter(bitVector, negVertexHash, &HashFunction::hash_extend, revNextCh, negHash0)))
										{
											outCount++;
										}
									}

									if (nextCh == posExtend)
									{
										++outCount;
									}
									else
									{
										uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
										uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);
										assert(posHash0 == posVertexHash[0]->hash(task.str.substr(pos, vertexLength) + nextCh));
										assert(negHash0 == negVertexHash[0]->hash(RevComp(task.str.substr(pos, vertexLength) + nextCh)));
										if ((posHash0 <= negHash0 && IsInBloomFilter(bitVector, posVertexHash, &HashFunction::hash_extend, nextCh, posHash0)) ||
											(posHash0 > negHash0 && IsInBloomFilter(bitVector, negVertexHash, &HashFunction::hash_prepend, revNextCh, negHash0)))
										{
											outCount++;
										}
									}
								}

								if (inCount > 1 || outCount > 1)
								{
									WriteCanonicalRecord(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue, task.str.begin() + pos, vertexLength, capacity, posExtend, posPrev, output);
								}
							}

							if (pos + vertexLength + 1 < task.str.size())
							{
								char negExtend = DnaString::Reverse(posExtend);
								char posPrev = task.str[pos];
								char negPrev = DnaString::Reverse(task.str[pos]);
								for (size_t i = 0; i < hashFunction.size(); i++)
								{
									posVertexHash[i]->update(posPrev, posExtend);
									negVertexHash[i]->reverse_update(negExtend, negPrev);
									assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
									assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaString::RevComp(task.str.substr(pos + 1, vertexLength))));
								}
							}
							else
							{
								break;
							}
						}
					}

					if (output.size() > 0)
					{
						boost::lock_guard<boost::mutex> guard(outMutex);
						totalRecords += output.size();
						outFile.write(reinterpret_cast<const char*>(&output[0]), output.size() * sizeof(output[0]));
						if (error != 0)
						{
							return;
						}

						if (!outFile)
						{
							error.reset(new StreamFastaParser::Exception("Can't write to the temporary file"));
							return;
						}
					}
				}
			}
		}

		const uint64_t BINS_COUNT = 1 << 24;
		const uint32_t MAX_COUNTER = UINT32_MAX >> 1;

		void FilterFillerWorker(uint64_t low,
			uint64_t high,
			const std::vector<HashFunctionPtr> & hashFunction,
			ConcurrentBitVector & filter,
			size_t edgeLength,
			TaskQueue & taskQueue)
		{
			std::vector<uint64_t> hvalue(hashFunction.size());
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{
					if (task.start == Task::GAME_OVER)
					{
						break;
					}

					if (task.str.size() < edgeLength)
					{
						continue;
					}

					size_t vertexLength = edgeLength - 1;
					std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
					std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
					InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength);
					for (size_t pos = 0; pos + edgeLength - 1 < task.str.size(); ++pos)
					{
						char nextCh = task.str[pos + edgeLength - 1];
						char revNextCh = DnaString::Reverse(nextCh);
						uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
						uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);
						uint64_t fistMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						char prevCh = task.str[pos];

						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							if (posHash0 < negHash0)
							{
								hvalue[i] = posVertexHash[i]->hash_extend(nextCh);
							}
							else
							{
								hvalue[i] = negVertexHash[i]->hash_prepend(revNextCh);
							}
						}

						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							posVertexHash[i]->update(prevCh, nextCh);
							assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
							negVertexHash[i]->reverse_update(DnaString::Reverse(nextCh), DnaString::Reverse(prevCh));
							assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(RevComp(task.str.substr(pos + 1, vertexLength))));
						}

						uint64_t secondMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						if (Within(fistMinHash0, low, high) || Within(secondMinHash0, low, high))
						{
							for (uint64_t hv : hvalue)
							{
								filter.SetConcurrently(hv);
							}
						}
					}
				}
			}
		}


		/*
		boost::mutex counterMutex;
		std::unordered_set<uint64_t> edgeCounter;
		std::vector<uint32_t> trueBin(BINS_COUNT, 0);*/

		void InitialFilterFillerWorker(uint64_t binSize,
			const std::vector<HashFunctionPtr> & hashFunction,
			ConcurrentBitVector & filter,
			size_t vertexLength,
			TaskQueue & taskQueue,
			std::atomic<uint32_t> * binCounter)
		{
			size_t edgeLength = vertexLength + 1;
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{
					if (task.start == Task::GAME_OVER)
					{
						break;
					}

					if (task.str.size() < edgeLength)
					{
						continue;
					}

					size_t vertexLength = edgeLength - 1;
					std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
					std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
					InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength);
					for (size_t pos = 0; pos + edgeLength - 1 < task.str.size(); ++pos)
					{
						char nextCh = task.str[pos + edgeLength - 1];
						char revNextCh = DnaString::Reverse(nextCh);
						uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
						uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);
						uint64_t fistMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						char prevCh = task.str[pos];
						bool wasSet = true;
						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							uint64_t hvalue = posHash0 <= negHash0 ? posVertexHash[i]->hash_extend(nextCh) : negVertexHash[i]->hash_prepend(revNextCh);
							if (!filter.Get(hvalue))
							{
								wasSet = false;
								filter.SetConcurrently(hvalue);
							}
						}
											
						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							posVertexHash[i]->update(prevCh, nextCh);
							assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
							negVertexHash[i]->reverse_update(DnaString::Reverse(nextCh), DnaString::Reverse(prevCh));
							assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(RevComp(task.str.substr(pos + 1, vertexLength))));
						}


						uint64_t secondMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						if (!wasSet)
						{
							uint64_t value[] = { fistMinHash0, secondMinHash0 };
							for (uint64_t v : value)
							{
								uint64_t bin = v / binSize;
								if (binCounter[bin] < MAX_COUNTER)
								{
									binCounter[bin].fetch_add(1);
								}
							}
						}
						/*
						{
							DnaString edge(task.str.substr(pos, edgeLength));
							if (!(posHash0 <= negHash0))
							{
								edge = edge.RevComp();
							}

							boost::lock_guard<boost::mutex> g(counterMutex);
							if (edgeCounter.count(edge.GetBody()) == 0)
							{
								edgeCounter.insert(edge.GetBody());
								uint64_t value[] = { fistMinHash0, secondMinHash0 };
								for (uint64_t v : value)
								{
									uint64_t bin = v / binSize;
									if (trueBin[bin] < MAX_COUNTER)
									{
										trueBin[bin] += 1;
									}
								}
							}
						}*/
					}
				}
			}
		}

		void DistributeTasks(const std::vector<std::string> & fileName, size_t overlapSize, std::vector<TaskQueuePtr> & taskQueue)
		{
			for (size_t file = 0; file < fileName.size(); file++)
			{
				size_t record = 0;
				const std::string & nowFileName = fileName[file];
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
				{
					char ch;
					std::string buf;
					uint64_t prev = 0;
					uint64_t start = 0;
					bool over = false;
					do
					{
						over = !parser.GetChar(ch);
						if (!over)
						{
							start++;
							buf.push_back(ch);
						}

						if (buf.size() >= overlapSize && (buf.size() == Task::TASK_SIZE || over))
						{
							for (bool found = false; !found;)
							{
								for (TaskQueuePtr & q : taskQueue)
								{
									if (q->write_available() > 0)
									{
										std::string overlap;
										if (!over)
										{
											overlap.assign(buf.end() - overlapSize, buf.end());
										}

										q->push(Task(prev, over, std::move(buf)));
										prev = start - overlapSize + 1;
										buf.swap(overlap);
										found = true;
										break;
									}
								}
							}
						}

					} while (!over);
				}
			}

			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				while (taskQueue[i]->write_available() == 0)
				{
					boost::this_thread::sleep_for(boost::chrono::nanoseconds(1000000));
				}

				taskQueue[i]->push(Task(Task::GAME_OVER, true, std::string()));
			}
		}

		typedef std::vector<VertexEnumerator::CompressedString>::iterator RecordIterator;
		
		struct Candidate
		{
			char prev;
			char extend;
		};

		Candidate Dereference(size_t size, RecordIterator it)
		{
			Candidate ret;
			ret.extend = it->GetChar(size);
			ret.prev = it->GetChar(size + 1);
			return ret;
		}

		bool EqualVertex(size_t vertexSize_, RecordIterator it1, RecordIterator it2)
		{
			size_t remain = vertexSize_;
			for (size_t i = 0; remain > 0; i++)
			{
				size_t current = std::min(remain, DnaString::UNIT_CAPACITY);
				uint64_t apiece = it1->str[i];
				uint64_t bpiece = it2->str[i];
				if (current != DnaString::UNIT_CAPACITY)
				{
					uint64_t mask = (uint64_t(1) << (current * 2)) - 1;
					apiece &= mask;
					bpiece &= mask;
				}

				if (apiece != bpiece)
				{
					return false;
				}

				remain -= current;
			}

			return true;
		}

		bool IsSelfRevComp(size_t vertexSize, const VertexEnumerator::CompressedString & str)
		{
			for (size_t i = 0; i < vertexSize; i++)
			{
				if (str.GetChar(i) != DnaString::Reverse(str.GetChar(vertexSize - i - 1)))
				{
					return false;
				}
			}

			return true;
		}

		struct TrueBifurcations
		{
			size_t vertexSize;
			uint64_t falsePositives;
			std::vector<VertexEnumerator::CompressedString> * candidate;
			std::vector<VertexEnumerator::CompressedString> * out;
			boost::mutex * outMutex;
			TrueBifurcations(std::vector<VertexEnumerator::CompressedString> * candidate,
				std::vector<VertexEnumerator::CompressedString> * out,
				boost::mutex * outMutex,
				size_t vertexSize) :
				candidate(candidate), out(out), vertexSize(vertexSize), outMutex(outMutex), falsePositives(0) {}

			uint64_t operator()(const tbb::blocked_range<RecordIterator> & range, uint64_t init) const
			{				
				RecordIterator vectorBegin = candidate->begin();
				RecordIterator vectorEnd = candidate->end();
				RecordIterator base = range.begin();
				if (base > vectorBegin)
				{
					if (EqualVertex(vertexSize, base, base - 1))
					{
						for (++base; base < range.end() && EqualVertex(vertexSize, base, range.begin()); ++base);
					}
				}

				uint64_t falsePositives = init;
				while(base < range.end())
				{					
					bool bifurcation = false;
					RecordIterator next = base;
					Candidate baseCandidate(Dereference(vertexSize, base));
					bool selfRevComp = IsSelfRevComp(vertexSize, *base);
					for (; next < vectorEnd; ++next)
					{
						Candidate nextCandidate(Dereference(vertexSize, next));
					//	std::cout << DnaString(vertexSize, &*next.GetIterator()).ToString() << ' ' << nextCandidate.extend << nextCandidate.prev << std::endl;
						if (!(EqualVertex(vertexSize, base, next)))
						{
							break;
						}
						else if (!bifurcation)
						{
							bifurcation = baseCandidate.prev != nextCandidate.prev || baseCandidate.extend != nextCandidate.extend;
							if (selfRevComp) 
							{
								bifurcation = bifurcation ||
									baseCandidate.prev != DnaString::Reverse(nextCandidate.extend) ||
									baseCandidate.extend != DnaString::Reverse(nextCandidate.prev);
							}
						}
					}

					if (bifurcation)
					{
						boost::lock_guard<boost::mutex> guard(*outMutex);
						out->push_back(VertexEnumerator::CompressedString());
						size_t remain = vertexSize;
						for (size_t i = 0; remain > 0; i++)
						{
							uint64_t piece = base->str[i];
							size_t current = std::min(remain, DnaString::UNIT_CAPACITY);
							if (current != DnaString::UNIT_CAPACITY)
							{
								uint64_t mask = (uint64_t(1) << (current * 2)) - 1;
								piece &= mask;
							}

							base->str[i] = piece;
							remain -= current;
						}

					}
					else
					{
						falsePositives++;
					}

					base = next;
				}

				return falsePositives;
			}
		};


		class VertexLess
		{
		public:
			VertexLess(size_t vertexSize) : vertexSize_(vertexSize)
			{

			}

			bool operator() (const uint64_t * v1, const uint64_t * v2) const
			{
				size_t remain = vertexSize_;
				for (size_t i = 0; remain > 0; i++)
				{
					size_t current = std::min(remain, DnaString::UNIT_CAPACITY);
					uint64_t apiece = v1[i];
					uint64_t bpiece = v2[i];
					if (current != DnaString::UNIT_CAPACITY)
					{
						uint64_t mask = (uint64_t(1) << (current * 2)) - 1;
						apiece &= mask;
						bpiece &= mask;
					}

					if (apiece != bpiece)
					{
						return apiece < bpiece;
					}

					remain -= current;
				}

				return false;
			}

		private:
			size_t vertexSize_;
		};
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName,
		size_t vertexLength,
		size_t filterSize,
		size_t hashFunctions,
		size_t rounds,
		size_t threads,
		size_t aggregationThreads,
		const std::string & tmpFileName) :
		vertexSize_(vertexLength)
	{		
		uint64_t realSize = uint64_t(1) << filterSize;
		hashMask = realSize - 1;
		std::cout << "Threads = " << threads << std::endl;		
		std::cout << "Aggregation threads = " << aggregationThreads << std::endl;
		std::cout << "Hash functions = " << hashFunctions << std::endl;
		std::cout << "Filter size = " << realSize << std::endl;
		std::cout << "Files: " << std::endl;
		for (const std::string & fn : fileName)
		{
			std::cout << fn << std::endl;
		}

		std::cout << std::string(80, '-') << std::endl;
		const uint64_t BIN_SIZE = std::max(uint64_t(1), realSize / BINS_COUNT);
		std::vector<HashFunctionPtr> hashFunction(hashFunctions);
		for (HashFunctionPtr & ptr : hashFunction)
		{
			ptr = HashFunctionPtr(new HashFunction(vertexLength, filterSize));
		}

		size_t edgeLength = vertexLength + 1;
		std::atomic<uint32_t> * binCounter = new std::atomic<uint32_t>[BINS_COUNT];
		{
			std::fill(binCounter, binCounter + BINS_COUNT, 0);
			std::vector<boost::thread> workerThread(threads);
			ConcurrentBitVector bitVector(realSize);
			std::vector<TaskQueuePtr> taskQueue;
			for (size_t i = 0; i < workerThread.size(); i++)
			{
				taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
				workerThread[i] = boost::thread(InitialFilterFillerWorker,					
					BIN_SIZE,
					boost::cref(hashFunction),
					boost::ref(bitVector),
					vertexLength,
					boost::ref(*taskQueue[i]),
					binCounter);
			}

			DistributeTasks(fileName, edgeLength, taskQueue);
			for (size_t i = 0; i < workerThread.size(); i++)
			{
				workerThread[i].join();
			}
		}


		rounds = 1;
		double roundSize = 0;
		for (;; ++rounds)
		{
			roundSize = std::accumulate(binCounter, binCounter + BINS_COUNT, 0.f) / rounds;
			if (realSize / roundSize > 16)
			{
				break;
			}
		}

		/*
		long long diff = 0;
		size_t nonZero = 0;
		std::ofstream diffOut("diff.txt");
		for (size_t i = 0; i < BINS_COUNT; i++)
		{
			long long c1 = binCounter[i];
			long long c2 = trueBin[i];
			diff += std::abs(c1 - c2);
			if (c1 || c2)
			{
				diffOut << c1 << " " << c2 << std::endl;
				nonZero++;
			}
		}

		std::cout << "Diff = " << diff << std::endl;
		std::cout << "Avg diff = " << double(diff) / nonZero << std::endl;*/

		uint64_t low = 0;
		uint64_t high = 0;
		size_t lowBoundary = 0;
		uint64_t totalFpCount = 0;
		for (size_t round = 0; round < rounds; round++)
		{
			size_t totalRecords = 0;
			time_t mark = time(0);			
			uint64_t accumulated = binCounter[lowBoundary];
			for (++lowBoundary; lowBoundary < BINS_COUNT; ++lowBoundary)
			{				
				if (accumulated <= roundSize || round + 1 == rounds)
				{
					accumulated += binCounter[lowBoundary];
				}
				else
				{					
					break;
				}
			}

			std::cout << "Ratio = " << double(realSize) / accumulated << std::endl;
			uint64_t high = lowBoundary * BIN_SIZE;
			{
				std::vector<TaskQueuePtr> taskQueue;
				ConcurrentBitVector bitVector(realSize);
				std::vector<boost::thread> workerThread(threads);
				std::cout << "Round " << round << ", " << low << ":" << high << std::endl;
				std::cout << "Counting\tEnumeration\tAggregation" << std::endl;
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
					workerThread[i] = boost::thread(FilterFillerWorker,
						low,
						high,
						boost::cref(hashFunction),
						boost::ref(bitVector),
						edgeLength,
						boost::ref(*taskQueue[i]));
				}

				DistributeTasks(fileName, edgeLength, taskQueue);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i].join();
				}

				std::cout << time(0) - mark << "\t";
				mark = time(0);
				boost::mutex tmpFileMutex;
				std::ofstream tmpFile(tmpFileName.c_str(), std::ios_base::binary);
				if (!tmpFile)
				{
					throw StreamFastaParser::Exception("Can't open the temporary file");
				}

				std::unique_ptr<StreamFastaParser::Exception> error;
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i] = boost::thread(CandidateCheckingWorker,
						std::make_pair(low, high),						
						boost::cref(hashFunction),
						boost::cref(bitVector),
						vertexLength,
						boost::ref(totalRecords),
						boost::ref(*taskQueue[i]),
						boost::ref(tmpFile),
						boost::ref(tmpFileMutex),
						boost::ref(error));
				}

				if (error != 0)
				{
					throw StreamFastaParser::Exception(*error);
				}

				DistributeTasks(fileName, vertexLength + 1, taskQueue);
				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					workerThread[i].join();
				}

				std::cout << time(0) - mark << "\t";
			}

			mark = time(0);
			std::vector<CompressedString> candidate(totalRecords);
			std::ifstream tmpFile(tmpFileName.c_str(), std::ios_base::binary);
			if (!tmpFile)
			{
				throw StreamFastaParser::Exception("Can't open the temporary file");
			}

			if (totalRecords > 0)
			{
				tmpFile.read(reinterpret_cast<char*>(&candidate[0]), totalRecords * sizeof(candidate[0]));
				if (!tmpFile)
				{
					throw StreamFastaParser::Exception("The temporary file is corrupted");
				}
			}

			boost::mutex outMutex;
			RecordIterator begin = candidate.begin();
			RecordIterator end = candidate.end();
		//	tbb::parallel_sort(begin, end, Comparator<VertexLess>(VertexLess(vertexSize_))); 
			uint64_t falsePositives = tbb::parallel_reduce(tbb::blocked_range<RecordIterator>(begin, end),
				uint64_t(0),
				TrueBifurcations(&candidate, &bifurcation_, &outMutex, vertexSize_),
				std::plus<uint64_t>());
			//uint64_t falsePositives = TrueBifurcations(&candidate, &bifurcation_, &outMutex, vertexSize_)(tbb::blocked_range<RecordIterator>(begin, end), 0);
//			size_t falsePositives = 0;
			std::cout << time(0) - mark << std::endl;
			std::cout << "Vertex count = " << bifurcation_.size() << std::endl;
			std::cout << "FP count = " << falsePositives << std::endl;
			std::cout << "Records = " << candidate.size() << std::endl;
			std::cout << std::string(80, '-') << std::endl;
			totalFpCount += falsePositives;
			low = high + 1;
		}

		std::cout << "Total FPs = " << totalFpCount << std::endl;
		delete[] binCounter;
		tbb::parallel_sort(bifurcation_.begin(), bifurcation_.end());
	}

	size_t VertexEnumerator::GetVerticesCount() const
	{
		return bifurcation_.size();
	}

	size_t VertexEnumerator::GetId(const DnaString & vertex) const
	{
		CompressedString check[2] = { vertex.GetBody()[0], vertex.RevComp().GetBody()[0] };
		for (CompressedString & str : check)
		{
			std::vector<CompressedString>::const_iterator it = std::lower_bound(bifurcation_.begin(), bifurcation_.end(), str, VertexLess(vertexSize_));
			if (it != bifurcation_.end() && *it == str)
			{
				return it - bifurcation_.begin();
			}
		}
		
		return INVALID_VERTEX;
	}
}