#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#define MAX_CAPACITY 10

#include <deque>
#include <cstdio>
#include <numeric>
#include <sstream>
#include <unordered_map>

#include <tbb/tbb.h>
#include <tbb/mutex.h>
#include <tbb/compat/thread>
#include <tbb/spin_rw_mutex.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>

#include <junctionapi/junctionapi.h>

#include "streamfastaparser.h"
#include "bifurcationstorage.h"
#include "candidateoccurence.h"
#include "concurrentbitvector.h"

#if _MSC_VER && !__INTEL_COMPILER
	#define MSVS
	#include <intrin.h> 
#endif

namespace TwoPaCo
{
	class VertexEnumerator
	{
	public:
		virtual size_t GetVerticesCount() const = 0;
		virtual size_t GetTotalVerticesCount() const = 0;
		virtual std::pair<uint64_t, uint64_t> GetId(const std::string & vertex) const = 0;

		virtual ~VertexEnumerator()
		{

		}
	};

	std::unique_ptr<VertexEnumerator> CreateEnumerator(const std::vector<std::string> & fileName,
		size_t vertexLength,
		size_t filterSize,
		size_t hashFunctions,
		size_t rounds,
		size_t threads,
		const std::string & tmpFileName,
		const std::string & outFileName);

	template<size_t CAPACITY>
	class VertexEnumeratorImpl : public VertexEnumerator
	{
	private:
		static const size_t BUF_SIZE = 1 << 24;
		BifurcationStorage<CAPACITY> bifStorage_;
		typedef CompressedString<CAPACITY> DnaString;
		typedef CandidateOccurence<CAPACITY> Occurence;

		class FilterFillerWorker;
		class InitialFilterFillerWorker;

		class OccurenceHash
		{
		public:
			uint64_t operator()(const Occurence & occurence) const
			{
				return occurence.Hash();
			}
		};

		class DnaStringHash
		{
		public:
			uint64_t operator()(const DnaString & dnaString) const
			{
				return dnaString.Hash();
			}
		};

		class OccurenceEquality
		{
		public:
			bool operator()(const Occurence & occurence1, const Occurence & occurence2) const
			{
				return occurence1.EqualBase(occurence2);
			}
		};

		typedef std::unordered_map<DnaString, uint64_t, DnaStringHash> BifurcationMap;
		typedef tbb::concurrent_unordered_set<Occurence, OccurenceHash, OccurenceEquality> OccurenceSet;

	public:

		size_t GetVerticesCount() const
		{
			return bifStorage_.GetDistinctVerticesCount();
		}

		size_t GetTotalVerticesCount() const
		{
			return bifStorage_.GetTotalVerticesCount();
		}

		std::pair<uint64_t, uint64_t> GetId(const std::string & vertex) const
		{
			return bifStorage_.GetId(vertex.begin());
		}

		VertexEnumeratorImpl(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			const std::string & tmpDirName,
			const std::string & outFileNamePrefix) :
			vertexSize_(vertexLength)
		{
			uint64_t realSize = uint64_t(1) << filterSize;
			std::cout << "Threads = " << threads << std::endl;
			std::cout << "Vertex length = " << vertexLength << std::endl;
			std::cout << "Hash functions = " << hashFunctions << std::endl;
			std::cout << "Filter size = " << realSize << std::endl;
			std::cout << "Capacity = " << CAPACITY << std::endl;
			std::cout << "Files: " << std::endl;
			for (const std::string & fn : fileName)
			{
				std::cout << fn << std::endl;
			}

			std::ofstream logFile((tmpDirName + "/log.txt").c_str());
			if (!logFile)
			{
				throw StreamFastaParser::Exception("Can't open the log file");
			}

			tbb::mutex errorMutex;
			std::unique_ptr<std::runtime_error> error;
			std::vector<HashFunctionPtr> hashFunction(hashFunctions);
			for (HashFunctionPtr & ptr : hashFunction)
			{
				ptr = HashFunctionPtr(new HashFunction(vertexLength, filterSize - 4));
			}

			size_t edgeLength = vertexLength + 1;
			std::vector<TaskQueuePtr> taskQueue(threads);
			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				taskQueue[i].reset(new TaskQueue());
				taskQueue[i]->set_capacity(QUEUE_CAPACITY);
			}

			const uint64_t BIN_SIZE = max(uint64_t(1), realSize / BINS_COUNT);
			std::atomic<uint32_t> * binCounter = 0;

			
			if (rounds > 1)
			{
				std::cout << "Splitting the input kmers set...";
				std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
				binCounter = new std::atomic<uint32_t>[BINS_COUNT];
				std::fill(binCounter, binCounter + BINS_COUNT, 0);
				ConcurrentBitVector bitVector(realSize);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					InitialFilterFillerWorker worker(BIN_SIZE,
						hashFunction,
						bitVector,
						vertexLength,
						*taskQueue[i],
						binCounter);
					workerThread[i].reset(new tbb::tbb_thread(worker));
				}

				DistributeTasks(fileName, edgeLength, taskQueue, error, errorMutex, logFile);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i]->join();
				}
			}

			double roundSize = 0;
			if (rounds > 1)
			{
				roundSize = double(std::accumulate(binCounter, binCounter + BINS_COUNT, size_t(0))) / rounds;
			}

			std::cout << std::string(80, '-') << std::endl;
			uint64_t low = 0;
			uint64_t high = 0;
			size_t lowBoundary = 0;
			uint64_t totalFpCount = 0;
			uint64_t verticesCount = 0;
			std::ofstream bifurcationTempWrite((tmpDirName + "/bifurcations.bin").c_str(), ios::binary);
			if (!bifurcationTempWrite)
			{
				throw StreamFastaParser::Exception("Can't create a temp file");
			}

			time_t mark;
			for (size_t round = 0; round < rounds; round++)
			{
				std::atomic<uint64_t> marks;
				marks = 0;
				mark = time(0);
				if (rounds > 1)
				{
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

					high = lowBoundary * BIN_SIZE;
				}
				else
				{
					high = realSize;
				}

				{
					ConcurrentBitVector bitVector(realSize);
					std::cout << "Round " << round << ", " << low << ":" << high << std::endl;
					std::cout << "Pass\tFilling\tFiltering" << std::endl << "1\t";
					{
						std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							FilterFillerWorker worker(low,
								high,
								std::cref(hashFunction),
								std::ref(bitVector),
								edgeLength,
								std::ref(*taskQueue[i]));
							workerThread[i].reset(new tbb::tbb_thread(worker));
						}

						DistributeTasks(fileName, edgeLength, taskQueue, error, errorMutex, logFile);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							workerThread[i]->join();
						}
					}

					std::cout << time(0) - mark << "\t";
					mark = time(0);
					{
						std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							CandidateCheckingWorker worker(std::make_pair(low, high),
								hashFunction,
								bitVector,
								vertexLength,
								*taskQueue[i],
								tmpDirName,
								marks,
								error,
								errorMutex);

							workerThread[i].reset(new tbb::tbb_thread(worker));
						}

						DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile);
						for (size_t i = 0; i < taskQueue.size(); i++)
						{
							workerThread[i]->join();
						}

						if (error != 0)
						{
							throw *error;
						}
					}

					std::cout << time(0) - mark << "\t" << std::endl;
				}

				mark = time(0);
				tbb::spin_rw_mutex mutex;
				std::cout << "2\t";
				OccurenceSet occurenceSet(1 << 20);
				{
					std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
					for (size_t i = 0; i < workerThread.size(); i++)
					{
						CandidateFilteringWorker worker(hashFunction,
							vertexLength,
							*taskQueue[i],
							occurenceSet,
							mutex,
							tmpDirName,
							error,
							errorMutex);

						workerThread[i].reset(new tbb::tbb_thread(worker));
					}

					DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile);
					for (size_t i = 0; i < taskQueue.size(); i++)
					{
						workerThread[i]->join();
					}

					if (error != 0)
					{
						throw std::runtime_error(*error);
					}

					std::cout << time(0) - mark << "\t";
				}

				size_t falsePositives = 0;
				mark = time(0);
				size_t truePositives = TrueBifurcations(occurenceSet, bifurcationTempWrite, vertexSize_, falsePositives);
				std::cout << time(0) - mark << std::endl;
				std::cout << "True junctions count = " << truePositives << std::endl;
				std::cout << "False junctions count = " << falsePositives << std::endl;
				std::cout << "Hash table size = " << occurenceSet.size() << std::endl;
				std::cout << "Candidate marks count = " << marks << std::endl;
				std::cout << std::string(80, '-') << std::endl;
				totalFpCount += falsePositives;
				verticesCount += truePositives;
				low = high + 1;
			}

			if (rounds > 1)
			{
				delete[] binCounter;
			}

			mark = time(0);
			std::string bifurcationTempReadName = (tmpDirName + "/bifurcations.bin");
			bifurcationTempWrite.close();
			{
				std::ifstream bifurcationTempRead(bifurcationTempReadName.c_str(), ios::binary);
				if (!bifurcationTempRead)
				{
					throw StreamFastaParser::Exception("Can't open the temp file");
				}

				bifStorage_.Init(bifurcationTempRead, verticesCount, vertexLength, threads);
			}

			std::remove(bifurcationTempReadName.c_str());
			std::cout << "Reallocating bifurcations time: " << time(0) - mark << std::endl;

			mark = time(0);

			std::atomic<uint64_t> occurence;
			std::atomic<uint64_t> currentPiece;
			std::atomic<uint64_t> currentStubVertex;
			JunctionPositionWriter posWriter(outFileNamePrefix);
			occurence = currentPiece = 0;
			currentStubVertex = verticesCount * 2;
			{
				std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					EdgeConstructionWorker worker(vertexLength,
						*taskQueue[i],
						bifStorage_,
						posWriter,
						currentPiece,
						occurence,
						currentStubVertex,
						error,
						errorMutex);

					workerThread[i].reset(new tbb::tbb_thread(worker));
				}

				DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile);
				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					workerThread[i]->join();
				}
			}

			if (error != 0)
			{
				throw std::runtime_error(*error);
			}

			std::cout << "True marks count: " << occurence << std::endl;
			std::cout << "Edges construction time: " << time(0) - mark << std::endl;
			std::cout << std::string(80, '-') << std::endl;
		}

	private:

		static const size_t QUEUE_CAPACITY = 16;
		static const uint64_t BINS_COUNT = 1 << 24;

		enum StrandComparisonResult
		{
			positiveLess,
			negativeLess,
			tie
		};


		static uint64_t NormHash(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash)
		{
			return std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
		}

		static bool Within(uint64_t hvalue, uint64_t low, uint64_t high)
		{
			return hvalue >= low && hvalue <= high;
		}

		static void InitializeHashFunctions(const std::vector<HashFunctionPtr> & seed,
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
					char ch = DnaChar::ReverseChar(*it);
					negEdgeHash[i]->eat(DnaChar::ReverseChar(*it));
				}

				assert(negEdgeHash[i]->hashvalue == negEdgeHash[i]->hash(DnaChar::ReverseCompliment(fragment.substr(offset, length))));
			}
		}

		static std::string CandidateMaskFileName(const std::string & directory, size_t sequence, size_t pos)
		{
			std::stringstream ss;
			ss << directory << "/" << sequence << "_" << pos << ".tmp";
			return ss.str();
		}

		static void ReportError(tbb::mutex & errorMutex, std::unique_ptr<std::runtime_error> & error, const std::string & msg)
		{
			errorMutex.lock();
			if (error == 0)
			{
				error.reset(new StreamFastaParser::Exception(msg));
			}

			errorMutex.unlock();
		}

		static int PopCount(uint16_t mask)
		{
#ifdef MSVS
			return __popcnt16(mask);
#else
			int ret = 0;		
			for (size_t j = 0; j < sizeof(mask) * 8; j++)
			{
				if (mask & (1 << j))
				{
					++ret;
				}
			}

			return ret;
#endif
		}

		class CandidateCheckingWorker
		{
		public:
			CandidateCheckingWorker(std::pair<uint64_t, uint64_t> bound,
				const std::vector<HashFunctionPtr> & hashFunction,
				const ConcurrentBitVector & bitVector,
				size_t vertexLength,
				TaskQueue & taskQueue,
				const std::string & tmpDirectory,
				std::atomic<uint64_t> & marksCount,
				std::unique_ptr<std::runtime_error> & error,
				tbb::mutex & errorMutex) : bound(bound), hashFunction(hashFunction), bitVector(bitVector), vertexLength(vertexLength), taskQueue(taskQueue),
				tmpDirectory(tmpDirectory), marksCount(marksCount), error(error), errorMutex(errorMutex)
			{

			}

			void operator()()
			{
				uint64_t low = bound.first;
				uint64_t high = bound.second;
				ConcurrentBitVector candidateMask(Task::TASK_SIZE);
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
					{
						if (task.start == Task::GAME_OVER)
						{
							break;
						}

						if (task.str.size() < vertexLength)
						{
							continue;
						}

						size_t edgeLength = vertexLength + 1;
						if (task.str.size() >= vertexLength + 2)
						{
							candidateMask.Reset();
							std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
							std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
							InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
							size_t definiteCount = std::count_if(task.str.begin() + 1, task.str.begin() + vertexLength + 1, DnaChar::IsDefinite);
							for (size_t pos = 1;; ++pos)
							{
								char posPrev = task.str[pos - 1];
								char posExtend = task.str[pos + vertexLength];
								assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
								if (Within(min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue), low, high) && definiteCount == vertexLength)
								{
									int bitCount = (DnaChar::IsDefinite(posPrev) && DnaChar::IsDefinite(posExtend)) ? 0 : 2;
									if (bitCount == 0)
									{
										uint16_t mask = (1 << 16) - 1;
										StrandComparisonResult res = DetermineStrand(posVertexHash, negVertexHash);
										for (size_t i = 0; i < posVertexHash.size() && PopCount(mask) > 1; i++)
										{
											if (res == positiveLess || res == tie)
											{
												mask &= bitVector.GetElement(posVertexHash[i]->hashvalue);
											}
											else
											{
												mask &= bitVector.GetElement(negVertexHash[i]->hashvalue);
											}
										}

										bitCount = PopCount(mask);
									}

									if (bitCount > 1)
									{
										++marksCount;
										candidateMask.SetBitConcurrently(pos);
									}
								}

								if (pos + edgeLength < task.str.size())
								{
									char negExtend = DnaChar::ReverseChar(posExtend);
									char posPrev = task.str[pos];
									char negPrev = DnaChar::ReverseChar(task.str[pos]);
									definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
									for (size_t i = 0; i < hashFunction.size(); i++)
									{
										posVertexHash[i]->update(posPrev, posExtend);
										negVertexHash[i]->reverse_update(negExtend, negPrev);
										assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
										assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
									}
								}
								else
								{
									break;
								}
							}

							try
							{
								candidateMask.WriteToFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start));
							}
							catch (std::runtime_error & err)
							{
								ReportError(errorMutex, error, err.what());
							}
						}
					}
				}
			}

		private:
			std::pair<uint64_t, uint64_t> bound;
			const std::vector<HashFunctionPtr> & hashFunction;
			const ConcurrentBitVector & bitVector;
			size_t vertexLength;
			TaskQueue & taskQueue;
			const std::string & tmpDirectory;
			std::atomic<uint64_t> & marksCount;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};


		class CandidateFilteringWorker
		{
		public:
			CandidateFilteringWorker(const std::vector<HashFunctionPtr> & hashFunction,
				size_t vertexLength,
				TaskQueue & taskQueue,
				OccurenceSet & occurenceSet,
				tbb::spin_rw_mutex & mutex,
				const std::string & tmpDirectory,
				std::unique_ptr<std::runtime_error> & error,
				tbb::mutex & errorMutex) : hashFunction(hashFunction), vertexLength(vertexLength), taskQueue(taskQueue), occurenceSet(occurenceSet),
				mutex(mutex), tmpDirectory(tmpDirectory), error(error), errorMutex(errorMutex)
			{

			}

			void operator()()
			{
				ConcurrentBitVector candidateMask(Task::TASK_SIZE);
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
					{
						if (task.start == Task::GAME_OVER)
						{
							break;
						}

						if (task.str.size() < vertexLength)
						{
							continue;
						}

						std::vector<HashFunctionPtr> posVertexHash(1);
						std::vector<HashFunctionPtr> negVertexHash(1);
						size_t edgeLength = vertexLength + 1;
						if (task.str.size() >= vertexLength + 2)
						{
							InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
							{
								try
								{
									candidateMask.ReadFromFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start), true);
								}
								catch (std::runtime_error & err)
								{
									ReportError(errorMutex, error, err.what());
								}
							}

							for (size_t pos = 1;; ++pos)
							{
								char posPrev = task.str[pos - 1];
								char posExtend = task.str[pos + vertexLength];
								if (candidateMask.GetBit(pos))
								{
									Occurence now;
									now.Set(posVertexHash[0]->hashvalue,
										negVertexHash[0]->hashvalue,
										task.str.begin() + pos,
										vertexLength,
										posExtend,
										posPrev,
										false);

									size_t count = 0;
									size_t inUnknownCount = now.Prev() == 'N' ? 1 : 0;
									size_t outUnknownCount = now.Next() == 'N' ? 1 : 0;
									auto ret = occurenceSet.insert(now);
									auto it = ret.first;
									if (!it->IsBifurcation())
									{
										if (it->Next() != now.Next() || it->Prev() != now.Prev() || inUnknownCount > 1 || outUnknownCount > 1)
										{
											it->MakeBifurcation();
										}
									}
								}

								if (pos + edgeLength < task.str.size())
								{
									char negExtend = DnaChar::ReverseChar(posExtend);
									char posPrev = task.str[pos];
									char negPrev = DnaChar::ReverseChar(task.str[pos]);
									for (size_t i = 0; i < 1; i++)
									{
										posVertexHash[i]->update(posPrev, posExtend);
										negVertexHash[i]->reverse_update(negExtend, negPrev);
										assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
										assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
									}
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

		private:
			const std::vector<HashFunctionPtr> & hashFunction;
			size_t vertexLength;
			TaskQueue & taskQueue;
			OccurenceSet & occurenceSet;
			tbb::spin_rw_mutex & mutex;
			const std::string & tmpDirectory;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};

		struct EdgeResult
		{
			uint32_t pieceId;
			std::vector<JunctionPosition> junction;
		};

		static bool FlushEdgeResults(std::deque<EdgeResult> & result,
			JunctionPositionWriter & writer,
			std::atomic<uint64_t> & currentPiece)
		{
			if (result.size() > 0 && result.front().pieceId == currentPiece)
			{
				for (auto junction : result.front().junction)
				{
					writer.WriteJunction(junction);
				}

				++currentPiece;
				result.pop_front();
				return true;
			}

			return false;
		}

		class EdgeConstructionWorker
		{
		public:
			EdgeConstructionWorker(size_t vertexLength,
				TaskQueue & taskQueue,
				const BifurcationStorage<CAPACITY> & bifStorage,
				JunctionPositionWriter & writer,
				std::atomic<uint64_t> & currentPiece,
				std::atomic<uint64_t> & occurences,
				std::atomic<uint64_t> & currentStubVertexId,
				std::unique_ptr<std::runtime_error> & error,
				tbb::mutex & errorMutex) : vertexLength(vertexLength), taskQueue(taskQueue), bifStorage(bifStorage),
				writer(writer), currentPiece(currentPiece), occurences(occurences),
				currentStubVertexId(currentStubVertexId), error(error), errorMutex(errorMutex)
			{

			}

			void operator()()
			{
				try
				{
					DnaString bitBuf;
					std::deque<EdgeResult> result;
					while (true)
					{
						Task task;
						if (taskQueue.try_pop(task))
						{
							if (task.start == Task::GAME_OVER)
							{
								break;
							}

							if (task.str.size() < vertexLength)
							{
								continue;
							}

							const std::vector<HashFunctionPtr> & hashFunction = bifStorage.GetHashFunctions();
							std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
							std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
							size_t edgeLength = vertexLength + 1;
							if (task.str.size() >= vertexLength + 2)
							{
								EdgeResult currentResult;
								currentResult.pieceId = task.piece;
								InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
								size_t definiteCount = std::count_if(task.str.begin() + 1, task.str.begin() + vertexLength + 1, DnaChar::IsDefinite);
								for (size_t pos = 1;; ++pos)
								{
									while (result.size() > 0 && FlushEdgeResults(result, writer, currentPiece));
									std::pair<uint64_t, uint64_t> bifId(INVALID_VERTEX, INVALID_VERTEX);
									char posPrev = task.str[pos - 1];
									char posExtend = task.str[pos + vertexLength];
									assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
									if (definiteCount == vertexLength)
									{
										bifId = bifStorage.GetId(task.str.begin() + pos, posVertexHash, negVertexHash);
										if (bifId.first != INVALID_VERTEX)
										{
											occurences++;
											currentResult.junction.push_back(JunctionPosition(task.seqId, task.start + pos - 1, bifId.first, bifId.second));
										}
									}

									if (((task.start == 0 && pos == 1) || (task.isFinal && pos == task.str.size() - vertexLength - 1)) && bifId.first == INVALID_VERTEX)
									{
										occurences++;
										if (definiteCount == vertexLength && DnaChar::IsSelfReverseCompliment(task.str.begin() + pos, vertexLength))
										{
											currentResult.junction.push_back(JunctionPosition(task.seqId, task.start + pos - 1, currentStubVertexId, currentStubVertexId));
											currentStubVertexId += 1;
										}
										else
										{
											currentResult.junction.push_back(JunctionPosition(task.seqId, task.start + pos - 1, currentStubVertexId, currentStubVertexId + 1));
											currentStubVertexId += 2;
										}

									}

									if (pos + edgeLength < task.str.size())
									{
										char negExtend = DnaChar::ReverseChar(posExtend);
										char posPrev = task.str[pos];
										char negPrev = DnaChar::ReverseChar(task.str[pos]);
										definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
										for (size_t i = 0; i < hashFunction.size(); i++)
										{
											posVertexHash[i]->update(posPrev, posExtend);
											negVertexHash[i]->reverse_update(negExtend, negPrev);
											assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
											assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
										}
									}
									else
									{
										break;
									}
								}

								result.push_back(currentResult);
							}
						}
					}

					while (result.size() > 0 && FlushEdgeResults(result, writer, currentPiece));
				}
				catch (std::runtime_error & e)
				{
					errorMutex.lock();
					error.reset(new std::runtime_error(e));
					errorMutex.unlock();
				}
			}

		private:
			size_t vertexLength;
			TaskQueue & taskQueue;
			const BifurcationStorage<CAPACITY> & bifStorage;
			JunctionPositionWriter & writer;
			std::atomic<uint64_t> & currentPiece;
			std::atomic<uint64_t> & occurences;
			std::atomic<uint64_t> & currentStubVertexId;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};

		static StrandComparisonResult DetermineStrand(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash)
		{
			for (size_t i = 0; i < posVertexHash.size(); i++)
			{
				uint64_t posHash = posVertexHash[i]->hashvalue;
				uint64_t negHash = negVertexHash[i]->hashvalue;
				if (posHash != negHash)
				{
					return posHash < negHash ? positiveLess : negativeLess;
				}
			}

			return tie;
		}

		static uint16_t EncodeSet(char prev, char next)
		{
			return 1 << (DnaChar::MakeUpChar(prev) * 4 + DnaChar::MakeUpChar(next));
		}

		static void DecodeSet(uint16_t set, char & prev, char & next)
		{
			prev = DnaChar::UnMakeUpChar(set >> 2);
			next = DnaChar::UnMakeUpChar(set & 12);
		}
		
		class InitialFilterFillerWorker
		{
		public:
			InitialFilterFillerWorker(uint64_t binSize,
				const std::vector<HashFunctionPtr> & hashFunction,
				ConcurrentBitVector & filter,
				size_t vertexLength,
				TaskQueue & taskQueue,
				std::atomic<uint32_t> * binCounter) : binSize(binSize), hashFunction(hashFunction), filter(filter),
				vertexLength(vertexLength), taskQueue(taskQueue), binCounter(binCounter)
			{

			}

			void operator()()
			{
				size_t edgeLength = vertexLength + 1;
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
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
							uint64_t hvalue;
							bool wasSet = true;
							char prevCh = task.str[pos];
							char nextCh = task.str[pos + edgeLength - 1];
							char revNextCh = DnaChar::ReverseChar(nextCh);
							uint64_t firstMinHash0 = min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
							uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
							uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);

							for (size_t i = 0; i < hashFunction.size(); i++)
							{
								if (posHash0 < negHash0 || (posHash0 == negHash0 && DnaChar::LessSelfReverseComplement(task.str.begin() + pos, vertexLength)))
								{
									hvalue = posVertexHash[i]->hash_extend(nextCh);
								}
								else
								{
									hvalue = negVertexHash[i]->hash_prepend(revNextCh);
								}

								if (!filter.GetBit(hvalue))
								{
									wasSet = false;
									filter.SetBitConcurrently(hvalue);
								}
							}

							for (size_t i = 0; i < hashFunction.size(); i++)
							{
								posVertexHash[i]->update(prevCh, nextCh);
								assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
								negVertexHash[i]->reverse_update(DnaChar::ReverseChar(nextCh), DnaChar::ReverseChar(prevCh));
								assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
							}

							uint64_t secondMinHash0 = min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
							if (!wasSet)
							{
								uint64_t value[] = { firstMinHash0, secondMinHash0 };
								for (uint64_t v : value)
								{
									uint64_t bin = v / binSize;
									if (binCounter[bin] < MAX_COUNTER)
									{
										binCounter[bin].fetch_add(1);
									}
								}
							}
						}
					}
				}
			}

		private:
			uint64_t binSize;
			const std::vector<HashFunctionPtr> & hashFunction;
			ConcurrentBitVector & filter;
			size_t vertexLength;
			TaskQueue & taskQueue;
			std::atomic<uint32_t> * binCounter;
		};
	
		class FilterFillerWorker
		{
		public:
			FilterFillerWorker(uint64_t low,
				uint64_t high,
				const std::vector<HashFunctionPtr> & hashFunction,
				ConcurrentBitVector & filter,
				size_t edgeLength,
				TaskQueue & taskQueue) : low(low), high(high), hashFunction(hashFunction), filter(filter), edgeLength(edgeLength), taskQueue(taskQueue)
			{

			}
			
			void operator()()
			{
				while (true)
				{
					Task task;
					const char DUMMY_CHAR = DnaChar::LITERAL[0];
					const char REV_DUMMY_CHAR = DnaChar::ReverseChar(DUMMY_CHAR);
					if (taskQueue.try_pop(task))
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
						if (task.str.size() >= vertexLength + 2)
						{
							std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
							std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
							InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
							size_t definiteCount = std::count_if(task.str.begin() + 1, task.str.begin() + vertexLength + 1, DnaChar::IsDefinite);
							for (size_t pos = 1;; ++pos)
							{
								char posPrev = task.str[pos - 1];
								char posExtend = task.str[pos + vertexLength];
								assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
								if (Within(min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue), low, high) && definiteCount == vertexLength)
								{
									StrandComparisonResult res = DetermineStrand(posVertexHash, negVertexHash);
									uint16_t posMask = EncodeSet(posPrev, posExtend);
									uint16_t negMask = EncodeSet(DnaChar::ReverseChar(posExtend), DnaChar::ReverseChar(posPrev));
									if (posPrev == 'N' || posExtend == 'N' || (tie && DnaChar::IsSelfReverseCompliment(task.str.begin() + pos, vertexLength)))
									{
										posMask = negMask = 1 | 2;
									}

									if (res == positiveLess || res == tie)
									{
										for (size_t i = 0; i < posVertexHash.size(); i++)
										{
											filter.OrElementCouncerrently(posVertexHash[i]->hashvalue, posMask);
										}
									}

									if (res == negativeLess || res == tie)
									{
										for (size_t i = 0; i < negVertexHash.size(); i++)
										{
											filter.OrElementCouncerrently(negVertexHash[i]->hashvalue, negMask);
										}
									}

								}

								if (pos + edgeLength < task.str.size())
								{
									char negExtend = DnaChar::ReverseChar(posExtend);
									char posPrev = task.str[pos];
									char negPrev = DnaChar::ReverseChar(task.str[pos]);
									definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
									for (size_t i = 0; i < hashFunction.size(); i++)
									{
										posVertexHash[i]->update(posPrev, posExtend);
										negVertexHash[i]->reverse_update(negExtend, negPrev);
										assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
										assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
									}
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

		private:
			uint64_t low;
			uint64_t high;
			const std::vector<HashFunctionPtr> & hashFunction;
			ConcurrentBitVector & filter;
			size_t edgeLength;
			TaskQueue & taskQueue;
		};

		static void DistributeTasks(const std::vector<std::string> & fileName,
			size_t overlapSize,
			std::vector<TaskQueuePtr> & taskQueue,
			std::unique_ptr<std::runtime_error> & error,
			tbb::mutex & errorMutex,
			std::ostream & logFile)
		{
			size_t record = 0;
			size_t nowQueue = 0;
			uint32_t pieceCount = 0;
			logFile << "Starting a new stage" << std::endl;
			for (size_t file = 0; file < fileName.size(); file++)
			{
				logFile << "Reading " << fileName[file] << std::endl;
				const std::string & nowFileName = fileName[file];
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
				{
					{
						errorMutex.lock();
						if (error != 0)
						{
							throw *error;
						}

						errorMutex.unlock();
					}

					std::stringstream ss;
					logFile << "Processing sequence " << parser.GetCurrentHeader() << " " << ss.str() << std::endl;

					char ch;
					uint64_t prev = 0;
					uint64_t start = 0;
					std::string buf = "N";
					bool over = false;
					do
					{
						over = !parser.GetChar(ch);
						if (!over)
						{
							start++;
							buf.push_back(DnaChar::IsDefinite(ch) ? ch : 'N');
						}

						if (buf.size() >= overlapSize && (buf.size() == Task::TASK_SIZE || over))
						{
							for (bool found = false; !found; nowQueue = nowQueue + 1 < taskQueue.size() ? nowQueue + 1 : 0)
							{
								TaskQueuePtr & q = taskQueue[nowQueue];
								if (q->capacity() - q->size() > 0)
								{
									std::string overlap;
									if (!over)
									{
										overlap.assign(buf.end() - overlapSize, buf.end());
									}
									else
									{
										buf.push_back('N');
									}

									q->push(Task(record, prev, pieceCount++, over, std::move(buf)));
									logFile << "Passed chunk " << prev << " to worker " << nowQueue << std::endl;
									prev = start - overlapSize + 1;
									buf.swap(overlap);
									found = true;
								}
							}

						}

					} while (!over);
				}
			}

			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				TaskQueuePtr & q = taskQueue[nowQueue];
				while (!taskQueue[i]->try_push(Task(0, Task::GAME_OVER, 0, true, std::string())))
				{
					
				}
			}
		}

		uint64_t TrueBifurcations(const OccurenceSet & occurenceSet, std::ofstream & out, size_t vertexSize, size_t & falsePositives) const
		{
			uint64_t truePositives = falsePositives = 0;
			std::vector<Occurence> store;
			for (auto it = occurenceSet.begin(); it != occurenceSet.end();)
			{
				Occurence base = *it;
				size_t inUnknownCount = 0;
				size_t outUnknownCount = 0;
				bool bifurcation = it->IsBifurcation();
				bool selfReverseCompliment = base.IsSelfReverseCompliment(vertexSize);

				auto jt = it;
				for (; jt != occurenceSet.end(); ++jt)
				{
					Occurence next = *jt;
					if (!base.EqualBase(next))
					{
						break;
					}

					inUnknownCount += DnaChar::IsDefinite(next.Prev()) ? 0 : 1;
					outUnknownCount += DnaChar::IsDefinite(next.Next()) ? 0 : 1;
					if (!bifurcation)
					{
						bifurcation = jt->IsBifurcation() || base.Prev() != next.Prev() || base.Next() != next.Next();
						if (selfReverseCompliment)
						{
							inUnknownCount += DnaChar::IsDefinite(next.Next()) ? 0 : 1;
							outUnknownCount += DnaChar::IsDefinite(next.Prev()) ? 0 : 1;
							bifurcation = bifurcation ||
								base.Prev() != DnaChar::ReverseChar(next.Next()) ||
								base.Next() != DnaChar::ReverseChar(next.Prev());
						}
					}
				}

				if (bifurcation || inUnknownCount > 1 || outUnknownCount > 1)
				{
					++truePositives;
					it->GetBase().WriteToFile(out);
					if (!out)
					{
						throw StreamFastaParser::Exception("Can't write to a temporary file");
					}
				}
				else
				{
					falsePositives++;
				}

				it = jt;
			}

			return truePositives;
		}

		size_t vertexSize_;
	};
}

#endif	
