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

#include <junctionapi/junctionapi.h>

#include "streamfastaparser.h"
#include "bifurcationstorage.h"
#include "candidateoccurence.h"

namespace TwoPaCo
{
	class VertexEnumerator
	{
	public:		
		virtual const VertexRollingHashSeed & GetHashSeed() const = 0;
		virtual std::unique_ptr<ConcurrentBitVector> ReloadBloomFilter() const = 0;
		virtual void ReloadJunctionCandidateMask(std::vector<std::vector<bool> > & ret) const = 0;
		virtual bool GetEdges(std::string::const_iterator it, const VertexRollingHash & hash, std::string & inEdges, std::string & outEdges) const = 0;
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
		const std::string & outFileName,
		std::ostream & logStream);

	template<size_t CAPACITY>
	class VertexEnumeratorImpl : public VertexEnumerator
	{
	private:
		typedef CompressedString<CAPACITY> DnaString;
		typedef CandidateOccurence<CAPACITY> Occurence;
		
		std::string filterDumpFile_;
		VertexRollingHashSeed hashFunctionSeed_;		
		static const size_t BUF_SIZE = 1 << 24;
		
		class FilterFillerWorker;

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
		
		typedef tbb::concurrent_unordered_set<Occurence, OccurenceHash, OccurenceEquality> OccurenceSet;		
		OccurenceSet occurenceSet_;
		std::string temporaryDirectory_;
		std::vector<std::vector<size_t> >  taskFragmentSize_;

	public:

		~VertexEnumeratorImpl<CAPACITY>()
		{
			std::remove(filterDumpFile_.c_str());
			for (size_t i = 0; i < taskFragmentSize_.size(); i++)
			{
				size_t pos = 0;
				for (size_t fragmentSize : taskFragmentSize_[i])
				{
					std::string fileName = CandidateMaskFileName(temporaryDirectory_, i, pos, 0);
					std::remove(fileName.c_str());
					pos += fragmentSize;
				}
			}
		}

		void ReloadJunctionCandidateMask(std::vector<std::vector<bool> > & ret) const
		{			
			for (size_t i = 0; i < taskFragmentSize_.size(); i++)
			{
				size_t pos = 0;
				ret.push_back(std::vector<bool>());
				for (size_t fragmentSize : taskFragmentSize_[i])
				{
					ConcurrentBitVector v(fragmentSize);
					std::string fileName = CandidateMaskFileName(temporaryDirectory_, i, pos, 0);
					v.ReadFromFile(fileName, false);
					for (size_t k = 0; k < fragmentSize; k++)
					{
						ret.back().push_back(v.GetBit(k));
					}

					pos += fragmentSize;
				}
			}
		}

		bool GetEdges(std::string::const_iterator pos, const VertexRollingHash & hash, std::string & inEdges, std::string & outEdges) const
		{
			Occurence query;
			inEdges.clear();
			outEdges.clear();
			bool reverse = !query.Init(hash, pos, hashFunctionSeed_.VertexLength(), 'N', 'N');			
			auto it = occurenceSet_.find(query);
			if (it != occurenceSet_.end())
			{
				for (char ch : DnaChar::LITERAL)
				{
					if (it->GetOutgoingEdge(ch))
					{
						if (reverse)
						{
							inEdges.push_back(DnaChar::ReverseChar(ch));
						}
						else
						{
							outEdges.push_back(ch);
						}
					}

					if (it->GetIngoingEdge(ch))
					{
						if (reverse)
						{
							outEdges.push_back(DnaChar::ReverseChar(ch));
						}
						else
						{
							inEdges.push_back(ch);
						}
					}
				}				

				return true;
			}

			return false;
		}	

		const VertexRollingHashSeed & GetHashSeed() const
		{
			return hashFunctionSeed_;
		}

		std::unique_ptr<ConcurrentBitVector> ReloadBloomFilter() const
		{
			uint64_t realSize = uint64_t(1) << hashFunctionSeed_.BitsNumber();
			std::unique_ptr<ConcurrentBitVector> ret(new ConcurrentBitVector(realSize));
			ret->ReadFromFile(filterDumpFile_, false);
			return ret;
		}

		VertexEnumeratorImpl(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			const std::string & tmpDirName,
			const std::string & outFileNamePrefix,
			std::ostream & logStream) :
			vertexSize_(vertexLength),
			hashFunctionSeed_(hashFunctions, vertexLength, filterSize),
			filterDumpFile_(tmpDirName + "/filter.bin"),
			occurenceSet_(1 << 24)
		{
			temporaryDirectory_ = tmpDirName;
			uint64_t realSize = uint64_t(1) << filterSize;
			logStream << "Threads = " << threads << std::endl;
			logStream << "Vertex length = " << vertexLength << std::endl;
			logStream << "Hash functions = " << hashFunctions << std::endl;
			logStream << "Filter size = " << realSize << std::endl;
			logStream << "Capacity = " << CAPACITY << std::endl;
			logStream << "Files: " << std::endl;
			for (const std::string & fn : fileName)
			{
				logStream << fn << std::endl;
			}
#ifdef LOGGING
			std::ofstream logFile((tmpDirName + "/log.txt").c_str());
			if (!logFile)
			{
				throw StreamFastaParser::Exception("Can't open the log file");
			}
#else
			std::ostream & logFile = std::cerr;
#endif

			tbb::mutex errorMutex;
			std::unique_ptr<std::runtime_error> error;
			
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
				logStream << "Splitting the input kmers set..." << std::endl;
				std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
				binCounter = new std::atomic<uint32_t>[BINS_COUNT];
				std::fill(binCounter, binCounter + BINS_COUNT, 0);
				ConcurrentBitVector bitVector(realSize);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					InitialFilterFillerWorker worker(BIN_SIZE,
						hashFunctionSeed_,
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


			logStream << std::string(80, '-') << std::endl;
			uint64_t low = 0;
			uint64_t high = realSize;
			uint64_t lowBoundary = 0;
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
					logStream << "Round " << round << ", " << low << ":" << high << std::endl;
					logStream << "Pass\tFilling\tFiltering" << std::endl << "1\t";
					{
						std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							FilterFillerWorker worker(low,
								high,
								std::cref(hashFunctionSeed_),
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

					bitVector.WriteToFile(filterDumpFile_);
					logStream << time(0) - mark << "\t";
					mark = time(0);
					{
						std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							CandidateCheckingWorker worker(std::make_pair(low, high),
								hashFunctionSeed_,
								bitVector,
								vertexLength,
								*taskQueue[i],
								tmpDirName,
								marks,
								round,
								error,
								errorMutex);

							workerThread[i].reset(new tbb::tbb_thread(worker));
						}

						DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile, taskFragmentSize_, true);
						for (size_t i = 0; i < taskQueue.size(); i++)
						{
							workerThread[i]->join();
						}

						if (error != 0)
						{
							throw *error;
						}
					}

					logStream << time(0) - mark << "\t" << std::endl;
				}

				mark = time(0);
				tbb::spin_rw_mutex mutex;
				logStream << "2\t";
				{
					std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
					for (size_t i = 0; i < workerThread.size(); i++)
					{
						CandidateFinalFilteringWorker worker(hashFunctionSeed_,
							vertexLength,
							*taskQueue[i],
							occurenceSet_,
							mutex,
							tmpDirName,
							round,
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

					logStream << time(0) - mark << "\t";
				}

				logStream << "Hash table size = " << occurenceSet_.size() << std::endl;
				logStream << "Candidate marks count = " << marks << std::endl;
				logStream << std::string(80, '-') << std::endl;
				low = high + 1;
			}

			if (rounds > 1)
			{
				delete[] binCounter;
			}

			mark = time(0);						
			logStream << std::string(80, '-') << std::endl;
		}

	private:

		static const size_t QUEUE_CAPACITY = 16;
		static const uint64_t BINS_COUNT = 1 << 24;		

		static bool Within(uint64_t hvalue, uint64_t low, uint64_t high)
		{
			return hvalue >= low && hvalue <= high;
		}

		static std::string CandidateMaskFileName(const std::string & directory, size_t sequence, size_t pos, size_t round)
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

		class InitialFilterFillerWorker
		{
		public:
			InitialFilterFillerWorker(uint64_t binSize,
				const VertexRollingHashSeed & hashFunction,
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

						std::vector<uint64_t> hashValue;
						size_t vertexLength = edgeLength - 1;						
						VertexRollingHash hash(hashFunction, task.str.begin(), hashFunction.HashFunctionsNumber());
						for (size_t pos = 0; pos + edgeLength - 1 < task.str.size(); ++pos)
						{							
							hashValue.clear();
							bool wasSet = true;
							char prevCh = task.str[pos];
							char nextCh = task.str[pos + edgeLength - 1];
							uint64_t startVertexHash = hash.GetVertexHash();							
							GetOutgoingEdgeHash(hash, nextCh, hashValue);
							for (auto hvalue : hashValue)
							{
								if (!filter.GetBit(hvalue))
								{
									wasSet = false;
									filter.SetBitConcurrently(hvalue);
								}
							}

							hash.Update(prevCh, nextCh);
							assert(hash.Assert(task.str.begin() + pos + 1));

							uint64_t endVertexHash = hash.GetVertexHash();
							if (!wasSet)
							{
								uint64_t value[] = { startVertexHash, endVertexHash };
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
			const VertexRollingHashSeed & hashFunction;
			ConcurrentBitVector & filter;
			size_t vertexLength;
			TaskQueue & taskQueue;
			std::atomic<uint32_t> * binCounter;
		};


		class CandidateCheckingWorker
		{
		public:
			CandidateCheckingWorker(std::pair<uint64_t, uint64_t> bound,
				const VertexRollingHashSeed & hashFunction,
				const ConcurrentBitVector & bitVector,
				size_t vertexLength,
				TaskQueue & taskQueue,
				const std::string & tmpDirectory,
				std::atomic<uint64_t> & marksCount,
				size_t round,
				std::unique_ptr<std::runtime_error> & error,				
				tbb::mutex & errorMutex) : bound(bound), hashFunction(hashFunction), bitVector(bitVector), vertexLength(vertexLength), taskQueue(taskQueue),
				tmpDirectory(tmpDirectory), marksCount(marksCount), error(error), errorMutex(errorMutex), round(round)
			{

			}

			void operator()()
			{
				uint64_t low = bound.first;
				uint64_t high = bound.second;
				std::vector<uint64_t> temp;
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
							VertexRollingHash hash(hashFunction, task.str.begin(), hashFunction.HashFunctionsNumber());
							size_t definiteCount = std::count_if(task.str.begin(), task.str.begin() + vertexLength, DnaChar::IsDefinite);
							for (size_t pos = 0;; ++pos)
							{								
								assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
								if (Within(hash.GetVertexHash(), low, high) && definiteCount == vertexLength)
								{
									size_t inCount = 0;
									size_t outCount = 0;
									for (int i = 0; i < DnaChar::LITERAL.size() && inCount < 2 && outCount < 2; i++)
									{
										char nextCh = DnaChar::LITERAL[i];
										if ((pos > 0 && nextCh == task.str[pos - 1]) || IsIngoingEdgeInBloomFilter(hash, bitVector, nextCh))
										{
											++inCount;
										}

										if ((pos + vertexLength < task.str.size() && task.str[pos + vertexLength] == nextCh) || IsOutgoingEdgeInBloomFilter(hash, bitVector, nextCh))
										{
											++outCount;
										}
									}

									if (inCount > 1 || outCount > 1 || inCount == 0 || outCount == 0)
									{
										++marksCount;
										candidateMask.SetBitConcurrently(pos);
									}
								}

								if (pos + vertexLength < task.str.size())
								{
									char posPrev = task.str[pos];
									definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
									hash.Update(posPrev, task.str[pos + vertexLength]);
									assert(hash.Assert(task.str.begin() + pos + 1));
								}
								else
								{
									break;
								}
							}

							try
							{
								candidateMask.WriteToFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start, round));
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
			const VertexRollingHashSeed & hashFunction;
			const ConcurrentBitVector & bitVector;
			size_t vertexLength;
			TaskQueue & taskQueue;
			const std::string & tmpDirectory;
			std::atomic<uint64_t> & marksCount;
			size_t round;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};



		class CandidateFinalFilteringWorker
		{
		public:
			CandidateFinalFilteringWorker(const VertexRollingHashSeed & hashFunction,
				size_t vertexLength,
				TaskQueue & taskQueue,
				OccurenceSet & occurenceSet,
				tbb::spin_rw_mutex & mutex,
				const std::string & tmpDirectory,
				size_t round,
				std::unique_ptr<std::runtime_error> & error,
				tbb::mutex & errorMutex) : hashFunction(hashFunction), vertexLength(vertexLength), taskQueue(taskQueue), occurenceSet(occurenceSet),
				mutex(mutex), tmpDirectory(tmpDirectory), round(round), error(error), errorMutex(errorMutex)
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

						size_t edgeLength = vertexLength + 1;
						if (task.str.size() >= vertexLength + 2)
						{
							VertexRollingHash hash(hashFunction, task.str.begin(), 1);
							{
								try
								{
									candidateMask.ReadFromFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start, round), false);
								}
								catch (std::runtime_error & err)
								{
									ReportError(errorMutex, error, err.what());
								}
							}

							for (size_t pos = 0;; ++pos)
							{																
								if (candidateMask.GetBit(pos))
								{
									Occurence now;									
									bool direct = now.Init(hash,
										task.str.begin() + pos,
										vertexLength,
										pos > 0 ? task.str[pos - 1] : 'N',
										pos + vertexLength < task.str.size() ? task.str[pos + vertexLength] : 'N');
									auto ret = occurenceSet.insert(now);
									typename OccurenceSet::iterator it = ret.first;
									if (!ret.second)
									{										
										it->Merge(now);
									}
								}

								if (pos + vertexLength < task.str.size())
								{
									char posPrev = task.str[pos];
									hash.Update(posPrev, task.str[pos + vertexLength]);
									assert(hash.Assert(task.str.begin() + pos + 1));
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
			const VertexRollingHashSeed & hashFunction;
			size_t vertexLength;
			TaskQueue & taskQueue;
			OccurenceSet & occurenceSet;
			tbb::spin_rw_mutex & mutex;
			const std::string & tmpDirectory;
			size_t round;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};
		
		class FilterFillerWorker
		{
		public:
			FilterFillerWorker(uint64_t low,
				uint64_t high,
				const VertexRollingHashSeed & hashFunction,
				ConcurrentBitVector & filter,
				size_t edgeLength,
				TaskQueue & taskQueue) : low(low), high(high), hashFunction(hashFunction), filter(filter), edgeLength(edgeLength), taskQueue(taskQueue)
			{

			}

			void operator()()
			{
				std::vector<uint64_t> setup;
				std::vector<uint64_t> hashValue;
				const char DUMMY_CHAR = DnaChar::LITERAL[0];
				const char REV_DUMMY_CHAR = DnaChar::ReverseChar(DUMMY_CHAR);
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

						uint64_t fistMinHash0;
						uint64_t secondMinHash0;
						size_t vertexLength = edgeLength - 1;
						size_t definiteCount = std::count_if(task.str.begin(), task.str.begin() + vertexLength, DnaChar::IsDefinite);
						VertexRollingHash hash(hashFunction, task.str.begin(), hashFunction.HashFunctionsNumber());						
						for (size_t pos = 0;; ++pos)
						{
							hashValue.clear();
							char prevCh = task.str[pos];
							char nextCh = task.str[pos + edgeLength - 1];
							assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
							if (definiteCount == vertexLength)
							{
								fistMinHash0 = hash.GetVertexHash();
								if (DnaChar::IsDefinite(nextCh))
								{
									GetOutgoingEdgeHash(hash, nextCh, hashValue);
								}
								else
								{
									GetOutgoingEdgeHash(hash, DUMMY_CHAR, hashValue);
									GetOutgoingEdgeHash(hash, REV_DUMMY_CHAR, hashValue);
								}

								if (pos > 0 && !DnaChar::IsDefinite(task.str[pos - 1]))
								{
									GetIngoingEdgeHash(hash, DUMMY_CHAR, hashValue);
									GetIngoingEdgeHash(hash, REV_DUMMY_CHAR, hashValue);
								}
							}

							hash.Update(prevCh, nextCh);
							assert(hash.Assert(task.str.begin() + pos + 1));
							if (definiteCount == vertexLength)
							{
								secondMinHash0 = hash.GetVertexHash();
								if (Within(fistMinHash0, low, high) || Within(secondMinHash0, low, high))
								{
									for (uint64_t value : hashValue)
									{
										setup.push_back(value);
									}
								}
							}

							if (pos + vertexLength < task.str.size() - 1)
							{
								definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(prevCh) ? 1 : 0);
							}
							else
							{
								break;
							}
						}
					}

					for (uint64_t hashValue : setup)
					{
						if (!filter.GetBit(hashValue))
						{
							filter.SetBitConcurrently(hashValue);
						}
					}

					setup.clear();
				}
			}

		private:
			uint64_t low;
			uint64_t high;
			const VertexRollingHashSeed & hashFunction;
			ConcurrentBitVector & filter;
			size_t edgeLength;
			TaskQueue & taskQueue;
		};


		static void DistributeTasks(const std::vector<std::string> & fileName,
			size_t overlapSize,
			std::vector<TaskQueuePtr> & taskQueue,
			std::unique_ptr<std::runtime_error> & error,
			tbb::mutex & errorMutex,
			std::ostream & logFile,
			std::vector<std::vector<size_t> > & taskFragmentSize = std::vector<std::vector<size_t> >(),
			bool writeTaskFragmentSize = false)
		{
			size_t record = 0;
			size_t nowQueue = 0;
			uint32_t pieceCount = 0;
#ifdef LOGGING
			logFile << "Starting a new stage" << std::endl;
#endif
			for (size_t file = 0; file < fileName.size(); file++)
			{
#ifdef LOGGING
				logFile << "Reading " << fileName[file] << std::endl;
#endif
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
#ifdef LOGGING
					logFile << "Processing sequence " << parser.GetCurrentHeader() << " " << ss.str() << std::endl;
#endif
					char ch;
					uint64_t prev = 0;
					uint64_t start = 0;
					std::string buf = "";
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

									if (writeTaskFragmentSize)
									{
										while(record >= taskFragmentSize.size())
										{
											taskFragmentSize.push_back(std::vector<size_t>());
										}

										taskFragmentSize[record].push_back(buf.size());
									}

									q->push(Task(record, prev, pieceCount++, over, std::move(buf)));
#ifdef LOGGING
									logFile << "Passed chunk " << prev << " to worker " << nowQueue << std::endl;
#endif
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
		

		size_t vertexSize_;
		DISALLOW_COPY_AND_ASSIGN(VertexEnumeratorImpl<CAPACITY>);
	};	
}

#endif	
