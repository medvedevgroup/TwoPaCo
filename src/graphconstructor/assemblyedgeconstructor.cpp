#include "assemblyedgeconstructor.h"

namespace TwoPaCo
{
	namespace
	{
		template<size_t CAPACITY>
		std::unique_ptr<AssemblyEdgeConstructor> CreateAssemblyEdgeConstructorImpl(const std::vector<std::string> & fileName,
			ConcurrentBitVector & bloomFilter,
			const VertexRollingHashSeed & seed,
			const std::string & edgeOutputFile,
			const VertexEnumerator & vertexEnumerator)
		{
			size_t neededCapacity = CalculateNeededCapacity(seed.VertexLength());
			if (CAPACITY == neededCapacity)
			{
				auto ptr = new AssemblyEdgeConstructorImpl<CAPACITY>(fileName,
					bloomFilter,
					seed,
					edgeOutputFile,
					dynamic_cast<const VertexEnumeratorImpl<CAPACITY> &>(vertexEnumerator));
				return std::unique_ptr<AssemblyEdgeConstructor>(ptr);
			}

			return CreateAssemblyEdgeConstructorImpl<CAPACITY + 1>(fileName,
				bloomFilter,
				seed,
				edgeOutputFile,
				vertexEnumerator);
		}

		template<>
		std::unique_ptr<AssemblyEdgeConstructor> CreateAssemblyEdgeConstructorImpl<MAX_CAPACITY>(const std::vector<std::string> & fileName,
			ConcurrentBitVector & bloomFilter,
			const VertexRollingHashSeed & seed,
			const std::string & edgeOutputFile,
			const VertexEnumerator & vertexEnumerator)
		{
			return 0;
		}
	}

	std::unique_ptr<AssemblyEdgeConstructor> CreateAssemblyEdgeConstructor(const std::vector<std::string> & fileName,
		ConcurrentBitVector & bloomFilter,
		const VertexRollingHashSeed & seed,
		const std::string & edgeOutputFile,
		const VertexEnumerator & vertexEnumerator)
	{
		return CreateAssemblyEdgeConstructorImpl<1>(fileName,
			bloomFilter,			
			seed,
			edgeOutputFile,
			vertexEnumerator);
	}
}