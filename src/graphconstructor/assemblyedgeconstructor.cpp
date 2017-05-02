#include "assemblyedgeconstructor.h"

namespace TwoPaCo
{
	namespace
	{
		template<size_t CAPACITY>
		std::unique_ptr<AssemblyEdgeConstructor> CreateAssemblyEdgeConstructorImpl(const std::vector<std::string> & inputFileName, const std::string & marksFileName, const VertexEnumerator & vertexEnumerator)
		{
			size_t neededCapacity = CalculateNeededCapacity(vertexEnumerator.GetHashSeed().VertexLength());
			if (CAPACITY == neededCapacity)
			{
				auto ptr = new AssemblyEdgeConstructorImpl<CAPACITY>(inputFileName, marksFileName, vertexEnumerator);
				return std::unique_ptr<AssemblyEdgeConstructor>(ptr);
			}

			return CreateAssemblyEdgeConstructorImpl<CAPACITY + 1>(inputFileName, marksFileName, vertexEnumerator);
		}

		template<>
		std::unique_ptr<AssemblyEdgeConstructor> CreateAssemblyEdgeConstructorImpl<MAX_CAPACITY>(const std::vector<std::string> & inputFileName, const std::string & marksFileName, const VertexEnumerator & vertexEnumerator)
		{
			return 0;
		}
	}

	std::unique_ptr<AssemblyEdgeConstructor> CreateAssemblyEdgeConstructor(const std::vector<std::string> & inputFileName, const std::string & marksFileName, const VertexEnumerator & vertexEnumerator)
	{
		return CreateAssemblyEdgeConstructorImpl<1>(inputFileName, marksFileName, vertexEnumerator);
	}
}