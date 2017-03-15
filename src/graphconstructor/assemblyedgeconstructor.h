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
	class AssemblyEdgeConstructorImpl: public AssemblyEdgeConstructor
	{
	public:
		AssemblyEdgeConstructorImpl(const std::vector<std::string> & fileName,
			ConcurrentBitVector & bloomFilter,
			const VertexRollingHashSeed & seed,
			const std::string & edgeOutputFile,
			const VertexEnumeratorImpl<CAPACITY> & vertexEnumerator) : 
				seed_(seed),
				bloomFilter_(bloomFilter),
				vertexEnumerator_(vertexEnumerator)
		{

		}
		
	private:
		VertexRollingHashSeed seed_;
		ConcurrentBitVector & bloomFilter_;
		const VertexEnumeratorImpl<CAPACITY> & vertexEnumerator_;
		DISALLOW_COPY_AND_ASSIGN(AssemblyEdgeConstructorImpl<CAPACITY>);
	};


}

#endif