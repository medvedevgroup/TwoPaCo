#ifndef _STREAM_FASTA_PARSER_H_
#define _STREAM_FASTA_PARSER_H_

#include <fstream>
#include <stdexcept>

#include <boost/thread.hpp>

#include "dnastring.h"

namespace Sibelia
{
	class StreamFastaParser
	{
	public:
		class Exception: public std::runtime_error
		{
		public:
			Exception(const std::string & msg);
		};

		static const char UNKNOWN_CHAR = 4;
		
		bool ReadRecord();
		~StreamFastaParser();
		bool GetChar(char & ch, bool makeUp = true);		
		static char MakeUp(char ch);
		static char Reverse(char ch);
		static char UnMakeUp(char ch);
		std::string GetErrorMessage() const;
		std::string GetCurrentHeader() const;
		StreamFastaParser(const std::string & fileName);
	private:		
		static const std::string VALID_CHARS;
		static const size_t BUF_SIZE = 1 << 20;

		bool Peek(char & ch);
		bool GetCh(char & ch);		

		std::ifstream stream_;
		std::string errorMessage_;
		std::string currentHeader_;
		char * buffer_;
		size_t bufferSize_;
		size_t bufferPos_;
		std::vector<char> isValid_;
	};
}

#endif