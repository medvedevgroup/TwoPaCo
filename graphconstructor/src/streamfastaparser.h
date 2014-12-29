#ifndef _STREAM_FASTA_PARSER_H_
#define _STREAM_FASTA_PARSER_H_

#include <string>
#include <fstream>
#include <exception>

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
		
		bool ReadRecord();
		bool GetChar(char & ch);
		std::string GetErrorMessage() const;
		std::string GetCurrentHeader() const;		
		StreamFastaParser(const std::string & fileName);
	private:
		static const std::string VALID_CHARS;

		std::ifstream stream_;
		std::string errorMessage_;
		std::string currentHeader_;		
	};
}

#endif