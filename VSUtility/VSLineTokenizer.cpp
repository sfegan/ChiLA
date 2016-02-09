//-*-mode:c++; mode:font-lock;-*-

/*! \file VSLineTokenizer.cpp

  Class to break lines into tokens

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       08/15/2005
  \note
*/

#include <cctype>
#include <vsassert>

#include <VSLineTokenizer.hpp>

using namespace VERITAS;

// ----------------------------------------------------------------------------
// VSToken
// ----------------------------------------------------------------------------

std::string VSToken::lower() const
{
  std::string temp;
  for(std::string::const_iterator ichar = fToken.begin();
      ichar!=fToken.end(); ichar++)temp += tolower(*ichar);
  return temp;
}

std::string VSToken::escaped() const
{
  std::string temp;
  temp += "'";
  for(std::string::const_iterator ichar = fToken.begin();
      ichar!=fToken.end(); ichar++)
    if(*ichar == '\'')temp += "'\\''";
    else temp += *ichar;
  temp += "'";  
  return temp;
}

// ----------------------------------------------------------------------------
// VSCharSupplier
// ----------------------------------------------------------------------------

VSCharSupplier::~VSCharSupplier()
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
// VSCSString
// ----------------------------------------------------------------------------

VSCSString::~VSCSString()
{
  // nothing to see here
}

bool VSCSString::hasChar()
{
  return fIString != fString.end();
}

char VSCSString::getChar()
{
  if(fIString != fString.end())return *fIString;
  else return 0;
}

void VSCSString::nextChar()
{
  if(fIString != fString.end())fIString++;
}

std::string VSCSString::remainder()
{
  if(fIString != fString.end())return fString.substr(fIString-fString.begin());
  else return std::string();
}

// ----------------------------------------------------------------------------
// VSCSStream
// ----------------------------------------------------------------------------

VSCSStream::~VSCSStream()
{
  // nothing to see here
}

bool VSCSStream::hasChar()
{
  if(fStream.eof())return false;
  else 
    {
      char c = fStream.peek();
      return c!=EOF;
    }
}

char VSCSStream::getChar()
{
  if(!fStream.eof())
    {
      char c = fStream.peek();
      return c;
    }
  else return EOF;
}

void VSCSStream::nextChar()
{
  if(!fStream.eof())
    {
      char c;
      fStream.get(c);
    }
}

// ----------------------------------------------------------------------------
// VSLineTokenize
// ----------------------------------------------------------------------------

void VSLineTokenizer::tokenize(VSCharSupplier* supplier, VSTokenList& tokens)
{
  tokens.clear();

  ParserState state = PS_WHITESPACE;
  bool literal_next = false;
  bool terminate_parser = false;

  std::string token;
  while((!terminate_parser)&&(supplier->hasChar()))
    {
      char c = supplier->getChar();

      if(literal_next)
	{
	  literal_next=false;
	  switch(c)
	    {
	    case '\'':
	    case '\"':
	    case '\\':
	    case ' ':
	    case '\t':
	    case '\n':
	      token += c;
	      supplier->nextChar();
	      break;

	    case 'n':
	      token += '\n';
	      supplier->nextChar();
	      break;

	    default:
	      if(fStrict)
		{
		  vsassert(0);
		}
	      else
		{
		  token += c;
		  supplier->nextChar();
		}
	    }
	}
      else if(c=='\n')
	{
	  if((fStrict)&&((state==PS_QUOTED_S)||(state==PS_QUOTED_D)))
	    {
	      vsassert(0);
	    }
	  terminate_parser=true;
	  supplier->nextChar();
	}
      else
	{
	  switch(state)
	    {
	    case PS_WHITESPACE:
	      switch(c)
		{
		case ' ':
		case '\t':
		  supplier->nextChar();
		  break;
		case '#':
		  state=PS_COMMENT;
		  break;
		default:
		  state=PS_NON_QUOTED;
		  break;
		}
	      break;

	    case PS_COMMENT:
	      supplier->nextChar();
	      break;

	    case PS_NON_QUOTED:
	      switch(c)
		{
		case ' ':
		case '\t':
		  tokens.push_back(token);
		  token.clear();
		  state=PS_WHITESPACE;
		  break;
		case '#':
		  tokens.push_back(token);
		  token.clear();
		  state=PS_COMMENT;
		  break;		  
		case '\\':
		  literal_next = true;
		  supplier->nextChar();
		  break;
		case '\'':
		  state=PS_QUOTED_S;
		  supplier->nextChar();
		  break;
		case '"':
		  state=PS_QUOTED_D;
		  supplier->nextChar();
		  break;
		default:
		  token += c;
		  supplier->nextChar();
		  break;
		}
	      break;

	    case PS_QUOTED_S:
	      switch(c)
		{
		case '\'':
		  state = PS_NON_QUOTED;
		  supplier->nextChar();
		  break;
		default:
		  token += c;
		  supplier->nextChar();
		  break;
		};
	      break;
	      
	    case PS_QUOTED_D:
	      switch(c)
		{
		case '\"':
		  state = PS_NON_QUOTED;
		  supplier->nextChar();
		  break;
		case '\\':
		  literal_next = true;
		  supplier->nextChar();
		  break;
		default:
		  token += c;
		  supplier->nextChar();
		  break;
		};
	      break;
	    }
	}
    }

  if((literal_next)&&(fStrict))
    {
      vsassert(0);
    }
      
  switch(state)
    {
    case PS_WHITESPACE:
    case PS_COMMENT:
      break;

    case PS_NON_QUOTED:
      tokens.push_back(token);
      token.clear();
      break;

    case PS_QUOTED_S:
    case PS_QUOTED_D:
      if(fStrict)
	{
	  vsassert(0);
	}
      else
	{
	  tokens.push_back(token);
	  token.clear();
	}
      break;
    }
}

void VSLineTokenizer::tokenize(std::string& string, VSTokenList& tokens)
{
  VSCSString supplier(string);
  tokenize(&supplier,tokens);
  string=supplier.remainder();
}

void VSLineTokenizer::tokenize(std::istream& stream, VSTokenList& tokens)
{
  VSCSStream supplier(stream);
  tokenize(&supplier,tokens);
}

#ifdef TEST_MAIN_2

int main()
{
  VSLineTokenizer tokenizer;
  std::string line;
  while(std::getline(std::cin,line))
    {
      VSTokenList tokens;
      tokenizer.tokenize(line,tokens);

      for(unsigned i=0;i<tokens.size();i++)
	std::cout << i+1 << ' ' <<  tokens[i].string() << std::endl;
    }
}

#endif

#ifdef TEST_MAIN

int main()
{
  VSLineTokenizer tokenizer(true);
  while(!std::cin.eof())
    {
      VSTokenList tokens;
      tokenizer.tokenize(std::cin,tokens);

      for(unsigned i=0;i<tokens.size();i++)
	std::cout << i+1 << ' ' <<  tokens[i].string() << std::endl;
    }
}

#endif
