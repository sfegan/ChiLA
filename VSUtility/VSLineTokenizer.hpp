//-*-mode:c++; mode:font-lock;-*-

/*! \file VSLineTokenizer.hpp

  Class to break lines into tokens

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n    

  \version    0.1
  \date       08/15/2005
  \note
*/

#ifndef VSLINETOKENIZER_HPP
#define VSLINETOKENIZER_HPP

#include <istream>
#include <string>
#include <vector>
#include <memory>
#include <stdint.h>

#include <VSDataConverter.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSToken
  {
  public:
    VSToken(const std::string& str = ""): fToken(str)
    { /* nothing to see here */ }
    template<typename T> bool convertTo(T& x) const
    { return VSDataConverter::fromString(x,fToken); }
    template<typename T> bool to(T& x) const { return convertTo(x); }
    std::string string() const { return fToken; } 
    std::string lower() const;
    std::string escaped() const;
    bool empty() { return fToken.empty(); }
    void clear() { fToken.clear(); }
    void set(const std::string& s) { fToken=s; }
  private:
    std::string fToken;
  };
    
  class VSTokenList
  {
  public:
    typedef std::vector<VSToken> token_list;

    VSTokenList(): fTokens() { /* nothing to see here */ }

    // Const accessors
    token_list::size_type size() const { return fTokens.size(); }
    bool empty() const { return fTokens.empty(); }
    token_list::const_iterator begin() const { return fTokens.begin(); }
    token_list::const_iterator end() const { return fTokens.end(); }
    const token_list::value_type& operator[](token_list::size_type i) const 
    { return fTokens[i]; }

    // Non-const accessors
    void clear() { fTokens.clear(); }
    void push_back(const token_list::value_type& token) 
    { fTokens.push_back(token); }
    token_list::iterator begin() { return fTokens.begin(); }
    token_list::iterator end() { return fTokens.end(); }
    token_list::value_type& operator[](token_list::size_type i)
    { return fTokens[i]; }
  private:
    token_list fTokens;
  };

  class VSCharSupplier
  {
  public:
    virtual ~VSCharSupplier();
    virtual bool hasChar() = 0;
    virtual char getChar() = 0;
    virtual void nextChar() = 0;
  };

  class VSCSString: public VSCharSupplier
  {
  public:
    VSCSString(const std::string& str): 
      fString(str), fIString(fString.begin()) { }
    virtual ~VSCSString();
    virtual bool hasChar();
    virtual char getChar();
    virtual void nextChar();
    std::string remainder();
  private:
    std::string fString;
    std::string::iterator fIString;
  };

  class VSCSStream: public VSCharSupplier
  {
  public:
    VSCSStream(std::istream& stream): fStream(stream) { }
    virtual ~VSCSStream();
    virtual bool hasChar();
    virtual char getChar();
    virtual void nextChar();
  private:
    std::istream& fStream;
  };

  class VSParserError
  {
  public:
  };

  class VSLineTokenizer
  {
  public:
    enum ParserState 
      { PS_WHITESPACE, PS_COMMENT, PS_NON_QUOTED, PS_QUOTED_S, PS_QUOTED_D };
    
    VSLineTokenizer(bool strict=false): fStrict(strict) { }

    void tokenize(VSCharSupplier* supplier, VSTokenList& tokens);
    void tokenize(std::string& string, VSTokenList& tokens);
    void tokenize(std::istream& stream, VSTokenList& tokens);

  private:
    bool fStrict;
  };

}

#endif // VSLINETOKENIZER_HPP
