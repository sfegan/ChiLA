//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleGraph.hpp

  Simple graph class.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/20/2007
*/

#ifndef VSSIMPLEGRAPH_HPP
#define VSSIMPLEGRAPH_HPP

#include <vector>

#ifndef NOHDF5
#include <VSOctaveIO.hpp>
#endif

namespace VERITAS
{
  
  template<typename TX, typename TY>
  class VSSimpleGraph
  {
  public:
    class iterator
    {
    public:
      iterator(const iterator& i): fGraph(i.fGraph), fVertex(i.fVertex) { }
      TX x() const { return fGraph->x(fVertex); }
      TY y() const { return fGraph->y(fVertex); }
      TX xerr() const { return fGraph->xerr(fVertex); }
      TY yerr() const { return fGraph->yerr(fVertex); }
      int vertex() const { return fVertex; }
      const iterator* operator->() { return this; }
      iterator& operator++() { fVertex++; return *this; }
      iterator operator++(int) { iterator i=*this; fVertex++; return i; }
      iterator& operator--() { fVertex--; return *this; }
      iterator operator--(int) { iterator i=*this; fVertex--; return i; }
      iterator& operator=(const iterator& i) 
      { fGraph=i.fGraph; fVertex=i.fVertex; return *this; }
      iterator& operator+=(int b) { fVertex+=b; return *this; }
      iterator& operator-=(int b) { fVertex-=b; return *this; }
      iterator operator+(int b) const { iterator i(*this); i+=b; return i; }
      iterator operator-(int b) const { iterator i(*this); i-=b; return i; }
      int operator-(const iterator& i) { return fVertex-i.fVertex; }
      bool operator!=(const iterator& o) const { return fVertex != o.fVertex; }
      bool operator==(const iterator& o) const { return fVertex == o.fVertex; }
      bool operator<(const iterator& o) const { return fVertex < o.fVertex; }
      bool operator<=(const iterator& o) const { return fVertex <= o.fVertex; }
      bool operator>(const iterator& o) const { return fVertex > o.fVertex; }
      bool operator>=(const iterator& o) const { return fVertex >= o.fVertex; }
    private:
      friend class VSSimpleGraph;
      iterator(const VSSimpleGraph* graph, int vertex): 
	fGraph(graph), fVertex(vertex) { }
      const VSSimpleGraph* fGraph;
      int                  fVertex;
    };

    class Vertex
    {
    public:
      Vertex(TX x, TX xerr, TY y, TY yerr):
	fX(x),fXErr(xerr),fY(y),fYErr(yerr) {}
      Vertex(TX x, TY y):
	fX(x),fXErr(0),fY(y),fYErr(0) {}
    private:
      friend class VSSimpleGraph;
      TX fX;
      TX fXErr;
      TY fY;
      TY fYErr;
    };

    VSSimpleGraph(const std::string& name = ""):
      fVertices(), fName(name) {}

    VSSimpleGraph(const std::vector< std::pair<TX,TY> >& xy,
		  const std::string& name = ""):
      fVertices(), fName(name) 
    {
      for(typename std::vector< std::pair<TX,TY> >::const_iterator itr = 
	    xy.begin(); itr != xy.end(); ++itr)
 	addVertex(itr->first,itr->second);
    }

#ifndef NOHDF5
    bool load(VSOctaveH5ReaderStruct* reader);
    bool save(VSOctaveH5WriterStruct* writer) const;    
#endif

    const std::string& name() const { return fName; }

    void clear() { fVertices.clear(); }

    inline TX x(int vertex) const;
    inline TX xerr(int vertex) const;
    inline TY y(int vertex) const;
    inline TY yerr(int vertex) const;

    inline iterator begin() const;
    inline iterator end() const;

    inline void addVertex(TX x, TY y);
    inline void addVertex(TX x, TX xerr, TY y, TY yerr);
    inline unsigned nvertex() { return fVertices.size(); }

  private:
    std::vector< Vertex > fVertices;
    std::string           fName;
  };

#ifndef NOHDF5
  template<typename TX, typename TY> 
  bool VSSimpleGraph<TX,TY>::load(VSOctaveH5ReaderStruct* reader)
  {
    if(!reader) return false;
    reader->readString("name",fName);
    std::vector<TX> x;
    std::vector<TX> xerr;
    std::vector<TY> y;
    std::vector<TY> yerr;

    reader->readVector("x",x);
    reader->readVector("x_error",xerr);
    reader->readVector("y",y);
    reader->readVector("y_error",yerr);

    if(x.size() != y.size())return false;
    clear();

    const unsigned nvertex = x.size();
    for(unsigned ivertex=0;ivertex<nvertex;ivertex++)
      fVertices.push_back(Vertex(x[ivertex],xerr[ivertex],
				 y[ivertex],yerr[ivertex]));

    return true;
  }

  template<typename TX, typename TY> 
  bool VSSimpleGraph<TX,TY>::save(VSOctaveH5WriterStruct* writer) const
  {
    std::vector<TX> x;
    std::vector<TX> xerr;
    std::vector<TY> y;
    std::vector<TY> yerr;
    const unsigned nvertex = fVertices.size();
    x.resize(nvertex);
    xerr.resize(nvertex);
    y.resize(nvertex);
    yerr.resize(nvertex);
    for(iterator itr=begin(); itr!=end();itr++)
      {
	x[itr.vertex()] = itr->x();
	xerr[itr.vertex()] = itr->xerr();
	y[itr.vertex()] = itr->y();
	yerr[itr.vertex()] = itr->yerr();
      }

    writer->writeString("name",fName);
    writer->writeVector("x",x);
    writer->writeVector("x_error",xerr);
    writer->writeVector("y",y);
    writer->writeVector("y_error",yerr);
    return true;
  }
#endif

  template<typename TX, typename TY>
  inline TX VSSimpleGraph<TX,TY>::x(int vertex) const
  {
    return fVertices[vertex].fX;
  }

  template<typename TX, typename TY>
  inline TX VSSimpleGraph<TX,TY>::xerr(int vertex) const
  {
    return fVertices[vertex].fXErr;
  }

  template<typename TX, typename TY>
  inline TY VSSimpleGraph<TX,TY>::y(int vertex) const
  {
    return fVertices[vertex].fY;
  }

  template<typename TX, typename TY>
  inline TY VSSimpleGraph<TX,TY>::yerr(int vertex) const
  {
    return fVertices[vertex].fYErr;
  }

  template<typename TX, typename TY>
  inline typename VSSimpleGraph<TX,TY>::iterator 
  VSSimpleGraph<TX,TY>::begin() const
  {
    return iterator(this,0);
  }

  template<typename TX, typename TY>
  inline typename VSSimpleGraph<TX,TY>::iterator 
  VSSimpleGraph<TX,TY>::end() const
  {
    return iterator(this,fVertices.size());
  }

  template<typename TX, typename TY>
  inline void VSSimpleGraph<TX,TY>::addVertex(TX x, TY y) 
  { 
    fVertices.push_back( Vertex(x,y) ); 
  }

  template<typename TX, typename TY>
  inline void VSSimpleGraph<TX,TY>::addVertex(TX x, TX xerr, TY y, TY yerr) 
  { 
    fVertices.push_back( Vertex(x,xerr,y,yerr) ); 
  }

}

#endif // VSSIMPLEGRAPH_HPP
