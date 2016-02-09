
//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRGRAPH1D_HPP
#define VSRGRAPH1D_HPP

#include <vector>
#include <memory>

#include <VSSimpleGraph.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRCanvasElement.hpp"

class VSRGraph1D : public VSRCanvasElement
{
public:
  VSRGraph1D();
  ~VSRGraph1D();

  // --------------------------------------------------------------------------
  // Virtual Methods
  // --------------------------------------------------------------------------
  virtual void draw();
  virtual VSRGraph1D* clone() const;

  virtual double getLoX() const { return m_loX; }
  virtual double getHiX() const { return m_hiX; }
  virtual double getLoY() const { return m_loY; }
  virtual double getHiY() const { return m_hiY; }
  virtual double getPosLoX() const { return m_posLoX; }
  virtual double getPosHiX() const { return m_posHiX; }
  virtual double getPosLoY() const { return m_posLoY; }
  virtual double getPosHiY() const { return m_posHiY; }

  virtual void addLegendEntry(TLegend *leg, const std::string& label) const
  {
    leg->AddEntry(m_graph.get(),label.c_str(),m_options.legendStyle.c_str());
  }

  virtual bool is1D() const { return true; }
  virtual bool is2D() const { return false; }

  virtual VSRGraph1D& operator*= (double x);

  // --------------------------------------------------------------------------
  // Set Methods
  // --------------------------------------------------------------------------
  void set(double x, double y);
  void set(double x, double y, double yerr);

  // --------------------------------------------------------------------------
  // Get Methods
  // --------------------------------------------------------------------------
  double getX(unsigned ipoint) const;
  double getErrX(unsigned ipoint) const;
  double getY(unsigned ipoint) const;
  double getErrY(unsigned ipoint) const;
  unsigned getN() const { return m_nPoints; }

  void print();
  void dump();


  TGraphErrors* getGraph() { return m_graph.get(); }

  double eval(double x) { return m_graph->Eval(x); }

  void clear();

  // --------------------------------------------------------------------------
  // Factory methods
  //---------------------------------------------------------------------------
  template< typename T >
  static VSRGraph1D* create(VERITAS::VSSimpleGraph<T,T>* graph);

  // --------------------------------------------------------------------------
  // Copy Constructor
  //---------------------------------------------------------------------------
  VSRGraph1D(const VSRGraph1D& hist);

  // --------------------------------------------------------------------------
  // Assignment Operator
  //---------------------------------------------------------------------------
  VSRGraph1D& operator=(const VSRGraph1D & hist);

  // --------------------------------------------------------------------------
  // Arithmetic Operator
  //---------------------------------------------------------------------------
  VSRGraph1D& operator/=(const VSRGraph1D & graph);

private:

  double                      m_loX;
  double                      m_hiX;
  double                      m_loY;
  double                      m_hiY;
  double                      m_posLoX;
  double                      m_posHiX;
  double                      m_posLoY;
  double                      m_posHiY;
  unsigned                    m_nPoints;

  std::auto_ptr<TGraphErrors> m_graph;
};

template< typename T >
VSRGraph1D* VSRGraph1D::create(VERITAS::VSSimpleGraph<T,T>* graph)
{
  VSRGraph1D* o = new VSRGraph1D;
  
  for(typename VERITAS::VSSimpleGraph<T,T>::iterator itr = graph->begin(); 
      itr != graph->end(); ++itr)
    {
      if(o->m_nPoints == 0 || itr->x() - itr->xerr() < o->m_loX)
	o->m_loX = itr->x() - itr->xerr();
      
      if(o->m_nPoints == 0 || itr->x() + itr->xerr() > o->m_hiX)
	o->m_hiX = itr->x() + itr->xerr();

      if(o->m_nPoints == 0 || itr->y() - itr->yerr() < o->m_loY)
	o->m_loY = itr->y() - itr->yerr();
      
      if(o->m_nPoints == 0 || itr->y() + itr->yerr() > o->m_hiY)
	o->m_hiY = itr->y() + itr->yerr();

      o->set(itr->x(),itr->y(),itr->yerr());
    }  

  return o;
}



#endif // VSRGRAPH1D_HPP
