//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRHISTOGRAM_HPP
#define VSRHISTOGRAM_HPP

#include <vector>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TCanvas.h>
#include <TPad.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRCanvasElement.hpp"

class VSRHistogram : public VSRCanvasElement
{

public:
  VSRHistogram();
  virtual ~VSRHistogram();

  virtual void draw() {}
  virtual VSRHistogram* clone() const = 0;

  virtual void addLegendEntry(TLegend *leg, 
			      const std::string& label) const = 0;
  
  virtual double getLoX() const = 0;
  virtual double getHiX() const = 0;
  virtual double getLoY() const = 0;
  virtual double getHiY() const = 0;
  virtual double getPosLoX() const = 0;
  virtual double getPosHiX() const = 0;
  virtual double getPosLoY() const = 0;
  virtual double getPosHiY() const = 0;

  virtual bool is1D() const = 0;
  virtual bool is2D() const = 0;

  virtual void multiply(double x) { }

  virtual void normalize() {}
  virtual void fill(double x) {}
  virtual void fill(double x, double y) {}
  virtual void set(unsigned i, double x) {}
  virtual void initialize(unsigned n, double lo, double hi) {}
};


#endif // VSRHISTOGRAM_HPP
