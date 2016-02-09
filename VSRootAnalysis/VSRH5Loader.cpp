#include "VSRH5Loader.hpp"


using namespace VERITAS;

TMatrixT<double> VSRH5Loader::createMatrix(const VSNSpace& nspace)
{
  TMatrixT<double> m(nspace.space().axes[0].nbin,
		     nspace.space().axes[1].nbin);
  
  VSNSpace::Cell c(2);
  for(c.i[0]=0; c.i[0]<nspace.space().axes[0].nbin; c.i[0]++)
    {
      for(c.i[1]=0; c.i[1]<nspace.space().axes[1].nbin; c.i[1]++)
	{
	  VSNSpace::Weight w = 0;
	  nspace.getWeight(c,w);
	  
	  if(std::isfinite(w))
	    m[c.i[0]][c.i[1]] = w;
	}
    }
  
  return m;
}

TVectorT<double> VSRH5Loader::createVector(const VSNSpace& nspace)
{
  TVectorT<double> v(nspace.space().axes[0].nbin);
  
  VSNSpace::Cell c(1);
  for(c.i[0]=0; c.i[0]<nspace.space().axes[0].nbin; c.i[0]++)
    {
      VSNSpace::Weight w = -1; // initialize any value
      nspace.getWeight(c,w);
      
      if(std::isfinite(w))
	v[c.i[0]] = w;
    }

  return v;
}


VERITAS::VSNSpace VSRH5Loader::createNSpace(const TVectorT<double>& v,
					    double lo, double hi)
{
  const unsigned nrow = v.GetNrows();

  VSNSpace::Space space(1);
  space.axes[0] = VSNSpace::Axis(lo, hi, nrow, "");

  VERITAS::VSNSpace nspace(space);

  for(unsigned irow = 0; irow < nrow; irow++)
    nspace[irow]=v[irow];

  return nspace;
}
