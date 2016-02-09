
#include "VSRFunctor.hpp"

using namespace VERITAS;

VSRFunctorPL::VSRFunctorPL(double constant,double index):
  m_constant(constant), m_index(index)
{

}

VSRFunctorDiscrete::VSRFunctorDiscrete(VERITAS::VSNSpace& nspace, 
				       double xmin, double xmax):
  m_gr(), m_func(), m_xmin(xmin), m_xmax(xmax)
{
  m_gr = new TGraph;
  
  VSNSpace::Cell c(1);
  for(c.i[0] = 0; c.i[0] < nspace.space().axes[0].nbin; c.i[0]++)
    {
      VSNSpace::Point p(1);
      nspace.space().midPointOfCellUnchecked(c,p);	
      m_gr->SetPoint(c.i[0],p.x[0], nspace[c]);	
    }
  
  m_func = new TF1("",this,xmin,xmax);
}

VSRFunctorDiscrete::VSRFunctorDiscrete(TVectorT<double>& vec, 
					 double min, double max)
{
  m_gr = new TGraph;

  const unsigned nrow = vec.GetNrows();
  double bin_size = (max-min)/(double)nrow;

  for(unsigned irow = 0; irow < nrow; irow++)
    {
      m_gr->SetPoint(irow,min+bin_size/2.+irow*bin_size,vec[irow]);
    }
  
  m_func = new TF1("",this,min,max);
}

VSRFunctorDiscrete::~VSRFunctorDiscrete() 
{ 
  delete m_gr; 
  delete m_func;
}


VSRFunctorDiscrete::VSRFunctorDiscrete(const VSRFunctorDiscrete& func)
{
  *this = func;
}

VSRFunctorDiscrete& 
VSRFunctorDiscrete::operator=(const VSRFunctorDiscrete& func) 
{
  m_gr = (TGraph*)func.m_gr->Clone();
  m_func = (TF1*)func.m_func->Clone();
  return *this;
}
