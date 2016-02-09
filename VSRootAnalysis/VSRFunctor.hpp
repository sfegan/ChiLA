//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRFUNCTOR_HPP
#define VSRFUNCTOR_HPP

#include <TF1.h>
#include <TGraph.h>
#include <TVector.h>

#include <VSNSpace.hpp>


class VSRFunctor
{
public:
  VSRFunctor() {}
  virtual ~VSRFunctor() {}

  virtual double operator() (const std::vector< double >& x) const = 0;
  virtual double operator() (const double* x, const double* p) const = 0;
  virtual double eval(const std::vector< double >& x) const = 0;
  virtual unsigned getNDim() const = 0;
  virtual VSRFunctor* clone() const = 0;

  // --------------------------------------------------------------------------
  // Copy Constructor
  //---------------------------------------------------------------------------
  //  VSRFunctor(const VSRFunctor& func);

  // --------------------------------------------------------------------------
  // Assignment Operator
  //---------------------------------------------------------------------------
  //  VSRFunctor& operator=(const VSRFunctor& func);
};

class VSRFunctorParameterized : public VSRFunctor
{
public:
  VSRFunctorParameterized() {}
  virtual ~VSRFunctorParameterized() {}


  virtual double operator() (const std::vector< double >& x) const = 0;
  virtual double eval(const std::vector< double >& x) const = 0;


  virtual double operator() (const std::vector< double >& x,
			     const std::vector< double >& p) const = 0;
  virtual double operator() (const double* x, const double* p) const = 0;

  virtual double eval(const std::vector< double >& x,
		      const std::vector< double >& p) const = 0;

  virtual unsigned getNDim() const = 0;
  virtual unsigned getNParam() const = 0;
  virtual double getParam(unsigned iparam) const = 0;

  virtual void setParam(unsigned iparam, double value) = 0;

  virtual VSRFunctorParameterized* clone() const = 0;
};


// ----------------------------------------------------------------------------
// Power-law function with 2 free parameters.
// ----------------------------------------------------------------------------

class VSRFunctorPL : public VSRFunctorParameterized
{
public:
  VSRFunctorPL(double constant, double index);
  virtual ~VSRFunctorPL() { }

  double operator() (const std::vector< double >& x) const
  {  
    return m_constant*std::pow(x[0],m_index);
  }

  double eval(const std::vector< double >& x) const
  { 
    return m_constant*std::pow(x[0],m_index); 
  }

  double operator() (const std::vector< double >& x,
		     const std::vector< double >& p) const
  {  
    return p[0]*std::pow(x[0],p[1]);
  }

  double operator() (const double* x, const double* p) const
  {
    return m_constant*std::pow(x[0],m_index);
  }

  double eval(const std::vector< double >& x,
	      const std::vector< double >& p) const
  { 
    return p[0]*std::pow(x[0],p[1]); 
  }

  unsigned getNDim() const { return 1; }
  unsigned getNParam() const { return 2; }

  double getParam(unsigned iparam) const 
  { 
    if(iparam == 0) return m_constant;
    else return m_index;
  }

  void setParam(unsigned iparam, double value)
  {
    if(iparam == 0) m_constant = value;
    else m_index = value;
  }

  virtual VSRFunctorPL* clone() const 
  {
    return new VSRFunctorPL(m_index,m_constant); 
  }

private:
  double    m_constant;
  double    m_index;
};

class VSRFunctorDiscrete : public VSRFunctor
{
public:
  typedef double (*FuncPtr) (double);

  VSRFunctorDiscrete(VERITAS::VSNSpace& nspace, double min, double max);
  VSRFunctorDiscrete(TVectorT<double>& vec, double min, double max);
  ~VSRFunctorDiscrete();

  double operator() (const std::vector< double >& x) const
  {
    return m_gr->Eval(x[0]); 
  }

  double operator() (const double* x, const double* p) const
  {
    return m_gr->Eval(x[0]); 
  }

  double eval(const std::vector< double >& x) const 
  { 
    return m_gr->Eval(x[0]); 
  }

  unsigned getNDim() const { return 1; }

  double integral(double x1, double x2)
  {
    return m_func->Integral(x1,x2);
  }

  void multiply(VSRFunctor* func)
  {
    std::vector< double > x(1);

    for(unsigned i = 0; i < (unsigned)m_gr->GetN(); i++)
      {
	x[0] = m_gr->GetX()[i];
	m_gr->GetY()[i] *= func->eval(x);
      }

    delete m_func;
    m_func = new TF1("",this,m_xmin,m_xmax);
  }

  virtual VSRFunctorDiscrete* clone() const
  {
    return new VSRFunctorDiscrete(*this); 
  }

  // --------------------------------------------------------------------------
  // Copy Constructor
  //---------------------------------------------------------------------------
  VSRFunctorDiscrete(const VSRFunctorDiscrete& func);

  // --------------------------------------------------------------------------
  // Assignment Operator
  //---------------------------------------------------------------------------
  VSRFunctorDiscrete& operator=(const VSRFunctorDiscrete& func);

private:
  TGraph* m_gr;
  TF1*    m_func;
  
  double m_xmin;
  double m_xmax;

};

template< typename T >
class VSRFunctorMF : public VSRFunctor
{
public:
  typedef double (T::* MFPtr) (double);

  VSRFunctorMF(T* obj, MFPtr mfp): m_obj(obj), m_mfp(mfp) { }
  virtual ~VSRFunctorMF() { }

  double operator() (const std::vector< double >& x) const
  {
    return (m_obj->*m_mfp)(x[0]);
  }
  
  double operator() (const double* x, const double* p) const
  {
    return (m_obj->*m_mfp)(x[0]);
  }

  double eval(const std::vector< double >& x) const
  {
    return (m_obj->*m_mfp)(x[0]);
  }

  unsigned getNDim() const { return 1; }

  virtual VSRFunctorMF* clone() const
  {
    return new VSRFunctorMF(*this); 
  }

private:

  T*            m_obj;
  MFPtr         m_mfp;
};

#endif // VSRFUNCTOR_HPP
