//-*-mode:c++; mode:font-lock;-*-

/*! \file VSALocalRegression.hpp

  Local Regression

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       05/01/2010
*/

#ifndef VSALOCALREGRESSION_HPP
#define VSALOCALREGRESSION_HPP

#include <vector>

#include <VSAFunction.hpp>
#include <VSAData.hpp>
#include <VSALinearLeastSquares.hpp>

namespace VERITAS
{
  namespace VSAMath
  {
    template<typename T>
    class LocalRegression1D
    {
    public:
      struct Sort
      {
	bool operator () (const std::pair<double,DataPoint<T> >& left, 
			  const std::pair<double,DataPoint<T> >& right)
	{
	  return left.first < right.first;
	}
      };

      LocalRegression1D(const Data<T>& data, double dx,
			unsigned norder, unsigned min_points = 30):
	m_data(data), m_dx(dx), m_norder(norder),
	m_min_points(min_points)
      { }

      double val(const T& x)
      {
	double v = 0, var = 0;
	val(x,v,var);
	return v;
      }

      void val(const T& x, double& val, double& var)
      {
	// --------------------------------------------------------------------
	// Sort all data points by distance from evaluation point x

	std::vector< std::pair<double,DataPoint<T> > > dxdp;
	for(typename Data<T>::iterator itr = m_data.begin(); 
	    itr != m_data.end(); ++itr)
	  {
	    const T xdiff = (x - itr->x);
	    double dx = fabs(xdiff/m_dx);
	    dxdp.push_back(std::make_pair(dx,*itr));
	  }

	std::sort(dxdp.begin(),dxdp.end(),Sort());
	double xmax = 0;
	unsigned n = 0;
	for(typename std::vector< std::pair<double,DataPoint<T> > >::iterator
	      itr = dxdp.begin(); itr != dxdp.end(); ++itr)
	  {
	    xmax = itr->first;
	    n++;

	    if( n > m_min_points ) break;
	  }

	Data<T> d;     
	for(typename std::vector< std::pair<double,DataPoint<T> > >::iterator
	      itr = dxdp.begin(); itr != dxdp.end(); ++itr)
	  {
	    double dx = itr->first;
	    if(xmax > 1) dx /= xmax;
	    
	    double wt = 0;

	    if(dx < 1) wt = std::pow(1-std::pow(dx,3),3);

	    if(wt != 0) 
	      {
		DataPoint<T> dp = itr->second;
		dp.sigma *= 1./sqrt(wt);
		d.insert(dp);
	      }
	  }

	if(d.size() == 0) return;

	VSAFunction::Poly fn(m_norder);

	VSAMath::Fitsvd<VSAFunction::Poly, double> fit(d,fn);
	fit.fit();
	VSAAlgebra::VecND param = fit.param();	
	VSAAlgebra::MatrixND cov = fit.cov();	
	
	VSAAlgebra::VecND dyda;
	fn.dyda(x,dyda);

	val = fn.val(x,param);

	var = 0;
	for(typename Data<T>::iterator itr = d.begin(); 
	    itr != d.end(); ++itr)
	  var += std::pow(itr->y-fn.val(itr->x,param),2)/(double)d.size();

      }

    private:

      Data<T>  m_data;
      double   m_dx;
      unsigned m_norder;
      unsigned m_min_points;
    };

    template<typename T>
    class LocalRegression2D
    {
    public:

      struct Sort
      {
	bool operator () (const std::pair<double,DataPoint<T> >& left, 
			  const std::pair<double,DataPoint<T> >& right)
	{
	  return left.first < right.first;
	}
      };

      LocalRegression2D(const Data<T>& data, double dx1, double dx2,
			unsigned norder, unsigned min_points = 30):
	m_data(data), m_dx1(dx1), m_dx2(dx2), m_norder(norder),
	m_min_points(min_points)
      { }
      
      double val(const T& x)
      {
	double v = 0, var = 0;
	val(x,v,var);
	return v;
      }

      void val(const T& x, double& val, double& var)
      {
	// --------------------------------------------------------------------
	// Sort all data points by distance from evaluation point x

	std::vector< std::pair<double,DataPoint<T> > > dxdp;
	for(typename Data<T>::iterator itr = m_data.begin(); 
	    itr != m_data.end(); ++itr)
	  {
	    const T xdiff = (x - itr->x);
	    double dx = sqrt(std::pow(xdiff.x()/m_dx1,2) + 
			     std::pow(xdiff.y()/m_dx2,2));

	    dxdp.push_back(std::make_pair(dx,*itr));
	  }

	std::sort(dxdp.begin(),dxdp.end(),Sort());
	double xmax = 0;
	unsigned n = 0;
	for(typename std::vector< std::pair<double,DataPoint<T> > >::iterator
	      itr = dxdp.begin(); itr != dxdp.end(); ++itr)
	  {
	    xmax = itr->first;
	    n++;

	    if( n > m_min_points ) break;

	    // std::cout << itr->first << " " 
	    // 	      << x[0] << " " << x[1] << " "
	    // 	      << itr->second.x[0] << " "
	    // 	      << itr->second.x[1] << std::endl;
	  }

	//	std::cout << "XMAX " << xmax << std::endl;

	Data<T> d;     
	for(typename std::vector< std::pair<double,DataPoint<T> > >::iterator
	      itr = dxdp.begin(); itr != dxdp.end(); ++itr)
	  {
	    double dx = itr->first;
	    if(xmax > 1) dx /= xmax;
	    
	    double wt = 0;

	    if(dx < 1) wt = std::pow(1-std::pow(dx,3),3);

	    if(wt != 0) 
	      {
		DataPoint<T> dp = itr->second;
		dp.sigma *= 1./sqrt(wt);
		d.insert(dp);
	      }
	  }

	// Data<T> d;     
	// for(typename Data<T>::iterator itr = m_data.begin(); 
	//     itr != m_data.end(); ++itr)
	//   {
	//     const T xdiff = (x - itr->x);
	//     double dx = sqrt(std::pow(xdiff.x()/m_dx1,2) + 
	// 		     std::pow(xdiff.y()/m_dx2,2));

	//     dxdp.push_back(std::make_pair(dx,*itr));

	//     double wt = 0;

	//     if(dx < 1) wt = std::pow(1-std::pow(dx,3),3);

	//     if(wt != 0) 
	//       {
	// 	DataPoint<T> dp = *itr;
	// 	dp.sigma *= 1./sqrt(wt);
	// 	d.insert(dp);

	// 	// std::cout << itr->x[0] << " " << itr->x[1] << " " << " "
	// 	// 	  << dx << " "
	// 	// 	  << wt << std::endl;
	//       }


	//   }

	if(d.size() == 0) return;

	VSAFunction::Poly2D<VSAAlgebra::Vec2D> fn(m_norder);

	VSAMath::Fitsvd<VSAFunction::Poly2D<VSAAlgebra::Vec2D>, 
	  VSAAlgebra::Vec2D> fit(d,fn);
	fit.fit();
	VSAAlgebra::VecND param = fit.param();	
	VSAAlgebra::MatrixND cov = fit.cov();	
	VSAAlgebra::VecND dyda;
	fn.dyda(x,dyda);

	val = fn.val(x,param);


	var = 0;
	for(typename Data<T>::iterator itr = d.begin(); 
	    itr != d.end(); ++itr)
	  var += std::pow(itr->y-fn.val(itr->x,param),2)/(double)d.size();

	//std::cout << val << " " << var << " " << dyda*(cov*dyda) << std::endl;

	//var = dyda*(cov*dyda);
      }

    private:

      Data<T>  m_data;
      double   m_dx1;
      double   m_dx2;
      unsigned m_norder;
      unsigned m_min_points;
    };

  }
}

#endif // VSALOCALREGRESSION_HPP
