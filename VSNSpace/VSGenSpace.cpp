/*! \file VSGenSpace.cpp

  GenSpace generalizeed multi-dimensional space

  Multi dimensional space whose cells are arbitrary types. Has less
  functionality than normal NSpace but can be useful in some
  circumstances.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/26/2007

*/

#include <VSGenSpace.hpp>

using namespace VERITAS;

VSNSpace::Weight 
OneSidedIntervalWeight::
interpolateArrayFraction(const std::vector<std::pair<VSNSpace::Weight,
			 VSNSpace::Weight> >& w, double sum, double frac)
{
  unsigned nw = w.size();
  if(nw==0)return 0;
  else if((nw==1)||(frac<0.0))return w.front().first;
  else if(frac>=1.0)return w.back().first;
  
  sum -=  w.back().second;
  double dp = 0;
  double p = 0;
  unsigned iindex=0;
  while(1)
    {
      dp = w[iindex].second/sum;
      if(p+dp > frac || iindex + 2 == nw) break;
      iindex++;
      p+=dp;
    }
  
  double y1 = w[iindex].first;
  double y2 = w[iindex+1].first;
  return y1 + (frac-p)*(y2-y1)/dp;
}
