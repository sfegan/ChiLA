/*! \file VSNSpace.cpp

  NSpace multi-dimensional cutting

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/16/2006

  Based on Utah version, made to agree with convetions of VS package

  This file includes the source code for the methods defined int
  NSpace.hpp.

  \author  Name            Jeter Hall \n
           Institution     U of U \n
           E-mail          jeter@physics.utah.edu

  \author  Name            V.V. Vassiliev \n
           Institution     U of U \n
           E-mail          vvv@physics.utah.edu

  \author  Name            T. Nagai \n
           Institution     U of U \n
           E-mail          tn68@physics.utah.edu

  OLD_date    June 1, 2003

  OLD_version 0.0

  OLD_note
 
*/

#include <cmath>
#include <algorithm>
#include <iostream>

#include <fast_alloc.hpp>
#include <VSAssert.hpp>
#include "VSNSpace.hpp"
#include "VSNSpaceOctaveH5IO.hpp"

using namespace VERITAS;

// ****************************************************************************
//
// VSNSpace::Volume
//
// ****************************************************************************

bool VSNSpace::Volume::
load(const Space& space, const std::vector<bool>& volume)
{
  if(space.size() != volume.size())return false;
  m_space = space;
  m_volume = volume;
  return true;
}

bool VSNSpace::Volume::
loadSparse(const Space& space, std::vector<Index>& index)
{
  unsigned n = space.size();
  m_space = space;
  m_volume.clear();
  m_volume.resize(n);
  for(std::vector<Index>::const_iterator iindex = index.begin();
      iindex != index.end(); iindex++)
    {
      if(*iindex >= n)return false;
      else m_volume[*iindex] = true;
    }
  return true;
}

void VSNSpace::Volume::
sparse(std::vector<Index>& index) const
{
  index.clear();
  Index iindex = 0;
  for(std::vector<bool>::const_iterator ivol = m_volume.begin();
      ivol != m_volume.end(); ivol++)
    {
      if(*ivol)index.push_back(iindex);
      iindex++;
    }
}

unsigned VSNSpace::Volume::volumeSize() const
{
  unsigned volume_size=0;
  for(std::vector<bool>::const_iterator ivol = m_volume.begin();
      ivol != m_volume.end(); ivol++)
    if(*ivol)volume_size++;
  return volume_size;
}

void VSNSpace::Volume::setInvert()
{
  for(std::vector<bool>::iterator ivol = m_volume.begin();
      ivol != m_volume.end(); ivol++)*ivol = !(*ivol);
}

bool VSNSpace::Volume::setIntersect(const Volume& o)
{
  if(!m_space.isSpaceEqual(o.m_space))return false;
  std::vector<bool>::const_iterator ovol = o.m_volume.begin();
  for(std::vector<bool>::iterator ivol = m_volume.begin();
      ivol != m_volume.end(); ivol++)
    *ivol = (*ivol && *(ovol++));
  return true;
}

bool VSNSpace::Volume::setUnion(const Volume& o)
{
  if(!m_space.isSpaceEqual(o.m_space))return false;
  std::vector<bool>::const_iterator ovol = o.m_volume.begin();
  for(std::vector<bool>::iterator ivol = m_volume.begin();
      ivol != m_volume.end(); ivol++)
    *ivol = (*ivol || *(ovol++));
  return true;
}

bool VSNSpace::Volume::setLess(const Volume& o)
{
  if(!m_space.isSpaceEqual(o.m_space))return false;
  std::vector<bool>::const_iterator ovol = o.m_volume.begin();
  for(std::vector<bool>::iterator ivol = m_volume.begin();
      ivol != m_volume.end(); ivol++)
    *ivol = (*ivol && !(*(ovol++)));
  return true;
}

bool VSNSpace::Volume::clearInsideRange(unsigned dim, 
					const Index lo, const Index hi)
{
  unsigned index = 0;
  for(std::vector<bool>::iterator ivol = m_volume.begin();
      ivol != m_volume.end(); ivol++)
    {
      Cell c;
      m_space.cellOfIndexUnchecked(index,c);
      if(c.i[dim] >= lo && c.i[dim] <= hi) *ivol = 0;
      index++;
    }
  return true;
}

unsigned VSNSpace::Volume::countCellsInVolume() const
{
  unsigned ncell = 0;
  for(std::vector<bool>::const_iterator ivol = m_volume.begin();
      ivol != m_volume.end(); ivol++)if(*ivol)ncell++;
  return ncell;
}

void VSNSpace::Volume::expandFromEdge()
{
  std::vector<bool> vol;

  for(unsigned icell = 0; icell < m_space.size(); icell++)
    {
      if(isIndexInVolumeUnchecked(icell)) vol.push_back(true);
      else if(isIndexAdjacent(icell)) vol.push_back(true);
      else vol.push_back(false);
    }

  m_volume = vol;
}

// ****************************************************************************
//
// VSNSpace
//
// ****************************************************************************

VSNSpace::VSNSpace(const Space& space, const std::string& comment,
		   const Weight zero_weight):
  m_comment(comment), m_space(space), m_weight()
{
  clear(zero_weight);
}

VSNSpace::VSNSpace(const Volume& volume, const std::string& comment,
		   const Weight zero_weight, const Weight one_weight):
  m_comment(comment), m_space(volume.m_space), m_weight()
{
  clear(zero_weight);
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight != nweight)
    {
      if(volume.m_volume[iindex])*iweight = one_weight;
      iindex++;
      iweight++;
    }
}

VSNSpace::VSNSpace(unsigned ndim, Coord lo, Coord hi, Index nbin):
  m_comment(), m_space(), m_weight()
{
  m_space = VSNSpace::Space(ndim);
  for(unsigned idim = 0; idim < ndim; idim++)
    m_space.axes[idim] = VSNSpace::Axis(lo, hi, nbin);
  clear(0.0);
}

VSNSpace::VSNSpace(Coord lo, Coord hi, Index nbin):
  m_comment(), m_space(), m_weight()
{
  const unsigned ndim = 1;

  m_space = VSNSpace::Space(ndim);
  m_space.axes[0] = VSNSpace::Axis(lo, hi, nbin);
  clear(0.0);
}

VSNSpace::VSNSpace(Coord lox, Coord hix, Index nbinx,
		   Coord loy, Coord hiy, Index nbiny):
  m_comment(), m_space(), m_weight()
{
  const unsigned ndim = 2;

  m_space = VSNSpace::Space(ndim);
  m_space.axes[0] = VSNSpace::Axis(lox, hix, nbinx);
  m_space.axes[1] = VSNSpace::Axis(loy, hiy, nbiny);
  clear(0.0);
}

VSNSpace::VSNSpace(const VSLimitedHist<double,double>& h):
  m_comment(), m_space(), m_weight()
{
  const unsigned ndim = 1;
  m_space = VSNSpace::Space(ndim);
  m_space.axes[0] = VSNSpace::Axis(h.loLimit(), h.hiLimit(), h.nBins());
  clear(0.0);

  for(VSLimitedHist<double,double>::iterator itr = h.begin(); itr !=
	h.end(); ++itr)
    {
      VSNSpace::Point p(1);
      p.x[0] = itr->center();
      setWeight(p,itr->count());
    }
}

VSNSpace::VSNSpace(const VSSimple2DHist<double,double>& h):
  m_comment(), m_space(), m_weight()
{
  const unsigned ndim = 2;
  m_space = VSNSpace::Space(ndim);
  m_space.axes[0] = VSNSpace::Axis(h.xLoLimit(), h.xHiLimit(), h.nXBins());
  m_space.axes[1] = VSNSpace::Axis(h.yLoLimit(), h.yHiLimit(), h.nYBins());
  clear(0.0);

  for(VSSimple2DHist<double,double>::iterator itr = h.begin(); itr !=
	h.end(); ++itr)
    {
      VSNSpace::Point p(2);
      p.x[0] = itr->x();
      p.x[1] = itr->y();
      setWeight(p,itr->count());
    }
}

bool VSNSpace::load(const Space& space, const std::vector<Weight>& weight,
		    const std::string& comment)
{
  if(space.size() != weight.size())return false;
  m_comment = comment;
  m_space = space;
  m_weight = weight;
  return true;
}

void VSNSpace::hist(const std::set< unsigned >& dims, 
		    const std::vector< double >& coords,
		    VSLimitedHist<double,double>& h)
{
  std::vector< unsigned > indices;

  std::vector<double>::const_iterator citr = coords.begin();

  for(std::set< unsigned >::const_iterator itr = dims.begin(); 
      itr != dims.end(); ++itr)
    {
      indices.push_back(m_space.axes[*itr].indexUnchecked(*citr));
      citr++;
    }

  hist(dims,indices,h);
}

void VSNSpace::hist(const std::set< unsigned >& dims, 
		    const std::vector< unsigned >& indices,
		    VSLimitedHist<double,double>& h)
{
  vsassert(m_space.ndim-dims.size()==1);

  unsigned idim = 0;
  for(idim = 0; idim < m_space.ndim; idim++)
    if(dims.find(idim) == dims.end()) break;

  vsassert(idim < m_space.ndim);

  
  h = VSLimitedHist<double,double>(m_space.axis(idim).bin_size,
				   m_space.axis(idim).lo_bound,
				   m_space.axis(idim).hi_bound);

  VSNSpace::Cell c = m_space.cell();

  std::vector< unsigned >::const_iterator iitr = indices.begin();
  for(std::set<unsigned>::const_iterator itr = dims.begin(); itr !=
	dims.end(); ++itr)
    {
      vsassert(*itr < m_space.ndim);

      if(iitr != indices.end()) c.i[*itr] = *iitr;
      else c.i[*itr] = 0;

      iitr++;
    }

  for(c.i[idim] = 0; c.i[idim] < m_space.axis(idim).nbin; c.i[idim]++)
    {
      double x = m_space.axis(idim).midCoordUnchecked(c.i[idim]);
      h.accumulate(x,(*this)[c]);
    }
}

// ----------------------------------------------------------------------------
//
// ACCESSORS
//
// ----------------------------------------------------------------------------

bool VSNSpace::interpolateWeight(const Point& p, Weight& weight) const
{
  assert(m_space.ndim < sizeof(unsigned)*8);
  if(p.ndim != m_space.ndim)return false;

  Coord* p_shift_array = FASTCALLOC(Coord, m_space.ndim);
  Point p_shift(m_space.ndim, p_shift_array);
  p_shift = p;

  for(unsigned idim=0;idim<m_space.ndim;idim++)
    p_shift.x[idim] -= m_space.axes[idim].bin_size * 0.5;

  Index* c_array = FASTCALLOC(Index, m_space.ndim);
  Cell c(m_space.ndim, c_array);

  Coord* x_array = FASTCALLOC(Coord, m_space.ndim);
  Point x(m_space.ndim, x_array);

  for(unsigned idim=0;idim<m_space.ndim;idim++)
    {
      if(p_shift.x[idim] <= m_space.axes[idim].lo_bound)
	c.i[idim] = 0;
      else 
	{
	  c.i[idim] = 
	    std::min(m_space.axes[idim].indexUnchecked(p_shift.x[idim]),
		     m_space.axes[idim].nbin-2);
	}

#warning temporary assert here
      vsassert(c.i[idim] < m_space.axes[idim].nbin - 1);
    }

  Point min_point;
  m_space.minPointOfCellUnchecked(c,min_point);
  x = p_shift;
  x -= min_point;
  for(unsigned idim=0;idim<m_space.ndim;idim++)
    x.x[idim] *= m_space.axes[idim].bin_factor;

  Index* cc_array = FASTCALLOC(Index, m_space.ndim);
  Cell cc(m_space.ndim, cc_array);

  unsigned icell = 0;
  unsigned icmax = 1<<m_space.ndim;
  weight = 0;
  while(icell < icmax)
    {
      double xfact;
      if(icell & 1)xfact = *x.x, *cc.i = *c.i + 1;
      else xfact = 1.0 - *x.x, *cc.i = *c.i;

      for(unsigned idim=1;idim<m_space.ndim;idim++)
	if(icell & (1<<idim))xfact *= x.x[idim], cc.i[idim] = c.i[idim] + 1;
	else xfact *= 1.0 - x.x[idim], cc.i[idim] = c.i[idim];

//      weight += xfact*getWeightUnchecked(m_space.indexOfCellUnchecked(cc));
      double w = 0;
      getWeight(m_space.indexOfCellUnchecked(cc),w);

      weight += xfact*w;
      icell++;
    }

  FASTFREE(cc_array);
  FASTFREE(x_array);
  FASTFREE(c_array);
  FASTFREE(p_shift_array);

  return true;
}

VSNSpace::Weight VSNSpace::totalWeight() const
{
  Weight total = 0;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(iweight!=nweight)total+=*(iweight++);
  return total;
}

VSNSpace::Weight VSNSpace::maxWeight() const
{
  Weight max_weight = m_weight.front();
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(++iweight!=nweight)if(*iweight > max_weight)max_weight=*iweight;
  return max_weight;
}

VSNSpace::Weight VSNSpace::maxWeight(Index& index) const
{
  Weight max_weight = m_weight.front();
  index = 0;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(++iweight!=nweight)if(*iweight > max_weight)
    {
      max_weight=*iweight;
      index = iweight - &(*m_weight.begin());
    }
  return max_weight;
}

VSNSpace::Weight VSNSpace::minWeight() const
{
  Weight min_weight = m_weight.front();
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(++iweight!=nweight)if(*iweight < min_weight)min_weight=*iweight;
  return min_weight;
}

VSNSpace::Weight VSNSpace::minWeight(Index& index) const
{
  Weight min_weight = m_weight.front();
  index = 0;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(++iweight!=nweight)if(*iweight < min_weight)
    {
      min_weight=*iweight;
      index = iweight - &(*m_weight.begin());
    }
  return min_weight;
}

bool VSNSpace::totalWeight(const Volume& volume, Weight& total) const
{
  if(!m_space.isSpaceEqual(volume.space()))return false;
  total = 0;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)
    {
      if(volume.m_volume[iindex])total+=*iweight;
      iindex++;
      iweight++;
    }
  return true;
}

bool VSNSpace::maxWeight(const Volume& volume, Weight& max) const
{
  bool found_one=false;
  if(!m_space.isSpaceEqual(volume.space()))return false;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)
    {
      if(volume.m_volume[iindex])
	if(found_one==false)max=*iweight,found_one=true;
	else if(*iweight>max)max=*iweight;
      iindex++;
      iweight++;
    }
  return found_one;  
}

bool VSNSpace::maxWeight(const Volume& volume, Weight& max, Index& index) const
{
  bool found_one=false;
  if(!m_space.isSpaceEqual(volume.space()))return false;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)
    {
      if(volume.m_volume[iindex])
	if(found_one==false)max=*iweight,index=iindex,found_one=true;
	else if(*iweight>max)max=*iweight,index=iindex;
      iindex++;
      iweight++;
    }
  return found_one;  
}

bool VSNSpace::minWeight(const Volume& volume, Weight& min) const
{
  bool found_one=false;
  if(!m_space.isSpaceEqual(volume.space()))return false;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)
    {
      if(volume.m_volume[iindex])
	if(found_one==false)min=*iweight,found_one=true;
	else if(*iweight<min)min=*iweight;
      iindex++;
      iweight++;
    }
  return found_one;  
}

bool VSNSpace::minWeight(const Volume& volume, Weight& min, Index& index) const
{
  bool found_one=false;
  if(!m_space.isSpaceEqual(volume.space()))return false;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)
    {
      if(volume.m_volume[iindex])
	if(found_one==false)min=*iweight,index=iindex,found_one=true;
	else if(*iweight<min)min=*iweight,index=iindex;
      iindex++;
      iweight++;
    }
  return found_one;  
}

#define LOGLIKE(x) (x)*(x)*::log(x)-gammaln((x)+1.0)

VSNSpace::Weight VSNSpace::logLikelihood() const
{
  Weight loglike = 0;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(iweight!=nweight)loglike += LOGLIKE(*iweight);
  return loglike;
}

bool VSNSpace::logLikelihood(const Volume& volume, Weight& loglike) const
{
  if(!m_space.isSpaceEqual(volume.space()))return false;
  loglike = 0;
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)
    {
      if(volume.m_volume[iindex])loglike += LOGLIKE(*iweight);
      iindex++;
      iweight++;
    }
  return true;
}

void VSNSpace::volumeAll(Volume& volume) const
{
  volume = Volume(m_space);
  std::vector<bool>::iterator iindex = volume.m_volume.begin();
  std::vector<bool>::iterator nindex = volume.m_volume.end();
  while(iindex!=nindex)*iindex = true;
}

bool VSNSpace::volumeSubSpace(const Space& space, Volume& volume) const
{
  if(!m_space.isSpaceCompatible(space))return false;
  volume = Volume(m_space);
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  Point p;
  while(iweight!=nweight)
    {
      m_space.midPointOfIndexUnchecked(iindex, p);
      if(space.isPointCompatible(p))volume.m_volume[iindex] = true;
      iindex++;
    }
  return true;
}

void VSNSpace::volumeAbove(Weight threshold, Volume& volume) const
{
  volume = Volume(m_space);
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)volume.m_volume[iindex++]=(*(iweight++) > threshold);
}

void VSNSpace::volumeBelow(Weight threshold, Volume& volume) const
{
  volume = Volume(m_space);
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)volume.m_volume[iindex++]=(*(iweight++) < threshold);
}

void VSNSpace::volumeAboveEqual(Weight threshold, Volume& volume) const
{
  volume = Volume(m_space);
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)volume.m_volume[iindex++]=(*(iweight++) >= threshold);
}

void VSNSpace::volumeBelowEqual(Weight threshold, Volume& volume) const
{
  volume = Volume(m_space);
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)volume.m_volume[iindex++]=(*(iweight++) <= threshold);
}

// ----------------------------------------------------------------------------
//
// ORDERING
//
// ----------------------------------------------------------------------------

bool VSNSpace::
forwardOrderingSimple(std::vector<Index>& ordering,
		      const VSNSpace& excess, const VSNSpace& variance,
		      Weight min_variance)
{
  if(!excess.m_space.isSpaceEqual(variance.m_space))return false;
  ordering.clear();
  ordering.reserve(variance.m_weight.size());

  const Weight* iexcess = &(*excess.m_weight.begin());
  const Weight* ivariance = &(*variance.m_weight.begin());
  const Weight* nvariance = &(*variance.m_weight.end());
  Index iindex = 0;

  VSNSpace::Volume vol(variance.m_space);
  while(ivariance!=nvariance)
    {
      if(*ivariance>=min_variance && *iexcess)
	ordering.push_back(iindex);

      ivariance++;
      iexcess++;
      iindex++;
    }

  Weight Q_num = 0;
  Weight Q_den = 0;

  std::vector<Index>::iterator zindex = ordering.begin();
  std::vector<Index>::iterator nindex = ordering.end();
  while(zindex != nindex)
    {
      std::vector<Index>::iterator jindex = zindex;
      std::vector<Index>::iterator jindex_max = ordering.end();
      Weight Q_max = 0;

      while(jindex != nindex)
	{
	  if(!vol.isIndexAdjacent(*jindex) && vol.volumeSize() != 0) 
 	    {
 	      jindex++;
 	      continue;
 	    }

	  Weight Q = 
	    (Q_num + excess.m_weight[*jindex])
	    /::sqrt(Q_den + variance.m_weight[*jindex]);

	  if(Q>Q_max && std::isfinite(Q))Q_max = Q, jindex_max = jindex;
	  jindex++;
	}

      if(jindex_max == ordering.end())
	{
	  ordering.erase(zindex,nindex);
	  break;
	}

      vsassert(jindex_max != ordering.end());

      Q_num += excess.m_weight[*jindex_max];
      Q_den += variance.m_weight[*jindex_max];
      vol.setIndexUnchecked(*jindex_max,true);

      iindex = *zindex;
      *zindex = *jindex_max;
      *jindex_max = iindex;
      zindex++;
    }

  return true;
}

bool VSNSpace::
cumulativeOrderingSimple(std::vector<Index>& ordering,
			 const VSNSpace& excess, const VSNSpace& variance,
			 Weight min_variance)
{
  // Find the sequence of elements of every length from 1:N that gives
  // the highest Q value [defined by sum(excess)/sqrt(sum(variance))]
  // or significance. Thanks for MDW for suggesting doing the search
  // from the end back to the start! In actuality this algorithm
  // imposes the requirement that every sequence of length "n" is a
  // subset of the one of length "n+1" which is not necessarily true.

  if(!excess.m_space.isSpaceEqual(variance.m_space))return false;
  ordering.clear();
  ordering.reserve(variance.m_weight.size());

  // Run through all elements finding those that have variance above
  // our threshold and computing the numerator and demominator of Q
  const Weight* iexcess = &(*excess.m_weight.begin());
  const Weight* ivariance = &(*variance.m_weight.begin());
  const Weight* nvariance = &(*variance.m_weight.end());
  Index iindex = 0;

  VSNSpace::Volume vol(variance.m_space);
  while(ivariance!=nvariance)
    {
      if(*ivariance>=min_variance && *iexcess)
	vol.setIndexUnchecked(iindex,true);
      ivariance++;
      iexcess++;
      iindex++;
    }

  Weight Q_num = 0;
  Weight Q_den = 0;
  iindex = 0;
  iexcess = &(*excess.m_weight.begin());
  ivariance = &(*variance.m_weight.begin());
  nvariance = &(*variance.m_weight.end());
  while(ivariance!=nvariance)
    {
      if(*ivariance>=min_variance && *iexcess && !vol.isIndexIsolated(iindex))
	{
	  vol.setIndexUnchecked(iindex,true);
	  ordering.push_back(iindex);
	  Q_num += *iexcess;
	  Q_den += *ivariance;
	}
      else vol.setIndexUnchecked(iindex,false);

      ivariance++;
      iexcess++;
      iindex++;
    }

  // Double loop over all elements to find ordering of elements that
  // maximimizes Q over all partial sums of elements

  // This is essentially the logic of an insertion sort - O(N^2)

  Index* zindex = &(*ordering.begin());
  Index* nindex = &(*ordering.end());
  while(zindex < (nindex-1))
    {
//       Index* jindex_max = zindex;
//       Weight Q_max = 
// 	(Q_num - excess.m_weight[*jindex_max])
// 	/::sqrt(Q_den - variance.m_weight[*jindex_max]);

//       Index* jindex = zindex;
//       jindex++;

      Weight Q_max = 0;
      Index* jindex_max = 0;
      Index* jindex = zindex;

      while(jindex != nindex)
	{
// 	  std::cout << *jindex << " " << vol.isIndexOnEdge(*jindex)
// 		    << " " << Q_max << " " << jindex_max << std::endl;

 	  if(!vol.isIndexOnEdge(*jindex)) 
 	    {
 	      jindex++;
 	      continue;
 	    }

	  Weight Q = 
	    (Q_num - excess.m_weight[*jindex])
	    /::sqrt(Q_den - variance.m_weight[*jindex]);

	  //	  std::cout << "Q: " << Q << std::endl;

	  if(Q>Q_max)Q_max = Q, jindex_max = jindex;
	  jindex++;
	}

      //      std::cout << jindex_max << std::endl;

      vsassert(jindex_max);

      vol.setIndexUnchecked(*jindex_max,false);
      Q_num -= excess.m_weight[*jindex_max];
      Q_den -= variance.m_weight[*jindex_max];

      // Swap index of element which maximizes Q and one at start of
      // index list

      nindex--;
      iindex = *nindex;
      *nindex = *jindex_max;
      *jindex_max = iindex;
    }
  
  return true;
}

// ----------------------------------------------------------------------------
//
// SETTERS AND OPERATORS
//
// ----------------------------------------------------------------------------

void VSNSpace::clear(const Weight weight)
{
  m_weight.clear();
  m_weight.resize(m_space.size(),weight);
}

//! Reset elements in space outside of volume
bool VSNSpace::
clearOutsideVolume(const Volume& volume, const Weight zero_weight)
{
  if(!m_space.isSpaceEqual(volume.space()))return false;
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  unsigned iindex = 0;
  while(iweight!=nweight)
    if(!volume.m_volume[iindex++])*(iweight++)=zero_weight;
    else iweight++;
  return true;
}  

void VSNSpace::clearKeepOnlyAbove(Weight threshold, const Weight zero_weight)
{
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(iweight!=nweight)
    if(*iweight <= threshold)*(iweight++)=zero_weight;
    else iweight++;
}

void VSNSpace::clearKeepOnlyBelow(Weight threshold, const Weight zero_weight)
{
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(iweight!=nweight)
    if(*iweight >= threshold)*(iweight++)=zero_weight;
    else iweight++;
}

VSNSpace& VSNSpace::sqrt()
{
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(iweight!=nweight)
    {
      *iweight = ::sqrt(*iweight);
      iweight++;
    }
  return *this;
}

VSNSpace& VSNSpace::log()
{
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(iweight!=nweight)
    {
      *iweight = ::log(*iweight);
      iweight++;
    }
  return *this;
}

VSNSpace& VSNSpace::normalize()
{
  Weight weight = totalWeight();
  *this *= 1/weight;
  return *this;
}

bool VSNSpace::normalize(const Volume& volume)
{
  Weight weight;
  if(!totalWeight(volume, weight))return false;
  *this *= 1/weight;
  return true;
}

bool VSNSpace::partiallyIntegrate()
{
  const Index nindex = m_space.size();
  Index*const c_array = FASTCALLOC(Index, m_space.ndim);
  for(unsigned idim=0;idim<m_space.ndim;idim++)
    for(Index iindex=0;iindex<nindex;iindex++)
      {
	Cell c(m_space.ndim, c_array);
	m_space.cellOfIndexUnchecked(iindex,c);
	if(c.i[idim]==0)continue;
	c.i[idim]--;
	Index ipartial = m_space.indexOfCellUnchecked(c);
	m_weight[iindex] += m_weight[ipartial];
      }
  FASTFREE(c_array);
  return true;
}

bool VSNSpace::marginalize(unsigned dim)
{
  Space marginal_space = m_space;
  if(!marginal_space.marginalize(dim))return false;
  std::vector<Weight> marginal_weight(marginal_space.size());
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  Index iindex = 0;
  while(iweight!=nweight)
    {
      if(*iweight)
	{
	  Cell c;
	  m_space.cellOfIndexUnchecked(iindex,c);
	  c.marginalize(dim);
	  Index iindex2 = marginal_space.indexOfCellUnchecked(c);
#warning temporary assert here
	  assert(iindex2 < marginal_weight.size());
	  marginal_weight[iindex2] += *iweight;
	}
      iweight++;
      iindex++;
    }
  m_space = marginal_space;
  m_weight = marginal_weight;
  return true;
}

bool VSNSpace::marginalize(unsigned dim, const Volume& volume)
{
  if(!clearOutsideVolume(volume))return false;
  return marginalize(dim);
}

bool VSNSpace::marginalize(unsigned dim, const VSNSpace& dim_weight)
{
  Space marginal_space = m_space;
  if(!marginal_space.marginalize(dim))return false;
  if(!dim_weight.m_space.ndim == 1)return false;
  if(!m_space.axes[dim].isAxisEqual(dim_weight.m_space.axes[0]))return false;
  std::vector<Weight> marginal_weight(marginal_space.size());
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  Index iindex = 0;
  Index*const c_array = FASTCALLOC(Index, m_space.ndim);
  while(iweight!=nweight)
    {
      if(*iweight)
	{
	  Cell c(m_space.ndim, c_array);
	  m_space.cellOfIndexUnchecked(iindex,c);
	  Weight dw = dim_weight.m_weight[c.i[dim]];
	  c.marginalize(dim);
	  Index iindex2 = marginal_space.indexOfCellUnchecked(c);
#warning temporary assert here
	  assert(iindex2 < marginal_weight.size());
	  marginal_weight[iindex2] += *iweight * dw;
	}
      iweight++;
      iindex++;
    }
  m_space = marginal_space;
  m_weight = marginal_weight;
  return true;
}

bool VSNSpace::slice(const std::vector< std::pair<unsigned,Coord> >& dim, 
		     VSNSpace& o)
{

  std::vector<std::pair<unsigned,Index> > slice_index;


  for(std::vector<std::pair<unsigned,Coord> >::const_iterator itr = dim.begin();
      itr != dim.end(); ++itr)
    {
      Axis& a = m_space.axes[itr->first];

      slice_index.
	push_back(std::make_pair(itr->first,a.indexUnchecked(itr->second)));
    }

  return slice(slice_index,o);
}

bool VSNSpace::slice(const std::vector< std::pair<unsigned,Index> >& dim, 
		     VSNSpace& o)
{
  if(m_space.ndim <= dim.size()) return false;

  Space slice_space;

  for(unsigned idim = 0; idim < m_space.ndim; idim++)
    {      
      std::vector<std::pair<unsigned,Index> >::const_iterator itr = dim.begin();
      for(; itr != dim.end(); ++itr)
	  if(itr->first == idim) break;

      if(itr == dim.end())
	slice_space.axes.push_back(m_space.axes[idim]);
    }

  slice_space.ndim = slice_space.axes.size();

  std::vector<Weight> slice_weight(slice_space.size());
  const Index nindex = slice_space.size();
  Cell c = slice_space.cell();
  for(Index iindex = 0; iindex<nindex; iindex++)
    {

      slice_space.cellOfIndexUnchecked(iindex,c);

      for(std::vector<std::pair<unsigned,Index> >::const_iterator itr = 
	    dim.begin(); itr != dim.end(); ++itr)
	c.addDimension(itr->first,itr->second);

      Index iindex2 = m_space.indexOfCellUnchecked(c);
      slice_weight[iindex] = m_weight[iindex2];
    }

  return true;
}

bool VSNSpace::slice(unsigned dim, Index index)
{
  Space slice_space = m_space;
  if(!slice_space.slice(dim))return false;
  if(m_space.axes[dim].nbin < index)return false;
  Index* c_array = FASTCALLOC(Index, m_space.ndim);
  std::vector<Weight> slice_weight(slice_space.size());
  const Index nindex = slice_space.size();
  for(Index iindex = 0; iindex<nindex; iindex++)
    {
      Cell c(m_space.ndim-1, c_array, m_space.ndim);
      slice_space.cellOfIndexUnchecked(iindex,c);
      c.addDimension(dim, index);
      Index iindex2 = m_space.indexOfCellUnchecked(c);
#warning temporary assert here
      assert(iindex2 < m_space.size());
      slice_weight[iindex] = m_weight[iindex2];
    }
  m_space = slice_space;
  m_weight = slice_weight;
  FASTFREE(c_array);
  return true;
}

bool VSNSpace::slice(unsigned dim, Coord coord)
{
  if(m_space.ndim <= dim) return false;
  return slice(dim,m_space.axes[dim].indexUnchecked(coord));
}

bool VSNSpace::project(unsigned dim)
{
  Space projected_space = m_space;
  if(!projected_space.project(dim))return false;
  std::vector<Weight> projected_weight(projected_space.size());
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  Index iindex = 0;
  while(iweight!=nweight)
    {
      if(*iweight)
	{
	  Cell c;
	  m_space.cellOfIndexUnchecked(iindex,c);
	  c.project(dim);
	  Index iindex2 = projected_space.indexOfCellUnchecked(c);
#warning temporary assert here
	  assert(iindex2 < projected_weight.size());
	  projected_weight[iindex2] += *iweight;
	}
      iweight++;
      iindex++;
    }
  m_space = projected_space;
  m_weight = projected_weight;
  return true;
}

bool VSNSpace::project(unsigned dim, const Volume& volume)
{
  if(!clearOutsideVolume(volume))return false;
  return project(dim);
}

bool VSNSpace::project(const std::set<unsigned>& dims)
{
  Space projected_space = m_space;
  if(!projected_space.project(dims))return false;
  std::vector<Weight> projected_weight(projected_space.size());
  const Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  Index iindex = 0;
  while(iweight!=nweight)
    {
      if(*iweight)
	{
	  Cell c;
	  m_space.cellOfIndexUnchecked(iindex,c);
	  c.project(dims);
	  Index iindex2 = projected_space.indexOfCellUnchecked(c);
#warning temporary assert here
	  assert(iindex2 < projected_weight.size());
	  projected_weight[iindex2] += *iweight;
	}
      iweight++;
      iindex++;
    }
  m_space = projected_space;
  m_weight = projected_weight;
  return true;

}

bool VSNSpace::project(const std::set<unsigned>& dims, const Volume& volume)
{
  if(!clearOutsideVolume(volume))return false;
  return project(dims);
}

VSNSpace& VSNSpace::operator *= (Weight w)
{
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  while(iweight!=nweight)*(iweight++) *= w;
  return *this;
}

VSNSpace VSNSpace::operator * (Weight w) const
{
  return VSNSpace(*this) *= w;
}

bool VSNSpace::operator += (const VSNSpace& s)
{
  if(!m_space.isSpaceEqual(s.m_space))return false;
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  const Weight* oweight = &(*s.m_weight.begin());
  while(iweight!=nweight)*(iweight++) += (*oweight++);
  return true;
}

bool VSNSpace::operator -= (const VSNSpace& s)
{
  if(!m_space.isSpaceEqual(s.m_space))return false;
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  const Weight* oweight = &(*s.m_weight.begin());
  while(iweight!=nweight)*(iweight++) -= (*oweight++);
  return true;
}

bool VSNSpace::operator *= (const VSNSpace& s)
{
  if(!m_space.isSpaceEqual(s.m_space))return false;
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  const Weight* oweight = &(*s.m_weight.begin());
  while(iweight!=nweight)*(iweight++) *= (*oweight++);
  return true;
}

bool VSNSpace::operator /= (const VSNSpace& s)
{
  if(!m_space.isSpaceEqual(s.m_space))return false;
  Weight* iweight = &(*m_weight.begin());
  const Weight* nweight = &(*m_weight.end());
  const Weight* oweight = &(*s.m_weight.begin());
  while(iweight!=nweight)*(iweight++) /= (*oweight++);
  return true;
}

#ifndef NOHDF5
bool VSNSpace::load(VSOctaveH5ReaderStruct* reader)
{
  VSNSpaceOctaveH5IO io;
  return io.readHistogram(reader,*this);
}

bool VSNSpace::save(VSOctaveH5WriterStruct* writer) const
{
  VSNSpaceOctaveH5IO io;
  return io.writeHistogram(writer,*this);
}    
#endif

// ----------------------------------------------------------------------------
// 
// PRIVATE AND STATIC FUNCTIONS
//
// ----------------------------------------------------------------------------

double VSNSpace::gammaln(double xx)
{
  static double cof[6]={76.18009172947146, -86.50532032941677,
                        24.01409824083091, -1.231739572450155,
			0.1208650973866179e-2, -0.5395239384953e-5};
  assert(xx>0);
  double x=xx;
  double y=xx;
  double tmp=x+5.5;
  tmp-=(x+0.5)*::log(tmp);
  double ser=1.000000000190015;
  for(unsigned j=0;j<=5;j++) ser+=cof[j]/++y;
  return -tmp+::log(2.5066282746310005*ser/x);
}

// ****************************************************************************
//
// VSNSpaceIO
//
// ****************************************************************************

VSNSpaceIO::~VSNSpaceIO()
{
  // nothing to see here
}

// ****************************************************************************
// 
// STREAM INSERTION
//
// ****************************************************************************

std::ostream& VERITAS::
operator<<(std::ostream& stream, const VSNSpace::Cell& c)
{
  for(unsigned i = 0; i < c.ndim; i++)
    stream << std::setw(10) << c.i[i];
  return stream;
}

std::ostream& VERITAS::
operator<<(std::ostream& stream, const VSNSpace::Axis& a)
{
  stream << a.name << ' '
	 << a.lo_bound << ' ' << a.hi_bound << ' ' << a.nbin;
  return stream;
}

std::ostream& VERITAS::
operator<<(std::ostream& stream, const VSNSpace::Space& s)
{
  for(std::vector<VSNSpace::Axis>::const_iterator iaxis = s.axes.begin();
      iaxis != s.axes.end(); iaxis++)stream << *iaxis << std::endl;
  return stream;
}

std::ostream& VERITAS::
operator<<(std::ostream& stream, const VSNSpace& n)
{
  vsassert(n.ndim() == 2);
  const unsigned nbin1 = n.axis(0).nbin;
  const unsigned nbin2 = n.axis(1).nbin;

  for(unsigned ibin1 = 0; ibin1 < nbin1; ibin1++)
    {
      for(unsigned ibin2 = 0; ibin2 < nbin2; ibin2++)
	{
	  VSNSpace::Cell c(2);

	  c.i[0] = ibin1;      
	  c.i[1] = ibin2;

	  std::cout << std::setw(18) << n[c];
	}
      std::cout << std::endl;
    }

  return stream;
}

// ****************************************************************************
// 
// TEST MAIN
//
// ****************************************************************************

#ifdef TEST_MAIN

#include <iostream>
#include <algorithm>
#include <iterator>

#include <RandomNumbers.hpp>
#include <VSNSpaceOctaveH5IO.hpp>

class Moment1
{
public:
  Moment1(unsigned dim): m_dim(dim) { }
  VSNSpace::Weight operator() (const VSNSpace::Point& p) const 
  { return p.x[m_dim]; }
private:
  unsigned m_dim;
};

class Moment2
{
public:
  Moment2(unsigned dim1, unsigned dim2): m_dim1(dim1), m_dim2(dim2) { }
  VSNSpace::Weight operator() (const VSNSpace::Point& p) const 
  { return p.x[m_dim1]*p.x[m_dim2]; }
private:
  unsigned m_dim1;
  unsigned m_dim2;
};

class Radius2
{
public:
  Radius2() { }
  VSNSpace::Weight operator() (const VSNSpace::Point& p) const 
  { VSNSpace::Weight w=0; for(unsigned i=0;i<p.ndim;i++)w+=p.x[i]*p.x[i];
    return w;}
};

int main(int argc, char**argv)
{
  RandomNumbers rng(RandomNumbers::defaultFilename().c_str());
  
  VSNSpace::Space space(3);
  space.axes[0] = VSNSpace::Axis(-10,10,20,"x");
  space.axes[1] = VSNSpace::Axis(-10,10,20,"y");
  space.axes[2] = VSNSpace::Axis(-10,10,20,"z");

  VSNSpace histo3(space);

  VSNSpace::Point point(3);
  for(unsigned i=0;i<1000000;i++)
    //for(unsigned i=0;i<10000;i++)
    {
      point.x[0] = rng.Normal();
      point.x[1] = rng.Normal();
      point.x[2] = rng.Normal();
      histo3.accumulate(point);
    }

  VSNSpaceOctaveH5IO io;
  io.writeHistogram("test.h5",histo3);
 
  double w = histo3.totalWeight();
  double r = histo3.integrate(Radius2());
  std::cerr << w << std::endl;
  std::cerr << r/w << std::endl;
#if 0
  double mx = histo3.integrate(Moment1(0));
  double my = histo3.integrate(Moment1(1));
  double mz = histo3.integrate(Moment1(2));
  double mxx = histo3.integrate(Moment2(0,0));
  double myy = histo3.integrate(Moment2(1,1));
  double mzz = histo3.integrate(Moment2(2,2));

  std::cerr << w << std::endl;
  std::cerr << mx/w << std::endl;
  std::cerr << my/w << std::endl;
  std::cerr << mz/w << std::endl;
  std::cerr << mxx/w-mx*mx/w/w << std::endl;
  std::cerr << myy/w-my*my/w/w << std::endl;
  std::cerr << mzz/w-mz*mz/w/w << std::endl;
#endif

  VSNSpace histo1;
  io.readHistogram("test.h5",histo1);

  VSNSpace::Volume vol;
  histo1.volumeAbove(30,vol);
  io.writeVolume("test2.h5",vol);
  io.writeVolume("test3.h5",vol,true);  
  
  std::cerr << histo1.space() << std::endl;
  histo1.marginalize(0);
  std::cerr << histo1.space() << std::endl;
  histo1.marginalize(0);
  std::cerr << histo1.space() << std::endl;

  for(unsigned i=0;i<histo1.weight().size();i++)
    std::cout << histo1.weight()[i] << std::endl;
}

#endif
