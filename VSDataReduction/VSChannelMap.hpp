//-*-mode:c++; mode:font-lock;-*-

/*! \file VSChannelMap.hpp

  Encapsulate knowledge of L2 channels

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/09/2007

  $Id: VSChannelMap.hpp,v 3.6 2009/10/29 22:24:57 matthew Exp $

*/

#ifndef VSCHANNELMAP_HPP
#define VSCHANNELMAP_HPP

#include<vector>

#include<VSAssert.hpp>
#include<VSTime.hpp>
#include<VSDataReductionTypes.hpp>

#define __NUMOF(x) (sizeof(x)/sizeof(*x))

namespace VERITAS
{

  class VSChannelMap
  {
  public:
    VSChannelMap(const VSTime& t): fTime(t) { }

    unsigned nscope() const { return 4; }

    pos_type posForScope(unsigned iscope) const
    {
      vsassert(iscope < nscope());

  /*

    Logbook URL         : http://veritas.sao.arizona.edu/private/elog/VBC/3677

    =================================

    We surveyed yesterday and measured the locations of the four telescopes. 
    The coordinates are given in a frame oriented with north as positive Y 
    and east as positive X. The origin is located at the mean of the distances 
    between telescopes. The units are in meters with an estimated uncertainty 
    of 0.1m.

             X        Y       Z    Distance to T1
    T1   -38.4    -23.6    -5.6
    T2    46.1    -49.9    -1.3      88.6
    T3    29.2     61.4     5.5     109.2
    T4   -36.9     12.1     1.3      36.4

    (Updated as of 5:40 PM PST)

  scope_pos.push_back(pos_type(-38.4, -23.6,  -5.6));
  scope_pos.push_back(pos_type( 46.1, -49.9,  -1.3));
  scope_pos.push_back(pos_type( 29.2,  61.4,   5.5));
  scope_pos.push_back(pos_type(-36.9,  12.1,   1.3));

  */

  /*

  Logbook URL         : http://veritas.sao.arizona.edu/private/elog/VBC/3679

  Jack and I spent another three, three and half hours on Tuesday
  checking these numbers again.  There are only very small differences
  between the values we derived and those above which Jen, Tim and Dan
  produced.

  The main differences are:

  1) I remeasured the bearing (angle) to the MMT before and after each
  telescope "shot" for a total of 9 times.  The std dev of the measured
  angle was ~7" (22cm at the MMT)

  2) We were able to just barely see Baboquivri (BBQ hereafter) and a
  third reference peak (Sardini)

  3) I calculated the expected bearing to the reference points (ie, MMT,
  BBQ and Sardini) using a "rhumbline" rather than a great-circle

  There are some caveats:
  
  a) I am not sure about the rhumbline vs g-c difference (it is very
  very small in any case)

  b) Sardini is a bit of a featureless wide bump so I couldn't be sure I
  always was sighting in at the same point.  For comparison, I took
  three bearings on BBQ after slewing the total-station left and right
  by about 60 degrees (to re/check the level) and they differed by 1"
  and 4" of arc.

  c) "b" notwithstanding, the angle I calculate from the MMT, through
  North and over to BBQ, compared to that which I measured is off by
  132"!  I haven't resolved this discrepancy but as you will see below,
  it doesn't make diddly-squat difference.

  d) To call what we did "surveying" is a pretty gross exaggeration.  We
  did rent a v. nice total-station for two days but good equipment can't
  make up for operator ignorance.

  When the weather clears up I can use our theodolite to check the
  overall angle.  Its been too cloudy to do any star shots.

  Anyway, here is the result (same convention as the one by Jen et al):

  Using the MMT as the angle (to North) reference):

           X(+E)   Y(+N)   Z
  <T1>    -37.99  -23.30  -5.20
  <T2A>    45.13  -49.14   -.95
  <T2B>    45.15  -49.14   -.93
  <T3>     28.93   60.71   4.51
  <T4A>   -36.12   11.72   1.63
  <T4B>   -36.04   11.73   2.31

  ("A" and "B" denote the two different measurements of T2 and T4 (ie
  T2->T1->T2->T4->T3->T4).

  or BBQ:

           X(+E)   Y(+N)   Z
  <T1>    -38.00  -23.32  -5.20
  <T2A>    45.13  -49.11   -.95
  <T2B>    45.16  -49.12   -.93
  <T3>     28.86   60.73   4.51
  <T4A>   -36.15   11.69   1.63
  <T4B>   -36.07   11.70   1.63

  As you can see the difference between using the MMT or BBQ is a couple
  of cm maximum which is just about the difference found when repeating
  the measurements on T2 and T4 (though I notice that the X value of T4
  "moved" by 8cm which I suspect was operator error).

  */

      if(fTime < VATime("2009-09-01 12:00:00"))
	{
	  static const double scope_pos[][3] =
	    { {-38.00, -23.32,  -5.20},
	      { 45.13, -49.11,  -0.95},
	      { 28.86,  60.73,   4.51},
	      {-36.15,  11.69,   1.63} };
	  return pos_type(scope_pos[iscope][0],
			  scope_pos[iscope][1],
			  scope_pos[iscope][2]);
	}
      else 
	{
	  static const double scope_pos[][3] =
	    { {135.48, -8.61,    7.23},
	      { 45.13, -49.11,  -0.95},
	      { 28.86,  60.73,   4.51},
	      {-36.15,  11.69,   1.63} };
	  return pos_type(scope_pos[iscope][0],
			  scope_pos[iscope][1],
			  scope_pos[iscope][2]);
	}
      
    }

    unsigned nchanPerBoard() const { return 10; }

    unsigned nchanForScope(unsigned iscope) const { return (iscope<4)?500:0; }

    unsigned ncratesForScope(unsigned iscope) const
    {
      if(iscope<4)return 4;
      else return 0;
    }

    unsigned nboardsForCrate(unsigned iscope, unsigned icrate) const
    {
      static const unsigned nboard[] = { 13, 12, 13, 12 };
      if(icrate<__NUMOF(nboard))return nboard[icrate];
      return 0;
    }

    unsigned crateForChannel(unsigned iscope, unsigned ichan) const
    {
      const unsigned ncrate = ncratesForScope(iscope);
      for(unsigned icrate=0;icrate<ncrate;icrate++)
	{
	  const unsigned nboard = nboardsForCrate(iscope,icrate);
	  const unsigned nchan = nboard*nchanPerBoard();
	  if(ichan<nchan)return icrate;
	  ichan -= nchan;
	}
      return ncrate;
    }

    unsigned l2ChannelForCrate(unsigned iscope, unsigned icrate) const
    {
      if(fTime < VATime("2007-09-11 12:00:00"))
	{
	  static const unsigned l2_chan[] = { 128, 249, 259, 498 };
	  if(icrate<ncratesForScope(iscope))return l2_chan[icrate];
	}
      else if(fTime >= VATime("2007-09-11 12:00:00") &&
	      fTime < VATime("2009-10-16 02:44:00"))
	{
	  static const unsigned l2_chan[][4] = { { 110, 249, 255, 404 },
						 { 128, 173, 259, 498 },
						 {  37, 159, 319, 499 },
						 {  99, 214, 259, 499 } };
	  if((iscope<nscope())&&(icrate<ncratesForScope(iscope)))
	    return l2_chan[iscope][icrate];
	}
      else
	{
	  static const unsigned l2_chan[][4] = { { 110, 249, 255, 404 },
						 { 128, 173, 259, 498 },
						 {  37, 159, 319, 499 },
						 {  99, 214, 333, 499 } };
	  if((iscope<nscope())&&(icrate<ncratesForScope(iscope)))
	    return l2_chan[iscope][icrate];
	}
      return nchanForScope(iscope);
    }

    unsigned l2ChannelForChannel(unsigned iscope, unsigned ichan) const
    {
      return l2ChannelForCrate(iscope,crateForChannel(iscope,ichan));
    }

    unsigned nfir() const { return 5; }
    
  private:
    VSTime fTime;    
  };
}

#undef __NUMOF

#endif

