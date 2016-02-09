//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCORSIKAEvent.cpp
  Simple CORSIKA Cerenkov File Dispatcher and Vistor

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/17/2005
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "VSCORSIKAEvent.hpp"

#include <initial.h>      /* This file includes others as required. */
#include <io_basic.h>     /* This file includes others as required. */
#include <mc_tel.h>

using namespace VERITAS;

/* ------------------- line_point_distance --------------------- */
/**
 *  Distance between a straight line and a point in space
 *
 *  From: sim_skeleton.c by Konrad Bernloehr
 *
 *  Copyright by Konrad Bernloehr (1997, 1999). All rights reserved.
 *  This file may be modified but all modified version must be
 *  declared as modified and by whom they were modified.
 *
 *  @param  x1,y1,z1  reference point on the line
 *  @param  cx,cy,cz  direction cosines of the line
 *  @param  x,y,z     point in space
 *
 *  @return distance
 *
*/

double line_point_distance (double x1, double y1, double z1, 
   double cx, double cy, double cz,
   double x, double y, double z)
{
   double a, a1, a2, a3, b;
   
   a1 = (y-y1)*cz - (z-z1)*cy;
   a2 = (z-z1)*cx - (x-x1)*cz;
   a3 = (x-x1)*cy - (y-y1)*cx;
   a  = a1*a1 + a2*a2 + a3*a3;
   b = cx*cx + cy*cy + cz*cz;
   if ( a<0. || b<= 0. )
      return -1;
   return sqrt(a/b);
}

VSCORSIKAEventVisitor::~VSCORSIKAEventVisitor()
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
registerDispatcher(VSCORSIKAEventDispatcherStop* dispatch_stop)
{
  fDispatchStop = dispatch_stop;
}

void VSCORSIKAEventVisitor::
visitFile(const char* filename)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
leaveFile()
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitRun(const VSCORSIKARunParameters& param)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitRunExtra(const VSCORSIKAExtraRunParameters& param)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
leaveRun()
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitInputConfigEntry(const char* line)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitArraySpec(const VSCORSIKAArraySpec& array)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitTelescopeSpec(const VSCORSIKATelescopeSpec& scope)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitEvent(const VSCORSIKAEvent& event, bool& veto)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
leaveEvent(bool veto)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitEventUse(const VSCORSIKAEventUse& use, bool& veto)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
leaveEventUse(bool veto)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope, bool& veto)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
leaveTelescopeEvent(bool veto)
{
  // nothing to see here
}

void VSCORSIKAEventVisitor::
visitPhotonBunch(const VSCORSIKAPhotonBunch& bunch)
{
  // nothing to see here
}

VSCORSIKAEventDispatcherStop::~VSCORSIKAEventDispatcherStop()
{
  // nothing to see here
}

VSCORSIKAEventDispatcher::
VSCORSIKAEventDispatcher(VSCORSIKAEventVisitor* visitor):
  fVisitor(visitor), fStopProcessingFlag(false), 
  fExtraParametersDispatched(false), fDispatchEmptyBunches(false)
{
  fVisitor->registerDispatcher(this);
}

VSCORSIKAEventDispatcher::~VSCORSIKAEventDispatcher()
{
  // nothing to see here
}

unsigned VSCORSIKAEventDispatcher::processFile(const char* filename)
{
  fStopProcessingFlag = false;
  fExtraParametersDispatched = false;

  /* I/O buffer for input needed */
  IO_BUFFER *iobuf = allocate_io_buffer(0);
  if (iobuf == NULL)
    {
      std::cerr << "Input I/O buffer not allocated" << std::endl;
      return 0;;
    }
  iobuf->max_length = 500000000;

  FILE* data_file_ptr = fopen(filename,"r");
  if (data_file_ptr == NULL)
    {
      perror(filename);
      return 0;
    }
  iobuf->input_file = data_file_ptr;  

  fVisitor->visitFile(filename);

  VSCORSIKARunParameters run_param;
  VSCORSIKAArraySpec array_spec;
  VSCORSIKAEvent event;
  
  int num_array_uses = VSCORSIKAEVENT_MAX_ARRAY;
  double array_use_toff = 0;
  double* array_use_xoff = new double[VSCORSIKAEVENT_MAX_ARRAY];
  double* array_use_yoff = new double[VSCORSIKAEVENT_MAX_ARRAY];
  
  int num_tel = VSCORSIKAEVENT_MAX_TEL;
  double* xtel = new double[VSCORSIKAEVENT_MAX_TEL]; 
  double* ytel = new double[VSCORSIKAEVENT_MAX_TEL]; 
  double* ztel = new double[VSCORSIKAEVENT_MAX_TEL]; 
  double* rtel = new double[VSCORSIKAEVENT_MAX_TEL]; 

  bunch* bunches = new bunch[VSCORSIKAEVENT_MAX_BUNCHES];
  
  bool vetoEvent;
  bool vetoEventUse;
  bool vetoScope;

  unsigned packet_num=0;
  while(!fStopProcessingFlag)
    {
      IO_ITEM_HEADER block_header;

      /* Find and read the next block of data. */
      /* In case of problems with the data, just give up. */
      if ( find_io_block(iobuf,&block_header) != 0 )break;
      if ( read_io_block(iobuf,&block_header) != 0 )break;
      packet_num++;

      switch ( block_header.type )
	{
	case IO_TYPE_MC_RUNH:
	  /* CORSIKA run header */
	  {
	    real runh[273];
	    read_tel_block(iobuf,IO_TYPE_MC_RUNH,runh,273);

            run_param.run_num          = (unsigned)roundf(runh[1]);
            run_param.run_date         = (unsigned)roundf(runh[2]);
            run_param.run_version      = runh[3];

	    unsigned nht               = (unsigned)roundf(runh[4]);
	    if (nht>0 && nht <= 10)
	      run_param.obs_height     = runh[4+nht]*0.01;
	    else 
	      run_param.obs_height     = -1.0;

	    run_param.sample_e_min     = runh[16]*0.001;
	    run_param.sample_e_max     = runh[17]*0.001;
	    run_param.sample_slope     = runh[15];

	    fVisitor->visitRun(run_param);
	  }
	  break;

	case IO_TYPE_MC_INPUTCFG:
	  /* CORSIKA inputs */
	  {
	    struct linked_string *xl = 
	      static_cast<struct linked_string *>
	      (calloc(1,sizeof(linked_string)));
	    struct linked_string *xln = 0;
	    xl->text=NULL;
	    xl->next=NULL;
 	    read_input_lines(iobuf,xl);
	
	    while(xl!=NULL)
	      {
		fVisitor->visitInputConfigEntry(xl->text);
		free(xl->text);
		xln = xl->next;
		free(xl);
		xl = xln;
	      }
	  }
	  break;

	case IO_TYPE_MC_TELPOS:
	  /* Telescope positions (relative positions in array) */
	  {
	    read_tel_pos(iobuf, VSCORSIKAEVENT_MAX_TEL,
			 &num_tel, xtel, ytel, ztel, rtel);

	    array_spec.array_max_tel   = VSCORSIKAEVENT_MAX_TEL;
	    array_spec.array_num_tel   = (unsigned)num_tel;
	    fVisitor->visitArraySpec(array_spec);

	    for(int itel=0; itel<num_tel; itel++)
	      {
		xtel[itel] *= 0.01;
		ytel[itel] *= 0.01;
		ztel[itel] *= 0.01;
		rtel[itel] *= 0.01;

		VSCORSIKATelescopeSpec scope_spec;
		scope_spec.scope_num   = unsigned(itel);
		scope_spec.scope_x     = xtel[itel];
		scope_spec.scope_y     = ytel[itel];
		scope_spec.scope_z     = ztel[itel];
		scope_spec.scope_r     = rtel[itel];
		fVisitor->visitTelescopeSpec(scope_spec);
	      }
	  }
	  break;

	case IO_TYPE_MC_EVTH:
	  /* CORSIKA event header */
	  {
	    real evth[273];
            read_tel_block(iobuf,IO_TYPE_MC_EVTH,evth,273);

	    if(!fExtraParametersDispatched)
	      {
		VSCORSIKAExtraRunParameters extra;
		extra.pri_start_alt    = evth[4];
		extra.magnetic_bx      = evth[70];
		extra.magnetic_bz      = evth[71];
		extra.sample_theta_lo  = evth[80];
		extra.sample_theta_hi  = evth[81];
		extra.sample_phi_lo    = evth[82];
		extra.sample_phi_hi    = evth[83];
		extra.sample_cone_lo   = evth[152];
		extra.sample_cone_hi   = evth[153];
		extra.cerenk_bunch     = (unsigned)roundf(evth[84]);
		extra.cerenk_event_use = (unsigned)roundf(evth[97]);
		extra.cerenk_lambda_lo = evth[95]; 
		extra.cerenk_lambda_hi = evth[96];
		extra.cerenk_atmo      = ((unsigned)roundf(evth[76]))&0x0010;
		fVisitor->visitRunExtra(extra);
		fExtraParametersDispatched=true;
	      }

            event.event_num            = (unsigned)roundf(evth[1]);
            event.pri_type             = (unsigned)roundf(evth[2]);
	    event.pri_energy           = evth[3]*0.001;
	    event.pri_interact_alt     = evth[6]*0.01;
	    event.pri_px               = evth[7]*0.001;
	    event.pri_py               = evth[8]*0.001;
	    event.pri_pz               = -evth[9]*0.001;

	    double az = fmod(3.0*M_PI-evth[11]+evth[92],2.0*M_PI)*(180.0/M_PI);
	    double el = (M_PI_2-evth[10])*(180.0/M_PI);

	    event.pri_azimuth          = az;
	    event.pri_elevation        = el;
	    
	    vetoEvent = false;
	    fVisitor->visitEvent(event,vetoEvent);

#if 0
	    float toffset = 
	      (event.pri_interact_alt-run_param.obs_height)/
	      (cos(evth[10])*0.299792458);
#endif
	  }
	  break;

	case IO_TYPE_MC_TELOFF:
	  /* Offsets of telescope array instances for the following event */
	  {
	    read_tel_offset(iobuf, VSCORSIKAEVENT_MAX_ARRAY,
			    &num_array_uses, &array_use_toff,
			    array_use_xoff, array_use_yoff);
	  }
	  break;

	case IO_TYPE_MC_TELARRAY:
	  /* Photons for a complete array (one of perhaps many instances) */
	  {
	    int event_use_num;
	    IO_ITEM_HEADER item_header;
            begin_read_tel_array(iobuf, &item_header, &event_use_num);

	    float cx = event.pri_px/event.pri_energy;
	    float cy = event.pri_py/event.pri_energy;
	    float cz = event.pri_pz/event.pri_energy;

	    VSCORSIKAEventUse event_use;
	    event_use.use_num          = event_use_num;
	    event_use.pri_xcore        = -array_use_xoff[event_use_num]*0.01;
	    event_use.pri_ycore        = -array_use_yoff[event_use_num]*0.01;
	    event_use.pri_impact       =
	      line_point_distance(event_use.pri_xcore,
				  event_use.pri_ycore,
				  0.0,
				  cx, cy, cz,
				  0.0, 0.0, 0.0);

	    vetoEventUse=vetoEvent;
	    if(!vetoEvent)
	      fVisitor->visitEventUse(event_use, vetoEventUse);

	    for(int itc=0; itc<num_tel; itc++)
	      {
		IO_ITEM_HEADER sub_item_header;
		sub_item_header.type = IO_TYPE_MC_PHOTONS;
		if(search_sub_item(iobuf,&item_header,&sub_item_header) < 0)
                  break;

		int ieventuse;
		int itel;
		double photons;
		int nbunches;
		if (read_tel_photons(iobuf,VSCORSIKAEVENT_MAX_BUNCHES,
				     &ieventuse, &itel, 
				     &photons, bunches, &nbunches) < 0)
		  {
		    fprintf(stderr,"Error reading %d photon bunches\n",
			    nbunches);
		    continue;
		  }
		
		assert(ieventuse == event_use_num);

		if((itel >= (int)array_spec.array_max_tel)||(itel < 0))
		  {
		    fprintf(stderr,
			    "Cannot process data for telescope #%d "
			    "because only %d are configured.\n",
			    itel+1,array_spec.array_max_tel);
		    continue;
		  }

		VSCORSIKATelescopeEvent scope_event;
		scope_event.scope_num  = itel;
		scope_event.rel_xcore  = event_use.pri_xcore-xtel[itel];
		scope_event.rel_ycore  = event_use.pri_ycore-ytel[itel];
		scope_event.rel_zcore  = -ztel[itel];
		scope_event.rel_impact = 
		  line_point_distance(event_use.pri_xcore,
				      event_use.pri_ycore,
				      0.0,
				      cx, cy, cz, 
				      xtel[itel], ytel[itel], ztel[itel]);
		scope_event.num_bunch  = nbunches;
		scope_event.num_ph     = photons;

		vetoScope = vetoEventUse;
		if(!vetoEventUse)
		  fVisitor->visitTelescopeEvent(scope_event, vetoScope);

		if(!vetoScope)
		  for(int ibunch=0; ibunch<nbunches; ibunch++)
		    {
		      float cx       = bunches[ibunch].cx;
		      float cy       = bunches[ibunch].cy;
		      float cz       = -1.0*sqrt(1.0-cx*cx-cy*cy);
		      
#if 0		    
		      float airmass  = -1./cz;
		      float dist     = 
			(bunches[ibunch].zem-run_param.obs_height)*airmass;
		      float tel_dist = ztel[itel]*airmass;
#endif		    
		      VSCORSIKAPhotonBunch bunch;
		      bunch.ph_count     = bunches[ibunch].photons;
		      if((!fDispatchEmptyBunches)&&(bunch.ph_count==0))
			continue;
#if 0
		      bunch.ph_rel_x     = bunches[ibunch].x + tel_dist*cx;
		      bunch.ph_rel_y     = bunches[ibunch].y + tel_dist*cy;
#else
		      bunch.ph_rel_x     = bunches[ibunch].x*0.01;
		      bunch.ph_rel_y     = bunches[ibunch].y*0.01;
#endif
		      bunch.ph_rel_impact =   
			line_point_distance(bunch.ph_rel_x,
					    bunch.ph_rel_y,
					    0.0,
					    cx, cy, cz, 
					    0, 0, 0);
		      
		      bunch.ph_cosine_x  = cx;
		      bunch.ph_cosine_y  = cy;
		      bunch.ph_cosine_z  = cz;
		      bunch.ph_lambda    = bunches[ibunch].lambda;
		      bunch.ph_time      = bunches[ibunch].ctime;
		      bunch.ph_height    = bunches[ibunch].zem*0.01;
		      fVisitor->visitPhotonBunch(bunch);
		    }
		
		if(!vetoEventUse)
		  fVisitor->leaveTelescopeEvent(vetoScope);
	      }

            end_read_tel_array(iobuf, &item_header);

	    if(!vetoEvent)
	      fVisitor->leaveEventUse(vetoEventUse);
	  }
	  break;

	case IO_TYPE_MC_EVTE:
	  /* CORSIKA event trailer */
	  fVisitor->leaveEvent(vetoEvent);
	  break;

	case IO_TYPE_MC_RUNE:
	  /* CORSIKA run trailer */
	  fVisitor->leaveRun();
	  break;

	default:
	  /* Unknown / any other material */
	  break;
	}
      
    }
  
  delete[] bunches;
  
  delete[] xtel;
  delete[] ytel;
  delete[] ztel;
  delete[] rtel;
  
  delete[] array_use_xoff;
  delete[] array_use_yoff;

  // Set this to false before we leave in case we are processing nested files
  fStopProcessingFlag = false;

  fVisitor->leaveFile();
  fclose(iobuf->input_file);
  iobuf->input_file = NULL;

  return packet_num;
}

void VSCORSIKAEventDispatcher::stopProcessingFile()
{
  fStopProcessingFlag = true;
}

std::ostream& operator <<(std::ostream& s, 
			  const VSCORSIKARunParameters& x)
{
  s << x.run_num << ' ' 
    << x.run_date << ' '
    << x.run_version << ' '
    << x.obs_height << ' '
    << x.sample_e_min << ' '
    << x.sample_e_max << ' ' 
    << x.sample_slope;
  return s;
}

std::ostream& operator <<(std::ostream& s, 
			  const VSCORSIKAExtraRunParameters& x)
{
  s << x.pri_start_alt << ' '
    << x.magnetic_bx << ' '
    << x.magnetic_bz << ' '
    << x.sample_theta_lo << ' ' 
    << x.sample_theta_hi << ' ' 
    << x.sample_phi_lo << ' '
    << x.sample_phi_hi << ' '
    << x.sample_cone_lo << ' '
    << x.sample_cone_hi << ' '
    << x.cerenk_bunch << ' ' 
    << x.cerenk_event_use << ' '
    << x.cerenk_lambda_lo << ' '
    << x.cerenk_lambda_hi << ' '
    << x.cerenk_atmo;
  return s;
}

std::ostream& operator <<(std::ostream& s, 
			  const VSCORSIKAArraySpec& x)
{
  s << x.array_max_tel << ' '
    << x.array_num_tel;
  return s;
}

std::ostream& operator <<(std::ostream& s, 
			  const VSCORSIKATelescopeSpec& x)
{
  s << x.scope_num << ' ' 
    << x.scope_x << ' ' 
    << x.scope_y << ' ' 
    << x.scope_z << ' ' 
    << x.scope_r;
  return s;
}

std::ostream& operator <<(std::ostream& s, 
			  const VSCORSIKAEvent& x)
{
  s << x.event_num << ' '
    << x.pri_type << ' '
    << x.pri_energy << ' '
    << x.pri_interact_alt << ' '
    << x.pri_px << ' '
    << x.pri_py << ' '
    << x.pri_pz << ' '
    << x.pri_azimuth << ' '
    << x.pri_elevation;
  return s;
}

std::ostream& operator <<(std::ostream& s, 
			  const VSCORSIKAEventUse& x)
{
  s << x.use_num << ' '
    << x.pri_xcore << ' '
    << x.pri_ycore << ' '
    << x.pri_impact;
  return s;
}

std::ostream& operator <<(std::ostream& s, 
			  const VSCORSIKATelescopeEvent& x)
{
  s << x.scope_num << ' '
    << x.rel_xcore << ' '
    << x.rel_ycore << ' '
    << x.rel_zcore << ' '
    << x.rel_impact << ' '
    << x.num_bunch << ' '
    << x.num_ph;
  return s;
}

std::ostream& operator <<(std::ostream& s, 
			  const VSCORSIKAPhotonBunch& x)
{
  s << x.ph_count << ' '
    << x.ph_rel_x << ' '
    << x.ph_rel_y << ' '
    << x.ph_rel_impact << ' '
    << x.ph_cosine_x << ' '
    << x.ph_cosine_y << ' '
    << x.ph_cosine_z << ' '
    << x.ph_lambda << ' '
    << x.ph_time << ' '
    << x.ph_height;
  return s;
}

void VSCORSIKAFileDumper::
visitFile(const char* filename)
{
  fStream << "FILE: " << filename << std::endl;
}

void VSCORSIKAFileDumper::
visitRun(const VSCORSIKARunParameters& param)
{
  fStream << "RUNPARAM: " << param << std::endl;
}

void VSCORSIKAFileDumper::
visitRunExtra(const VSCORSIKAExtraRunParameters& param)
{
  fStream << "EXTRAPARAM: " << param << std::endl;
}

void VSCORSIKAFileDumper::
visitInputConfigEntry(const char* line)
{
  fStream << "CONFIGLINE: " << line << std::endl;
}

void VSCORSIKAFileDumper::
visitArraySpec(const VSCORSIKAArraySpec& arrayspec)
{
  fStream << "ARRAYSPEC: " << arrayspec << std::endl;
}

void VSCORSIKAFileDumper::
visitTelescopeSpec(const VSCORSIKATelescopeSpec& scopespec)
{
  fStream << "SCOPESPEC: " << scopespec << std::endl;
}

void VSCORSIKAFileDumper::
visitEvent(const VSCORSIKAEvent& event, bool& veto)
{
  fStream << "EVENT: " << event << std::endl;
}

void VSCORSIKAFileDumper::
visitEventUse(const VSCORSIKAEventUse& use, bool& veto)
{
  fStream << "EVENTUSE: " << use << std::endl;
}

void VSCORSIKAFileDumper::
visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope, bool& veto)
{
  fStream << "SCOPEEVENT: " << scope << std::endl;
}

void VSCORSIKAFileDumper::
visitPhotonBunch(const VSCORSIKAPhotonBunch& bunch)
{
  std::cout << "PHOTON: " << bunch << std::endl;
}
