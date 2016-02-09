//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCentralizedDBAccessData.hpp

  Centralize all access to the database and provide asynchnonous data
  retrieval through thread

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       12/12/2006

  $Id: VSCentralizedDBAccessData.hpp,v 3.5 2008/04/14 22:29:59 sfegan Exp $

*/

#ifndef VSCENTRALIZEDDBACCESSDATA_HPP
#define VSCENTRALIZEDDBACCESSDATA_HPP

#include <VSTime.hpp>
#include <VSOctaveIO.hpp>
#include <CorrectionParameters.h>

namespace VERITAS
{
  namespace VSCentralizedDBAccessData
  {

    struct RunInfoDatum
    {
      RunInfoDatum():
	run_id(), run_type(), observing_mode(), run_status(),
	db_start_time(), db_end_time(), data_start_time(), data_end_time(),
	duration_msec(), weather(), config_mask(), pointing_mode(),
	trigger_config(), trigger_multiplicity(), trigger_coincidence(),
	/* offset_ra(), offset_dec(), */ offset_distance(), offset_angle(),
	source_id() { /* nothing to see here */ }

      uint32_t          run_id;
      std::string       run_type;
      std::string       observing_mode;
      std::string       run_status;
      VSTime            db_start_time;
      VSTime            db_end_time;
      VSTime            data_start_time;
      VSTime            data_end_time;
      uint32_t          duration_msec;
      std::string       weather;
      uint32_t          config_mask;
      std::string       pointing_mode;
      std::string       trigger_config;
      uint32_t          trigger_multiplicity;
      float             trigger_coincidence;
      //float             offset_ra;
      //float             offset_dec;
      float             offset_distance;
      float             offset_angle;
      std::string       source_id;            
    };

    struct CorrectionParametersDatum: public SEphem::CorrectionParameters
    {
      CorrectionParametersDatum(): 
	SEphem::CorrectionParameters(),
        db_start_time(), db_end_time(), db_end_time_is_null() 
      { /* nothing to see here */ }
	
      VSTime            db_start_time;
      VSTime            db_end_time;
      bool              db_end_time_is_null;
      std::string       comment;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	// Elements from this class
	H5_ADDSIMPLECOMPOSITE(c,CorrectionParametersDatum,db_start_time);
	H5_ADDSIMPLECOMPOSITE(c,CorrectionParametersDatum,db_end_time);
	H5_ADDMEMBER(c,CorrectionParametersDatum,db_end_time_is_null);
	H5_ADDMEMBER(c,CorrectionParametersDatum,comment);

	// Elements from base class (SEphem::CorrectionParameters)
	H5_ADDMEMBER(c,CorrectionParametersDatum,enable_offsets);
	H5_ADDMEMBER(c,CorrectionParametersDatum,enable_corrections);
	H5_ADDMEMBER(c,CorrectionParametersDatum,enable_vff);
	H5_ADDMEMBER(c,CorrectionParametersDatum,az_ratio);
	H5_ADDMEMBER(c,CorrectionParametersDatum,el_ratio);
	H5_ADDMEMBER(c,CorrectionParametersDatum,az_offset);
	H5_ADDMEMBER(c,CorrectionParametersDatum,el_offset);
	H5_ADDMEMBER(c,CorrectionParametersDatum,az_ns);
	H5_ADDMEMBER(c,CorrectionParametersDatum,az_ew);
	H5_ADDMEMBER(c,CorrectionParametersDatum,el_udew);
	H5_ADDMEMBER(c,CorrectionParametersDatum,fp_az);
	H5_ADDMEMBER(c,CorrectionParametersDatum,flex_el_A);
	H5_ADDMEMBER(c,CorrectionParametersDatum,flex_el_B);
	H5_ADDMEMBER(c,CorrectionParametersDatum,el_pos_vff_s);
	H5_ADDMEMBER(c,CorrectionParametersDatum,el_pos_vff_t);
	H5_ADDMEMBER(c,CorrectionParametersDatum,el_neg_vff_s);
	H5_ADDMEMBER(c,CorrectionParametersDatum,el_neg_vff_t);
	H5_ADDMEMBER(c,CorrectionParametersDatum,az_pos_vff_s);
	H5_ADDMEMBER(c,CorrectionParametersDatum,az_pos_vff_t);
	H5_ADDMEMBER(c,CorrectionParametersDatum,az_neg_vff_s);
	H5_ADDMEMBER(c,CorrectionParametersDatum,az_neg_vff_t);
      }
    };

    struct TrackingTargetDatum
    {
      TrackingTargetDatum():
        db_start_time(), db_end_time(), db_end_time_is_null(),
	target_type(),
	target_angle1_rad(), target_angle2_rad(), target_epoch(), 
	has_extended_target_info(), target_name(),
        tracking_mode(), tracking_param1(), tracking_param2(),
        pointing_mode(), pointing_array_mask(), pointing_param()
      { /* nothing to see here */ }

      VSTime            db_start_time;
      VSTime            db_end_time;
      bool              db_end_time_is_null;
      std::string       target_type;
      double            target_angle1_rad;
      double            target_angle2_rad;
      double            target_epoch;
      bool              has_extended_target_info;
      std::string       target_name;
      std::string       tracking_mode;
      double            tracking_param1;
      double            tracking_param2;
      std::string       pointing_mode;
      uint32_t          pointing_array_mask;
      double            pointing_param;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDSIMPLECOMPOSITE(c,TrackingTargetDatum,db_start_time);
	H5_ADDSIMPLECOMPOSITE(c,TrackingTargetDatum,db_end_time);
	H5_ADDMEMBER(c,TrackingTargetDatum,db_end_time_is_null);
        H5_ADDMEMBER(c,TrackingTargetDatum,target_type);
        H5_ADDMEMBER(c,TrackingTargetDatum,target_angle1_rad);
        H5_ADDMEMBER(c,TrackingTargetDatum,target_angle2_rad);
        H5_ADDMEMBER(c,TrackingTargetDatum,target_epoch);
        H5_ADDMEMBER(c,TrackingTargetDatum,has_extended_target_info);
        H5_ADDMEMBER(c,TrackingTargetDatum,target_name);
        H5_ADDMEMBER(c,TrackingTargetDatum,tracking_mode);
        H5_ADDMEMBER(c,TrackingTargetDatum,tracking_param1);
        H5_ADDMEMBER(c,TrackingTargetDatum,tracking_param2);
        H5_ADDMEMBER(c,TrackingTargetDatum,pointing_mode);
        H5_ADDMEMBER(c,TrackingTargetDatum,pointing_array_mask);
        H5_ADDMEMBER(c,TrackingTargetDatum,pointing_param);
      }
    };

    struct PowerDatum
    {
      PowerDatum(): db_start_time(), db_end_time(), ichan(), power_on()
      { /* nothing to see here */ }

      VSTime            db_start_time;
      VSTime            db_end_time;
      unsigned          ichan;
      bool              power_on;
    };

    struct HVStatusDatum
    {
      HVStatusDatum(): timestamp(), ichan(), voltage(), current()
      { /* nothing to see here */ }

      VSTime            timestamp;
      unsigned          ichan;
      float             voltage;
      float             current;
    };

    struct L1RateDatum
    {
      L1RateDatum(): timestamp(), ichan(), rate()
      { /* nothing to see here */ }

      VSTime            timestamp;
      unsigned          ichan;
      float             rate;
    };

    struct FIRDatum
    {
      FIRDatum(): timestamp(), ambient(), sky()
      { /* nothing to see here */ }

      VSTime            timestamp;
      float             ambient;
      float             sky;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDSIMPLECOMPOSITE(c,FIRDatum,timestamp);
	H5_ADDMEMBER(c,FIRDatum,ambient);
        H5_ADDMEMBER(c,FIRDatum,sky);
      }
    };

    struct L3ScopeDatum
    {
      L3ScopeDatum(): timestamp(), scaler_l2(), scaler_qi(), scaler_hm(),
		      scaler_np(), scaler_l2_less_l3(), scaler_l3(),
		      scaler_vdacq_busy(), scaler_ten_mhz()
      { /* nothing to see here */ }

      VSTime            timestamp;
      uint32_t          scaler_l2;
      uint32_t          scaler_qi;
      uint32_t          scaler_hm;
      uint32_t          scaler_np;
      uint32_t          scaler_l2_less_l3;
      uint32_t          scaler_l3;
      uint32_t          scaler_vdacq_busy;
      uint32_t          scaler_ten_mhz;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDSIMPLECOMPOSITE(c,L3ScopeDatum,timestamp);
	H5_ADDMEMBER(c,L3ScopeDatum,scaler_l2);
	H5_ADDMEMBER(c,L3ScopeDatum,scaler_qi);
	H5_ADDMEMBER(c,L3ScopeDatum,scaler_hm);
	H5_ADDMEMBER(c,L3ScopeDatum,scaler_np);
	H5_ADDMEMBER(c,L3ScopeDatum,scaler_l2_less_l3);
	H5_ADDMEMBER(c,L3ScopeDatum,scaler_l3);
	H5_ADDMEMBER(c,L3ScopeDatum,scaler_vdacq_busy);
	H5_ADDMEMBER(c,L3ScopeDatum,scaler_ten_mhz);
      }
    };

    struct FADCChannelSettingsDatum
    {
      FADCChannelSettingsDatum(): 
	fadc_id(), fadc_channel(), lookback_time(), 
	area_width(), area_offset(), hi_lo_offset(), reread_offset(), 
	ped_offset(), is_enabled(), zero_sup_thresh(), delay()
      { /* nothing to see here */ }

      uint32_t          fadc_id;
      uint32_t          fadc_channel;
      uint32_t          lookback_time;
      uint32_t          area_width;
      uint32_t          area_offset;
      uint32_t          hi_lo_offset;
      uint32_t          reread_offset;
      uint32_t          ped_offset;
      bool              is_enabled;
      uint32_t          zero_sup_thresh;
      uint32_t          delay;
    };

    struct FADCBoardSettingsDatum
    {
      FADCBoardSettingsDatum():
	fadc_id(), mode(), event_counter(), dwords_per_chan(), is_enabled(),
	hi_lo_window_size(), reread_window_size(), ped_window_size()
      { /* nothing to see here */ }
      
      enum Mode { M_FADC, M_FULL, M_QADC, M_WORD };

      uint32_t          fadc_id;
      Mode              mode;
      uint32_t          event_counter;
      uint32_t          dwords_per_chan;
      bool              is_enabled;
      uint32_t          hi_lo_window_size;
      uint32_t          reread_window_size;
      uint32_t          ped_window_size;
    };

    struct FADCChannelRelation
    {
      FADCChannelRelation():
	telescope_id(), fadc_crate(), fadc_slot(), fadc_channel(), pixel_id(),
	chan_type(), channel_id()
      { /* nothing to see here */ }

      enum Type { T_PMT, T_L2, T_LASER, T_EMPTY };

      uint32_t          telescope_id;
      uint32_t          fadc_crate;
      uint32_t          fadc_slot;
      uint32_t          fadc_channel;
      uint32_t          pixel_id;
      Type              chan_type;
      uint32_t          channel_id;
    };

    struct FADCSlotRelation
    {
      FADCSlotRelation():
	telescope_id(), fadc_crate(), fadc_slot(), fadc_id()
      { /* nothing to see here */ }
	
      uint32_t          telescope_id;
      uint32_t          fadc_crate;
      uint32_t          fadc_slot;
      uint32_t          fadc_id;
    };

    struct TargetTableCoord
    {
      TargetTableCoord(): name(), ra_rad(), dec_rad(), epoch() 
      { /* nothing to see here */ }

      std::string       name;
      double            ra_rad;
      double            dec_rad;
      double            epoch;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDMEMBER(c,TargetTableCoord,name);
	H5_ADDMEMBER(c,TargetTableCoord,ra_rad);
	H5_ADDMEMBER(c,TargetTableCoord,dec_rad);
	H5_ADDMEMBER(c,TargetTableCoord,epoch);
      }
    };

    struct VPMCentroids
    {
      VPMCentroids()
	: timestamp(),
	  x01(), y01(), b01(), x02(), y02(), b02(), x03(), y03(), b03(),
	  x04(), y04(), b04(), x05(), y05(), b05(), x06(), y06(), b06(),
	  x07(), y07(), b07(), x08(), y08(), b08(), x09(), y09(), b09(),
	  x10(), y10(), b10(), x11(), y11(), b11(), x12(), y12(), b12(),
	  x13(), y13(), b13(), x14(), y14(), b14(), x15(), y15(), b15(),
	  x16(), y16(), b16(), x17(), y17(), b17(), x18(), y18(), b18(),
	  x19(), y19(), b19(), x20(), y20(), b20(), x21(), y21(), b21(),
	  x22(), y22(), b22(), x23(), y23(), b23(), x24(), y24(), b24(),
	  x25(), y25(), b25(), x26(), y26(), b26(), x27(), y27(), b27(),
	  x28(), y28(), b28(), x29(), y29(), b29(), x30(), y30(), b30()
      { /* nothing to see here */ }

      VSTime            timestamp;
      uint16_t          x01, y01, b01, x02, y02, b02, x03, y03, b03;
      uint16_t          x04, y04, b04, x05, y05, b05, x06, y06, b06;
      uint16_t          x07, y07, b07, x08, y08, b08, x09, y09, b09;
      uint16_t          x10, y10, b10, x11, y11, b11, x12, y12, b12;
      uint16_t          x13, y13, b13, x14, y14, b14, x15, y15, b15;
      uint16_t          x16, y16, b16, x17, y17, b17, x18, y18, b18;
      uint16_t          x19, y19, b19, x20, y20, b20, x21, y21, b21;
      uint16_t          x22, y22, b22, x23, y23, b23, x24, y24, b24;
      uint16_t          x25, y25, b25, x26, y26, b26, x27, y27, b27;
      uint16_t          x28, y28, b28, x29, y29, b29, x30, y30, b30;

#define ADDXYB(x,y,b) H5_ADDMEMBER(c,VPMCentroids,x); \
	    H5_ADDMEMBER(c,VPMCentroids,y);	      \
	    H5_ADDMEMBER(c,VPMCentroids,b)	      \

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDSIMPLECOMPOSITE(c,VPMCentroids,timestamp);
	ADDXYB(x01, y01, b01);	ADDXYB(x02, y02, b02);	ADDXYB(x03, y03, b03);
	ADDXYB(x04, y04, b04);	ADDXYB(x05, y05, b05);	ADDXYB(x06, y06, b06);
	ADDXYB(x07, y07, b07);	ADDXYB(x08, y08, b08);	ADDXYB(x09, y09, b09);
	ADDXYB(x10, y10, b10);	ADDXYB(x11, y11, b11);	ADDXYB(x12, y12, b12);
	ADDXYB(x13, y13, b13);	ADDXYB(x14, y14, b14);	ADDXYB(x15, y15, b15);
	ADDXYB(x16, y16, b16);	ADDXYB(x17, y17, b17);	ADDXYB(x18, y18, b18);
	ADDXYB(x19, y19, b19);	ADDXYB(x20, y20, b20);	ADDXYB(x21, y21, b21);
	ADDXYB(x22, y22, b22);	ADDXYB(x23, y23, b23);	ADDXYB(x24, y24, b24);
	ADDXYB(x25, y25, b25);	ADDXYB(x26, y26, b26);	ADDXYB(x27, y27, b27);
	ADDXYB(x28, y28, b28);	ADDXYB(x29, y29, b29);	ADDXYB(x30, y30, b30);
      }

#undef ADDXYB
    };

    struct VPMLEDs
    {
      VPMLEDs()
	: timestamp(),
	  n1(), x1(), y1(), 
	  n2(), x2(), y2(), 
	  n3(), x3(), y3(), 
	  n4(), x4(), y4()
      { /* nothing to see here */ }

      VSTime            timestamp;
      uint16_t          n1, x1, y1;
      uint16_t          n2, x2, y2;
      uint16_t          n3, x3, y3;
      uint16_t          n4, x4, y4;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDSIMPLECOMPOSITE(c,VPMLEDs,timestamp);
	H5_ADDMEMBER(c,VPMLEDs,n1);
	H5_ADDMEMBER(c,VPMLEDs,x1);
	H5_ADDMEMBER(c,VPMLEDs,y1);
	H5_ADDMEMBER(c,VPMLEDs,n2);
	H5_ADDMEMBER(c,VPMLEDs,x2);
	H5_ADDMEMBER(c,VPMLEDs,y2);
	H5_ADDMEMBER(c,VPMLEDs,n3);
	H5_ADDMEMBER(c,VPMLEDs,x3);
	H5_ADDMEMBER(c,VPMLEDs,y3);
	H5_ADDMEMBER(c,VPMLEDs,n4);
	H5_ADDMEMBER(c,VPMLEDs,x4);
	H5_ADDMEMBER(c,VPMLEDs,y4);
      }
    };

  }
}

#endif // defined VSCENTRALIZEDDBACCESSDATA_HPP
