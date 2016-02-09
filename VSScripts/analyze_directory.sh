#!/bin/bash

ulimit -c unlimited

code_dir=/home/mdwood/raid/diagnostics/bin
toffset=$code_dir/laser
stage2=$code_dir/stage2
#dbopt="-VSDBHost=vraid.neutrino.hosted.ats.ucla.edu -VSDBUser=readonly"
dbopt="-VSDBHost=romulus.ucsc.edu -VSDBUser=readonly"
#dbopt="-VSDBHost=vela.physics.umass.edu -VSDBUser=readonly"
s1opt="$dbopt -s1_no_nsb_suppress=true"
toopt="$dbopt -no_nsb_suppress=true -no_pedestals_in_core -ped_suppress_lo=0.333 -ped_suppress_hi=3.0"
s2opt="$s1opt -s2_nthreads=3 -s2_permissive_laser -s1_no_pedestals_in_core -s2_method=1 -s2_weighting=size_ellipticity -s2_nimage_cut=2 -s2_nscope_cut=2 -s2_qc=3/250/1.5/2 -s2_cleaning=regional,4.5,3,10 -s2_window_width=5 -s2_wobble_off_region_radius=0.18 -s2_muon_nimage_cut=50 -s2_muon_radius_min_cut=0.7 -s2_muon_radius_max_cut=1.4 -s2_muon_width_max_cut=0.2 -s2_muon_ring_edge_dist_max_cut=1.8 -s2_muon_centroid_radius_ratio_max_cut=0.75 -s2_ped_suppress_mode=interval -s2_ped_suppress_fraction=0.85"
src_dir=$1
tar_dir=$2
laser=$3

if ( test "$4" != "" )
then
  overwrite=$4
else
  overwrite=1
fi

###############################################################################
# CREATE THE OUTPUT DIRECTORY
###############################################################################

if ( ! test -d $tar_dir )
then
  mkdir $tar_dir
fi

###############################################################################
# GET THE LASER OPTION / CALCULATE THE LASER GAINS
###############################################################################

if ( test "$laser" != "" ) && ( ! test -f $laser )
then
  laser=$src_dir/$laser.cvbf
fi

if ( test "$laser" == "" )
then
  lasopt="-s2_no_laser"
else
  if ( file $laser | grep "Hierarchical Data Format (version 5)" )
  then
    # assume laser has already been calculated
    lasopt="-laser=$laser"
  elif ( file $laser | grep "ASCII text" )
  then
    # assume laser has already been calculated
    lasopt="-laser=$laser"
  else
    # assume we have been given a VBF file so run toffset
    runno=`basename $laser .cvbf`
    runno=`basename $runno .vbf`
    $toffset $toopt $laser 2>&1 | tee $tar_dir/x${runno}_laser.log
    lasopt="-laser=x${runno}_s1.h5"
  fi
fi

###############################################################################
# OPERATE ON ALL FILES
###############################################################################

for vbf in $src_dir/*vbf
do
  runno=`basename $vbf .cvbf`
  runno=`basename $runno .vbf`
  s1file=$tar_dir/x${runno}_s1.h5
  s2file=$tar_dir/x${runno}_s2.h5;
  if( test "$overwrite" != 0 -o ! -f ${s2file} )
  then
    $stage2 $s2opt $vbf $lasopt -stage1=${s1file} -o ${s2file} 2>&1 | tee $tar_dir/x${runno}.log
  else
    echo Skipped ${runno} since ${s2file} exists
  fi
done
