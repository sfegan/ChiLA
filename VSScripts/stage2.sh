code_dir=/home/sfegan/vbf
stage2=$code_dir/stage2
s1opt="-VSDBHost=vela.physics.umass.edu -VSDBUser=readonly"
s2opt="$s1opt -s2_method=1 -s2_cleaning=regional,4,3,9 -s2_nthreads=3 -s2_nscope_cut=1 -s2_wobble_off_region_radius=0.18 -s2_muon_nimage_cut=50 -s2_muon_radius_min_cut=0.7 -s2_muon_radius_max_cut=1.4 -s2_muon_width_max_cut=0.2 -s2_muon_ring_edge_dist_max_cut=1.8 -s2_muon_centroid_radius_ratio_max_cut=0.75 -s2_ped_suppress_mode=interval -s2_ped_suppress_fraction=0.85 -s2_nimage_cut=4 -s2_nscope_cut=2"

if ( test "$2" == "" )
then
  lasopt="-s2_no_laser"
else
  lasopt="-laser=$2"
fi

$stage2 $s2opt $1 $lasopt

