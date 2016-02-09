#!/bin/bash
code_dir=/home/sfegan/vbf
laser=$code_dir/laser
s1opt="-VSDBHost=vela.physics.umass.edu -VSDBUser=readonly"
laseropt="$s1opt -ped_suppress_lo=0.333 -ped_suppress_hi=3.0"
$laser $laseropt $1

