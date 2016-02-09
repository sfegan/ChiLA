#!/bin/awk -f
#
# Simple script to extract live and elapsed times (L3 and GPS) from the
# diagnostics.dat file (which must be in the current directory). Input
# should be list of run numbers, one per line.
#
# $Id: gettime.awk,v 3.1 2007/06/18 18:17:31 sfegan Exp $
#
BEGIN{
  while((getline < "diagnostics.dat")>0)
    {
      mjd[$1]     = $5;
      l3elap[$1]  = $10;
      l3live[$1]  = $11;
      gpselap[$1] = $41; 
      gpslive[$1] = $42;
    };
};
NF>1{
  print $1,mjd[$1],l3elap[$1],l3live[$1],gpselap[$1],gpslive[$1];
  all_l3elap  = all_l3elap  + l3elap[$1];
  all_l3live  = all_l3live  + l3live[$1];
  all_gpselap = all_gpselap + gpselap[$1];
  all_gpslive = all_gpslive + gpslive[$1];
}
END{
  print "total",all_l3elap,all_l3live,all_gpselap,all_gpslive;
}
