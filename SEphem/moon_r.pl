#!/usr/bin/perl -w

#
# moon_r.pl
#
# Convert moon distance corefficients from moon1.dat to C code
#
# Original Author: Stephen Fegan
# $Author: sfegan $
# $Date: 2007/11/04 20:33:05 $
# $Revision: 1.1 $
# $Tag$
#

use strict;
use FileHandle;

my $moon = new FileHandle("<moon1.dat");
my $line = $moon->getline();
while ( defined $line )
  {
    chomp $line;
    $line =~ s/#.*//;
    $line =~ s/^\s*$//;
    if ( !$line )
      {
	$line = $moon->getline();
	next;
      }
    my ($D,$M,$MP,$F,$Sl,$Sr) = split /\t/,$line;
    $Sl =~ s/\s*//g;
    $Sr =~ s/\s*//g;

    if ( $Sr == 0 )
      {
	$line = $moon->getline();
	next;
      }

    print "  Sr ";
    if($Sr > 0) { print "+= "; }
    else { print "-= "; }

    if($M == 0) { }
    elsif(($M == 1)||($M == -1)) { print "E*"; }
    elsif(($M == 2)||($M == -2)) { print "E2*"; }
    elsif(($M == 3)||($M == -3)) { print "E3*"; }
    elsif(($M == 4)||($M == -4)) { print "E4*"; }
    else { die; }

    print abs($Sr),"*cos(";

    my $first = "";

    if($D>0) { print "D",abs($D); $first = "+"; }
    elsif($D<0) { print "-D",abs($D); $first = "+"; }
    if($M>0) { print $first,"M",abs($M); $first = "+"; }
    elsif($M<0) { print "-M",abs($M); $first = "+"; }
    if($MP>0) { print $first,"MP",abs($MP); $first = "+"; }
    elsif($MP<0) { print "-MP",abs($MP); $first = "+"; }
    if($F>0) { print $first,"F",abs($F); $first = "+"; }
    elsif($F<0) { print "-F",abs($F); $first = "+"; }

    print ");";

    print "\n";
    $line = $moon->getline();
  }
