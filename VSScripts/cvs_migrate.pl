#!/usr/bin/perl
#
# Program:     cvs_migrate.pl
# Author:      Matthew Wood <mdwood@astro.ucla.edu>
# Date:        07/24/08
#
# Description: Script for migrating a cvs module to a new repository.
#
# $Id: cvs_migrate.pl,v 1.2 2009/02/06 22:46:07 matthew Exp $
#

use warnings;
use strict;

my $hostname = "romulus.ucsc.edu";

my(@dirs,@cvs_dirs);

push(@dirs,@ARGV);

while(@dirs)
{
    my $dir = pop(@dirs);
    next if ! -d $dir;
    $dir =~ s/^(.*)([^\/])$/$1$2\//;

    if($dir =~ /CVS\/$/)
    {
	push @cvs_dirs, $dir;
    }
    else
    {
	push @dirs, <$dir*>;
    }   
}

foreach my $dir (@cvs_dirs)
{
    open(ROOT,"$dir/Root") or die $!;
    my $root = <ROOT>;
    close(ROOT);

    $root =~ s/^(.+):(.+)\@(.+):(.+)$/$1:$2\@$hostname:$4/;

    open(ROOT,">$dir/Root") or die $!;
    print ROOT "$root";
    close(ROOT);
}
