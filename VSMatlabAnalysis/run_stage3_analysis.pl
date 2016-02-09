#! /usr/bin/perl -w
#
#    This script runs the stage3_analysis.m file, without display in the
#    MATLAB environment, and creates an output .ps file.
#
#    NOTE: This will not display the correct .ps file, unless the created file
#          is the only file named '_analysis.ps' in the current directory 
#          (which it should be).
#
#    Timothy C. Arlen
#    UCLA
#    arlen@astro.ucla.edu
#

use strict;
use Getopt::Long;
use Pod::Usage;

# -----------------------------------------------------------------------------
# Parse the Options
# -----------------------------------------------------------------------------
my $opt_help = 0;
my $opt_s3file = "";
my $opt_run_convert = 0;

GetOptions( 'help!'             => \$opt_help,
	    's3file=s'          => \$opt_s3file,
	    'run_convert!'      => \$opt_run_convert
            ) or pod2usage(2);

pod2usage(1) if $opt_help;   # prints help file

my $matlab_cmd = "\"stage3_analysis(\'$opt_s3file\'); exit;\"";

print `DISPLAY="" /veritas/matlab/bin/matlab -logfile stage3_analysis.log  -r $matlab_cmd`;

my $psfile = "";
$psfile = `ls *_analysis.ps`;
chomp($psfile);
print "File created: $psfile\n";

my $pngfile = $psfile;
$pngfile =~ s/.ps/.png/;

#print "Png File: $pngfile\n";

if ($opt_run_convert) {   
    print `convert -density 300 -quality 90 -resize 1200 $psfile $pngfile`;
}

#///////////////////////////////
# Description of all options:   
#///////////////////////////////
=head1 NAME

run_stage3_analysis.pl - Scripts for running stage3 Matlab Analysis script: stage3_analysis.m.

=head1 SYNOPSIS

B<run_stage3_analysis.pl> [options]

Options:

=over 8

=item B<--help>

Print this help message.

=item B<--s3file>

stage3 file which is processed by stage3_analysis.m

=item B<--run_convert>

Set to true for running the convert program, to convert output .ps file to .png file (more useful most of the time).

=back
