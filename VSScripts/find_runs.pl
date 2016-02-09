#!/usr/bin/perl -w
#
# find_runs.pl -- Finds runs filtered by user inputs.
#


use strict;
use Getopt::Long;
use Pod::Usage;

# -----------------------------------------------------------------------------
# Global Parameters
# -----------------------------------------------------------------------------
my $diagnostics_file = "/veritas/chila_diagnostics/diagnostics.dat";
my $reduced_file = "/veritas/chila_reduced/";


# -----------------------------------------------------------------------------
# Parse the Options
# -----------------------------------------------------------------------------
my $opt_help = 0;
my $opt_source = "";
my $opt_run = "";
my $opt_el_lo = 0;
my $opt_el_hi = 90;
my $opt_date_lo = 00000000;
my $opt_date_hi = 99999999;
my $opt_stage2_mode = 0;
my $opt_tel_min = 3;
my $opt_tel_mask = 0;
my $opt_max_moon_el = 90;
my $opt_elapsed_time = 0;
my $opt_chisq = 100.0;
my $opt_max_fir1_rms = 0.0;

GetOptions( 'help!'            => \$opt_help,
	    's2_mode!'         => \$opt_stage2_mode,
	    'date_lo=s'        => \$opt_date_lo,  # =s forces string argument
	    'date_hi=s'        => \$opt_date_hi,
            'source=s'         => \$opt_source,
            'el_lo=i'          => \$opt_el_lo,    # =i forces integer argument
            'el_hi=i'          => \$opt_el_hi,
	    'tel_min=i'        => \$opt_tel_min,
	    'tel_mask!'        => \$opt_tel_mask,
	    'fir1_rms=f'       => \$opt_max_fir1_rms,
	    'max_moon_el=i'    => \$opt_max_moon_el,
	    'elapsed_time=f'   => \$opt_elapsed_time,
	    'chisq=f'          => \$opt_chisq     # =f forces float argument
            ) or pod2usage(2);

pod2usage(1) if $opt_help;   # prints help file

open(DIAGNOSTICS,$diagnostics_file);

my $sum_time = 0;

while(<DIAGNOSTICS>) 
{
    chomp $_;
    
    my %diagnostics = parse_diagnostics($_);
    
    next if $opt_source ne "" && $diagnostics{'source'} !~ /$opt_source/;

    if($opt_max_fir1_rms > 0.0 &&
       $diagnostics{'fir1_rms'} > $opt_max_fir1_rms) { next; }

    if($diagnostics{'mean_el'} < $opt_el_lo) { next; }
    elsif($diagnostics{'mean_el'} > $opt_el_hi) { next; }
    
    if($diagnostics{'date'} < $opt_date_lo) { next; }
    elsif($diagnostics{'date'} > $opt_date_hi) { next; }

    #---------------------------------------------------------------
    # NOTE: Use either $opt_tel_min or $opt_tel_mask but NOT BOTH
    #---------------------------------------------------------------
    if($diagnostics{'ntel'} < $opt_tel_min) { next; }

    # Telescope cut; if 3 telescopes and T1 or T4 is out, then use. Otherwise,
    #   skip.
    if ($opt_tel_mask) {
	if ($diagnostics{'ntel'} < 3) { next; }
	if ($diagnostics{'ntel'} == 3) {
	    
	    if ($diagnostics{'config_mask'} == 1101  || 
		$diagnostics{'config_mask'} == 1011 )  
	    { next; }
	
	}
    }


    # Moon elevation cut:
    if($diagnostics{'moon_el'} > $opt_max_moon_el) { next; }

    if($diagnostics{'gps_elapsed_time'} < $opt_elapsed_time ) { next; }

    if($diagnostics{'chisq'} > $opt_chisq ) { next; }

    if ($opt_stage2_mode) { stage2_mode(\%diagnostics); }
    else { vbf_mode(\%diagnostics); }
    
    $sum_time = $sum_time + $diagnostics{'gps_elapsed_time'};
}

print "TOTAL TIME: $sum_time\n";

sub vbf_mode
{
    my %diagnostics=%{$_[0]};
    
    print 
	"$diagnostics{'date'} " .
	"$diagnostics{'run'} " .
	"$diagnostics{'laser_run'} " .
	"$diagnostics{'source'} " .
	"$diagnostics{'time'} " .
	"$diagnostics{'mode'} " .	 
	"$diagnostics{'config_mask'} " .	 
	"$diagnostics{'mean_el'} " .	 
	"$diagnostics{'gps_elapsed_time'} " .
	sprintf("%6.1f",$diagnostics{'l3_rate'}) . 
	sprintf("%6.2f",$diagnostics{'chisq'}) . 
	sprintf("%6.1f",$diagnostics{'moon_el'}) . 
	sprintf("%6.2f",$diagnostics{'fir0_rms'}) . 
	sprintf("%6.2f",$diagnostics{'fir1_rms'}) . 
	sprintf("%6.2f",$diagnostics{'fir3_rms'}) . 
	sprintf("%8s",$diagnostics{'t1_laser_run'}) . 
	sprintf("%8s",$diagnostics{'t2_laser_run'}) . 
	sprintf("%8s",$diagnostics{'t3_laser_run'}) . 
	sprintf("%8s",$diagnostics{'t4_laser_run'}) . 
	"\n";
#	"$diagnostics{'moon_el'} \n";
}


sub stage2_mode
{

    my %diagnostics=%{$_[0]};
    
    my $filename = $reduced_file."d"."$diagnostics{'date'}"."/".
	"x"."$diagnostics{'run'}"."_s2.h5";
    print "$filename\n";
}


sub parse_diagnostics
{
    chomp $_[0];
    my @line = split ' ', $_[0];
    my %diagnostics = ( 'run' => 0 ,
			'laser_run' => 0 ,
			'mean_el' => 0,
			'date' => 0,
			'time' => 0,
			'source' => 'none',
			'mode' => 'none',
			'gps_elapsed_time' => 0,
			'ntel' => 0);

    
  
    if(scalar(@line) < 41)
    {
	return %diagnostics;
    }	

    my $version = $line[1];

    $diagnostics{'run'} = $line[0];
    $diagnostics{'laser_run'} = sprintf("%7s",$line[42]);
    $diagnostics{'date'} = $line[2];
    $diagnostics{'mean_el'} = $line[7];
    $diagnostics{'mean_az'} = $line[8];
    $diagnostics{'time'} = substr($line[3],0,8);
    $diagnostics{'source'} = $line[5];
    $diagnostics{'mode'} = sprintf("%14s",$line[6]);
    $diagnostics{'gps_elapsed_time'} = sprintf("%5.2f",$line[40]);
    $diagnostics{'l3_rate'} = sprintf("%5.2f",$line[12]);
    $diagnostics{'chisq'} = sprintf("%5.2f",$line[14]);
    
    $diagnostics{'moon_el'} = $line[43];

    $diagnostics{'t1_current'} = $line[50];
    $diagnostics{'t2_current'} = $line[64];
    $diagnostics{'t3_current'} = $line[78];
    $diagnostics{'t4_current'} = $line[92];

    if($version == 2)
    {
	$diagnostics{'fir0_rms'} = $line[106];
	$diagnostics{'fir1_rms'} = $line[111];
	$diagnostics{'fir3_rms'} = $line[121];

	$diagnostics{'t1_laser_run'} = $line[107];
	$diagnostics{'t2_laser_run'} = $line[112];
	$diagnostics{'t3_laser_run'} = $line[117];
	$diagnostics{'t4_laser_run'} = $line[122];
    }
    else
    {
	$diagnostics{'fir0_rms'} = 0;
	$diagnostics{'fir1_rms'} = 0;
	$diagnostics{'fir3_rms'} = 0;

	$diagnostics{'t1_laser_run'} = $diagnostics{'laser_run'};
	$diagnostics{'t2_laser_run'} = $diagnostics{'laser_run'};
	$diagnostics{'t3_laser_run'} = $diagnostics{'laser_run'};
	$diagnostics{'t4_laser_run'} = $diagnostics{'laser_run'};
    }

    my @cfd_threshold = ($line[17],$line[23],$line[29],$line[35]);

    $diagnostics{'t1_cfd_threshold'} = $line[17];
    $diagnostics{'t2_cfd_threshold'} = $line[23];
    $diagnostics{'t3_cfd_threshold'} = $line[29];
    $diagnostics{'t4_cfd_threshold'} = $line[35];

    if($diagnostics{'date'} =~ /(....)-(..)-(..)/)
    {
	$diagnostics{'date'} = "${1}${2}${3}";
    }

    my $config_mask;

    my $itel = 0;
    for(my $i = 0; $i < 4; $i++)
    {
	if($cfd_threshold[$i] > 0)
	{
	    $config_mask = $config_mask . "1";
	    $itel++;
	}
	else
	{
	    $config_mask = $config_mask . "0";
	}

    }

    $diagnostics{'ntel'} = $itel;
    
    $diagnostics{'config_mask'} = $config_mask;

    return %diagnostics;
}


=head1 NAME

find_runs.pl - Script for obtaining runs for either stage 1/2 analysis, or for obtaining the paths to processed stage2 files, but chila next day analysis.

=head1 SYNOPSIS

B<find_runs.pl> [options]

Options:

=over 8

=item B<--help>

Print this help message.

=item B<--s2_mode>

Prints out paths to all stage 2 files instead of data relevant to undergoing stage 1/2 analysis on vbf files.

=item B<--date_lo=yyyymmdd>

Lists runs if they occured after this date.

=item B<--date_hi=yyyymmdd>

Lists runs if they occured before this date.

=item B<--source=s>

Lists only runs for this source.

=item B<--el_lo=i>

Lists only runs with elevation above this threshold.

=item B<--el_hi=i>

Lists only runs with elevation below this threshold.

=item B<--tel_min=i>

Min number of telescopes in the data run to list

=item B<--tel_mask>

Masks 3 telescope runs unless it was T4 or T1 down

=item B<--max_moon_el=i>

Maximum moon elevation allowed

=item B<--elapsed_time=i>

Minimum elapsed time for run to be active for.

=item B<--chisq=i>

Maximum allowed chisq value of L3 rate for constant rate model.

=back




