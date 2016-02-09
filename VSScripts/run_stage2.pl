#!/usr/bin/perl -w
#
# Program:     run_stage2.pl
# Author:      Matthew Wood <mdwood@astro.ucla.edu>
# Date:        10/25/07
#
# Description: Script for running stage1/stage2 chila analysis.
#
# $Id: run_stage2.pl,v 1.21 2010/06/18 21:27:55 matthew Exp $
#
#$ -e $HOME/logs/stage2_logs
#$ -o $HOME/logs/stage2_logs

use strict;
use Getopt::Long;
use Pod::Usage;
use Fcntl qw(:flock);

$ENV{'LD_LIBRARY_PATH'}=
    "$ENV{'LD_LIBRARY_PATH'}:$ENV{'HOME'}/local/lib:$ENV{'HOME'}/local/lib/mysql";

# -----------------------------------------------------------------------------
# Default/Global Parameters
# -----------------------------------------------------------------------------
my $max_telescopes = 16;

my $chila_dir;
my $s2_bin;
my $laser_bin;

if(exists($ENV{'CHILADIR'}) && (-d $ENV{'CHILADIR'}))
{
    $chila_dir = $ENV{'CHILADIR'};
}
else
{
    $chila_dir = "$ENV{'HOME'}/ChiLA";
}

if(!(-e "$chila_dir/VSDataReduction/stage2") || 
   !(-e "$chila_dir/VSDataReduction/laser"))
{
    print "Could not locate stage2/laser executables in: $chila_dir\n";
    exit(1);
}
else
{
    $s2_bin = "$chila_dir/VSDataReduction/stage2";
    $laser_bin = "$chila_dir/VSDataReduction/laser";
}

# -----------------------------------------------------------------------------
# Default stage1/stage2 options
# -----------------------------------------------------------------------------
my %default_cfg_stg2;

$default_cfg_stg2{'s2_method'}="1";
$default_cfg_stg2{'s2_qc'}="3/250/1.5/2";
$default_cfg_stg2{'s2_cleaning'}="regional,4.5,3,10";
$default_cfg_stg2{'s2_weighting'}="size_ellipticity";
$default_cfg_stg2{'s2_window_width'}="5";
$default_cfg_stg2{'s2_nthreads'}="4";
$default_cfg_stg2{'s2_nscope_cut'}="2";
$default_cfg_stg2{'s2_nimage_cut'}="3";
$default_cfg_stg2{'s2_no_muon_analysis'}="true";
$default_cfg_stg2{'s2_ped_suppress_mode'}="interval";
$default_cfg_stg2{'s2_ped_suppress_fraction'}="0.85";
$default_cfg_stg2{'s2_msc_weight_power'}="1";
$default_cfg_stg2{'s2_sc_parameter_lookup'} = "";
# $default_cfg_stg2{'s2_tracking_recorrections_date'}="0123,2007-11-01";
$default_cfg_stg2{'s2_tracking_recorrections_date'}="";
$default_cfg_stg2{'s2_energy_weight_power'}="1";
$default_cfg_stg2{'s2_energy_lookup'} = "";
$default_cfg_stg2{'s2_demand_tracking_target'} = "";
$default_cfg_stg2{'s2_cam_rotation'} = "";
$default_cfg_stg2{'s2_npackets'}="0";
$default_cfg_stg2{'s2_no_set_scope_pos_from_sims'}="false";
$default_cfg_stg2{'s2_permissive_laser'}="false";
$default_cfg_stg2{'s2_pad_zero_suppressed_chan'}="false";
$default_cfg_stg2{'s2_lo_gain_zero_sample'}="3,3,3,3";
$default_cfg_stg2{'s1_no_nsb_suppress'}="true";

# -----------------------------------------------------------------------------
# Default DB options
# -----------------------------------------------------------------------------
my %default_cfg_db;

$default_cfg_db{'VSDBHost'}="romulus.ucsc.edu";
$default_cfg_db{'VSDBUser'}="readonly";
$default_cfg_db{'VSDBRdbms'}="MYSQL3X";
$default_cfg_db{'no_db'}="false";

# -----------------------------------------------------------------------------
# Default Laser options
# -----------------------------------------------------------------------------
my %default_cfg_lsr;

$default_cfg_lsr{'ped_suppress_lo'}="0.333";
$default_cfg_lsr{'ped_suppress_hi'}="3.0";

# -----------------------------------------------------------------------------
# Parse the Options
# -----------------------------------------------------------------------------
my %opt;
$opt{'help'} = 0;
$opt{'man'} = 0;
$opt{'sim'} = 0;
$opt{'laser'} = "";
$opt{'mask'} = "";
$opt{'scopes'} = "";
$opt{'multiplicity'} = 2;
$opt{'output_file'} = "";
$opt{'overwrite'} = 0;
$opt{'verbose'} = 0;
$opt{'s2_dir'} = "";

my %opt_hash = 
(
    'h|help!'        => \$opt{'help'},
    'man!'           => \$opt{'man'}, 
    'sim!'           => \$opt{'sim'},
    'v|verbose!'     => \$opt{'verbose'},
    'cfg=s'          => \$opt{'cfg'},
    'laser=s'        => \$opt{'laser'},
    'mask=s'         => \$opt{'mask'},
    'scopes=s'       => \$opt{'scopes'},
    'o=s'            => \$opt{'output_file'},
    'overwrite!'     => \$opt{'overwrite'},
    'multiplicity=s' => \$opt{'multiplicity'},
    's2_dir=s'       => \$opt{'s2_dir'}
);

my %cfg_lsr, my %cfg_stg2, my %cfg_stg3, my %cfg_db;

foreach my $setting ( keys(%default_cfg_lsr) )
{
    $opt_hash{"${setting}=s"} = \$cfg_lsr{$setting};
}

foreach my $setting ( keys(%default_cfg_stg2) )
{
    $opt_hash{"${setting}=s"} = \$cfg_stg2{$setting};
}

foreach my $setting ( keys(%default_cfg_db) )
{
    $opt_hash{"${setting}=s"} = \$cfg_db{$setting};
}

GetOptions( %opt_hash ) or pod2usage(2);

# -----------------------------------------------------------------------------
# Print help/man if asked or if no file
# -----------------------------------------------------------------------------
pod2usage(1) if $opt{'help'};
pod2usage(-verbose => 2) if $opt{'man'};
pod2usage("$0: No list or vbf file given.") if ((@ARGV == 0) && (-t STDIN));

my @arglist = @ARGV;

# -----------------------------------------------------------------------------
# Load stage1/stage2/stage3/laser options from cfg file
# -----------------------------------------------------------------------------
if($opt{'cfg'})
{
    open(CFG,$opt{'cfg'}) or die "Can't open $opt{'cfg'}: $!";
    while(<CFG>)
    {
	next if /^\#/ || /^$/;
	my ($opt,$value) = split ' ',$_;
	next if !$opt || !$value;

	if(($opt =~ /VSDB/ || $opt eq "no_db") && !defined($cfg_db{$opt}))
	{
	    $cfg_db{$opt} = $value;
	}
	elsif(($opt =~ /^s2_/ || $opt =~ /^s1_/) && !defined($cfg_stg2{$opt}))
	{
	    $cfg_stg2{$opt} = $value;	
	}
	elsif(!defined($cfg_lsr{$opt}))
	{
	    $cfg_lsr{$opt} = $value;
	}
    }
}

# -----------------------------------------------------------------------------
# Load default options
# -----------------------------------------------------------------------------
load_default_cfg(\%cfg_lsr,\%default_cfg_lsr);
load_default_cfg(\%cfg_stg2,\%default_cfg_stg2);
load_default_cfg(\%cfg_db,\%default_cfg_db);

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------
$opt{'s2_dir'} =~ s/^(.*)([^\/])$/$1$2\//;

if(!(-d $opt{'s2_dir'}) && $opt{'s2_dir'})
{
    mkdir $opt{'s2_dir'} or die $!;
}

my %vbf_file_hash = ();

# -----------------------------------------------------------------------------
# Determine DB settings
# -----------------------------------------------------------------------------
if($opt{'sim'}) 
{
    $cfg_db{'no_db'}="true";
    delete($cfg_db{'VSDBHost'});
    delete($cfg_db{'VSDBUser'});
}
 
# =============================================================================
# Main Analysis Loop
# =============================================================================
if(exists $ENV{'TMPDIR'})
{
    chdir($ENV{'TMPDIR'}) or die "Can't chdir $ENV{'TMPDIR'}: $!";
}

my @run_list = ();

foreach my $file (@arglist)
{
    parse_list($file,\@run_list);
}

# Run stage2 analysis ---------------------------------------------------------
foreach my $run (@run_list) 
{
    # -------------------------------------------------------------------------
    # Check for the existence of the laser file.  
    # -------------------------------------------------------------------------
    if(exists($run->{'laser'}))
    {
	if($run->{'laser'} =~ /([0-9]+)\.(vbf|cvbf)$/) 
	{       
	    $run->{'laser_s1'} =  "x$1_s1.h5";
	} 
	elsif($run->{'laser'} =~ /([^\/]+)\.(vbf|cvbf)$/) 
	{       
	    $run->{'laser_s1'} = "$1_s1.h5";
	} 
	else
	{
	    print "ERROR\n";
	    exit(1);
	}
    }

    my $has_laser = 1;
    $has_laser = run_laser($run) if exists($run->{'laser'});
    run_stage2($run) if $has_laser;
}

sub parse_list
{
    my $file = $_[0];

    if($file =~ /\.(cvbf|vbf)$/)
    {
	# Data Run ------------------------------------------------------------
	my $data = find_vbf_file($file,"");
	
	# Laser Run -----------------------------------------------------------
	if($opt{'laser'})
	{
	    my $laser = find_vbf_file($opt{'laser'},"");
	    push(@{$_[1]}, { 'data' => $data, 'laser' => $laser } );
	}
	else
	{
	    push(@{$_[1]}, { 'data' => $data } );
	}	
    }
    else
    {
	print "Opening list file: $file\n";
	
	open(LIST,"$file") or die print "Couldn't open file: $file\n";
	while(<LIST>) 
	{
	    next if /^\#/ || /^$/;

	    chomp;
	    my @line = split ' ',$_;
	    my $ncol = scalar(@line);

	    my ($data, $laser);

	    if($ncol == 1)
	    {
		$data = find_vbf_file($line[0],"");
	    }
	    elsif($ncol == 2)
	    {
		$data = find_vbf_file($line[1],$line[0]);
	    }
	    elsif($ncol >= 3)
	    {
		$data = find_vbf_file($line[1],$line[0]);
		$laser = find_vbf_file($line[2],$line[0]);
	    }
	    else
	    {
		print "ERROR\n";
		exit(1);
	    }

	    if($opt{'laser'}) { $laser = $opt{'laser'}; }

	    if($laser)
	    {
		push(@{$_[1]}, { 'data' => $data, 'laser' => $laser } );
	    }
	    else
	    {
		push(@{$_[1]}, { 'data' => $data } );
	    }
	}
    }
}

# =============================================================================
# Subroutine for stage1/stage2
# =============================================================================
sub run_stage2 
{
    my %cfg_run = %{$_[0]};
    my %cfg = %cfg_stg2;
    my $data_file = $cfg_run{'data'};
    my ($log_file,$log_settings);

    # Output File -------------------------------------------------------------
    my $stg2_output_file;

    if($cfg_run{'data'} =~ /([0-9]{3,})\.(vbf|cvbf)$/) 
    {
	$cfg{'stage1'} = "x$1_s1.h5";
	$cfg{'o'}      = "x$1_s2.h5";
	$stg2_output_file   = "$opt{'s2_dir'}$1_s2.h5";
	$log_file           = "$opt{'s2_dir'}x$1_s2.log";
    } 
    elsif($cfg_run{'data'} =~ /([^\/]+)\.(vbf|cvbf|vbf.gz|cvbf.gz)$/) 
    {
	$cfg{'stage1'} = "$1_s1.h5";
	$cfg{'o'}      = "$1_s2.h5";
	$stg2_output_file   = "$opt{'s2_dir'}$1_s2.h5";
	$log_file           = "$opt{'s2_dir'}$1_s2.log";
    } 
    else
    {
	print "ERROR parsing filename\n";
	exit(1);
    }
 
    if(-e $stg2_output_file && !$opt{'overwrite'})
    {
	print "stage2 output file already exists: $stg2_output_file\n";
	return;
    }
    elsif($opt{'overwrite'})
    {
	`rm $cfg{'stage1'}`;
	`rm $cfg{'o'}`;
	`rm $log_file`;
    }

    open(LOG,">>$log_file");
    return if !flock(LOG, LOCK_EX | LOCK_NB);    

    my $tmp_dir=`mktemp -d run_stage2.XXXXXXXXXX`;
    chomp $tmp_dir;

    $cfg{'o'} = "$tmp_dir/$cfg{'o'}";
    $cfg{'stage1'} = "$tmp_dir/$cfg{'stage1'}";

    my $tmp_file;

    if($data_file =~ /([^\/]+)\.(vbf\.gz|cvbf\.gz)$/) 
    {
	$tmp_file = `mktemp $tmp_dir/tmp.XXXXXXXXXX`;
	chomp $tmp_file;
	
	print "gunzip -c $data_file > $tmp_file\n";
	`gunzip -c $data_file > $tmp_file`;
	$data_file = "$tmp_file";
    }

    # -------------------------------------------------------------------------
    # Derive stage1/stage2 settings
    # -------------------------------------------------------------------------

    # Laser Options -----------------------------------------------------------
    if(exists($cfg_run{'laser_s1'}))
    {
	$cfg{'s2_no_laser'}="false";
	$cfg{'laser'}=$cfg_run{'laser_s1'};
    }
    else
    {
	$cfg{'s2_no_laser'}="true";
    }

    # Simulation Options ------------------------------------------------------
    if($opt{'sim'})
    {
	$cfg{'s2_pointing'}="l3direct";
	$cfg{'s2_lo_gain_zero_sample'}="0,0,0,0";
	$cfg{'s1_no_pedestals_in_core'}="true";	
	$cfg{'s2_pad_zero_suppressed_chan'}="true";
    }

    if($opt{'scopes'})
    {
	my $scope_suppress = join(',',get_suppressed_scopes($opt{'scopes'}));

	if($scope_suppress ne "")
	{
	    $cfg{'s2_scope_suppress'}=$scope_suppress;
	}

	$cfg{'s2_software_trigger_masks'}=
	    join(',',get_allowed_trigger_masks($opt{'scopes'},
					       $opt{'multiplicity'}));
    }

    if($opt{'verbose'})
    {
	$log_settings = "2>&1 | tee $log_file";
    }
    else
    {
	$log_settings = "&> $log_file";
    }

    # Merge DB options --------------------------------------------------------
    @cfg{ keys %cfg_db } = values %cfg_db;

    my $cfg = "";
    while(my ($opt,$val) = each(%cfg)) 
    {
	next if $val eq "";
	$cfg = $cfg . "-$opt=$val "; 
    }
    
    # -------------------------------------------------------------------------
    # Run Stage1/Stage2
    # -------------------------------------------------------------------------
    print "----------------------------------------------------------------\n";
    print "Analyzing Data Run\n";
    print "Data Path...........: $cfg_run{'data'}\n";    
    if(exists($cfg_run{'laser'})) 
    { print "Laser Path..........: $cfg_run{'laser'}\n"; }
    print "Log File............: $log_file\n";
    print "Stage1 File.........: $cfg{'stage1'}\n";
    print "Stage2 File.........: $cfg{'o'}\n";
    print "Command: $s2_bin $cfg $data_file\n";
    system("$s2_bin $cfg $data_file $log_settings");

    `mv $cfg{'stage1'} $opt{'s2_dir'}`;
    `mv $cfg{'o'} $opt{'s2_dir'}`;

    `rm -rf $tmp_dir`;

    close(LOG);
}

# =============================================================================
# Subroutine for laser
# =============================================================================
sub run_laser 
{
    my %cfg_run = %{$_[0]};

    my $vbf_file = $cfg_run{'laser'};    
    my $output_file = $cfg_run{'laser_s1'};
    my $log_file = "laser.log";
    my $log_settings;

    my %cfg = %cfg_lsr;

    if(!(-e $vbf_file))
    {
	print "ERROR Laser file $vbf_file does not exist.\n";
	return 0;
    }
    
    if($output_file =~ /(.+)_s1.h5/) 
    {
	$cfg{'o'}="$opt{'s2_dir'}$output_file";
	$log_file = "$opt{'s2_dir'}$1_laser.log";
    }
    else
    {
	print "ERROR Cannot parse laser output file: $output_file\n";
	return 0;
    }

    if(-e $log_file && !$opt{'overwrite'})
    {
	open(LOG,">>$log_file");
	my $locked = !flock(LOG, LOCK_EX | LOCK_NB);   
	close(LOG);
	if($locked) { return 0; } 
	else { return 1; }
    }
    elsif($opt{'overwrite'})
    {
	`rm $cfg{'o'}`;
	`rm $log_file`;
    }

    open(LOG,">>$log_file");
    return 0 if !flock(LOG, LOCK_EX | LOCK_NB);    

    if($opt{'verbose'})
    {
	$log_settings = "2>&1 | tee $log_file";
    }
    else
    {
	$log_settings = "&> $log_file";
    }

    # Merge DB options --------------------------------------------------------
    @cfg{ keys %cfg_db } = values %cfg_db;

    my $cfg = "";
    while(my ($opt,$val) = each(%cfg)) 
    {
	next if $val eq "";
	$cfg = $cfg . "-$opt=$val "; 
    }

    print "----------------------------------------------------------------\n";
    print "Analyzing Laser Run\n";
    print "Path................: $vbf_file\n";
    print "Log File............: $log_file\n";
    print "Stage1 File.........: $output_file\n";
    print "Command: $laser_bin $cfg $vbf_file\n";
    system("$laser_bin $cfg $vbf_file $log_settings");
    return 1;
}

sub find_vbf_file
{
    my $vbf_file = $_[0];
    my $date = $_[1];

    if(-e $vbf_file)
    {
	return $vbf_file;
    }
    elsif(-e "/veritas/data/d${date}/${vbf_file}.cvbf")
    {
	return "/veritas/data/d${date}/${vbf_file}.cvbf";
    }

    if(!scalar(keys(%vbf_file_hash)))
    {
	build_vbf_hash();
    }

    if(exists $vbf_file_hash{$vbf_file}) 
    {
	return $vbf_file_hash{$vbf_file};
    } 
    else 
    {
	return $vbf_file;
    }
}

sub build_vbf_hash
{
    # -------------------------------------------------------------------------
    # Build a hash of vbf file paths
    # -------------------------------------------------------------------------
    my @vbf_data_dir = ("/veritas/data");
#     my @vbf_data_dir = ("/net/gamma3/srv/raid2/data",
# 		    "/net/gamma3/srv/raid3/data",
# 		    "/net/gamma3/srv/raid6/data");

    foreach my $data_dir (@vbf_data_dir)
    {
	my @vbf_files = <$data_dir/d*/*cvbf>;
	
	foreach my $vbf_file (@vbf_files) 
	{    
	    $vbf_file =~ m/(\d+).cvbf/;
	    my $runno = $1;
	    
	    $vbf_file_hash{$runno} = $vbf_file;
	}    
    }
}

sub mask_to_scopes
{
    my @mask = split(//,$_[0]);
    my $nscope = scalar(@mask);
    my @scopes;

    for(my $iscope = 0; $iscope < $nscope; $iscope++) 
    {
	if($mask[$iscope] eq "1")
	{
	    push(@scopes,$iscope); 
	}
    }

    return join(',',@scopes);
}

sub get_suppressed_scopes 
{
    my @scopes = split(/,/,$_[0]);
    my @suppressed_scopes;

    for(my $iscope = 0; $iscope < 4; $iscope++) 
    {
	if( !(grep /^$iscope$/,@scopes) )
	{
	    push(@suppressed_scopes,$iscope);
	}
    }

    return @suppressed_scopes;
}

sub get_allowed_trigger_masks
{
    my $mask = scopes_to_mask($_[0]);
    my $multiplicity = $_[1];
    my @trigger_mask_list;

    for(my $pattern = 0; $pattern < 16; $pattern++) 
    {
	if(count_bits($mask & $pattern) >= $multiplicity)
	{
	    my $trigger_mask = "";

	    for(my $scope = 0; $scope < $max_telescopes; $scope++)
	    {
		if((1<<$scope) & $pattern)
		{
		    $trigger_mask = $trigger_mask . $scope;
		}
	    }

	    push(@trigger_mask_list,$trigger_mask);
	}
    }

    return @trigger_mask_list;
}

sub scopes_to_mask
{
    my @scopes = split(/,/,$_[0]);
    my $mask = 0;
    
    for(my $i = 0; $i < scalar(@scopes); $i++)
    {
 	$mask |= (0x1<<($scopes[$i]));
    }
      
    return $mask;
}

sub string_to_mask
{
    my @mask = split(//,$_[0]);
    my $nbits = scalar(@mask);

    my $mask = 0;

    for(my $i = 0; $i < $nbits; $i++)
    {
	if($mask[$i] eq "1")
	{
	    $mask |= (0x1<<($nbits-$i-1));
	}
    }
      
    return $mask;
}

sub count_bits 
{
    my $mask = $_[0];

    my $count = 0;
    for($count = 0; $mask; $count++)
    {
	$mask &= $mask-1;
    }

    return $count;
}

sub load_default_cfg
{
    my $cfg = shift;
    my $default_cfg = shift;

    foreach my $opt (keys(%$cfg))
    {
	if(!defined($cfg->{$opt}) && exists($default_cfg->{$opt}))
	{
	    $cfg->{$opt} = $default_cfg->{$opt};
	}
	elsif(!defined($cfg->{$opt}))
	{
	    delete($cfg->{$opt});
	}
    }
}


__END__

=head1 NAME

run_stage2.pl - Script for running chila stage1/stage2 data analysis.

=head1 DESCRIPTION

This script runs the stage1 and stage2 ChiLA analysis on a vbf file or list
of vbf files.  If using a file list, it should have the following format:

<YYYYMMDD> <data run #> <laser run #>

=head1 SYNOPSIS

run_stage2.pl [options] [vbf_file]

run_stage2.pl [options] [file_list]

Options:

=over 8

=item B<--help>

Print this help message.

=item B<--man>

Print the full man page.

=item B<--cfg>

Specify a configuration file from which command line options will be
read.

=item B<--sim>

Simulation file option.  Turns off DB access.

=item B<--laser>

Set the run number or path for laser file.

=item B<--mask>

Set the mask of telescopes to be analyzed (e.g. 1111).

=item B<--multiplicity>    

Apply an array trigger multiplicity in analysis.

=item B<-v>

Enable verbose output.  Print output of stage2 to STDOUT/STDERR in
addition to the log file.

=back
