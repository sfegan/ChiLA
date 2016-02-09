#!/usr/bin/perl -w
#
# Program:     vse_batch.pl
# Author:      Matthew Wood <mdwood@astro.ucla.edu>
# Date:        12/12/07
#
# Description: Script for running vse batch jobs.
#
# $Id: vse_run_sge.pl,v 1.9 2010/09/30 19:40:57 matthew Exp $
#
#$ -e $HOME/vse_logs
#$ -o $HOME/vse_logs

BEGIN 
{ 
    # Turn on autoflush
    $| = 1;

    if(!exists($ENV{'CHILADIR'}))
    {
	$ENV{'CHILADIR'} = "$ENV{'HOME'}/ChiLA";
    }
    elsif(!(-d $ENV{'CHILADIR'}))
    {
	print "CHILADIR is not a valid directory: $ENV{'CHILADIR'}\n";
	exit(1);
    }

    push @INC, "$ENV{'CHILADIR'}/VSScripts";
}

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use threads;
use threads::shared;
use ThreadPool;
use Fcntl qw(:flock);

# -----------------------------------------------------------------------------
# Parse the Options
# -----------------------------------------------------------------------------
my %opt = ();
$opt{'help'} = 0;
$opt{'man'} = 0;
$opt{'tag'} = "";
$opt{'nthread'} = 1;
$opt{'vsebin'} = "";
$opt{'database'} = "";
$opt{'work_dir'} = "/work/mdwood";
$opt{'target_dir'} = "/u/work/VERITAS/vbf";
$opt{'njob'} = 1;
$opt{'az'} = "000,090,180,270";
$opt{'zn'} = "00.0,12.0,19.0,26.0,34.0,41.0,49.0,56.0";
$opt{'nsb'} = "0.07,0.12";
$opt{'config'} = "";
$opt{'nevents'} = "4.0";

GetOptions( 'h|help!'        => \$opt{'help'},
            'man!'           => \$opt{'man'},
            'mode=s'         => \$opt{'mode'},
            'tag=s'          => \$opt{'tag'},
	    'nthread=i'      => \$opt{'nthread'},
	    'vsebin=s'       => \$opt{'vsebin'},
	    'work_dir=s'     => \$opt{'work_dir'},
	    'target_dir=s'   => \$opt{'target_dir'},
	    'database=s'     => \$opt{'database'},
	    'njob=i'         => \$opt{'njob'},
	    'az=s'           => \$opt{'az'},
	    'zn=s'           => \$opt{'zn'},
	    'nsb=s'          => \$opt{'nsb'},
	    'config=s'       => \$opt{'config'},
            ) or pod2usage(2);


pod2usage(1) if $opt{'help'};
pod2usage(-verbose => 2) if $opt{'man'};
#pod2usage("$0: No job list specified.") if ((@ARGV < 1) && (-t STDIN));

if(!$opt{'vsebin'} || ! -e $opt{'vsebin'})
{
    print STDERR "Error: Could not find vse binary: $opt{'vsebin'}";
    exit(1);
}

$ENV{'LD_LIBRARY_PATH'}=
    "$ENV{'HOME'}/local/lib:$ENV{'HOME'}/local/lib/mysql:" .
    "$ENV{'HOME'}/veritas/lib";

my $WORKDIR = File::Spec->rel2abs( $opt{'work_dir'} );
my $vsebin = $opt{'vsebin'};
my $vbf_compress = "/u/home9/vvv/mdwood/ChiLA/VSMiscAnalysis/vbf_compress";
$opt{'work_dir'} =~ s/\/$//;
$opt{'target_dir'} =~ s/\/$//;

#
# ARCHITECTURE AND NODE VARIABLES
#

my $ARCH=`uname -m`;
chomp($ARCH);
my $NODE=`hostname -s`;
chomp($NODE);
my $IP=`hostname -i`;
chomp($IP);
my $HOSTNAME1=`hostname -a`;
chomp($HOSTNAME1);
my $HOSTNAME2=`hostname`;
chomp($HOSTNAME2);
my $HOSTNAME;
my $JOB_ID="";
my $TASK_ID="";

if(exists($ENV{'JOB_ID'}))
{
    $JOB_ID = $ENV{'JOB_ID'};
}

if(exists($ENV{'SGE_TASK_ID'}))
{
    $TASK_ID = $ENV{'SGE_TASK_ID'};
}

if($HOSTNAME1 =~ /hoffman2/ || $HOSTNAME2 =~ /hoffman2/)
{
    $WORKDIR=$ENV{'TMPDIR'};
    $HOSTNAME = $HOSTNAME1;
}
else
{
    $HOSTNAME = $HOSTNAME2;
}

print "ARCH:           $ARCH\n";
print "NODE:           $NODE\n";
print "HOSTNAME:       $HOSTNAME\n";
print "IP:             $IP\n";      
print "JOB_ID:         $JOB_ID\n";
print "TASK_ID:        $TASK_ID\n";
print "WORKDIR:        $WORKDIR\n";
print "OUTPUTDIR:      $opt{'target_dir'}\n";

if(! -d $WORKDIR) 
{ 
    mkdir($WORKDIR, 0755) or die "Can't mkdir $WORKDIR: $!";; 
}

my $sim_table = "";

if(scalar(@ARGV) > 0)
{
    $sim_table = $ARGV[0];
}

my @jobs;
my %cfg;
my @cfg;

if(-e $sim_table)
{
    load_table(\@cfg,$sim_table);
}
else
{
    generate_jobs(\@cfg);
}

$cfg{'cfg'} = \@cfg;

for( my $n = 0; $n < $opt{'nthread'}; $n++)
{
    push(@jobs,\%cfg);        
}


ThreadPool::run_jobs(\&run_vse,\@jobs,$opt{'nthread'});

sub get_base_config
{
    my %config;

    my $config_file = $_[0];

    open(BASE_CONFIG,"$config_file") or die "Can't open $config_file: $!";
    
    while(<BASE_CONFIG>)
    {
	chomp;
	next if /^\#/ || /^$/;
	my @line = split ' ',$_;
	$config{$line[0]} = $line[1];
    }

    close(BASE_CONFIG);

    return %config;
}

sub create_filename
{
    my %config = %{$_[0]};

    my $filename;

    my $nevent;

    if($config{'evgNEvents'} >  100E6)
    {
	$nevent = sprintf("%06.2f",$config{'evgNEvents'}/1E6);
    }
    else
    {
	$nevent = sprintf("%05.2f",$config{'evgNEvents'}/1E6);
    }

    if($config{'evgMode'} eq "PL_SPECTRUM")
    {
	my $spectrum_index = $config{'evgSpectrumIndex'};
	$spectrum_index =~ s/\-//;

	$filename = 
	    $config{'evgDBDatabase'} . "_" .
	    "PL" . $spectrum_index . "_" .
	    "Zn" . $config{'evgTargetZenith'} . "_";

	if($config{'evgTargetZenithOnly'} ne "true")
	{
	    $filename = $filename . "Az" . $config{'evgTargetAzimuth'} . "_";
	}

	$filename = $filename . 
	    $opt{'tag'} . "_" .
	    "NSB" . $config{'calNSBRate'} . "_" .
	    $nevent . "M";
    }
    elsif($config{'evgMode'} eq "HE_WEIGHTED")
    {
	$filename = 
	    $config{'evgDBDatabase'} . "_" .
	    "HE_WTD_" .
	    "Zn" . $config{'evgTargetZenith'} . "_";

	if($config{'evgTargetZenithOnly'} ne "true")
	{
	    $filename = $filename . "Az" . $config{'evgTargetAzimuth'} . "_";
	}
	
	$filename = $filename . 
	    $opt{'tag'} . "_" .
	    "NSB" . $config{'calNSBRate'} . "_" .
	    $nevent . "M";
    }
    elsif($config{'evgMode'} eq "ONE_TABLE")
    {
	$filename = 
	    $config{'evgDBDatabase'} . "_" .
	    $config{'evgDBTable'} . "_" .
	    $opt{'tag'} . "_" .
	    "NSB" . $config{'calNSBRate'} . "_" .
	    $nevent . "M";
    }

    return $filename;
}

sub generate_jobs
{
    my @az = split ',',$opt{'az'};
    my @zn = split ',',$opt{'zn'};
    my @nsb = split ',',$opt{'nsb'};
    my @db = split ',',$opt{'database'};

    my $njob = scalar(@az)*scalar(@zn)*scalar(@nsb)*scalar(@db);

    for(my $i = 0; $i < $njob; $i++)
    {
	my %cfg;

	my $j = $i;

	my $izn = $j % scalar(@zn);
	$j = $j/scalar(@zn);

	my $iaz = $j % scalar(@az);
	$j = $j/scalar(@az);

	my $insb = $j % scalar(@nsb);
	$j = $j/scalar(@nsb);

	my $idb = $j % scalar(@db);

	my %options = get_base_config($opt{'config'});

	$options{'evgDBDatabase'} = $db[$idb];
	$options{'evgMode'} = "HE_WEIGHTED";	
	$options{'evgTargetZenith'} = $zn[$izn];
	$options{'evgTargetAzimuth'} = $az[$iaz];
	$options{'calNSBRate'} = $nsb[$insb];
	$options{'evgNEvents'} = 1E6*$opt{'nevents'};
	$options{'evgTargetZenithOnly'} = "false";

	my $filename = create_filename(\%options);

	$cfg{'options'} = \%options;
	$cfg{'vbf_file'} = "$filename.vbf";
	$cfg{'config_file'} = "$opt{'target_dir'}/$filename.cfg";
	$cfg{'log_file'} = "$opt{'target_dir'}/$filename.log";

	next if $zn[$izn] == 0 && $az[$iaz] != 0;

	push @{$_[0]}, \%cfg;
    }
}

sub load_table
{
    my $table_file = $_[1]; 

    open(TABLE,$sim_table) or die "Can't open $sim_table: $!";
    while(<TABLE>)
    {
	chomp;
	next if /^\#/ || /^$/;
	my @line = split ' ',$_; 
	my $base_config;
	my %cfg;

	if(! -e $line[0])
	{
	    print STDERR "Could not find config file: $line[0]\n";
	    next;
	}

	$base_config = shift(@line);
	
	if($base_config =~ /([^\/]+).cfg/ && $opt{'tag'} eq "")
	{
	    $opt{'tag'} = $1;
	}

	my %options = get_base_config($base_config);

	if($line[1] eq "HE_WEIGHTED")
	{
	    ($options{'evgDBDatabase'},
	     $options{'evgMode'},
	     $options{'evgSpectrumIndex'},
	     $options{'evgTargetZenith'},
	     $options{'evgTargetAzimuth'},
	     $options{'calNSBRate'},
	     $options{'evgNEvents'}) = @line;
	    
	    $options{'evgSpectrumIndex'}*=-1;
	}
	elsif($line[1] eq "PL_SPECTRUM")
	{
	    ($options{'evgDBDatabase'},
	     $options{'evgMode'},
	     $options{'evgSpectrumIndex'},
	     $options{'evgTargetZenith'},
	     $options{'evgTargetAzimuth'},
	     $options{'calNSBRate'},
	     $options{'evgNEvents'}) = @line;

	    $options{'evgSpectrumIndex'}*=-1;
	}
	else
	{
	    die "Error parsing line: \n $_\n";
	}

	if($options{'evgTargetAzimuth'} eq "*")
	{
	    $options{'evgTargetZenithOnly'} = "true";
	}
	else
	{
	    $options{'evgTargetZenithOnly'} = "false";
	}

	if($opt{'database'})
	{
	    $options{'evgDBDatabase'} = $opt{'database'};
	}	

	$options{'evgNEvents'}*=1E6;

	my $filename = create_filename(\%options);		    
	
	$cfg{'options'} = \%options;
	$cfg{'vbf_file'} = "$filename.vbf";
	$cfg{'config_file'} = "$opt{'target_dir'}/$filename.cfg";
	$cfg{'log_file'} = "$opt{'target_dir'}/$filename.log";

	push @{$_[0]}, \%cfg;
    }
    close(TABLE);
}

sub run_vse
{
    my %arg = %{$_[0]};
    my $njob = $opt{'njob'};

    # MAKE TEMPORARY DIRECTORY ------------------------------------------------
    my $TMPDIR=`mktemp -d $WORKDIR/vs_run.XXXXXXXXXX`;
    chomp($TMPDIR);

    chdir($TMPDIR) or die "Can't chdir $TMPDIR: $!";

    while($njob != 0 && scalar(@cfg) > 0)
    {
	my %cfg = %{pop(@{$arg{'cfg'}})};

	next if -e $cfg{'config_file'};

	open(CFG,">>$cfg{'config_file'}") or 
	    die "Can't open $cfg{'config_file'}: $!";
	next if !flock(CFG, LOCK_EX | LOCK_NB);

	my $vbf_file = "$TMPDIR/$cfg{'vbf_file'}";
	my $cfg_file = "$cfg{'config_file'}";

	$cfg{'options'}->{'vbfOutputFile'} = $vbf_file;

	foreach my $opt ( sort(keys(%{$cfg{'options'}})) )
	{
	    print CFG sprintf("%-30s",$opt) . "$cfg{'options'}->{$opt}\n";
	}

	print `date`;
	print `pwd`;
	print "$vsebin -config=$cfg_file 2>&1 | tee $cfg{'log_file'}\n";
	print `$vsebin -config=$cfg_file 2>&1 | tee $cfg{'log_file'}`;
	print "gzip $vbf_file\n";
	print `gzip $vbf_file`;
	print "mv ${vbf_file}.gz $opt{'target_dir'}\n";
	print `mv ${vbf_file}.gz $opt{'target_dir'}`;
	print `date`;

	close(CFG);

	$njob = $njob - 1;
    }

    # CLEAN UP TEMPORARIES ----------------------------------------------------
    chdir($WORKDIR) or die "Can't chdir $WORKDIR: $!";
    print `rm -rf $TMPDIR`;
}


__END__

=head1 NAME

vse_run_sge.pl - Script for running vse detector simulation.

=head1 DESCRIPTION

This script automates running a set of simulations with the vse
detector code either locally or on the cluster.  To run locally use
the following syntax:

vse_run_sge.pl --work_dir=[work_dir] --target_dir=[target_dir]
--nthread=X --njob=-1 --vsebin=[vsebin] [list]

The B<nthread> option can be used to start X threads each of which
will run simulation jobs in parallel.  Setting B<njob> option to -1
causes each script instance (thread) to run indefinitely until all
simulation jobs are complete.  The B<work_dir> option specifies the
location where temporary files will be written.  

To submit this script as an sge job use the following syntax:

qsub -t 1-N vse_run_sge.pl [list] --vsebin=[vsebin] --njob=1

This will submit N sge tasks each of which will run one job before
exiting.  Note that the number of tasks can be larger than the
number of simulation jobs in the list.

The list of simulation jobs to be generated is provided to the script
as an input text file with each line corresponding to a single
simulation job.  An instance of the script will open the list file and
find the first configuration which is not yet completed by checking
for the existence of the completed vbf file in the destination
directory (set with B<--target_dir>).  The columns of the input list
define the configuration to be used for each file and should have the
following format:

B<Sim Table Format:>

<cfgFile> <evgDBDatabase> <evgMode> <evgSpectrumIndex> <evgTargetZenith> <evgTargetAzimuth> <calNSBRate> <evgNEvents>

B<cfgFile:> Configuration file which can be used to override default
vse options which are not already specified by the other columns in
the table.

B<evgDBDatabase:> Sets the simulation database which will be
substituted for the evgDBDatabase option.

B<evgMode:> Sets the mode controlling the generated energy distribution.

B<evgSpectrumIndex:> Set the spectral index used when running with
the evgMode=PL_SPECTRUM.

B<evgTargetZenith:> Set the target zenith angle.  vse will select the database
pointing that is closest to this value.

B<evgTargetAzimuth:> Set the target azimuth angle.  vse will select
the database pointing that is closest to this value.

B<calNSBRate:> Set the NSB photon rate.

B<evgNEvents:> Set the maximum number of events to simulate in millions.



=head1 SYNOPSIS

vse_run_sge.pl [options] [list] 

Options:

=over 8

=item B<--target_dir>

Set the output directory for completed vbf files.

=item B<--vsebin>

Set the path to the vse binary.

=item B<--njob>=N

Set the number of jobs to run before exiting.  If this is set to -1 then
the script will run indefinitely until all workunits are complete.  

=item B<--target_dir>

Set the output directory for completed vbf files.

=item B<--work_dir>

Set the working directory where temporary files will be written.  This
option should not be used when running on the cluster.

=item B<--nthread>=N

Set the number of parallel threads to initiate.  Note that each thread
has an independent job counter.   This option should not be used when running
on the cluster.

=back
