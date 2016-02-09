#!/usr/bin/perl -w
#
# \author     Matthew Wood                \n
#             UCLA                        \n
#             mdwood@astro.ucla.edu       \n        
#
# \version    0.1
# \date       06/24/2008
#$ -e $HOME/logs
#$ -o $HOME/logs

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

my $opt_help = 0;
my $opt_man = 0;
my $opt_target_dir = "./";
my $opt_work_dir = "./";
my $opt_code_dir = "$ENV{'HOME'}/VS";
my $opt_njob = -1;
my $opt_log = 0;
my $opt_event_use = 1;
my $opt_vsdb_host = "";
my $opt_vsdb_user = "";

GetOptions( 'h|help!'        => \$opt_help,
	    'man!'           => \$opt_man,
	    'log!'           => \$opt_log,
            'target_dir=s'   => \$opt_target_dir,
            'work_dir=s'     => \$opt_work_dir,
            'code_dir=s'     => \$opt_code_dir,
            'vsdb_host=s'    => \$opt_vsdb_host,
            'vsdb_user=s'    => \$opt_vsdb_user,
            'njob=i'         => \$opt_njob,
	    'event_use=i'    => \$opt_event_use
            ) or pod2usage(2);

# -----------------------------------------------------------------------------
# Print help/man if asked or if no file
# -----------------------------------------------------------------------------
pod2usage(1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;
pod2usage("$0: No database specified.") if ((@ARGV == 0) && (-t STDIN));

$ENV{'LD_LIBRARY_PATH'}=
    "$ENV{'LD_LIBRARY_PATH'}:" .
    "$ENV{'HOME'}/local/lib:$ENV{'HOME'}/local/lib/mysql";

my @DATABASES=@ARGV;
my $TARGET_DIR_BASE = File::Spec->rel2abs( $opt_target_dir );
my $WORKDIR = File::Spec->rel2abs( $opt_work_dir );
my $CODE_DIR = File::Spec->rel2abs( $opt_code_dir );

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

print "ARCH:           $ARCH\n";
print "NODE:           $NODE\n";
print "HOSTNAME1:      $HOSTNAME1\n";
print "HOSTNAME2:      $HOSTNAME2\n";
print "IP:             $IP\n";      

#
# CODE AND TEMPORARY LOCATION VARIABLES
#

if($HOSTNAME1 =~ /hoffman2/ || $HOSTNAME2 =~ /hoffman2/)
{
    $ENV{'VSDB_HOST'}="login3";
    $ENV{'VSDB_USER'}="simulation";
    $WORKDIR=$ENV{'TMPDIR'};
    $TARGET_DIR_BASE = "/u/home9/vvv/mdwood/simulations";
    $HOSTNAME = $HOSTNAME1;
}
else
{
    $ENV{'VSDB_HOST'}="gamma5.astro.ucla.edu";
    $ENV{'VSDB_USER'}="simulation";
    $HOSTNAME = $HOSTNAME2;
}

if($opt_vsdb_user ne "")
{
    $ENV{'VSDB_USER'}=$opt_vsdb_user;
}

if($opt_vsdb_host ne "")
{
    $ENV{'VSDB_USER'}=$opt_vsdb_user;
}

if(!(-d $CODE_DIR))
{
    print "Error. Could not find code directory: $CODE_DIR\n";
    exit(1);
}

#
# MISC VARIABLES
#

my $CORSIKA="$CODE_DIR/corsika";
my $MAKE_STEERING="$CODE_DIR/make_steering";
my $MANAGE_SIMDB="$CODE_DIR/manage_simdb";
my $RAYTRACE="$CODE_DIR/raytrace_corsika_to_simdb";
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

#
# MAKE TEMPORARY DIRECTORY BASE
#

if(!(-d $WORKDIR))
{
    mkdir $WORKDIR, 0770 or die "Can't mkdir $WORKDIR: $!";
}

if(! -d $TARGET_DIR_BASE)
{
    mkdir $TARGET_DIR_BASE, 0770 or die "Can't mkdir $TARGET_DIR_BASE: $!";
}

foreach my $DATABASE (@DATABASES)
{
    my $TARGET_DIR = "$TARGET_DIR_BASE/$DATABASE";

    if($TARGET_DIR ne "" && !(-d $TARGET_DIR))
    {
	mkdir $TARGET_DIR, 0770 or die "Can't mkdir $TARGET_DIR: $!";
    }
}

if($opt_log)
{
    my $log_file = `mktemp $WORKDIR/vs_log.XXXXXX`;
    chomp($log_file);

    open STDERR, ">$log_file";
    open STDOUT, ">$log_file";
}




print "JOB_ID:         $JOB_ID\n";
print "TASK_ID:        $TASK_ID\n";
print "CODE_DIR:       $CODE_DIR\n";
print "DATABASES:      @DATABASES\n";
print "WORKDIR:        $WORKDIR\n";
print "TARGET_DIR:     $TARGET_DIR_BASE\n";

print `date`;

# -----------------------------------------------------------------------------
# MAIN LOOP
# -----------------------------------------------------------------------------

my $workunit_table=`$MANAGE_SIMDB -get_workunit $DATABASES[0]`;
chomp($workunit_table);

my $njob = $opt_njob;
while( $workunit_table ne "" && $njob != 0 )
{
    # MAKE TEMPORARY DIRECTORY ------------------------------------------------
    my $RUNDIR=`mktemp -d $WORKDIR/vs_run.XXXXXXXXXX`;
    chomp($RUNDIR);
    chdir($RUNDIR) or die "Can't chdir $RUNDIR: $!";

    my @workunit_run_id;

    foreach my $DATABASE (@DATABASES)
    {
	# REGISTER WORKUNIT START ---------------------------------------------
	my $START_WORKUNIT_ARGS = 
	    "-start_workunit_run $DATABASE $workunit_table " .
	    "$HOSTNAME \"$JOB_ID\"";

	print "$MANAGE_SIMDB $START_WORKUNIT_ARGS\n";
	my $workunit_run_id= `$MANAGE_SIMDB $START_WORKUNIT_ARGS`;
	chomp($workunit_run_id);

	if($workunit_run_id !~ /^[0-9]+$/)
	{
	    print "Error: Invalid workunit run id: $workunit_run_id " .
		" for database $DATABASE.\n";
	    exit(1);
	}

	push @workunit_run_id, $workunit_run_id;
    }    

    # MAKE STEERING FILE ------------------------------------------------------
    my $MAKE_STEERING_ARGS = "-INIT -db $DATABASES[0] -table $workunit_table ";
    if($opt_event_use > 1)
    {
	$MAKE_STEERING_ARGS = 
	    $MAKE_STEERING_ARGS . " -cerenkov_event_use=$opt_event_use";
    }
    print "$MAKE_STEERING $MAKE_STEERING_ARGS > steering.dat\n";
    print `$MAKE_STEERING $MAKE_STEERING_ARGS > steering.dat`;

    # CHECK SIZE OF STEERING FILE ---------------------------------------------
    my $size = -s "steering.dat";
    if($size == 0)
    {
	print "Error: Empty steering file.\n";
	exit(1);
    }

    # RUN CORSIKA -------------------------------------------------------------
    print `$CORSIKA < steering.dat`; 

    for(my $i = 0; $i < scalar(@DATABASES); $i++)
    {
	raytrace($DATABASES[$i],$workunit_run_id[$i]);
    }

    # CLEAN UP TEMPORARIES ----------------------------------------------------
    chdir($WORKDIR) or die "Can't chdir $WORKDIR: $!";
    print `rm -rf $RUNDIR`;

    $workunit_table=`$MANAGE_SIMDB -get_workunit $DATABASES[0]`;
    chomp($workunit_table);
    $njob = $njob-1;
}

print `date`;

sub raytrace 
{
    my $DATABASE = $_[0];
    my $workunit_run_id = $_[1];

    # OUTPUT FILE -------------------------------------------------------------
    my $telfil = `awk '\$1=="TELFIL"{print \$2}' steering.dat`;
    chomp($telfil);

    $telfil =~ s/^.//;
    $telfil =~ s/.$//;

    # PROCESS OPTICS ----------------------------------------------------------
    my $RAYTRACE_OPT = 
	"--table=$workunit_table --workunit_runid=$workunit_run_id $telfil";

    print "$RAYTRACE $DATABASE $RAYTRACE_OPT\n";
    print `$RAYTRACE $DATABASE $RAYTRACE_OPT`;

    # IF OUTPUT IS TO AN HDF FILE THEN COPY WORKUNIT TO DIRECTORY -------------
    my $data_storage=
	`$MANAGE_SIMDB -get_parameter $DATABASE DataStorage mode`;
    chomp($data_storage);

    my $TARGET_DIR = "$TARGET_DIR_BASE/$DATABASE";

    if( $data_storage eq "hdf5" )
    {
	if($TARGET_DIR eq "")
	{
	    $TARGET_DIR = 
		`$MANAGE_SIMDB -get_parameter $DATABASE DataStorage directory`;
	    chomp($TARGET_DIR);
	}

	my $target_file=`$MANAGE_SIMDB -get_hdf_name $workunit_run_id`;
	chomp($target_file);
	print `chmod g+w $target_file`;
	print `mv $target_file $TARGET_DIR`;
	print `$MANAGE_SIMDB -mark_all_workunit_events_as_complete $DATABASE $workunit_table $workunit_run_id`;
    }
    
    # REGISTER COMPLETION OF WORKUNIT -----------------------------------------
    print "$MANAGE_SIMDB -finish_workunit_run $DATABASE $workunit_run_id\n";
    print `$MANAGE_SIMDB -finish_workunit_run $DATABASE $workunit_run_id`;
}


=head1 NAME

vs_run_sge.pl - Script for running corsika/vsoptics simulation chain.  

=head1 DESCRIPTION

This script can be used to submit batch jobs to the sge queue system
or to run the simulation chain locally.  To run locally use the following
syntax:

vs_run_sge.pl --work_dir=[work_dir] --target_dir=[target_dir] [database]

The B<work_dir> option specifies the location where temporary files will
be written.  After generating each workunit the script will copy the
output simulation files to the directory <target_dir>/<database>.  To
submit as an sge job use the following syntax:

qsub -t 1-N vs_run_sge.pl [database] -njob=1

This will submit N sge tasks each of which will complete 1 simulation
workunit before exiting.  The number of workunits per task can be
increased by changing the 'njob' option.  The number of workunits per
task should be set so as not to exceed the maximum allowed execution
time of the queue.

=head1 SYNOPSIS

B<vs_run_sge.pl> [options] [database]

Options:

=over 8

=item B<--help>

Print this help message.

=item B<--man>

Print the full manual page for this script.

=item B<--log>

Write STDOUT/STDERR to a log file.  This option should not be used
when running on the cluster as STDOUT/STDERR are already logged by the
batch system

=item B<--njob>=N

Set the number of jobs to run before exiting.  If this is set to -1 then
the script will run indefinitely until all workunits are complete.  When
running on the cluster 'njob' should generally be ~1-10 depending on the
time required to complete a single workunit.

=item B<--event_use>=N

Set the event reuse number used by corsika.  This specifies the number
of times the core position for each shower will be rethrown.  

=item B<--target_dir>

Set the output directory for completed workunit files.

=item B<--work_dir>

Set the working directory where temporary files will be written.

=item B<--code_dir>

Set the location of the directory containing the simulation binaries.

=item B<--vsdb_host>

Set the simulation database hostname.

=item B<--vsdb_user>

Set the simulation database username.


=back




