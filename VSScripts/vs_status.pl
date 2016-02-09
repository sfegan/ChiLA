#!/usr/bin/perl
# 
# Example of how to run this:
# 
#   > vs_status.pl base_g_wob050_2009b -host=gamma5
#

use strict;
use FileHandle;
use Getopt::Long;
use Pod::Usage;

my $help = 0;
my $host = "gamma5.astro.ucla.edu";
my $user = "mdwood";
my $database = "base_g_wob100_2008b_survey";
my $password = "";

GetOptions( 'help|h!'          => \$help,
	    'host=s'         => \$host,
	    'user=s'         => \$user
	    ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage("$0: No database specified.") if ((@ARGV == 0) && (-t STDIN));

$database = $ARGV[0];

my $mysql_options = 
    " mysql -h $host -u $user --password=$password $database -ABN -e ";

my $query = 
    "\"SELECT TableID, EnergyGeV, ZenithMinRad, TargetEventCount ".
    "FROM VSD_TableDirectory\"";

my $mysql_handle = new FileHandle($mysql_options.$query.' |');
die $! if(!$mysql_handle);

my $total_workunits = 0;

while(<$mysql_handle>)
{
    chomp;
    my ($table_id, $energy_gev, $zenith_min_rad, $table_event_count) = 
	split /\s+/;

    my $zenith_deg = $zenith_min_rad*180/3.14159;

#    next if $energy_gev < 100;
#    next if $zenith_deg > 45;
#    print "$table_id $energy_gev $table_event_count\n";

    $total_workunits = $total_workunits + 10;
}


$query = 
    "\"SELECT TableID, DATE(StartTime), DATE(FinishTime), " .
    "TIME_TO_SEC(TIMEDIFF(FinishTime,StartTime)), WorkunitRunComplete ".
    "FROM VSD_WorkunitRun ORDER BY FinishTime\"";

$mysql_handle = new FileHandle($mysql_options.$query.' |');
die $! if(!$mysql_handle);

my %ncomplete;
my %nincomplete;

my $complete_workunits = 0;

my %complete_day;
my %incomplete_day;
my %duration_day;

while(<$mysql_handle>)
{
    chomp;
    my ($table_id, $start_time,
	$end_time,$duration,$workunit_complete) = split /\s+/;

    if($workunit_complete == 0)
    {
	$incomplete_day{$start_time} = $incomplete_day{$start_time}+1;
	$nincomplete{$table_id}=$nincomplete{$table_id}+1;
    }

    next if $workunit_complete != 1;

    $ncomplete{$table_id}=$ncomplete{$table_id}+1;

    if($ncomplete{$table_id} <= 10)
    {
	$complete_workunits = $complete_workunits+1;
    }

    $complete_day{$end_time} = $complete_day{$end_time}+1;
    $duration_day{$end_time} += $duration;

#    print "$table_id $end_time\n";
}

my $percent_complete = $complete_workunits/$total_workunits;

print sprintf("%15s","DATE") . 
    sprintf("%15s","COMPLETE") .     
    sprintf("%15s","INCOMPLETE") .     
    sprintf("%25s","WORKUNIT DURATION [HR]") . 
    "\n";


foreach my $day (sort(keys(%complete_day)))
{
    my $avg_duration = $duration_day{$day}/$complete_day{$day};

    print sprintf("%15s",$day) .
	sprintf("%15i",$complete_day{$day}) . 
	sprintf("%15i",$incomplete_day{$day}) . 
	sprintf("%25.3f",$avg_duration/3600.) . 
	"\n";
}

foreach my $table_id (sort(keys(%ncomplete)))
{
#    print "$table_id $ncomplete{$table_id}\n";
}

my $remaining_workunits = $total_workunits - $complete_workunits;

print "\n";

print "TOTAL WORKUNITS:    " . sprintf("%10i\n",$total_workunits);
print "COMPLETE/REMAINING: " .
sprintf("%10i/%10i\n",$complete_workunits,$remaining_workunits);
print "% COMPLETE:         " . sprintf("%10.5f\n",$percent_complete);
