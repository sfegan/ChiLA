#!/usr/bin/perl

use DBI;

my $host = "vraid.neutrino.hosted.ats.ucla.edu";
my $database = $ARGV[0];

my $dbh = DBI->connect("DBI:mysql:$database:$host", 
		       "mdwood", "*Nizamu")
    or die "Could not connect to DB: !\n";

my $query_tables = 
    "SELECT TableID, TableName FROM VSD_TableDirectory;";

my %tablename;

my $stmt = $dbh->prepare($query_tables)
    or die "Couldn't prepare statement: " . $dbh->errstr;
$stmt->execute();

while(my @data = $stmt->fetchrow_array())
{
    $tablename{$data[0]} = $data[1];
}

undef $stmt;

my $query_workunits = 
    "SELECT WorkunitRunID,TableID,WorkunitRunComplete,Hostname " .
    "FROM VSD_WorkunitRun";

$stmt = $dbh->prepare($query_workunits)
    or die "Couldn't prepare statement: " . $dbh->errstr;
$stmt->execute();

my %workunits;

while(my ($workunit_run_id,$table_id,$workunit_complete,$hostname) = 
      $stmt->fetchrow_array())
{
    $workunits{$workunit_run_id}->{'table_id'} = $table_id;
    $workunits{$workunit_run_id}->{'workunit_complete'} = $workunit_complete;
}

undef $stmt;

foreach my $workunit_run_id (sort { $a <=> $b } (keys(%workunits)))
{
    my $table_id = $workunits{$workunit_run_id}->{'table_id'};
    my $workunit_complete =
	$workunits{$workunit_run_id}->{'workunit_complete'};
    my $table_name = $tablename{$table_id};
    my $nevents = 0;

    if(exists($tablename{$table_id}))
    {

	my $query_nevents = 
	    "SELECT COUNT(*) from VSD_$tablename{$table_id}_Events where " .
	    "WorkunitRunID=$workunit_run_id;";
	
	$stmt = $dbh->prepare($query_nevents)
	    or die "Couldn't prepare statement: " . $dbh->errstr;
	$stmt->execute();
	$nevents = $stmt->fetchrow_array();
    }


    print sprintf("%10s %40s %10s %10s %5s\n",
		  $workunit_run_id,$table_name,$table_id,$nevents,
		  $workunit_complete);

    if($nevents == 0)
    {
	my $query_update = 
	    "UPDATE VSD_WorkunitRun set WorkunitRunComplete=0 " .
	    "where WorkunitRunID=$workunit_run_id;";
	$dbh->do($query_update);

	print "CLEARING WORKUNIT\n";
    }

}
