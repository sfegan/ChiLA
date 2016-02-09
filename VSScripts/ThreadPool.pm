#!/usr/bin/perl -w
#
# Program:     ThreadPool.pm
# Author:      Matthew Wood <mdwood@astro.ucla.edu>
# Date:        03/20/07
#
# Description: Perl module for handling thread pools.
#
# $Id: ThreadPool.pm,v 1.4 2010/05/19 21:43:42 matthew Exp $
#

package ThreadPool;

use strict;
use Exporter;
use Time::HiRes qw (sleep);

our ($VERSION,@ISA,@EXPORT,@EXPORT_OK,%EXPORT_TAGS);


$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ( );
@EXPORT_OK   = ( );

# Variables for threading -----------------------------------------------------

my %thread_status : shared;
my $nthread : shared;

$nthread = 0;

sub clear_threads
{
    @_ = ();

    lock(%thread_status);
    lock($nthread);

    foreach my $tid ( keys(%thread_status) )
    {
	if($thread_status{$tid})
	{
	    my $th = threads->object( $tid );
	    $th->join();
	    delete $thread_status{$tid};
	    $nthread = scalar(keys(%thread_status));
	}
    }
}

sub run_jobs
{
    my $func = $_[0];
    my @jobs = @{$_[1]};
    my $max_thread = $_[2]; 
   
    @_ = ();

    my $njob = @jobs;
    while( scalar(@jobs) || $nthread )
    {
	if( $nthread < $max_thread && scalar(@jobs))
	{
	    my %job = %{pop(@jobs)};
#	    create_thread($func,\%{pop(@jobs)});
	    lock(%thread_status);
	    lock($nthread);
	    my $th = threads->new(\&run_job,$func,\%job);
	    my $tid = $th->tid();
	    $thread_status{ $tid } = 0;
	    $nthread = scalar(keys(%thread_status));
	}
	else
	{
	    Time::HiRes::sleep(0.01);
	}

	clear_threads();
    }
}

sub create_thread
{
    my $func = $_[0];
    my %job = %{$_[1]};

    lock(%thread_status);
    lock($nthread);
    my $th = threads->new(\&run_job,$func,\%job);
    my $tid = $th->tid();
#    print "$tid $nthread\n";
    $thread_status{ $tid } = 0;
    $nthread = scalar(keys(%thread_status));
}

sub run_job
{
    my $func = $_[0];
    my %args = %{$_[1]};

    @_ = ();

    $func->(\%args);

    lock(%thread_status);
    $thread_status{ threads->self->tid() } = 1;
}

1;
