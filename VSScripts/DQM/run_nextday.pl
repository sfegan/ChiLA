#!/usr/bin/perl
#
# Script for running nextday analysis.
#
# \author     Matthew Wood                \n
#             UCLA                        \n
#             mdwood@astro.ucla.edu       \n        
#
# \version $Id: run_nextday.pl,v 1.5 2010/05/19 21:40:56 matthew Exp $
#
#

BEGIN 
{ 
    if(!exists($ENV{'CHILADIR'}))
    {
	print "CHILADIR environment variable is not defined.\n";
	exit(1);
    }
    elsif(!(-d $ENV{'CHILADIR'}))
    {
	print "CHILADIR is not a valid directory: $ENV{'CHILADIR'}\n";
	exit(1);
    }

    push @INC, "$ENV{'CHILADIR'}/VSScripts";
}



use strict;
use Time::localtime;
use POSIX qw(strftime);
use Getopt::Long;
use Pod::Usage;
use threads;
use threads::shared;
use ThreadPool;
use Fcntl qw(:flock);
use Time::localtime;
use File::stat;
use File::Basename;
use File::Path;
use XML::Simple;
use List::Util qw(shuffle);

# -----------------------------------------------------------------------------
# Default/Global Parameters
# -----------------------------------------------------------------------------

my $chila_dir;

if(exists($ENV{'CHILADIR'}))
{
    $chila_dir = "$ENV{'CHILADIR'}";
}
else
{
    $chila_dir = "$ENV{'HOME'}";
}

$ENV{'MATLABPATH'}="$ENV{'HOME'}/matlab:$chila_dir/VSMatlabAnalysis";
$ENV{'LD_LIBRARY_PATH'}="/usr/local/hdf5/lib:/usr/local/veritas/lib";

my $s3_bin = "/home/mdwood/raid/nextday/stage3";
my $octaveio = "/home/mdwood/raid/nextday/octaveio";

my %opt = ();
$opt{'help'} = 0;
$opt{'man'} = 0;
$opt{'date'} = '';
$opt{'nthread'} = 8;
$opt{'overwrite'} = 0;
$opt{'cfg'} = '';
$opt{'nextday_dir'} = "/veritas/chila_nextday";
$opt{'reduced_dir'} = "/net/gamma3/srv/raid5/chila_analysis/analysis";
$opt{'analyze_dates'}=0;
$opt{'lck'}=0;

GetOptions( 'help!'           => \$opt{'help'},
            'man!'            => \$opt{'man'}, 
            'lck!'            => \$opt{'lck'}, 
            'analyze_dates!'  => \$opt{'analyze_dates'}, 
	    'date=s'          => \$opt{'date'},
	    'nthread=i'       => \$opt{'nthread'},
	    'nextday_dir=s'   => \$opt{'nextday_dir'},
	    'reduced_dir=s'   => \$opt{'reduced_dir'}
            ) or pod2usage(2);

# -----------------------------------------------------------------------------
# Print help/man if asked or if no file
# -----------------------------------------------------------------------------
pod2usage(1) if $opt{'help'};
pod2usage(-verbose => 2) if $opt{'man'};
pod2usage("$0: No source file.") if ((@ARGV < 2) && (-t STDIN));

print "-" x 80 . "\n";
print "Starting Nextday Analysis ";
print `date`;

my $lck : shared;
my $s3_reduced_dir = File::Spec->rel2abs( $opt{'reduced_dir'} );
my $s3_results_dir = File::Spec->rel2abs( $opt{'nextday_dir'} );

my $source_file = $ARGV[0];
my @cfg;

for(my $i = 1; $i < scalar(@ARGV); $i++)
{
    push @cfg, $ARGV[$i];
}

if($opt{'lck'})
{
    my $lck_file = "/tmp/nextday.lck";
    open(NEXTDAY,">>$lck_file");
    exit(0) if !flock(NEXTDAY, LOCK_EX | LOCK_NB);
}

# -----------------------------------------------------------------------------
# Load list of runs to be excluded from analysis
# -----------------------------------------------------------------------------
my %excluded_runs;

open(EXCLUDED,"/home/mdwood/raid/nextday/excluded_runs.txt");

while(<EXCLUDED>)
{
    chomp;
    my ($date, $run, $source) = split ' ',$_;
    $excluded_runs{$run} = 1;
}

close(EXCLUDED);

print "Loading diagnostics...\n";
# -----------------------------------------------------------------------------
# Define diagnostics hash as a shared variable so that it doesn't get
# duplicated in memory when new threads are created
# -----------------------------------------------------------------------------
my %diagnostics : shared;
load_diagnostics(\%diagnostics);

my %sources = ();
my @jobs;

print "Loading source list...\n";
open(SOURCES,"$ARGV[0]") or die "Can't open $ARGV[0]: $!";
flock(SOURCES, LOCK_EX);

while(<SOURCES>)
{
    chomp $_;
    my @line = split ' ',$_;
    
    next if /^\#/ || /^$/;

    my ($source_name, $date, $run) = @line;
    my %data;

    $source_name =~ tr/\//_/;

    $data{'source'} = $source_name;
    $data{'date'} = $date;
    $data{'run'} = $run;
    push @{$sources{$source_name}}, \%data; 
}
close(SOURCES);

foreach my $source ( sort(keys(%sources)) )
{
    my @runs;
    my %dates;

    foreach my $run (@{$sources{$source}})
    {
	next if exists($excluded_runs{$run->{'run'}});
	push @runs, "$run->{'date'} $run->{'run'}";
	push @{$dates{$run->{'date'}}}, "$run->{'date'} $run->{'run'}";	
    }

    next if (scalar(@runs)==0);

    
    if($opt{'analyze_dates'})
    {
	foreach my $date (keys(%dates))
	{
	    my @tmp = @{$dates{$date}};
	    
	    foreach my $cfg (@cfg)
	    {
		my $cfg_name = basename $cfg, '.txt';

		push(@jobs, { 'source' => $source,
			      'cfg' => $cfg,
			      'cfg_name' => $cfg_name,
			      'date' => $date,
			      'runlist' => \@tmp } );
	    }
	}
    }

    foreach my $cfg (@cfg)
    {
	my $cfg_name = basename $cfg, '.txt';
	push(@jobs, { 'source' => $source,
		      'cfg' => $cfg,
		      'cfg_name' => $cfg_name,
		      'runlist' => \@runs } );
    }
}

# Randomize jobs --------------------------------------------------------------
@jobs = shuffle @jobs;

if(scalar(@jobs))
{
    print "Starting threads...\n";
    ThreadPool::run_jobs(\&run_analysis,\@jobs,$opt{'nthread'});
}

print "Finished Nextday Analysis ";
print `date`;


sub find_s2_files
{    
    my @files = <$_[0]/*s2.h5>;
    my @s2_files;

    foreach my $file (@files)
    {
	my $source = `$octaveio -l -cat $file observation.name`;	
	chomp $source;
	my $livetime = `$octaveio -l -cat $file diagnostics.gps_livetime_sec`;
	chomp $livetime;

	if($source !~ /drift/ && $source !~ /dark/ && $source ne "unknown" &&
	   $livetime > 400 && $file =~ /([^\/]+)_s2.h5/)
	{
	    my @file_list = ($file);
	    push(@s2_files,{ 'name' => $1, 'list' => \@file_list });
	}
    }

    return @s2_files;
}

sub run_analysis 
{ 
    my $source = $_[0]->{'source'};
    my $cfg_name = $_[0]->{'cfg_name'};
    my $results_dir;
    my $reduced_dir;

    if(exists($_[0]->{'date'}))
    {
	my $date = $_[0]->{'date'};

	$reduced_dir = "$s3_reduced_dir/d$date";
	$results_dir = "$s3_results_dir/d$date";
    }
    else
    {	
	$reduced_dir = "$s3_reduced_dir";
	$results_dir = "$s3_results_dir";
    }

    mkpath($results_dir, 1, 0755) if ! -d $results_dir;
    mkpath($reduced_dir, 1, 0755) if ! -d $reduced_dir;

    my $lck_file = "$reduced_dir/${source}_${cfg_name}.lck";
    open(LCK,">>$lck_file");
    return 0 if !flock(LCK, LOCK_EX | LOCK_NB);

    my $update = run_stage3(@_);
    $update |= run_matlab(@_);
    
    # -------------------------------------------------------------------------
    # Perform a lock here so that multiple threads don't try to update
    # analysis pages at the same time
    # -------------------------------------------------------------------------
    if($update)
    {
	lock($lck);
	generate_analysis_page($results_dir,$reduced_dir);
	generate_source_page($results_dir,$reduced_dir,
			     $source,\@{$_[0]->{'runlist'}});
    }

    print `rm $lck_file`;
    close(LCK);
}

sub generate_source_page
{
    my $results_dir = $_[0];
    my $reduced_dir = $_[1];
    my $source = $_[2];
    my @runlist = @{$_[3]};
    my $htmlfile = "$results_dir/$source/index.html";

    open(HTML,"+> $htmlfile") or die "Can't open $htmlfile: $!";
    flock(HTML, LOCK_EX);
    seek(HTML, 0, 0); truncate(HTML, 0);

    print HTML 
       "<html>\n" .
       "<head>\n<title>$source</title>\n" .
       "<meta http-equiv=\"content-type\" " .
       "content=\"text/html; charset=utf-8\">" .
       "<style type=\"text/css\">\n" .
	"a { text-decoration:none }\n" .
	"td.sigma5 { color: red; font-weight: bolder; }\n" .
	"td.sigma3 { color: orange; font-weight: bolder; }\n" .
	"</style>\n" . 
       "</head>\n" . 	    
       "<body>\n";

    
    my $table_style = 
	"border=\"2\" cellpadding=\"4\" cellspacing=\"0\" style=\"margin-top:1em; margin-bottom:1em; background:\#f9f9f9; border:1px \#aaa solid; border-collapse:collapse;\"";

    print HTML "<table $table_style>\n";

    my $src_ra_hms = "";
    my $src_dec_dms = "";

    open(ANALYSIS,"$results_dir/analysis.dat");
    flock(ANALYSIS, LOCK_EX);

    while(<ANALYSIS>)
    {	
	chomp;
	my ($name,$analysis,$ra_hms,$dec_dms,$livetime,$significance,
	   $excess_rate,$excess_rate_err) = split ' ',$_;

	if($name eq $source)
	{
	    $src_ra_hms = $ra_hms;
	    $src_dec_dms = $dec_dms;
	    last;
	}
    }

    close(ANALYSIS);    

    print HTML "<tr style=\"background:\#dfdfdf;\">\n";
    print HTML "<td><b>Source Name</b></td>\n";
    print HTML "<td><b>RA [hms]</b></td>\n";
    print HTML "<td><b>Dec [dms]</b></td>\n";
    print HTML "</tr>\n";


    print HTML "<tr>\n";
    print HTML "<td>$source</td>\n";
    print HTML "<td>$src_ra_hms</td>\n";
    print HTML "<td>$src_dec_dms</td>\n";
    print HTML "</tr>\n";    


    print HTML "</table>\n";

    print HTML "<table $table_style>\n";

    print HTML "<tr style=\"background:\#dfdfdf;\">\n";
    print HTML "<td><b>Analysis/Cuts</b></td>\n";
    print HTML "<td><b>Livetime [min]</b></td>\n";
    print HTML "<td><b>Significance</b></td>\n";
    print HTML "<td><b>Excess [min<sup>-1</sup>]</b></td>\n";
    print HTML "<td><b>Analysis</b></td>\n";
    print HTML "<td></td>\n";
    print HTML "<td><b>Processed</b></td>\n";
    print HTML "</tr>\n";

    open(ANALYSIS,"$results_dir/analysis.dat");
    flock(ANALYSIS, LOCK_EX);

    while(<ANALYSIS>)
    {	
	chomp;
	my ($name,$analysis,$ra_hms,$dec_dms,$livetime,$significance,
	   $excess_rate,$excess_rate_err) = split ' ',$_;

	my $htmlfile = "${source}_${analysis}_analysis.html";
	my $pdffile = "${source}_${analysis}_analysis.pdf";
	my $pngfile0 = "${source}_${analysis}_analysis-0.png";
	my $pngthumb0 = "${source}_${analysis}_analysis_thumb-0.png";
	my $pngfile1 = "${source}_${analysis}_analysis-1.png";
	my $pngthumb1 = "${source}_${analysis}_analysis_thumb-1.png";
	my $s3_file = "$reduced_dir/${source}_${analysis}_s3.h5";

	next if $name ne $source || !(-e $s3_file);

	my $time = stat($s3_file)->mtime;
	my $time_string = strftime "%Y%m%d %H:%M:%S", @{localtime($time)};

	print HTML "<tr>\n";
	print HTML "<td>$analysis</td>\n";
	print HTML "<td>$livetime</td>\n";
	if($significance < 3)
	{
	    print HTML "<td>$significance</td>\n";
	}
	elsif($significance > 3 && $significance < 5)
	{
	    print HTML "<td class=\"sigma3\">$significance</td>\n";
	}
	else
	{
	    print HTML "<td class=\"sigma5\">$significance</td>\n";
	}

	print HTML "<td>$excess_rate ± $excess_rate_err</td>\n";
	print HTML "<td>";
	print HTML "<a href=\"$htmlfile\">html</a> ";
	print HTML "<a href=\"$pdffile\">pdf</a> ";
	print HTML "</td>\n";
	print HTML "<td><a href=\"$pngfile0\"><img src=\"$pngthumb0\"></a>" .
	    "<a href=\"$pngfile1\"><img src=\"$pngthumb1\"></a>" .
	    "</td>\n";
	
	print HTML "<td>$time_string</td>\n";
	print HTML "</tr>\n";
    }

    close(ANALYSIS);

    print HTML "</table>\n";

    print HTML "<table $table_style>\n";
   
    print HTML "<tr style=\"background:\#dfdfdf;\">\n";

    print HTML "<td><b>Date</b></td>\n";
    print HTML "<td><b>Run</b></td>\n";
    print HTML "<td><b>UTC</b></td>\n";
    print HTML "<td><b>Duration</b></td>\n";
    print HTML "<td><b>Mode</b></td>\n";
    print HTML "<td><b>Scopes</b></td>\n";
    print HTML "<td><b>El</b></td>\n";
    print HTML "<td><b>Az</b></td>\n";
    print HTML "<td><b>L3 Rate</b></td>\n";
    print HTML "<td><b>L3 RChisq</b></td>\n";
    print HTML "<td><b>FIR RMS</b></td>\n";
    print HTML "<td><b>T1 FIR RMS</b></td>\n";
    print HTML "<td><b>T3 FIR RMS</b></td>\n";
    print HTML "<td><b>Moon El</b></td>\n";
    print HTML "<td><b>Diagnostics</b></td>\n";
    print HTML "<td><b>Analysis Selection</b></td>\n";
    print HTML "</tr>\n";

    my %selected_runs;

    foreach my $line(@runlist)
    {
	my ($date, $run) = split ' ',$line;
	$selected_runs{$run} = 1;
    }

    foreach my $run (sort(keys(%diagnostics)))
    {
	next if $diagnostics{$run}->{'source'} ne $source;
	my $date = $diagnostics{$run}->{'date'};

	print HTML "<tr>\n";
	print HTML "<td>$diagnostics{$run}->{'date'}</td>\n";
	print HTML "<td>$run</td>\n";
	print HTML "<td>$diagnostics{$run}->{'time'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'gps_elapsed_time'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'mode'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'config_mask'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'mean_el'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'mean_az'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'l3_rate'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'chisq'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'fir0_rms'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'fir1_rms'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'fir3_rms'}</td>\n";
	print HTML "<td>$diagnostics{$run}->{'moon_el'}</td>\n";

	my $url = "http://gamma1.astro.ucla.edu/diagnostics/d${date}/";

	print HTML 
	    "<td>" . 
	    "<a href=\"${url}x${run}_diagnostics.pdf\">pdf</a> " .
	    "<a href=\"${url}x${run}_diagnostics.html\">html</a> " .
	    "</td>\n";

	if(exists($selected_runs{$run}))
	{
	    print HTML  "<td>1</td>\n";
	}
	else
	{
	    print HTML  "<td>0</td>\n";
	}

	print HTML "</tr>\n";
    }
    
    print HTML "</table>\n";

    print HTML "<h3>Analysis Run List</h3>\n";
    print HTML "<pre>\n";
    foreach my $run(@runlist)
    {
	print HTML "$run\n";
    }
    print HTML "</pre>\n";

    print HTML "</body>\n</html>\n";
    close(HTML); 
}

sub generate_analysis_page
{
    my $results_dir = $_[0];
    my $reduced_dir = $_[1];

    my $datfile = "$results_dir/analysis.dat";

    open(ANALYSIS,"+> $datfile") || die;
    flock(ANALYSIS, LOCK_EX);
    seek(ANALYSIS, 0, 0); truncate(ANALYSIS, 0);
    `cat $results_dir/*/*_analysis.dat > $results_dir/analysis.dat`;

    open(ANALYSIS_XML,">$results_dir/analysis.xml");
    my $xml = new XML::Simple (NoAttr=>1, RootName=>'source');
    my $htmlfile = "$results_dir/index.html";
    
    open(HTML,">$htmlfile");
    print HTML 
       "<html>\n" .
       "<head>\n<title>ChiLA Next Day Analysis Page</title>\n" .
       "<meta http-equiv=\"content-type\" " .
       "content=\"text/html; charset=utf-8\">" .
       "<style type=\"text/css\">\n" .
	"a { text-decoration:none }\n" .
	"td.sigma5 { color: red; font-weight: bolder; }\n" .
	"td.sigma3 { color: orange; font-weight: bolder; }\n" .
	"</style>\n" . 
       "</head>\n" . 	    
       "<body>\n";

    print HTML
       "<table border=\"2\" cellpadding=\"4\" cellspacing=\"0\" style=\"margin-top:1em; margin-bottom:1em; background:\#f9f9f9; border:1px \#aaa solid; border-collapse:collapse;\">\n";
    
    print HTML "<tr style=\"background:\#dfdfdf;\">\n";
    print HTML "<td><b>Source</b></td>\n";
    print HTML "<td><b>Analysis/Cuts</b></td>\n";
    print HTML "<td><b>RA [hms]</b></td>\n";
    print HTML "<td><b>Dec [dms]</b></td>\n";
    print HTML "<td><b>Livetime [min]</b></td>\n";
    print HTML "<td><b>Significance</b></td>\n";
    print HTML "<td><b>Excess [min<sup>-1</sup>]</b></td>\n";
    print HTML "<td><b>Analysis</b></td>\n";
    print HTML "<td><b>Processed</b></td>\n";
    print HTML "</tr>\n";
   
   my $nrun = 0;

    print ANALYSIS_XML "<chila>\n";

   while(<ANALYSIS>)
   {
       chomp;
       my @line = split ' ',$_;

       my ($source,$analysis,$ra_hms,$dec_dms,$livetime,$significance,
	   $excess_rate,$excess_rate_err) = @line;

       my $pngfile = "${source}/${source}_${analysis}_analysis-0.png";
       my $pngthumb = "${source}/${source}_${analysis}_analysis_thumb-0.png";
       my $htmlfile = "${source}/${source}_${analysis}_analysis.html";
       my $pdffile = "${source}/${source}_${analysis}_analysis.pdf";
       my $s3_file = "$reduced_dir/${source}_${analysis}_s3.h5";
       next if ! -e $s3_file;

       my $time = stat($s3_file)->mtime;
       my $time_string = strftime "%Y%m%d %H:%M:%S", @{localtime($time)};

       my %src_xml = (
		      'source_name'      => $source,
		      'analysis'         => $analysis,
		      'ra_hms'           => $ra_hms,
		      'dec_dms'          => $dec_dms,
		      'livetime'         => $livetime,
		      'significance'     => $significance,
		      'excess_rate'      => $excess_rate,
		      'excess_rate_err'  => $excess_rate_err,
		      'processed_date'   => $time_string
		      );
       print ANALYSIS_XML $xml->XMLout(\%src_xml);
       
       if($nrun % 2 == 0)
       {
	   print HTML "<tr style=\"background:\#f9f9f9;\">\n";
       }
       else
       {
	   print HTML "<tr style=\"background:\#f4f4f4;\">\n";
       }
       
       print HTML "<td><a href=\"$source\">$source</a></td>\n";
       print HTML "<td>$analysis</td>\n";
       print HTML "<td>$ra_hms</td>\n";
       print HTML "<td>$dec_dms</td>\n";
       print HTML "<td>$livetime</td>\n";

       if($significance < 3)
       {
	   print HTML "<td>$significance</td>\n";
       }
       elsif($significance > 3 && $significance < 5)
       {
	   print HTML "<td class=\"sigma3\">$significance</td>\n";
       }
       else
       {
	   print HTML "<td class=\"sigma5\">$significance</td>\n";
       }

       print HTML "<td>$excess_rate ± $excess_rate_err</td>\n";
       print HTML "<td>";
       print HTML "<a href=\"$htmlfile\">html</a> ";
       print HTML "<a href=\"$pdffile\">pdf</a> ";
       print HTML "</td>\n";
       print HTML "<td>$time_string</td>\n";

       my $vegas_url = $source;

       $vegas_url =~ tr/_/ /;
       $vegas_url =~ tr/a-z/A-Z/;

       print HTML "<td><a href=\"http://veritas.sao.arizona.edu/private/nextday/detail.pl?$vegas_url\">vegas</a></td>\n";
       print HTML "</tr>\n";
       $nrun = $nrun + 1;
   }
   
    print ANALYSIS_XML "</chila>\n";

    print HTML "</table>\n";
    print HTML "</body>\n</html>\n";
    close(HTML); 
    close(ANALYSIS_XML);
    close(ANALYSIS);
}

sub run_stage3
{
    my $source = $_[0]->{'source'};
    my $cfg_name = $_[0]->{'cfg_name'};
    my $cfg = $_[0]->{'cfg'};

    my %cfg_stg3 = load_cfg($cfg);
    my @runlist = @{$_[0]->{'runlist'}};

    my $dir;

    if(exists($_[0]->{'date'}))
    {
	$dir = "$s3_reduced_dir/d$_[0]->{'date'}";
    }
    else
    {
	$dir = "$s3_reduced_dir";
    }

    my $log_file = "$dir/${source}_${cfg_name}_s3.log";
    my $s3_file = "$dir/${source}_${cfg_name}_s3.h5";
    my $runlist = "$dir/${source}_${cfg_name}_s3.txt";
    my @s2_files;

    foreach my $daterun (@runlist)
    {
	my ($date, $run) = split ' ',$daterun;

	my $s2_file = 
	    "/net/gamma3/srv/raid5/chila_analysis/d${date}/x${run}_s2.h5";

	if(! -e $s2_file)
	{
	    print "Error: Could not find file: $s2_file\n";
	    return 0;
	}

	push @s2_files, $s2_file;
    }    

    my $run_stage3 = 0;

    if(-e $s3_file)
    {
	my $files = `$octaveio -v $s3_file files`;
	my @files;
	my @tmp = split '\n',$files;
	foreach my $tmp(@tmp)
	{
	    my @line = split ' ', $tmp;
	    if(scalar(@line)==2)
	    {
		push @files, $line[1];
	    }
	}	

	# ---------------------------------------------------------------------
	# Check if the file list differs or if any one stage2 file is
	# newer than the stage3 file
	# ---------------------------------------------------------------------
	if(scalar(@files) != scalar(@s2_files))
	{
	    $run_stage3 = 1;
	}
	else
	{
	    my $s3_time = stat($s3_file)->mtime;

	    foreach my $s2_file (@s2_files)
	    {
		if( ! grep { $_ eq $s2_file } @files )
		{
		    $run_stage3 = 1;
		    last;
		}
		
		my $s2_time =  stat($s2_file)->mtime;
		if($s2_time > $s3_time)
		{
		    $run_stage3 = 1;
		    last;
		}
	    }
	}
    }
    else
    {
	$run_stage3 = 1;
    }

    return 0 if $run_stage3 == 0;

    open(RUNLIST,">$runlist");
    foreach my $s2_file (@s2_files)
    {
	print RUNLIST "$s2_file\n";
    }

    my $cfg = "";      
    foreach my $setting (keys(%cfg_stg3))
    {
 	$cfg = $cfg . "-$setting=$cfg_stg3{$setting} ";
    }

    print "RUNNING STAGE3 == SRC: $source CFG: $cfg_name LOG: $log_file\n";
#    print "$s3_bin $runlist $cfg -o $s3_file &> $log_file\n";
    `$s3_bin $runlist $cfg -o $s3_file &> $log_file`;

    close(RUNLIST);
    return 1;
}

sub run_matlab
{
    my $source = $_[0]->{'source'};
    my $cfg_name = $_[0]->{'cfg_name'};
    my $output_dir = "$s3_results_dir/$source";
    my $reduced_dir;

    if(exists($_[0]->{'date'}))
    {
	$output_dir = "$s3_results_dir/d$_[0]->{'date'}/$source";
	$reduced_dir = "$s3_reduced_dir/d$_[0]->{'date'}";
    }
    else
    {
	$output_dir = "$s3_results_dir/$source";
	$reduced_dir = "$s3_reduced_dir";
    }

    if(! -d "$output_dir")
    {
	mkdir "$output_dir", 0755 or die "Can't mkdir $output_dir: !"; 
    }

    my $convert_bin = "/usr/local/ImageMagick/bin/convert";


    my $matlab_options = "-nodisplay -nodesktop ";
    my $convert_options = "-density 300 -quality 90";

    my $file = "$reduced_dir/${source}_${cfg_name}_s3.h5";
    my $base = "${source}_${cfg_name}";
    my $log="$output_dir/${base}_analysis.log";
    my $cmd_log="$output_dir/${base}_analysis_cmd.log";
    my $analysis="$output_dir/${base}_analysis.dat";
    my $htmlfile="$output_dir/${base}_analysis.html";
    my $psfile="$output_dir/${base}_analysis.ps";
    my $pdffile="$output_dir/${base}_analysis.pdf";
    my $pngfile="$output_dir/${base}_analysis";
    my $pngthumb="$output_dir/${base}_analysis_thumb";

    return 0 if ! -e $file;

    my $time1 = stat($file)->mtime;
    my $time2 = $time1;
    if(-f $pdffile) { $time2 = stat($pdffile)->mtime; }

    return 0 if -f $pdffile && !$opt{'overwrite'} && $time1 < $time2;
    
    my $matlab_command = 
	"\"try; stage3_analysis(\'$file\',\'$psfile\',\'$analysis\'," .
	"\'$source\',\'$cfg_name\'); " .
	"catch; " .
	"msg=lasterror; disp(msg.message); disp(msg.identifier); " .
	"end; exit;\"";

    my $cmd = 
	"DISPLAY=\"\" /veritas/matlab/bin/matlab -logfile $log " .
	"$matlab_options -r $matlab_command";

#    print "/veritas/matlab/bin/matlab $matlab_options -r $matlab_command\n";

    print "RUNNING MATLAB == FILE: $file ANALYSIS: $analysis SRC: $source " .
	" CFG: $cfg_name\n";
    `$cmd`;
    `sed '1i\\$cmd' $log > $cmd_log`;
    `mv $cmd_log $log`;

    `ps2pdf $psfile $pdffile`;

#    if(! -d "$s3_results_dir/png") { mkdir '$s3_results_dir/png', 0755; }
#    print "$convert_bin $convert_options -resize 1200 $pdffile $pngfile.png\n";
    `$convert_bin $convert_options -resize 1200 $pdffile $pngfile.png`;
#    print "$convert_bin $convert_options -resize 350 $pdffile $pngthumb.png\n";
    `$convert_bin $convert_options -resize 350  $pdffile $pngthumb.png`; 
    `rm -f $psfile`;
        
    `cp $reduced_dir/${source}_${cfg_name}_s3.txt $output_dir`;

    open(HTML,">$htmlfile");
    print HTML 
	"<html>\n" .
	"<head>\n<title>${base}_analysis</title>\n" .
	"<meta http-equiv=\"content-type\" " .
	"content=\"text/html; charset=utf-8\">\n</head>\n" . 	    
	"<body>\n";
    for(my $i = 0; $i < 2; $i++)
    {
	my $tmp = fileparse($pngfile);

	print HTML 
	    "<a href=\"$tmp-$i.png\">" .
	    "<img src=\"${tmp}_thumb-$i.png\"></a>\n";
    }
    print HTML "</body>\n</html>\n";
    close(HTML);

    return 1;
}

sub load_cfg
{
    my $cfg = $_[0];
    my %cfg;
    
    open(CFG_FILE,$cfg) or die "Can't open $cfg: $!";
    while(<CFG_FILE>)
    {
	next if /^\#/ || /^$/;
	
	chomp;
	my @line = split ' ',$_;
	if(scalar(@line) != 2)
	{
	    print "Error parsing the following line: $_\n";
	    exit(1);
	}
	
	$cfg{$line[0]} = $line[1];	
    }
    
    close(CFG_FILE);    

    return %cfg;
}

sub load_diagnostics
{
    my $diagnostics = shift;

    my $diagnostics_file = "/veritas/chila_diagnostics/diagnostics.dat";
    open(DIAGNOSTICS,$diagnostics_file) or 
	die "Can't open $diagnostics_file: $!";
    while(<DIAGNOSTICS>)
    {
	my %data : shared;
	%data = parse_diagnostics($_);
	$diagnostics->{$data{'run'}} = \%data;
    }
    close(DIAGNOSTICS);
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
			'gps_elapsed_time' => 0);

    
   

    if(scalar(@line) < 41)
    {
	return %diagnostics;
    }	

    if($line[6] =~ /\@(.+)/)
    {
	$diagnostics{'wobble_direction'} = $1;
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
    }
    else
    {
	$diagnostics{'fir0_rms'} = 0;
	$diagnostics{'fir1_rms'} = 0;
	$diagnostics{'fir3_rms'} = 0;
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

    for(my $i = 3; $i >= 0; $i--)
    {
	if($cfd_threshold[$i] > 0)
	{
	    $config_mask = $config_mask . "1";
	}
	else
	{
	    $config_mask = $config_mask . "0";
	}

    }
    
    

    if($version >= 2)
    {
	my $scope_mask = $line[102];
	my $scopes = "";
	for(my $i = 0; $i < 4; $i++)
	{
	    my $scope_id = $i+1;
	    my $tmp = (1 << $i);
	    
	    if($tmp & $scope_mask)
	    {
		$scopes = $scopes . "$scope_id";
	    }
	    else
	    {
		$scopes = $scopes . "-";
	    }
	}

	$diagnostics{'config_mask'} = $scopes;
    }
    else
    {
	my $config_mask;
	for(my $i = 0; $i < 4; $i++)
	{
	    if($cfd_threshold[$i] > 0)
	    {
		my $scope_id = $i+1;
		$config_mask = $config_mask . $scope_id;
	    }
	    else
	    {
		$config_mask = $config_mask . "-";
	    }
	}
	$diagnostics{'config_mask'} = $config_mask;
    }


    

    return %diagnostics;
}



__END__

=head1 NAME

run_nextday.pl - Script for running next-day analysis.

=head1 SYNOPSIS

run_nextday.pl [options]

Options:

  -help                    Print this help message.
  -man                     Print the man page.

