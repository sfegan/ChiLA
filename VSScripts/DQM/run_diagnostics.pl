#!/usr/bin/perl -w
#
# \author     Matthew Wood                \n
#             UCLA                        \n
#             mdwood@astro.ucla.edu       \n        
#
# \version    0.1
# \date       06/24/2008

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::stat;
use Time::localtime;
use POSIX qw(strftime);
use DBI;
use Fcntl qw(:flock);

my $opt_help = 0;
my $opt_man = 0;
my $opt_utc_date = 0;
my $opt_overwrite = 0;
my $opt_overwrite_matlab = 0;
my $opt_nday = 4;
my $opt_cfg_file = "/veritas/chila_diagnostics/diagnostics_options.txt";
my $opt_vbf_dir_base = "/veritas/data";
my $opt_local_dir = "/net/gamma3/srv/raid4/chila_analysis";
my $opt_diagnostics_dir = "/veritas/chila_diagnostics";
my $opt_code_dir = "/veritas/chila_diagnostics/bin";
my $opt_laser_table = "/veritas/chila_diagnostics/laser.dat";
my $opt_verbose = 0;

GetOptions( 'h|help!'           => \$opt_help,
	    'man!'              => \$opt_man,
	    'overwrite!'        => \$opt_overwrite,
	    'overwrite_matlab!' => \$opt_overwrite_matlab,
	    'v!'                => \$opt_verbose,
	    'utc=s'             => \$opt_utc_date,
	    'nday=i'            => \$opt_nday,
	    'cfg=s'             => \$opt_cfg_file,
	    'vbf_dir=s'         => \$opt_vbf_dir_base,
	    'local_dir=s'       => \$opt_local_dir,
	    'diagnostics_dir=s' => \$opt_diagnostics_dir,
	    'code_dir=s'        => \$opt_code_dir,
	    'laser_table=s'     => \$opt_laser_table
            ) or pod2usage(2);

# -----------------------------------------------------------------------------
# Print help/man if asked or if no file
# -----------------------------------------------------------------------------
pod2usage(1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;

my $vbf_dir_base = $opt_vbf_dir_base;
my $red_dir_base = "/veritas/chila_reduced";
my $dia_dir_base = $opt_diagnostics_dir;
#my $code_dir = "$ENV{'HOME'}/raid/diagnostics/bin";
#my $local_dir = "$ENV{'HOME'}/raid/reduced";
my $code_dir = $opt_code_dir;
my $local_dir = $opt_local_dir;
my @utc_dates;
my %opt_s2;
my %opt_lsr;

$ENV{'LD_LIBRARY_PATH'}="/usr/local/hdf5/lib:/usr/local/veritas/lib";

if(!$opt_local_dir)
{
    print "No destination directory for reduced files specified.\n";
    exit(1);
}
elsif(!$opt_code_dir)
{
    print "Directory for chila executables undefined.  " .
	"Set this with the --code_dir option.\n";
    exit(1);
}

print "Starting Diagnostics ";
print `date`;

# -----------------------------------------------------------------------------
# Load stage2 options
# -----------------------------------------------------------------------------

if($opt_cfg_file)
{
    open(OPT_FILE,$opt_cfg_file) or die 
	"Can't open configuration file $opt_cfg_file: $!";
    while(<OPT_FILE>)
    {
	next if /^\#/ || /^$/;

	chomp;
	my @line = split ' ',$_;
	if(scalar(@line) != 2)
	{
	    print "Error parsing the following line: $_\n";
	    exit(1);
	}
	
	if($line[0] =~ /^s2_/ || $line[0] =~ /^s1_/ || $line[0] =~ 'hilo')
	{
	    $opt_s2{$line[0]} = $line[1];	
	}
	elsif($line[0] =~ /^VSDB/)
	{
	    $opt_s2{$line[0]} = $line[1];	
	    $opt_lsr{$line[0]} = $line[1];
	}
	else
	{
	    $opt_lsr{$line[0]} = $line[1];
	}
    }
}

# -----------------------------------------------------------------------------

my %run_date;

my @vbf_data_dir = ("$vbf_dir_base");

foreach my $dir (@vbf_data_dir)
{
    my @vbf_files = <$dir/d*/*cvbf>;
    
    foreach my $vbf_file (@vbf_files) 
    {    
	if($vbf_file =~ m/d(........)\/(\d+).cvbf/)
	{
	    $run_date{$2} = $1;
	}
    }    
}

# -----------------------------------------------------------------------------
# Process the list of dates specified with the -utc option
# -----------------------------------------------------------------------------
if($opt_utc_date)
{    
    my @dates = split ',',$opt_utc_date;

    foreach my $date (@dates)
    {
	my @vbf_dir = <${vbf_dir_base}/d$date>;
	foreach my $dir(@vbf_dir)
	{
	    push @utc_dates,$1 if($dir =~ /d(........)$/);
	}
    }
}
# -----------------------------------------------------------------------------
# Process all dates which are within -nday of the current one
# -----------------------------------------------------------------------------
else
{
    for(my $iday = 0; $iday < $opt_nday; $iday++)
    {
	my @utc_time = gmtime(time()-86400*$iday);
	my $date = sprintf("%04s%02s%02s",
			   $utc_time[5]+1900,$utc_time[4]+1,$utc_time[3]);

#	print "$iday $date\n";
	push @utc_dates, $date;
    }
}

foreach my $utc_date (@utc_dates)
{
    my $red_dir = "$red_dir_base/d$utc_date";
    my $dia_dir = "$dia_dir_base/d$utc_date";
    my $vbf_dir = "$vbf_dir_base/d$utc_date";

    my @vbf_files = <$vbf_dir/*>;
    my $num_vbf_files = scalar(@vbf_files);
    print "-" x 80 . "\n";
    print "PROCESSING UTC DATE:    $utc_date NFILES $num_vbf_files\n";

    # Skip this date if there are no VBF files --------------------------------
    next if $num_vbf_files == 0;

    if(! -e $red_dir)
    {
	if(! -d "${local_dir}/d$utc_date")
	{
	    mkdir "${local_dir}/d$utc_date", 0775 or 
		die "Can't mkdir $red_dir: !";
	}

	# Create a softlink if the local directory and reduced data
	# directory differ
	if(! -e $red_dir)
	{
	    `ln -s ${local_dir}/d$utc_date $red_dir`;
	}
    }

    print "VBF DIR:     $vbf_dir\n";
    print "REDUCED DIR: $red_dir\n";

    open(LCKFILE,">>$red_dir/diagnostics.lck") or 
	die "Can't open $red_dir/diagnostics.lck: $!";
    if(!flock(LCKFILE, LOCK_EX | LOCK_NB))
    {
	close(LCKFILE);
	next;
    }	

    # -------------------------------------------------------------------------
    # Run laser and stage1/stage2 analysis
    # -------------------------------------------------------------------------
    analyze_directory($utc_date);

    my @reduced_files = <$red_dir/*.h5>;
    my $nfiles = scalar(@reduced_files);

    if($nfiles > 0)
    {
	if(!-d $dia_dir)
	{
	    mkdir $dia_dir, 0775 or die "Can't mkdir $dia_dir: $!";
	}

	# ---------------------------------------------------------------------
	# Run matlab analysis
	# ---------------------------------------------------------------------
	run_matlab($red_dir,$dia_dir);
	if($opt_verbose)
	{
	    print "cat $red_dir/*diagnostics.dat > $dia_dir/diagnostics.dat\n";
	}
	`cat $red_dir/*diagnostics.dat > $dia_dir/diagnostics.dat`;
	chdir($dia_dir);

	# ---------------------------------------------------------------------
	# Generate diagnostics webpage
	# ---------------------------------------------------------------------
	generate_diagnostics_page($utc_date);
    }
    `rm $red_dir/diagnostics.lck`;
    close(LCKFILE);
}

my $dat_tmp = "$dia_dir_base/diagnostics_tmp.dat";

open(DATFILE,">$dat_tmp") or die "Can't open $dat_tmp: $!";
flock(DATFILE, LOCK_EX);

if($opt_verbose)
{
    print "cat $dia_dir_base/d????????/diagnostics.dat > $dat_tmp\n";
}
`cat $dia_dir_base/d????????/diagnostics.dat > $dat_tmp`;

generate_main_page();


`mv $dat_tmp $dia_dir_base/diagnostics.dat`;

flock(DATFILE, LOCK_UN);
close(DATFILE);

print "Finished Diagnostics ";
print `date`;

sub analyze_directory
{
    my $utc_date = $_[0];
	
    my $red_dir = "$red_dir_base/d$utc_date";
    my $dia_dir = "$dia_dir_base/d$utc_date";
    my $vbf_dir = "$vbf_dir_base/d$utc_date";

    chdir $red_dir or die "$!\n";

    my @vbf_files = <$vbf_dir/*vbf>;
    foreach my $vbf_file (@vbf_files)
    {
	my $run = fileparse($vbf_file,".cvbf");
	my $run_type = get_runtype($run);
	my %opt = %opt_s2;

	next if $run_type ne "laser" && $run_type ne "observing" && 
	    $run_type ne "flasher";

	my @laser;
	if($run_type eq "laser" || $run_type eq "flasher")
 	{
 	    run_laser($run);
 	    $opt{'laser'}="x${run}_s1.h5"; 
 	    $opt{'s2_no_laser'}="false"; 
 	    $laser[0]=$run;
 	    $laser[1]=$run;
 	    $laser[2]=$run;
 	    $laser[3]=$run;
 	}	
	elsif(find_laser($utc_date,$vbf_file,\@laser))
	{
	    my @lsr_file;
	    for(my $iscope = 0; $iscope < scalar(@laser); $iscope++)
	    {
		if($laser[$iscope] > 0)
		{
		    run_laser($laser[$iscope]);
		    push @lsr_file, 
		    "$iscope/$red_dir_base/d$run_date{$laser[$iscope]}/" .
			"x$laser[$iscope]_s1.h5"; 
		}
	    }	    
	    
	    $opt{'laser'}=join ',',@lsr_file; 
	    $opt{'s2_no_laser'}="false";
	}   
	else
	{
	    delete($opt{'laser'});
	    $opt{'s2_no_laser'}="true";
	}
    	   
	my $s1_file = "x${run}_s1.h5";
	my $s2_file = "x${run}_s2.h5";
	my $s2_tmp_file = "x${run}_tmp_s2.h5";
	my @s2_lsr_run = get_s2_laser($s2_file);
	my @s2_nchan = get_s2_mask($s2_file);

	my $lsr_match = 1;
	for(my $iscope = 0; $iscope < scalar(@laser); $iscope++)
	{
	    next if !defined $laser[$iscope];
	    next if !defined $s2_nchan[$iscope];

	    if(! defined $s2_lsr_run[$iscope])
	    {
		$lsr_match = 0;
		last;
	    }
	    elsif($laser[$iscope] != $s2_lsr_run[$iscope] && 
		  $s2_nchan[$iscope] > 0)
	    {
		$lsr_match = 0;
		last;
	    }
	}

	next if(-f $s2_file && !$opt_overwrite && $lsr_match);


	if(! -f $s2_file)
	{
	    print "STAGE2 FILE $s2_file does not exist.\n";
	}
	elsif($lsr_match == 0)
	{
	    print "LASER CALIBRATION for $s2_file has changed.\n";
	    print "@laser\n";
	    print "@s2_lsr_run\n";
	}

	$opt{'stage1'}="$red_dir/$s1_file";
	$opt{'o'}="$red_dir/$s2_tmp_file";

	my $opt_s2 = "";
	foreach my $opt (sort keys(%opt))
	{
	    $opt_s2 = $opt_s2 . "-$opt=$opt{$opt} ";
	}

	my $opt_log = "2>&1 | tee $red_dir/x${run}.log";
	
	print "RUNNING STAGE1/STAGE2 == FILE: ${vbf_file} " .
	    "LOG: $red_dir/x${run}.log\n";
	if($opt_verbose)
	{
	    print "$code_dir/stage2 $opt_s2 $vbf_file $opt_log\n";
	}
	`$code_dir/stage2 $opt_s2 $vbf_file $opt_log`;

	# ---------------------------------------------------------------------
	# Obtain a lock on the stage2 file before overwriting it
	# ---------------------------------------------------------------------
	open(S2FILE,">>$s2_file") or die "Can't open $s2_file: $!";
	flock(S2FILE, LOCK_EX);	
	if($opt_verbose)
	{
	    print "mv $s2_tmp_file $s2_file\n";
	}
	`mv $s2_tmp_file $s2_file`;	
	flock(S2FILE, LOCK_UN);
	close(S2FILE);
    }
}

sub run_laser
{
    my ($run) = @_;
    my $date = $run_date{$run};

    my $vbf_file = "$vbf_dir_base/d${date}/$run.cvbf";
    my $lsr_file = "$red_dir_base/d${date}/x${run}_s1.h5"; 
    my $lsr_log = "$red_dir_base/d${date}/x${run}_laser.log"; 

    if(! -e "$red_dir_base/d${date}")
    {
	mkdir "${local_dir}/d${date}", 0775 or 
	    die "Can't mkdir ${local_dir}/d${date}: !";

	`ln -s ${local_dir}/d${date} $red_dir_base/d${date}`;
    }
	
    if(-f $lsr_file && !$opt_overwrite)
    {
#	print "Skipping $lsr_file because it already exists.\n";
    }
    else
    {
	my $opt_lsr = "-o $lsr_file ";
	foreach my $opt (sort keys(%opt_lsr))
	{
	    $opt_lsr = $opt_lsr . "-$opt=$opt_lsr{$opt} ";
	}
	
	my $opt_log = "2>&1 | tee $lsr_log";

	print "RUNNING LASER == FILE: ${vbf_file} LOG: $lsr_log\n";

	if($opt_verbose)
	{
	    print "$code_dir/laser $opt_lsr $vbf_file $opt_log\n";
	}
	`$code_dir/laser $opt_lsr $vbf_file $opt_log`;
    }
}

sub run_matlab
{
    my ($red_dir,$dia_dir) = @_;

    chdir $dia_dir or die "$!\n";

    my $matlab_options = "-nodisplay -nodesktop ";
    my $convert_bin = "/usr/local/ImageMagick/bin/convert";
    my $convert_options = " -density 300 -quality 90 ";

    my @files = <$red_dir/x*_s2.h5>;

    foreach my $file (@files)
    {
	my $base = fileparse($file,"_s2.h5");
	my $diagnostics="${base}_diagnostics.dat";
	my $diagnostics_log="${base}_diagnostics.log";
	my $diagnostics_cmd_log="${base}_diagnostics_cmd.log";
	my $htmlfile="${base}_diagnostics.html";
	my $psfile="${base}_diagnostics.ps";
	my $pdffile="${base}_diagnostics.pdf";
	my $pngfile="png/${base}_diagnostics";
	my $pngthumb="png/${base}_diagnostics_thumb";

	my $time1 = stat($file)->mtime;
	my $time2 = $time1;
	if(-f $pdffile) { $time2 = stat($pdffile)->mtime; }

	next if -f $pdffile && !$opt_overwrite_matlab && $time1 < $time2;

	my $matlab_command = 
	    "\"try; one_diagnostics(\'$file\',\'$diagnostics\'); catch; " .
	    "msg=lasterror; disp(msg.message); disp(msg.identifier); " .
	    "end; exit;\"";
	my $cmd = "DISPLAY=\"\" /veritas/matlab/bin/matlab " .
	    "$matlab_options -logfile $diagnostics_log -r $matlab_command";

	# open(LOG,">$diagnostics_cmd_log");
	# print LOG "DISPLAY=\"\" /veritas/matlab/bin/matlab $matlab_options " .
	#     "-logfile $diagnostics_log -r $matlab_command\n";
	# close(LOG);
	print "RUNNING MATLAB == FILE: $file LOG: $diagnostics_log\n";
	`$cmd`;
	`sed '1i\\$cmd' $diagnostics_log > $diagnostics_cmd_log`;
	`mv $diagnostics_cmd_log $diagnostics_log`;

	if(! -e $diagnostics)
	{
	    print "ERROR: Diagnostics failed.\n";
	    next;
	}

	`ps2pdf $psfile`;
	if(! -d "png") { mkdir 'png', 0775; }
	print "$convert_bin $convert_options -resize 1200 $pdffile $pngfile.png\n";
	print "$convert_bin $convert_options -resize 350  $pdffile $pngthumb.png\n";
	`$convert_bin $convert_options -resize 1200 $pdffile $pngfile.png`;
	`$convert_bin $convert_options -resize 350  $pdffile $pngthumb.png`; 

	print `pwd`;
	`rm -f $psfile`;
	`mv $diagnostics $red_dir`;

	open(HTML,">$htmlfile");
	print HTML 
	    "<html>\n" .
	    "<head>\n<title>${base}_diagnostics</title>\n" .
	    "<meta http-equiv=\"content-type\" " .
	    "content=\"text/html; charset=utf-8\">\n</head>\n" . 	    
	    "<body>\n";
	for(my $i = 0; $i < 10; $i++)
	{
	    if(-e "$pngfile-$i.png")
	    {
		print HTML 
		    "<a href=\"$pngfile-$i.png\">" .
		    "<img src=\"$pngthumb-$i.png\"></a>\n";
	    }
	    elsif(-e "$pngfile.png.$i")
	    {
		print HTML 
		    "<a href=\"$pngfile.png.$i\">" .
		    "<img src=\"$pngthumb.png.$i\"></a>\n";
	    }
	}
	print HTML "</body>\n</html>\n";
	close(HTML);
    }
}

sub generate_main_page
{
    my $htmlfile = "$dia_dir_base/index.html";
    my @dirs = <$dia_dir_base/d*>;
    
    @dirs = sort { $b cmp $a } @dirs;

    open(HTML,">$htmlfile");
    print HTML 
	"<html>\n" .
	"<head>\n<title>ChiLA Diagnostics Page</title>\n" .
	"<meta http-equiv=\"content-type\" " .
	"content=\"text/html; charset=utf-8\">\n" .
	"<style type=\"text/css\">\n" .
	"a { text-decoration:none }\n" .
	"</style>\n" . 
	"</head>\n" . 	    
	"<body>\n";
    
    print HTML "<p><a href=diagnostics.dat>diagnostics.dat</a></p>\n"; 
    print HTML "<p><a href=\"http://gamma1.astro.ucla.edu/nextday/\">" .
	"Cumulative Next Day Analysis</a></p>\n"; 

    print HTML
       "<table border=\"2\" cellpadding=\"4\" cellspacing=\"0\" style=\"margin-top:1em; margin-bottom:1em; background:\#f9f9f9; border:1px \#aaa solid; border-collapse:collapse;\">\n";

#    print HTML "<tr style=\"background:\#dfdfdf;\">\n";
#    print HTML "<td><b>Date</b></td>\n";
#    print HTML "<td><b>Diagnostics</b></td>\n";
#    print HTML "<td><b>Next Day Analysis</b></td>\n";
#    print HTML "</tr>\n";

    my $ndate = 0;
    foreach my $dir (@dirs)
    {
	next if $dir !~ /d(........)$/;
	my $date = $1;


	if($ndate % 2 == 0)
	{
	    print HTML "<tr style=\"background:\#fafafa;\">\n";
	}
	else
	{
	   print HTML "<tr style=\"background:\#f4f4f4;\">\n";
       }

	print HTML "<td><b>$date</b></td>\n"; 
	print HTML "<td><a href=\"d$date\">diagnostics</a></td>\n"; 

	if(-e "/veritas/chila_nextday/d$date/index.html")
	{
	    print HTML 
		"<td><a href=\"http://gamma1.astro.ucla.edu/nextday/d$1\">" .
		"nextday</a></td>\n"; 
	}
	else
	{
	    print HTML "<td></td>\n"; 
	}
       
	print HTML "</tr>\n";
	$ndate = $ndate + 1;
    }
   
    print HTML "</table>\n";    
    print HTML "</body>\n</html>\n";
    close(HTML); 
}

sub generate_diagnostics_page
{
   my $date = $_[0];
   my $htmlfile = "index.html";

   open(DIAGNOSTICS,"diagnostics.dat");
   open(HTML,">$htmlfile");
   print HTML 
       "<html>\n" .
       "<head>\n<title>d${date}_diagnostics</title>\n" .
       "<meta http-equiv=\"content-type\" " .
       "content=\"text/html; charset=utf-8\">\n</head>\n" . 	    
       "<body>\n";
   
   print HTML
       "<table border=\"2\" cellpadding=\"4\" cellspacing=\"0\" style=\"margin-top:1em; margin-bottom:1em; background:\#f9f9f9; border:1px \#aaa solid; border-collapse:collapse;\">\n";
   
   print HTML "<tr style=\"background:\#dfdfdf;\">\n";
   print HTML "<td><b>Run</b></td>\n";
   print HTML "<td><b>UTC</b></td>\n";
   print HTML "<td><b>Duration</b></td>\n";
   print HTML "<td><b>Source</b></td>\n";
   print HTML "<td><b>Mode</b></td>\n";
   print HTML "<td><b>Scopes</b></td>\n";
   print HTML "<td><b>Laser</b></td>\n";
   print HTML "<td><b>El</b></td>\n";
   print HTML "<td><b>Az</b></td>\n";
   print HTML "<td><b>L3 Rate</b></td>\n";
   print HTML "<td><b>L3 RChisq</b></td>\n";
   print HTML "<td><b>FIR RMS</b></td>\n";
   print HTML "<td><b>T1 FIR RMS</b></td>\n";
   print HTML "<td><b>T3 FIR RMS</b></td>\n";
   print HTML "<td><b>Moon El</b></td>\n";
   print HTML "<td><b>Diagnostics</b></td>\n";
   print HTML "<td><b>Processed</b></td>\n";
   print HTML "</tr>\n";
   
   my $nrun = 0;

   while(<DIAGNOSTICS>)
   {
       chomp;
       my @line = split ' ',$_;
       
       my $run = $line[0];
       my $version = $line[1];
       my $utc_time = substr($line[3],0,8);
       my $laser = $line[42];
       my $source = $line[5];
       my $mode = $line[6];
       my $elevation = $line[7];
       my $azimuth = $line[8];
       my $gps_time = $line[40];
       my $l3_rate = $line[12];
       my $l3_chisq = $line[14];
       my $moon_elevation = $line[43];
       my $fir0_rms = 0;
       my $fir1_rms = 0;
       my $fir3_rms = 0;
       my $scopes = "1234";
       my $s2_file = "$red_dir_base/d${date}/x${run}_s2.h5";

       my $time = stat($s2_file)->mtime;
       my $time_string = strftime "%Y%m%d %H:%M:%S", @{localtime($time)};

       if($version == 2)
       {
	   $fir0_rms = $line[106];
	   $fir1_rms = $line[111];
	   $fir3_rms = $line[121];
	   my $scope_mask = $line[102];
	   $scopes = "";
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
       }
       else
       {

       }


       if($nrun % 2 == 0)
       {
	   print HTML "<tr style=\"background:\#f9f9f9;\">\n";
       }
       else
       {
	   print HTML "<tr style=\"background:\#f4f4f4;\">\n";
       }
       
       my $run_type = get_runtype($run);

       if($run_type eq "laser")
       {
	   $source = "laser";
       }

       if($run_type eq "flasher")
       {
	   $source = "flasher";
       }

       print HTML "<td>$run</td>\n";
       print HTML "<td>$utc_time</td>\n";
       print HTML "<td>$gps_time</td>\n";
       print HTML "<td>$source</td>\n";
       print HTML "<td>$mode</td>\n";
       print HTML "<td>$scopes</td>\n";
       print HTML "<td>$laser</td>\n";
       print HTML "<td>$elevation</td>\n";
       print HTML "<td>$azimuth</td>\n";
       print HTML "<td>$l3_rate</td>\n";
       print HTML "<td>$l3_chisq</td>\n";
       print HTML "<td>$fir0_rms</td>\n";
       print HTML "<td>$fir1_rms</td>\n";
       print HTML "<td>$fir3_rms</td>\n";
       print HTML "<td>$moon_elevation</td>\n";
       print HTML 
	   "<td>" . 
	   "<a href=x${run}_diagnostics.pdf>pdf</a> " .
	   "<a href=x${run}_diagnostics.html>html</a> " .
	   "<a href=diagnostics.dat>txt</a> " .
	   "</td>\n";
       print HTML "<td>$time_string</td>\n";
       print HTML "</tr>\n";
       $nrun = $nrun + 1;
   }
   
   print HTML "</table>\n";
   
   open(LOGSHEETS,"$dia_dir_base/logsheets.txt");
   while(<LOGSHEETS>)
   {
       chomp;
       my @line = split ' ',$_;

       if($line[0] eq $date)
       {
	   

	   my $log = 
	       `lynx -width=1000 -dump -auth=veritas:Al3thia http://veritas.sao.arizona.edu/private/elog/$line[1] | sed -n '/Message ID/,\$p'`;

	   

	   print HTML "<h3>Observing Log</h3>\n";
	   print HTML "<pre>\n";
	   print HTML "$log";
	   print HTML "</pre>\n";
	   last;
       }
   }

   close(LOGSHEETS);

   print HTML "</body>\n</html>\n";
   close(HTML); 
}

sub get_s2_mask
{
    my $s2_file = $_[0];

    my $mask = `$code_dir/octaveio -l -cat $s2_file stage1.run_info.nchan`;
    my @mask = split '\n', $mask;
    return @mask;
}

sub get_s2_laser
{
    my $s2_file = $_[0];
    my @laser;

    if(! -e $s2_file)
    {
	return @laser;
    }

    my $laser = `$code_dir/octaveio -l $s2_file stage1`;

    if($laser =~ /stage1\.laser/)
    {
	$laser = `$code_dir/octaveio -l -cat $s2_file stage1.laser.m_runno`;
	chomp $laser;

	$laser[0] = $laser;
	$laser[1] = $laser;
	$laser[2] = $laser;
	$laser[3] = $laser;	
    }

    for(my $iscope = 0; $iscope < 4; $iscope++)
    {
	my $path = "stage1.laser.scope{$iscope}.runno";
	$laser = `$code_dir/octaveio -l -cat $s2_file $path`;
	chomp $laser;

	if($laser =~ /^[\d]+$/ && $laser ne "0")
	{
	   $laser[$iscope] = $laser; 
	}
    }

    return @laser;
}

sub find_laser
{
    my $utc = shift;
    my $data_run = shift;
    my $laser = shift;

    $data_run =~ s/(.+)\/([^\/]+)$/$2/;
    $data_run =~ s/(\d+)\.cvbf/$1/;

    # -------------------------------------------------------------------------
    open(LASER,"$dia_dir_base/laser.dat");
    while(<LASER>)
    {
	chomp;
	next if /^\#/ || /^$/;

	my @line = split / +/,$_;
	my $d = $line[0];

	next if $d ne $utc;

	for(my $i = 1; $i < scalar(@line); $i++)
	{
	    last if $line[$i] eq '#';

	    my @tmp = split /\//,$line[$i];

	    if(scalar(@tmp) == 1)
	    {
		for(my $iscope = 0; $iscope < 4; $iscope++)
		{
		    if(!defined @$laser[$iscope])
		    {
			@$laser[$iscope] = $tmp[0];
		    }
		}
	    }
	    else
	    {
		@$laser[$tmp[0]] = $tmp[1];
	    }	       
	}
    }
    close(LASER);

    if(scalar(@$laser)>0)
    {
	return 1;
    }
    else
    {
	return find_laser_db($data_run,$laser);
    }
}

sub find_laser_db
{
    my $run = shift;
    my $laser = shift;
    my $dbh =
	DBI->connect("DBI:mysql:VERITAS:romulus.ucsc.edu", "readonly", "")
	or die "Could not connect to DB: !\n";

#     my $query= 
# 	"select t1.group_id, t2.excluded_telescopes from " .
# 	"tblRun_Group as t1, " .
# 	"tblRun_GroupComment as t2 where t1.group_id=t2.group_id " .
# 	"and t2.group_type=\"laser\" and t1.run_id=$run;";
    my $query = qq{        
	SELECT info.run_id, grp_cmt.excluded_telescopes 
	    FROM tblRun_Info AS info, tblRun_Group AS grp, 
	    tblRun_GroupComment AS grp_cmt, 
	    (SELECT group_id FROM tblRun_Group WHERE run_id=$run) AS run_grp 
	    WHERE grp_cmt.group_id = run_grp.group_id 
	    AND grp_cmt.group_type='laser' 
	    AND grp_cmt.group_id=grp.group_id 
	    AND grp.run_id=info.run_id 
	    AND (info.run_type='laser' OR info.run_type='flasher')
	};

    my @group_id;
    my $group_count = 0;
    my @laser_id;
    my %laser_excluded;

    my $stmt = $dbh->prepare($query);
    $stmt->execute();

    my @row;
    while(@row = $stmt->fetchrow_array())
    {
	push @laser_id, $row[0];
	$laser_excluded{$row[0]} = $row[1]; 
    }

    my $max_scope = 4;
    for( my $iscope=0; $iscope < $max_scope; $iscope++)
    {
      	my $scope_bit = 2**($iscope);
	
      	# Check all laser run exclusions associated with 
        # this data run against telescopes in this run.

      	foreach my $laser_run (@laser_id)
      	{
	   # Determine telescopes not excluded from this laser group.
	   # Use '0+' to ensure masks are interpreted as integer values.
	   if( ( ~(0+$laser_excluded{$laser_run}) & $scope_bit ) == 
	       $scope_bit )
	   {
	     @{$laser}[$iscope] = $laser_run;
	   }
      	}
      }

#     $query = 
# 	"select t2.run_id, t2.db_start_time " .
# 	"from tblRun_Group as t1, tblRun_Info as t2 " .
#  	"where t1.run_id=t2.run_id and t2.run_type=\"laser\" and " .
# 	"t1.group_id=$group_id;";

#     $stmt = $dbh->prepare($query);
#     $stmt->execute();
#     if($stmt->rows) 
#     {
# 	my @row = $stmt->fetchrow_array();
# 	$row[1] =~ s/\-//g;	
# 	my ($date,$time) = split ' ',$row[1];
# 	push @laser, $row[0];
#     } 

    if(scalar(@{$laser})) 
    {
	return 1;
    }
    else
    {
	return 0;
    }
}

sub get_runtype
{
    my $run = $_[0];
    my $dbh =
	DBI->connect("DBI:mysql:VERITAS:romulus.ucsc.edu", "readonly", "")
	or die "Could not connect to DB: !\n";

    my $query= 
	"select run_type from tblRun_Info where run_id=$run";

    my $stmt = $dbh->prepare($query);
    $stmt->execute();
    if($stmt->rows) 
    {
	my @run_type = $stmt->fetchrow_array();
	return $run_type[0];
    }
    
    return "none";
}



__END__

=head1 NAME

Script for running chila stage1/stage2 data reduction and diagnostics.
The script will search for VBF files in the subdirectories of
B<vbf_dir> (default: /veritas/data) that are named dYYYYMMDD.  Each
VBF file will be processed with the appropriate calibration laser run
defined in B<laser_table> or the database if no laser run is defined
in B<laser_table>.  If no laser run is specified in either the laser
table or the database, runs will be processed with all gains set to 1.
The UTC dates to be processed can be specified using either the B<utc>
option to define a list of dates (e.g. --utc=20101001,20101002) or the
B<nday> option that will process all dates that are less than N days
old.  The default stage1/stage2 command-line options can be overriden
by specifying a two-column configuration file with the B<cfg> option
which has the format:

option      value

Reduced stage1/stage2 files for each UTC date are written to a
corresponding subdirectory in B<local_dir> named dYYYYMMDD and
diagnostics plots are written to B<diagnostics_dir> (default:
/veritas/chila_diagnostics).  A soft link is created in
/veritas/chila_reduced for every reduced data directory if the
B<local_dir> destination differs from that directory.

=head1 SYNOPSIS

run_diagnostics.pl [options] 

Options:

=over 8

=item B<--help> 

Print this help message.

=item B<--man> 

Print the man page.

=item B<--utc>

Run diagnostics for all data from one or more UTC dates.  Date should
be in the format YYYYMMDD.  To run diagnostics for multiple days, UTC
dates can be specified as a comma-separated list and/or with wildcards
(e.g. 20090101,200902\* will run Jan. 1st and all dates in Feb).

=item B<--nday> (default: 4)

Run diagnostics for all data that is less than ndays old. 

=item B<--overwrite>

When processing diagnostics overwrite existing stage1/stage2 files.

=item B<--vbf_dir> (default: /veritas/data)

Directory that will be searched for VBF files.  Data files are assumed
to be organized according to UTC date in subdirectories with the
naming convention dYYYYMMDD.

=item B<--local_dir> 

Directory to which stage1/stage2 reduced data files will be written.
Data files will be organized by UTC date with directory names
dYYYYMMDD.

=item B<--diagnostics_dir> (default: /veritas/chila_diagnostics)

Directory to which diagnostics plots/data will be written.  Diagnostics
will be organized by UTC date in subdirectories named dYYYYMMDD.

=item B<--cfg>

Set the configuration file defining the command-line options for
stage1/stage2.

=item B<--laser_table> (default: /veritas/chila_diagnostics/laser.dat)

Set the table defining the calibration laser run to be used to process all
runs on the corresponding date.  The format of this table should be

YYYYMMDD <RUN>

Separate laser runs for each telescope can be specified by prepending
the run number with the telescope ID and a slash.  For instance to
process telescopes 2-4 with <RUN1> and telescope 1 with <RUN2> use:

YYYYMMDD <RUN1> 0/<RUN2>



=back
