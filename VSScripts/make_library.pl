#!/usr/bin/perl -w

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
use Getopt::Long;
use Pod::Usage;
use threads;
use threads::shared;
use ThreadPool;

my %opt = ();
$opt{'help'} = 0;
$opt{'man'} = 0;
$opt{'verbose'} = 0;
$opt{'use_size_tables'} = "true";
$opt{'s3_simple_cuts'} = "";
$opt{'s3_energy_lookup'} = "";
$opt{'s3_sc_parameter_lookup'} = "";
$opt{'mode'} = "sp";
$opt{'o'} = "";
$opt{'nthread'} = "8";
$opt{'overwrite'} = 0;

GetOptions( 'help!'              => \$opt{'help'},
            'man!'               => \$opt{'man'}, 
            's3_simple_cuts=s'   => \$opt{'s3_simple_cuts'},
            's3_energy_lookup=s' => \$opt{'s3_energy_lookup'},
            's3_sc_parameter_lookup=s' => \$opt{'s3_sc_parameter_lookup'},
            'use_size_tables=s'  => \$opt{'use_size_tables'},
            'o=s'                => \$opt{'o'},
	    'mode=s'             => \$opt{'mode'},
            'v!'                 => \$opt{'verbose'},
	    'overwrite!'         => \$opt{'verbose'}
            ) or pod2usage(2);

# -----------------------------------------------------------------------------
# Print help/man if asked or if no file
# -----------------------------------------------------------------------------
pod2usage(1) if $opt{'help'};
pod2usage(-verbose => 2) if $opt{'man'};
pod2usage("$0: No list given.") if ((@ARGV == 0) && (-t STDIN));

my $runlist = $ARGV[0];

open(RUNLIST,$runlist);

my $egy_table_opt =
    "-use_size_tables=$opt{'use_size_tables'} " .
    "-median=true " .
    "-theta_max=0.15 -no_scope_tables=true ";

#    "-s3_energy_lookup=$opt{'s3_energy_lookup'} " .

my $s3_cfg = 
    "-s3_simple_cuts=$opt{'s3_simple_cuts'} " .
    "-s3_sc_parameter_lookup=$opt{'s3_sc_parameter_lookup'} " .
    "-s3_energy_weight_power=1 -s3_msc_weight_power=1 ";

my @s2_files;

while(<RUNLIST>)
{
    next if /^\#/;
    next if /^$/;

    chomp;

    if(-e $_)
    {
	push(@s2_files,$_);
    }
    else
    {
	print "No such stage2 file: $_\n";
    }
}

if($opt{'mode'} eq "sp")
{
    make_sp_library(\@s2_files);
}
elsif($opt{'mode'} eq "egy")
{
    make_egy_library(\@s2_files);
}
elsif($opt{'mode'} eq "krn")
{
    make_krn_library(\@s2_files);
}

#make_egy_tables x_zn20.1_s2.h5 -o egy_table_md_noscope_size_th0.15_recimpact_noreccut.h5 -median=true -no_reconstructed_impact=false -no_reconstruction_cut=true -theta_max=0.15 -use_size_tables=true

sub make_krn_library
{
    my @jobs;

    my @s2_files = @{$_[0]};
    open(KRNLIST,">krn_table_list.txt");

    foreach my $s2_file (@s2_files)
    {
	my $krn_file, my $s3_file, my $egy_file;

	if($s2_file =~ /([^\/]+)_s2.h5/)
	{
	    $s3_file = "${1}_s3.h5";
	    $krn_file = "${1}_krn.h5";
	    $egy_file = "${1}_egy.h5";
	}

	my ($zn,$az,$ped_rms) = get_run_data($s2_file);

	print KRNLIST "$zn $az $ped_rms $krn_file\n";

	push(@jobs, { 's2_file' => $s2_file,
		      's3_cfg' => $s3_cfg,
		      's3_file' => $s3_file,
		      'krn_file' => $krn_file,
		      'egy_file' => $egy_file,
		      'overwrite' => $opt{'overwrite'} } );
			  

	# Run stage3 ----------------------------------------------------------
#	print "stage3 $s2_file $s3_cfg -o $s3_file\n";
#	`stage3 $s2_file $s3_cfg -o $s3_file`;

	# Extract the kernel --------------------------------------------------
#	print "extract_egy_kernel $s3_file -o $krn_file\n";
#	`extract_egy_kernel $s3_file -o $krn_file`;

    }

    close(KRNLIST);

    ThreadPool::run_jobs(\&run_stage3,\@jobs,$opt{'nthread'});


    print "combine_sp_tables krn_table_list.txt -o $opt{'o'}\n";
    print `combine_sp_tables krn_table_list.txt -o $opt{'o'}`;
}

sub make_egy_library
{
    my @s2_files = @{$_[0]};

    open(LIST,">egy_table_list.txt");


    foreach my $s2_file (@s2_files)
    {
	my $egy_table;
	my $s3_file;
	my $krn_file;
		
	if($s2_file =~ /([^\/]+)_s2.h5/)
	{
	    $egy_table = "${1}_egy.h5";
	}

	my ($zn,$az,$ped_rms) = get_run_data($s2_file);

	print LIST "$zn $az $ped_rms $egy_table\n";

	# Create the energy lookup table --------------------------------------
	if(!(-e $egy_table) || $opt{'overwrite'})
	{
	    print "make_egy_tables $s2_file $egy_table_opt -o $egy_table\n";
	    `make_egy_tables $s2_file $egy_table_opt -o $egy_table`;
	}
    }

    close(LIST);

    print "combine_sp_tables egy_table_list.txt -o $opt{'o'}\n";
    print `combine_sp_tables egy_table_list.txt -o $opt{'o'}`;
}

sub make_sp_library
{
    my @s2_files = @{$_[0]};

    open(LIST,">table_list.txt");

    foreach my $s2_file (@s2_files)
    {
	my $sp_table;

	if($s2_file =~ /([^\/]+)_s2.h5/)
	{
	    $sp_table = "${1}_scp.h5";
	}
	else
	{
	    print "Error parsing $s2_file\n";
	    exit(1);
	}

	my $cfg = "-median=true -no_scope_tables=true ";

	print "make_sp_tables $cfg $s2_file -o $sp_table\n";
	`make_sp_tables $cfg $s2_file -o $sp_table`;
	
	my ($zn,$az,$ped_rms) = get_run_data($s2_file);	

	print LIST "$zn $az $ped_rms $sp_table\n";
    }    

    close(LIST);

    print "combine_sp_tables table_list.txt -o $opt{'o'}\n";
    print `combine_sp_tables table_list.txt -o $opt{'o'}`;
}


sub get_run_data
{
    my $file = $_[0];
    my $zn = 
	`octaveio -cat $file sim_header.table.zenith_min_deg | head -n 1`; 
    my $az = 
	`octaveio -cat $file sim_header.table.azimuth_min_deg | head -n 1`;
    my $ped_rms =
	`octaveio -cat $file array_info.mean_scaled_dev | head -n 1`;

    chomp($zn);
    chomp($az);
    chomp($ped_rms);

    return ($zn,$az,$ped_rms);
}

sub run_stage3
{
    my %args = %{$_[0]};

    @_ = ();

    my $overwrite = $args{'overwrite'};


    # Create the energy lookup table --------------------------------------
    if(!(-e $args{'egy_file'}) || $overwrite)
    {
	print "make_egy_tables $args{'s2_file'} " .
	    "$egy_table_opt -o $args{'egy_file'}\n";
	`make_egy_tables $args{'s2_file'} $egy_table_opt -o $args{'egy_file'}`;
    }	

    $args{'s3_cfg'} = $args{'s3_cfg'} . 
	"-s3_energy_lookup=$args{'egy_file'} " .
	"-s3_sim_energy_weight=powerlaw,-2.5 ";

    # Run stage3 ----------------------------------------------------------
    if(!(-e $args{'s3_file'}) || $overwrite)
    {
	print "stage3 $args{'s2_file'} $args{'s3_cfg'} -o $args{'s3_file'}\n";
	`stage3 $args{'s2_file'} $args{'s3_cfg'} -o $args{'s3_file'}`;
    }    

    # Extract the kernel --------------------------------------------------
    print "extract_egy_kernel $args{'s3_file'} -o $args{'krn_file'}\n";
    `extract_egy_kernel $args{'s3_file'} -o $args{'krn_file'}`;
}

