#!/usr/bin/perl

my $output_dir = $ARGV[0];


@h5files = <*.h5>;

foreach $h5file (@h5files)
{
    if(!(-e "$output_dir/$h5file"))
    {
	print "$h5file\n";
	print `h5repack -i $h5file -o $output_dir/$h5file -f GZIP=1`;
    }
}
