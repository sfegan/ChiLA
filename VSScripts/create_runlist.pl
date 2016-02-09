#!/usr/bin/perl -w

my $list = $ARGV[0];


# my %h5_file_hash = ();
# my @h5_data_dir = ("/home/mdwood/reduced");


# foreach my $data_dir (@h5_data_dir)
# {
#     print "Searching for data files in $data_dir.\n";
    
#     my @h5_files = <$data_dir/d*/x*_s2.h5>;
    
#     foreach my $h5_file (@h5_files) 
#     {    
# 	if($h5_file =~ /x(\d+)_s2.h5/)
# 	{
# 	    $h5_file_hash{$1} = $h5_file;
# 	}
#     }
# }


open(LIST,"$list") or die print "Couldn't open list: $list\n";
while(<LIST>) 
{
    next if /^\#/ || /^$/;
    
    chomp;
    my @line = split ' ',$_;
    my $ncol = scalar(@line);

    my ($date, $run) = @line;
    if(-e "/home/mdwood/reduced/d${date}/x${run}_s2.h5")
    {
	print "/home/mdwood/reduced/d${date}/x${run}_s2.h5\n";
    }
#     else
#     {
# 	print "$h5_file_hash{$run}\n";
#     }
}
