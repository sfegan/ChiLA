#!/usr/bin/perl -w
#
# \author     Matthew Wood                \n
#             UCLA                        \n
#             mdwood@astro.ucla.edu       \n        
#
# \version    0.1
# \date       01/20/2009

use File::Basename;
use Cwd 'abs_path';

my @dir;
my @lnk_dir;

push @dir, </net/gamma3/srv/raid4/chila_analysis/*>;
push @dir, </net/gamma3/srv/raid5/chila_analysis/*>;
push @lnk_dir, </veritas/chila_reduced/*>;

foreach my $dir (@lnk_dir)
{
    if(-l "$dir")
    {
	`rm $dir`;
    }
}

foreach my $dir (@dir)
{
    my $date = basename($dir);

    my @tmp = grep(/$date/,@lnk_dir);
    my @files = <$dir/*>;

    
    if(scalar(@files) == 0)
    {
	print "NO FILES $dir\n";
	`rm -rf $dir`;
    }
    elsif(scalar(@tmp) == 0)
    {
	print "ln -s $dir /veritas/chila_reduced/$date\n";
	`ln -s $dir /veritas/chila_reduced/$date`;
    } 

}

