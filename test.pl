#!/usr/bin/env perl
#
# original filename: test.pl
#
# Steven Bickerton
#  Dept. of Astrophysical Sciences, Princeton University
#  bick@astro.princeton.edu
#  Created: Wed Aug 20, 2008  19:19:58 DST
#  Host: bender.astro.Princeton.EDU
#  Working Directory: /Users/bick/usr/src/Canalysis
#


use strict;
use warnings;
use File::Basename;

my $exe = basename($0);
my $usage = "Usage: $exe exefiles\n";

my @binaries = @ARGV;
die $usage unless @binaries;

my %execinfo = ( 
    "fftSmooth"   => ["test_timeseries.dat", "1:2 2.0", 41629, 0],
    "rmsSlide"    => ["test_timeseries.dat", "2.0 1 2", 11485, ""],
    "smoothC"     => ["test_timeseries.dat", "1:2 1 2.0 smoothC.par", 53320, ""],
    "lombscargle" => ["test_sinx.dat", "", 16979, ""],
    "getPeaks"    => ["test_sinx.dat", "0.9", 4849, ""], 
    "pca"         => ["test_line4D.dat", "", 20342, ""],
    "normC"       => ["test_timeseries.dat", "1:2 1 2.0 normC.par", 41678, 0],
    "fitsplitN2"  => ["test_timeseries.fits", "10 13", 51348, "test_timeseries.00.fits"],
    "fitsTSdump"  => ["test_timeseries.fits", "100 10", 28616, ""],
    "fitsTSsample"=> ["test_timeseries.fits", "10", 29562, ""],
    "ran_poisson" => ["", "100 10 2", 0, ""],
    "convolve"    => ["test_timeseries.dat", "test_kernel.dat", 53693, ""],
    "text2fits"   => ["test_kernel.dat", "", 3930, "test_kernel.dat.fits"],
    "text2fitsTS" => ["test_kernel.dat", "", 63551, "test_kernel.dat.fits"],
    "fits2text"   => ["test_timeseries.fits", "", 60185, ""],
    "bin"         => ["test_timeseries.dat", "2 0.1 1", 28387, ""],
    "bin2d"       => ["test_data4col.dat", "2:4 2:2", 41816, ""],
    "cstats"      => ["test_data4col.dat", "", 55023, ""],
    "pdm"         => ["test_sinx.dat", "", 9151, ""],
    "notchFilter" => ["test_sinx.dat", "", 9532, ""],
    "erf"         => ["", "2.0", 29290, ""],
    "nsigma"      => ["", "43.9558", 40668, ""],
    "interp"      => ["test_kernel.dat", "0:0.1:1.4", 31117, ""]
    );


sub test_command($$$$$) {
    my ($cmd, $args, $out, $size, $special_out_file) = @_;

    # run it
    printf STDOUT "Testing %-16s %3s", "$cmd", "...  ";
    ( -f $out ) and unlink $out;
    ( length($special_out_file) > 1 and -f $special_out_file ) and 
	unlink $special_out_file;
    if ($special_out_file) {
	system("./$cmd $args");
	$out = $special_out_file;
    } else {
	system("./$cmd $args > $out");
    }

    # see if file was written
    if ( -s $out ) {
	my ($testsize) = (split /\s+/, `sum $out`);
	#unlink $out;
	if ( ($testsize == $size) or ($size == 0 && $testsize > 0) ) {
	    printf STDOUT "Success\n";
	} else {
	    printf STDOUT "FAILED - Outfile wrong size\n";
	}
	
    } else {
	printf STDOUT "FAILED - No outfile written\n";
    }
}

foreach my $binary ( sort @binaries ) {

    if ( $execinfo{$binary} ) {
	my ($infile, $args, $size, $out) = @{$execinfo{$binary}};
	test_command($binary, "$infile $args", $binary.".testout", $size, $out);
    } else {
	printf STDOUT "Don't know about $binary ... skipping\n";
    }

}


exit 0;
