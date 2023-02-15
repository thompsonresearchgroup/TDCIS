#! /usr/bin/perl

use strict;
use warnings;

my %result = ();

open(INPUT,$ARGV[0]) or die "Can't open file $ARGV[0]\n";

my $time = 0;
my $dipole = '';
my $field = '';

print "#      Time (au)    |                                            Field (au)                                              |                                               Dipole (au)\n";
print "#                   |      real x          aimag x         real y          aimag y         real z          aimag z       |         real x          aimag x         real y          aimag y         real z          aimag z\n";

while(<INPUT>){
        if (/^ Time step:\s+(\d+\.\d+)\s+au/){
                $time = $1;
        }
        if (/Applied field vector/){
		for(my $i = 0; $i <= 2; $i++){
			@_ = split(' ', <INPUT>);
			if (scalar @_ == 2){
				$field .= sprintf("\t% .10f\t 0.0000000000",$_[1]);
			}elsif(scalar @_ == 3){
				$field .= sprintf("\t% .10f\t% .10f",$_[1],substr($_[2],0,-1));	
			}else{
				die('Applied field data in unknown format');
			}
        	}
        }
        if (/Total dipole/){
		for(my $i = 0; $i <= 2; $i++){
			@_ = split(' ', <INPUT>);
			if (scalar @_ == 2){
				$dipole .= sprintf("\t% .10f\t 0.0000000000",$_[1]);
			}elsif(scalar @_ == 3){
				$dipole .= sprintf("\t% .10f\t% .10f",$_[1],substr($_[2],0,-1));	
			}else{
				die('Dipole data in unknown format');
			}
        	}
		printf("%15.4f\t$field\t$dipole\n",$time);
		$field = '';
		$dipole = '';
        }
}
#
#if ($finish ne 1){
#        print "$file did not finish\n";
#}
#close INPUT;
#
#open(OUTPUT,'>td_data.dat') or die "Can't open td_data.dat\n";
#
#foreach my $name (sort keys %result) {
#    printf OUTPUT "%4.2f\t%28s\n", $name, $result{$name};
#}

