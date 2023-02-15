#! /usr/bin/perl

use strict;
use warnings;

use Parallel::ForkManager;
use Cwd;

my $nOcc = 0;
my $nVirt = 0;
my $chkfile = ""; 
my $g16root = $ENV{g16root};
my $i = 0;
my $j = 0;
my $MAX_PROCESSES = 4;
my $nproc=1;
my $mem="1GB";
my $sub = 1;

foreach my $temp (@ARGV){
      chomp($temp);
      if($temp =~ /^-nproc=(\d+)$/i){
        $nproc=$1;
      }elsif($temp =~ /^-mem=(\d+[a-zA-Z][a-zA-Z])$/i){
        $mem=$1;
      }elsif($temp =~ /^-nosub/i){
        $sub = 0;
      }elsif($temp =~ /^-nocc=(\d+)$/i){
        $nOcc = $1;
      }elsif($temp =~ /^-nvirt=(\d+)$/i){
        $nVirt = $1;
      }elsif($temp =~ /^-chk=(\S+)$/i){
        $chkfile=$1;
      }
}
die "Script requires nocc, nvirt and chk options to be set.\n" if ($nOcc eq 0 || $nVirt eq 0 || $chkfile eq "");

my $curdir = getcwd;

for (my $i=1; $i<=$nOcc; $i++){
	for (my $j=$nOcc+1; $j<=$nOcc+$nVirt; $j++){
        	my $gjf = "nocis_alter_$i"."_$j";
        	open  (GJF,">$gjf.com") or die "Could not create $gjf.\n";
        	print GJF "%oldchk=$curdir/$chkfile\n";
        	print GJF "%chk=$curdir/$gjf.chk\n";
        	print GJF "#P uhf chkbas guess=(read,alter) geom=check output=matrixelement\n";
        	print GJF "# IOp(5/200=3) scf=(conver=10,novaracc,maxcycles=250,fermi) nosymm\n\n";
        	print GJF "NOCIS orbital swap. Occ = $i, Virt = $j\n\n";
        	print GJF "0 1\n\n";
        	printf GJF "$i $j\n\n\n";
        	print GJF "$curdir/$gjf.mat\n\n";
        	close GJF;
	}
}

if ($sub eq 1){
        my @files = <*.com>;

        my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

        foreach my $file (@files){
                $pm->start and next;

                my $log_file = $file;
                unless($log_file =~ s/\.(?:gjf|com)/.log/){
                        $log_file .= ".log";
                }
                my $sys_cmd = "$g16root/g16/g16 -m=$mem -p=$nproc < $file > $log_file";
                system("$sys_cmd");
                print "Job $file submitted.\n";

                $pm->finish;
        }
        $pm->wait_all_children;
}
