#!/usr/local/perl

open (FILE,"$ARGV[0]");
while (<FILE>){
    if ($_=~/^NAME\s+(.+)/){
	$id=$1;
	print "$_";
	print "ACC\t$id\n";
    }else{
	print "$_";
    }
}
close(FILE);
