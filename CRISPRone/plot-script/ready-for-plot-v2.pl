#!/usr/local/perl

$output_plot=$ARGV[1];
$output_spacer=$ARGV[2];
$all_spacer_file=$output_spacer."/all-spacer.fa";

open (FILE,"$ARGV[0]"); #k6a6cfIi.plot.spacer
while (<FILE>){
    if ($_=~/^\#\#(.+)/){
	$n_spacer=1;
	$new_file=$output_plot."/".$1.".plot.txt";
	$new_spacer=$output_spacer."/".$1.".plot.spacer";
	$contig_name=$1;
    }else{
	open (NEW,">>$new_file");
	@tmp=split(/\s+/,$_);
	if ($tmp[4]=~/^repeat:/){
	    @repeat_spacer=split(/\-/,$tmp[5]);
	    $spacer_uniq_id=$tmp[1]."-".$tmp[2];
	    $spacer_id="s".$n_spacer;
	    $tmp_order="";
	    for ($i=0;$i<=@repeat_spacer-1;$i++){
		if ($repeat_spacer[$i]=~/^r\d+$/){
		    $tmp_order=$tmp_order.$repeat_spacer[$i]."-";
		    next;
		}
		$spacer_id="s".$n_spacer.":unk";
		$tmp_order=$tmp_order.$spacer_id."-";
		$len_spacer=length($repeat_spacer[$i]);
		open (SPNEW,">>$new_spacer");
		open (ALL,">>$all_spacer_file");
		print SPNEW ">$contig_name\_s$n_spacer\[$spacer_uniq_id\]\:$len_spacer\n";
		print SPNEW "$repeat_spacer[$i]\n";
		print ALL ">$contig_name\_s$n_spacer\[$spacer_uniq_id\]\:$len_spacer\n";
		print ALL "$repeat_spacer[$i]\n";
		close(ALL);
		close(SPNEW);
		$n_spacer++;
	    }
	    chop $tmp_order;
	    print NEW "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp_order\n";
	}else{
	    print NEW $_;
	}
	close(NEW);
    }
}
close(FILE);
