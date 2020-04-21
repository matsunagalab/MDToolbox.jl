#!/usr/bin/perl
#
use strict;
use warnings;

my $opt={};
$opt->{nreplica}=40;
$opt->{k}=100;
$opt->{mdstep}=4500000;
$opt->{mdout}=5000;
$opt->{template}="pathcv.in";
$opt->{tbcpath}="/swork/chig/genesis/genesis_gat_rpath_fast_smd/bin/pathcv_analysis";

&readin($opt);

for (my $nrep = 1;$nrep<=$opt->{nreplica};$nrep++) {
	$opt->{intcl}=sprintf("pathcv.%d.in",$nrep);
	$opt->{num}=$nrep;
	&readout($opt);
	my $command=sprintf("%s %s",$opt->{tbcpath}, $opt->{intcl});
	system($command);
	unlink($opt->{intcl});
}

sub readout {
	my ($opt)=@_;
	open(TCL,">".$opt->{intcl}) || die "Cannot open $opt->{intcl}\n";
	for (my $i = 0;$i<scalar(@{$opt->{lines}});$i++) {
		my $str=$opt->{lines}->[$i];

		foreach my $key(keys(%{$opt})){
			$str=~s/\<$key\>/$opt->{$key}/g;
		}
		print TCL $str;
	}
	close(TCL);
}

sub readin {
	my ($opt)=@_;
	@{$opt->{lines}}=();
	open(INFILE,$opt->{template}) || die "Cannot open $opt->{template}\n";
	while(<INFILE>) {
		push(@{$opt->{lines}},$_);
	}
	close(INFILE);
}
