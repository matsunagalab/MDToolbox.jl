#!/usr/bin/perl
#
use strict;
use warnings;
my $opt={};
$opt->{nimage}=40;
$opt->{header}="pr";
$opt->{path_suf}="rpath";
$opt->{pdbfile} = "../0_system/sys_selected_13_no.pdb";

readpdb($opt);
readpath($opt);
outpdb($opt);

sub outpdb {
	my ($opt)=@_;
	for (my $ifile =0;$ifile<scalar(@{$opt->{data}});$ifile++) {
		my $ofile=sprintf("%s.%d.pdb",$opt->{header},$ifile+1);
		open(OUTFILE,">".$ofile) || die "Cannot open $ofile\n";
		my $ic=0;
		for (my $j = 0;$j<scalar(@{$opt->{data}->[$ifile]});$j+=3) {
			printf OUTFILE "%s%8.3f%8.3f%8.3f%s",
				   $opt->{atom}->[$ic],
				   $opt->{data}->[$ifile]->[$j],
				   $opt->{data}->[$ifile]->[$j+1],
				   $opt->{data}->[$ifile]->[$j+2],
				   $opt->{segid}->[$ic];
				   $ic++;
		}
		printf OUTFILE "END\n";
		close(OUTFILE);
	}
}
sub readpath {
	my ($opt)=@_;
	$opt->{path}=sprintf("%s.tot.%s_last",$opt->{header},$opt->{path_suf});
	open(OUTFILE,">".$opt->{path}) || die "Cannot open $opt->{path}\n";
	for (my $i = 1;$i<=$opt->{nimage};$i++) {
		my $file=sprintf("%s.%d.%s",$opt->{header},$i,$opt->{path_suf});
		open(INFILE,$file)|| die "Cannot open $file\n";
		my $line="";
		while(<INFILE>){
			if (/^\s+[0-9]./) {
				$line=$_;
			}
		}
		close(INFILE);
		printf OUTFILE "%s",$line;
		my $file_last=sprintf("%s.%d.%s_last",$opt->{header},$i,$opt->{path_suf});
	    open(OUTFILE_IN,">".$file_last) || die "Cannot open $file_last\n";
		printf OUTFILE_IN "%s",$line;
		close(OUTFILE_IN);
	}
	close(OUTFILE);

	open(INFILE,$opt->{path}) || die "Cannot open $opt->{path}\n";
	@{$opt->{data}}=();
	my $ii=0;
	while(<INFILE>){
		if (/^\s+[0-9]./) {
			@{$opt->{data}->[$ii]}=();
			my (@dum)=split(' ');
			for (my $i = 1;$i<scalar(@dum);$i++) {
				$opt->{data}->[$ii]->[$i-1]=$dum[$i];
			}
			$ii++;
		}
	}
	close(INFILE);
}

sub readpdb {
	my ($opt)=@_;
	open(INFILE,$opt->{pdbfile}) || die "Cannot open $opt->{pdbfile}\n";
	@{$opt->{atom}}=();
	@{$opt->{segid}}=();
	while(<INFILE>){
		if (/^ATOM/) {
			my $top=substr($_,0,30);
			push(@{$opt->{atom}},$top);
			my $last=substr($_,54);
			push(@{$opt->{segid}},$last);
		}
	}
	close(INFILE);
}
