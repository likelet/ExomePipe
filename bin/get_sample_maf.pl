#! /cygdrive/c/Perl64/bin/perl  -w
use strict;
use Data::Dumper;
#use List::MoreUtils qw/arrayify uniq/;

open MERGE,shift;
my $outdir=$ARGV[1]

my %patient_snp;
my %sample_maf;
while(<MERGE>){
	chomp;
	next if /Hugo_Symbol/;
	my @line = split /\t/,$_,-1;
	my $snp = join("_",@line[0,1,2,3,6,8]);
	my $sample_id = $line[-1];
	push @{$sample_maf{$sample_id}},$_;
	$sample_id =~ s/-/\t/;
	my ($patient,$type) = split /\t/,$sample_id;
	$snp .= "\t$type";
	my $VAF = $line[12-1];
	push @{$patient_snp{$patient}},$snp;
	
}
close MERGE;

##split the total snp_merge.maf into patients
for my $sample (sort keys %sample_maf){
	open SAMPLE_OUT,">$outdir/$sample.imputed.maf";
	print SAMPLE_OUT join("\n",@{$sample_maf{$sample}});
	print SAMPLE_OUT "\n";
}
