#!/usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl annovar2maf.pl Mut_multianno.lst > Merge_mut.maf
#
#  DESCRIPTION: Change annovar output into MAF format
#
#  INPUT FILES: List file containing annovar multi-anno result and sample ID (tab separated)
#
# REQUIREMENTS:
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 12/11/2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 1;
use strict;
use File::Basename;

print "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tRef_allele_depth\tAlt_allele_depth\tVAF\tCDS_Change\tProtein_Change\tTumor_Sample_Barcode\tcosmic70\tavsnp147\tALL.sites.2015_08\tesp6500siv2_all\tExAC_ALL\n";

open IN, shift or die $!;

while (<IN>){
	chomp;
	my $file = (split /\t/, $_)[0];
	my $sample_id = (split /\t/, $_)[1];
	open IN2, $file or die $!;
	my %index;
	while (<IN2>){
		chomp;
		if (/^Chr/){
			my @ID = split /\t/;
			for (0..$#ID){
				$index{$ID[$_]} = $_;
			}
			next;
		}
		my @F = split /\t/;
		my @G = split /:/, $F[$index{"Otherinfo"}+13];
		my @H = split /,/, $G[1];
		my $dp_ref = $H[0];
		my $dp_alt = $H[1];
		next if $dp_alt + $dp_ref == 0;
		my $vaf = int($dp_alt * 100 / ($dp_alt + $dp_ref)) / 100;
		$F[$index{"ExonicFunc.refGene"}] = "Splice_Site" if $F[$index{"Func.refGene"}] =~ "splicing";
		$F[$index{"ExonicFunc.refGene"}] = "Intron" if $F[$index{"Func.refGene"}] =~ "intronic";
		$F[$index{"ExonicFunc.refGene"}] = "IGR" if $F[$index{"Func.refGene"}] =~ "intergenic";
		$F[$index{"ExonicFunc.refGene"}] = "3'Flank" if $F[$index{"Func.refGene"}] =~ "downstream";
		$F[$index{"ExonicFunc.refGene"}] = "5'Flank" if $F[$index{"Func.refGene"}] =~ "upstream";
		$F[$index{"ExonicFunc.refGene"}] = "3'UTR" if $F[$index{"Func.refGene"}] =~ "UTR3";
		$F[$index{"ExonicFunc.refGene"}] = "5'UTR" if $F[$index{"Func.refGene"}] =~ "UTR5";
		$F[$index{"ExonicFunc.refGene"}] = "RNA" if $F[$index{"Func.refGene"}] =~ /ncRNA_/;
		next if $F[$index{"ExonicFunc.refGene"}] eq "unknown";
		$F[$index{"ExonicFunc.refGene"}] =~ s/^nonsynonymous SNV/Missense_Mutation/;
		$F[$index{"ExonicFunc.refGene"}] =~ s/^synonymous SNV/Silent/;
		$F[$index{"ExonicFunc.refGene"}] =~ s/^nonframeshift insertion/In_Frame_Ins/;
		$F[$index{"ExonicFunc.refGene"}] =~ s/^nonframeshift deletion/In_Frame_Del/;
		$F[$index{"ExonicFunc.refGene"}] =~ s/^frameshift deletion/Frame_Shift_Del/;
		$F[$index{"ExonicFunc.refGene"}] =~ s/^frameshift insertion/Frame_Shift_Ins/;
		$F[$index{"ExonicFunc.refGene"}] =~ s/^stopgain/Nonsense_Mutation/;
		$F[$index{"ExonicFunc.refGene"}] =~ s/^stoploss/Nonstop_Mutation/;
		$F[$index{"AAChange.refGene"}] = "." unless $F[$index{"AAChange.refGene"}];
		my @J = split /,/, $F[$index{"AAChange.refGene"}];
		my ($cds_change, $aa_change);
		if ($J[0] =~ /(c\.\S+)\:/){
			$cds_change = $1;
		}else{
			$cds_change = "NA";
		}
		if ($J[0] =~ /(p\.\S+)/){
			$aa_change = $1;
		}else{
			$aa_change = "NA";
		}
		my $start_pos = $F[$index{"Start"}];
		my $end_pos = $start_pos + length($F[$index{"Ref"}]) - 1;
		my $mut_type;
		if (length($F[$index{"Ref"}]) > length($F[$index{"Alt"}])){
			$mut_type = "DEL";
		}elsif(length($F[$index{"Ref"}]) < length($F[$index{"Alt"}])){
			$mut_type = "INS";
		}else{
			$mut_type = "SNP";
		}
		my $mut_allele1;
		#$mut_allele1 = $F[$index{"Ref"}] if $G[0] =~ /0\/1/;
		$mut_allele1 = $F[$index{"Alt"}] ;
		#if $G[0] =~ /1\/1/;
		$F[$index{"Chr"}] =~ s/chr//;
		$F[$index{"cosmic70"}] = "." unless $F[$index{"cosmic70"}];
		$F[$index{"avsnp147"}] = "." unless $F[$index{"avsnp147"}];
		$F[$index{"ALL.sites.2015_08"}] = "." unless $F[$index{"ALL.sites.2015_08"}];
		$F[$index{"esp6500siv2_all"}] = "." unless $F[$index{"esp6500siv2_all"}];
		$F[$index{"ExAC_ALL"}] = "." unless $F[$index{"ExAC_ALL"}];
		print join "\t", $F[$index{"Gene.refGene"}], $F[$index{"Chr"}], $start_pos, $end_pos, $F[$index{"ExonicFunc.refGene"}], $mut_type, $F[$index{"Ref"}], $mut_allele1, $F[$index{"Alt"}], $dp_ref, $dp_alt, $vaf, $cds_change, $aa_change, $sample_id, $F[$index{"cosmic70"}], $F[$index{"avsnp147"}], $F[$index{"ALL.sites.2015_08"}], $F[$index{"esp6500siv2_all"}], $F[$index{"ExAC_ALL"}];
		print "\n";
	}
}