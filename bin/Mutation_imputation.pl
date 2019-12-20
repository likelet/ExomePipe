#!/usr/bin/perl -w
=head1 #===============================================================================
#        USAGE: perl Mutation_imputation.pl <Maf_file> <fasta_file> [bamdir] [outdir]
#
#  DESCRIPTION: Imputate mutation sites of samples which come from the same patient (Muti-sampling)
#
#  INPUT FILES: Merged mutation maf file
#
# REQUIREMENTS:
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: //2017
#     REVISION: ---
#===============================================================================
=cut
die `pod2text $0` unless @ARGV == 4;
use strict;
use File::Basename;
use Data::Dumper;

my $maf = shift or die $!;
my $ref_fa = shift or die $!;
my $bamdir = shift or die $!;
my $outdir = shift or die $!;
my (%info, %mut, %type, %ref, %alt, %sample, %id);
my $file_name = basename($maf);
$file_name =~ /(\S+)\.maf/;

open IN, $maf or die $!;
open OUT, ">$outdir/$1.imputed.maf" or die $!;
open OUT2, ">$outdir/$1.imputed_raw.maf" or die $!;
my $header = <IN>;
print OUT $header;
print OUT2 $header;
chomp $header;
my @ID = split /\t/, $header;
for (0..$#ID){
	$id{$ID[$_]} = $_;
}

while(<IN>){
	chomp;
	my @F = split /\t/;
	#$F[$id{"Tumor_Sample_Barcode"}] =~ /(P\d+)/;
	$F[$id{"Tumor_Sample_Barcode"}] =~ /(31\d+-.+)$/;
	$info{$1}{$F[$id{"Tumor_Sample_Barcode"}]}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $_;
	$mut{$1}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $_;
	$type{$1}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $F[$id{"Variant_Type"}];
	$ref{$1}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $F[$id{"Reference_Allele"}];
	$alt{$1}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $F[$id{"Tumor_Seq_Allele2"}];
	$sample{$1}{$F[$id{"Tumor_Sample_Barcode"}]} = ();
}
print(Dumper%sample);
`mkdir -p $outdir/Binary_mat`;

my (%binary, %count);
for my $k1(sort keys %sample){
	my $m = keys %{$sample{$k1}};
	my $n = 0;
	for my $k2(sort keys %{$sample{$k1}}){
		$n ++;
		for my $k3(1..22,"X","Y"){
			next unless exists $mut{$k1}{$k3};
			for my $k4(sort {$a<=>$b} keys %{$mut{$k1}{$k3}}){
				my @F = split /\t/, $mut{$k1}{$k3}{$k4};
#				if(exists $info{$k1}{$k2}{$k3}{$k4}){
#					print OUT $info{$k1}{$k2}{$k3}{$k4}."\n"; 
#					print OUT2 $info{$k1}{$k2}{$k3}{$k4}."\n";
#					$binary{$k1}{"$k3:$k4-$F[0]"} .= "\t1";
#					$count{$k1}{"$k3:$k4-$F[0]"} += 1 * 10**($m-$n); ## Score this mut_site of this sample
#				}else{
					my $line = `samtools mpileup -A -x -B -q 20 -Q 20 -r chr$k3:$k4-$k4 -f $ref_fa $bamdir/$k2\_sort_dedup_realigned_recal.bam`;
					$binary{$k1}{"$k3:$k4-$F[0]"} .= "\tNA" and next unless $line;  ## There is no mapped reads.
					my @G = split /\t/, $line;
					my ($mutant, $alt_read_count);
					if ($F[$id{"Reference_Allele"}] eq "SNP"){
						$mutant = $F[$id{"Tumor_Seq_Allele2"}];
						$alt_read_count = () = $G[4] =~ /($mutant)/ig;
					}else{
						$mutant = $F[$id{"Tumor_Seq_Allele2"}];
						$mutant =~ s/^$F[$id{"Reference_Allele"}]//;
						$alt_read_count = () = $G[4] =~ /($mutant)/ig;
					}
					my $ref_read_count = () = $G[4] =~ /(\.|,)/ig;
					$binary{$k1}{"$k3:$k4-$F[0]"} .= "\tNA" and next if $alt_read_count + $ref_read_count == 0;  ## Special type SNV
					my $vaf = int($alt_read_count * 10000 / ($alt_read_count + $ref_read_count)) / 10000;
					if ($alt_read_count < 3 and $G[3] <= 10){  ## Reads depth is not enough or no difference between cancer and ref
						$binary{$k1}{"$k3:$k4-$F[0]"} .= "\tNA";
						print OUT2 join "\t", $F[$id{"Hugo_Symbol"}], $F[$id{"Chromosome"}], $F[$id{"Start_Position"}], $F[$id{"End_Position"}], $F[$id{"Variant_Classification"}], $F[$id{"Variant_Type"}], $F[$id{"Reference_Allele"}], $F[$id{"Tumor_Seq_Allele1"}], $F[$id{"Tumor_Seq_Allele2"}], $ref_read_count, $alt_read_count, $vaf, $F[$id{"CDS_Change"}], $F[$id{"Protein_Change"}], $k2;
						print OUT2 "\n";
						next;
					}elsif($alt_read_count < 3 and $G[3] > 10){  ## No mutations
						$binary{$k1}{"$k3:$k4-$F[0]"} .= "\t0";
						print OUT2 join "\t", $F[$id{"Hugo_Symbol"}], $F[$id{"Chromosome"}], $F[$id{"Start_Position"}], $F[$id{"End_Position"}], $F[$id{"Variant_Classification"}], $F[$id{"Variant_Type"}], $F[$id{"Reference_Allele"}], $F[$id{"Tumor_Seq_Allele1"}], $F[$id{"Tumor_Seq_Allele2"}], $ref_read_count, $alt_read_count, $vaf, $F[$id{"CDS_Change"}], $F[$id{"Protein_Change"}], $k2;
                                                print OUT2 "\n";
						next;
					}elsif($alt_read_count >= 3){  ## Inputed mutations
						$binary{$k1}{"$k3:$k4-$F[0]"} .= "\t1";
						$count{$k1}{"$k3:$k4-$F[0]"} += 1 * 10**($m-$n);  ## Score this mut_site of this sample
					}
					print OUT join "\t", $F[$id{"Hugo_Symbol"}], $F[$id{"Chromosome"}], $F[$id{"Start_Position"}], $F[$id{"End_Position"}], $F[$id{"Variant_Classification"}], $F[$id{"Variant_Type"}], $F[$id{"Reference_Allele"}], $F[$id{"Tumor_Seq_Allele1"}], $F[$id{"Tumor_Seq_Allele2"}], $ref_read_count, $alt_read_count, $vaf, $F[$id{"CDS_Change"}], $F[$id{"Protein_Change"}], $k2;
					print OUT "\n";
					print OUT2 join "\t", $F[$id{"Hugo_Symbol"}], $F[$id{"Chromosome"}], $F[$id{"Start_Position"}], $F[$id{"End_Position"}], $F[$id{"Variant_Classification"}], $F[$id{"Variant_Type"}], $F[$id{"Reference_Allele"}], $F[$id{"Tumor_Seq_Allele1"}], $F[$id{"Tumor_Seq_Allele2"}], $ref_read_count, $alt_read_count, $vaf, $F[$id{"CDS_Change"}], $F[$id{"Protein_Change"}], $k2;
                                        print OUT2 "\n";
#				}
			}
		}
	}
}

##  Binary matrix output
for my $k1(sort keys %sample){
	open BM, ">$outdir/Binary_mat/$k1\_binary_mat.txt" or die $!;
	print BM "\t";
	my $line = join "\t", sort keys %{$sample{$k1}};
	$line =~ s/P\d+//g;
	print BM "$line\tN\n";
	for my $k2(sort {$count{$k1}{$a} <=> $count{$k1}{$b}} keys %{$count{$k1}}){
		print BM "$k2$binary{$k1}{$k2}\t0\n";
	}
	close BM;
}


=t

	for my $k3(1..22,"X","Y"){
		next unless exists $mut{$k1}{$k3};
		for my $k4(sort {$a<=>$b} keys %{$mut{$k1}{$k3}}){
			my @F = split /\t/, $mut{$k1}{$k3}{$k4};
			my $mark = "$k3:$k4-$F[0]";
			print BM "$mark$binary{$k1}{$mark}\t0\n";
		}
	}
	close BM;
}























