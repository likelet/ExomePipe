#! /bin/perl

my $sample= $ARGV[1];
open MAF1,shift;
open MAF2,">${sample}.snv.txt";
print MAF2 "Chromosome\tStart_Position\tReference_Allele\tTumor_Allele\tRef_allele_depth\tAlt_allele_depth\tmutation_id\n";
while (<MAF1>) {
	chomp;
	my @lines = split/\t/;
	next if $lines[6-1] ne "SNP";
	my $info = (join "\t",@lines[1,2,7,8,9,10])."\t$lines[0]:$lines[1]:$lines[2]";	
	print MAF2 "$info\n";
}
close MAF1;
close MAF2;



my $patientID = $ARGV[1];
open MAF1,shift;
open MAF2,">${sample}.snv.txt";
print MAF2 "Chromosome\tStart_Position\tReference_Allele\tTumor_Allele\tRef_allele_depth\tAlt_allele_depth\tmutation_id\n";
while (<MAF1>) {
	chomp;
	my @lines = split/\t/;
	next if $lines[6-1] ne "SNP";
	my $info = (join "\t",@lines[1,2,7,8,9,10])."\t$lines[0]:$lines[1]:$lines[2]";	
	print MAF2 "$info\n";
}
close MAF1;
close MAF2;