#!usr/bin/perl -w
#change results of controlfreec(.ratio) to inputfiles of GISTIC
die "perl $0 <ratio.txt> <chrNameLength.txt> (samplename)" unless @ARGV >= 1;

#open IN, shift or die $!;
open IN, $ARGV[0] or die $!;

#define log2
sub log2{
my $n = shift;
my $result = log($n)/log(2);
return $result;
}

readline IN;	#跳过首行
my $head = <IN>;
my @tmp;
while(<IN>){
	chomp;
	@tmp = split;
	if ($tmp[3] >= 0){
		last;
	}
}
my $marker = 1;
my $start = $tmp[1];
my @GIS = qw/sample 0 0 0 0 0/;
$GIS[0] = $ARGV[2];

#print "ID\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n";

while(<IN>){
	chomp;
	my @freec = split;
	next if $freec[3] < 0;		#如果medianratio<0 则舍去
	if ($freec[3] == $tmp[3]){	#如果两行的medianratio相等
		if ($freec[0] eq $tmp[0]){	#如果在同一条染色体上
		$GIS[1] = $tmp[0];
		$GIS[2] = $start;
		$marker = $marker + 1;
		@tmp = @freec;
		}else{
			open IN2, $ARGV[1] or die $!;   #如果不在同一条染色体上
			while (<IN2>){
				 chomp;
				 my @c = split; 
				 if ($c[0] eq "chr$tmp[0]"){
					$GIS[3]= $c[1];
				}
			}
			close IN2;
			$GIS[1] = $tmp[0];
			$GIS[2] = $start;
			$GIS[4] = $marker * 2;
			if ($tmp[3] == 0){
				$GIS[5] = -2;
			}else{
				$GIS[5] = &log2($tmp[3]);
			}
			print "$GIS[0]\t$GIS[1]\t$GIS[2]\t$GIS[3]\t$GIS[4]\t$GIS[5]\n";
			$start = $freec[1];
			$marker = 1;
			@tmp = @freec;
		}
	}else{
		if ($freec[0] eq $tmp[0]){	#如果在同一条染色体上
			$GIS[3] = $freec[1]-1;
		}else{				#如果不在同一条染色体上
			open IN2, $ARGV[1] or die $!;
			while (<IN2>){
				chomp;
				my @c = split;
				if ($c[0] eq "chr$tmp[0]"){
					$GIS[3]= $c[1];
				}
			}
			close IN2;
		}
		$GIS[1] = $tmp[0];
		$GIS[2] = $start;
		$GIS[4] = $marker * 2;
		if ($tmp[3] == 0){
			$GIS[5] = -2;
		}else{
			$GIS[5] = &log2($tmp[3]);
		}
		print "$GIS[0]\t$GIS[1]\t$GIS[2]\t$GIS[3]\t$GIS[4]\t$GIS[5]\n";
		$start = $freec[1];
		$marker = 1;
		@tmp = @freec;
	}
}
open IN, $ARGV[1] or die $!;
while (<IN>){
chomp;
my @C = split /\t/;
if ($C[0] eq "chr$tmp[0]"){
	$GIS[3] = $C[1];
	}
}
$GIS[1] = $tmp[0];
$GIS[2] = $start;
$GIS[4] = $marker * 2;
if ($tmp[3]>0){
	$GIS[5] = &log2($tmp[3]);
	print "$GIS[0]\t$GIS[1]\t$GIS[2]\t$GIS[3]\t$GIS[4]\t$GIS[5]\n";
}elsif($tmp[3]==0){
	$GIS[5] = -2;
	print "$GIS[0]\t$GIS[1]\t$GIS[2]\t$GIS[3]\t$GIS[4]\t$GIS[5]\n";
}