#!usr/bin/perl


my ($filein, $fileout) = @ARGV;

open(DATA, "$filein");

open(OUTPUT, ">$fileout.wig");

my $chr="chr";
my $start=-1;
my $sig=0;
my $count=1;
my $startC=0;
my $posn=1;

while(<DATA>)
{
	my @line=split/\s+/, $_;
	my $temp_chr=$line[0];
	my $temp_start=$line[1];
	my $temp_sig=$line[3];
	my $posn=$line[1];
	if ($count < 2){
		print OUTPUT "track type=wiggle_0\nvariableStep chrom=$temp_chr span=1\n";
		$count++;
		print "$count\n";
	}

if (($temp_chr !~ /liftover/) && ($temp_chr !~ /#/) && ($temp_start !~ /e/)){
 
	if ($posn<$startC) {
		print OUTPUT "variableStep chrom=$temp_chr span=1\n";
		$startC=$posn;
		}
		
		if ($posn>$startC) {	print OUTPUT "$temp_start\t$temp_sig\n";
		$sig=$temp_sig;
		$startC=$posn;
	}
	}
}
close(OUTPUT);
close(DATA);