use strict;
open(CONF,"FitnessShift.Config.txt") or die $!;
my $in_line = <CONF>;
$in_line =~ m/\=/;
my $in_file = "$'";
open (IN,$in_file) or die $!; 

my $n=1;

my $idd = <IN>;
my $seqq = <IN>;
chomp($seqq);
my $length = length($seqq);

seek(IN,0,0);

while ($n <= $length) {
   
    open (OUT, ">>$n.site") or die $!;
    
    seek (IN,0,0);
    
    while (<IN>)  {
	
	my $id = $_;
	$id =~ s/\>//;
	my $seq = <IN>;	

	print OUT $id;
	
	my $letter = substr ($seq,$n-1,1);
	print OUT $letter."\n";
    }
    $n = $n+1;
}


    
