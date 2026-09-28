use Getopt::Std;
getopt('fm');
my $min_length=0;
if($opt_m){$min_length=$opt_m;}
my $fasta;

if($opt_f){$fasta=$opt_f;}else{$fasta = $ARGV[0];}
open my $in,"<", $fasta or die "can't open file";
my $seqname;
while(<$in>){chomp;

if(/^>/){$seqname = $_;}

else{


if(length($_)>$min_length){

my $seq = $_;
my $firstN = substr($seq,0,1);
my $seqlength=length($seq);
$seqname =~ s/ +/:/;
$seqname = join("_",$seqname, $firstN, $seqlength);

print $seqname, "\n";
                            print $_, "\n";}

}

}



