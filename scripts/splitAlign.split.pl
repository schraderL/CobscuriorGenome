use English     qw<$RS>;
use File::Slurp qw<write_file>;
my $id     = 0;
my $source = $ARGV[0]; #align file
my $source2 = $ARGV[1]; #bed intersect file
my %hash;
my %overlap;
my %nonoverlap;

#usage perl ../splitAlign3.pl Cobs.fa.align.new Cobs.fa.align.new.TEisland.bed

{   local $RS = "\n";
open my $in2, $source2 or die "can t read $source2: $!\n";
    while ( <$in2> ) {
      chomp;
      my @line = split('\t', $_);
      $overlap{$line[6]} = '';
      #print $line[6];
    }
    close $in2;    # an explicit close is not necessary
    }

#print %TEislands;
{   local $RS = "\n\n\n";
    open my $in, $source or die "can t read $source: $!\n";
    while ( <$in> ) {
        $_ =~ /((scaffold[0-9]+) ([0-9]+) ([0-9]+) .*?(m_b.*))\n/;
        my @line = split(' ', $1);
        #my @pos = split(' ', $2);
        chomp; # removes the line "\n/DOCUMENTS/\n"
        #write_file( 'file' . ( ++$id ) . '.txt', $_ );
        #print $2,"\t",$3,"\t",$4,"\t",$5,"\n";
        if (exists $overlap{$5}){
          $overlap{$5}=$_."\n";
        }else{
          $nonoverlap{$5}=$_."\n";
        }
    }
    # being scoped by the surrounding brackets (my "local block"),
    close $in;    # an explicit close is not necessary
}
open(FH1, '>', "overlap.align") or die $!;
open(FH2, '>', "nonoverlap.align") or die $!;

print FH1 values %overlap;
print FH2 values %nonoverlap;
close(FH1);
close(FH2);
