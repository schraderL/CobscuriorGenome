use English     qw<$RS>;
use File::Slurp qw<write_file>;
my $id     = 0;
my $source = $ARGV[0];
my %hash;
{   local $RS = "\n\n\n";
    open my $in, $source or die "can t read $source: $!\n";
    while ( <$in> ) {
        $_ =~ /((scaffold[0-9]+) ([0-9]+) ([0-9]+) .*?(m_b.*))\n/;
        my @line = split(' ', $1);
        #my @pos = split(' ', $2);
        chomp; # removes the line "\n/DOCUMENTS/\n"
        #write_file( 'file' . ( ++$id ) . '.txt', $_ );
        print $2,"\t",$3,"\t",$4,"\t",$5,"\n";
        #$hash{$5} = $_;
    }
    # being scoped by the surrounding brackets (my "local block"),
    close $in;    # an explicit close is not necessary
}
