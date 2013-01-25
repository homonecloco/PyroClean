#!/bin/perl -w
use Bio::SeqIO;
use File::Basename;


my($filename, $directories, $suffix) = fileparse($ARGV[0]);

my %SEQS;
my $in  = Bio::SeqIO->new(-file => $ARGV[0],
                       -format => 'Fasta');



my $fh;
while ( my $seq = $in->next_seq() ){
  my $s = $seq->seq;
  
     
  unless ($SEQS{$s}){
    #$fh = *FILEHANDLE; 
      $SEQS{$s} = 0;
    }
    $SEQS{$s}++;
    #$fh = $IDS{$short};
    
    
}
  #print $name;

#my $file = "U_".$filename;
#open FILEHANDLE, ">$file" or die $!; 
#open FILEHANDLE ">$file" or die $!;

my $count =0;
foreach my $key ( keys %SEQS )
{
  $count++;
#print "key: $key, value: $IDS{$key}\n";

  print   ">Seq".$count."_".$SEQS{$key}."\n".$key."\n";
#close $IDS{$key};
}


#close(FILEHANDLE);

#$out = Bio::SeqIO->new(-file => ">outputfilename",
#                       -format => 'EMBL');