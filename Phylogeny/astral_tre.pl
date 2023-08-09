#!/usr/bin/env perl
use Bio::Seq;
use Bio::SeqIO;
use strict;
use warnings;
defined($ARGV[0]) or die "Usage: $0 Speicies_id Sequences_id folders_containing_orthogrps and minmum_taxa\n";
chomp(@ARGV);
open Spe,"<$ARGV[0]";
open SPE, ">SpeTable.txt";
my %speID;

while(<Spe>)
{
chomp();
my @info=split(/\s/,);
$info[0] =~ s/://;
$speID{$info[0]} = $info[1];
}
close Spe;

open Seq, "<$ARGV[1]";

my %seqID;

while(<Seq>)
{
chomp();
my @info=split(/\s+/,);
$info[0] =~ s/://;
$info[1] =~ s/\|/_/g;
$info[0] =~ s/\s//g;
$info[1] =~ s/\s//g;
$seqID{$info[1]} = $info[0];

}
my %speTable;

open ALN, ">aln.sh";
open TRIM, ">trimm.sh";
open RAXML, ">raxml.sh";
my @ortho=`ls $ARGV[2] | grep fa`;
my $min = $ARGV[3];
foreach my $id (@ortho)
{
my %length;
chomp($id);
my %taxa;
my $seqIO = Bio::SeqIO->new(-file => "$ARGV[2]/$id" ,
                         -format => 'fasta');



while ( my $seq = $seqIO->next_seq() ) {
    # the output handle is reset for every file
    my $Sid=$seq->primary_id();
    my $len=$seq->length();
    ($len<40) and next;
    $Sid=~s/\|/_/g;
    my $transID = $seqID{$Sid};
    $transID =~ s/(\d+?)_.*/$1/;
    my $spe=$speID{$transID};
    $taxa{$spe} = 1;

    $length{$Sid}=$len;
}
my @len=sort values(%length);
my $count = int(scalar(@len) / 2);
my $len_mean = $len[$count];
#print "$len_mean\n";
my %taxa2;
my %selected;
for my $Sid (keys %length)
{

($length{$Sid} < ($len_mean * 0.6)) and next;
$selected{$Sid}=1;
my $transID = $seqID{$Sid};
$transID =~ s/(\d+?)_.*/$1/;
my $spe=$speID{$transID};
$taxa2{$spe} = 1;
}

if (scalar(keys(%taxa2))>=$min)
{
my $seqIO2= Bio::SeqIO->new(-file => "$ARGV[2]/$id" ,
                         -format => 'fasta');

my $out=`basename $id`;
chomp($out);
my $seqOUT = Bio::SeqIO->new(-file => ">$out" ,
                         -format => 'fasta');

while ( my $seq = $seqIO2->next_seq() ) {
    # the output handle is reset for every file
    my $Sid=$seq->primary_id();
    $Sid=~s/\|/_/g;
    #print "Hi1:$Sid\t$out\n";
    (exists($selected{$Sid})) or next;
    #print "Hi2\n";
    my $transID = $seqID{$Sid};
    $seq->id($transID);
    $seqOUT->write_seq($seq);
    my @split=split(/_/,$transID);
    push @{$speTable{$split[0]}}, $transID;
}
print ALN "muscle -align $id -output $id.aln.fa\n";
print TRIM "trimal -automated1 -in $id.aln.fa -out $id.trimmed.fa -colnumbering >$id.aln.columns\n";
my $random=int(1000*(rand(1000)+0.001));
print RAXML "raxmlHPC-PTHREADS-AVX2 -s $id.trimmed.fa -n $id -m PROTCATGTR -e 0.1 -f d -p $random -T 2\n";
#print RAXML "FastTree -wag -out $id.fasttree.tre $id.trimmed.fa\n";
}
}

foreach my $spe (sort keys %speTable)
{
foreach my $gene (sort @{$speTable{$spe}})
{
print SPE "$gene\t$spe\n";
}

}
print "All files prepared\n";
system "ParaFly -c aln.sh -CPU 100";
system "ParaFly -c trimm.sh -CPU 100";
system "ParaFly -c raxml.sh -CPU 100";
system "cat RAxML_bestTree.* >geneTre.tre";
#system "cat *.fasttree.tre >geneTre.tre";
system "java -D\"java.library.path=/nfs_genome/BioInformatics/A-pro/ASTRAL-MP/lib\" -jar /nfs_genome/BioInformatics/A-pro/ASTRAL-MP/astral.1.1.6.jar -i geneTre.tre -a SpeTable.txt -o Species.tre";
system "rename_tree.pl Species.tre $ARGV[0]";
