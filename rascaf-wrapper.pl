#!/bin/perl

# Wrapper for Rascaf. 

use strict ;
use warnings ;

use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename;

my @rascafARGV ;
my @rascafJoinARGV ;
my @bamFiles ;
my @bamFilePrefix ;
my $WD = dirname( abs_path( $0 ) ) ;

my $usage = "perl rascaf-wrapper.pl [OPTIONS]:\n".
	    "\t-b STRING: bam file. Can use multiple -b to specify multiple alignment files[required]\n".
	    "\t-f STRING: path to the raw assembly fasta file\n".
	    "\t-o STRING : prefix of the output file (default: rascaf_scaffold)".
	    "\t-ms INT: minimum support for connecting two contigs(default: 2)".
	    "\t-ml INT: minimum exonic length if no intron (default: 200)".
	    "\t-k INT: the size of a kmer(<=32. default: 21)" ;

my %oneParaOptions ;
$oneParaOptions{ "-f" } = 1 ;
$oneParaOptions{ "-ms" } = 1 ;
$oneParaOptions{ "-ml" } = 1 ;
$oneParaOptions{ "-k" } = 1 ;

for ( my $i = 0 ; $i < scalar( @ARGV ) ; ++$i )
{
	if ( $ARGV[ $i ] eq "-b" )
	{
		my $tmp = $ARGV[$i + 1] ;
		push @bamFiles, $tmp ;
		
		$tmp = basename( $tmp ) ;
		$tmp =~ s/\.bam$// ;
		push @bamFilePrefix, $tmp ;

		++$i ;

	}
	elsif ( $ARGV[ $i ] eq "-o" )
	{
		push @rascafARGV, $ARGV[$i] ;
		push @rascafARGV, $ARGV[$i + 1] ;

		push @rascafJoinARGV, $ARGV[$i] ;
		push @rascafJoinARGV, $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( defined $oneParaOptions{ $ARGV[ $i ] } )
	{
		push @rascafARGV, $ARGV[$i] ;
		push @rascafARGV, $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "Unknown option ", $ARGV[$i], "\n" ;
	}
}
if ( scalar( @bamFiles ) <= 0 )
{
	die "Must speicfy bam files through -b option.\n" ;
}


my $joinList = "" ;
for ( my $i = 0 ; $i < scalar( @bamFiles ) ; ++$i )
{
	my $cmd = "$WD/rascaf -b ".$bamFiles[$i]." -o ".$bamFilePrefix[$i]."_$i @rascafARGV" ;
	print STDERR $cmd, "\n" ;
	die "rascaf failed.\n" if ( system( $cmd ) != 0 ) ;
	$joinList .= "-r ".$bamFilePrefix[$i]."_$i".".out" ;
}

my $cmd = "$WD/rascaf-join $joinList @rascafJoinARGV" ;
print STDERR $cmd, "\n" ;
die "rascaf-join failed.\n" if ( system( $cmd ) != 0 ) ;

