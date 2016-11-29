#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: perl FilterRascafConnection.pl evaluate_blast.out rascaf_cs.fa rascaf.out > rascaf_filter.out\n" if ( @ARGV == 0 ) ;

# Store the evaluation result
open FP1, $ARGV[0] ;
my %csValidation ;
while ( <FP1> )
{
	chomp ;
	my @cols = split ;
	$csValidation{ $cols[0] } = $cols[1] ;
}
close FP1 ;

# Extract the blocks from rascaf_cs.fa file
open FP1, $ARGV[1] ;
my %blocksLength ;
my %blocksValidation ;
while ( <FP1> )
{
	next if ( !/^>/ ) ;
	chomp ;
	my @cols = split ;
	#print substr( $cols[0], 1), "\n" ;
	my $validation = $csValidation{ substr( $cols[0], 1 ) } ;

	$blocksLength{ $cols[1] } = $cols[3] ;
	$blocksValidation{ $cols[1]} = $validation ;
}
close FP1 ;

# Filter out the blocks that can not validated. Only keep the blocks with validation==1.
open FP1, $ARGV[2] ;
my $start = 0 ;
while ( <FP1> )
{
	chomp ;
	if ( /WARNING/ ) 
	{
		$start = 0 ;
		print $_, "\n" ;
		next ;
	}

	if ( /command/ )
	{
		$start = 1 ;
		print $_, "\n" ;
		next ;
	}
	elsif ( $start == 0 ) 
	{
		print $_, "\n" ;
		next ;
	}


	if ( !/^[0-9]+:/ ) 
	{
		print $_, "\n" ;
		next ;
	}
	my @cols = split ;

	my $line = $_ ;
	my $cnt ;
	my $i ;
	my $j ;
	my @blocksId ;
	my @connections ;

	if ( $line =~ /^([0-9]+):/ )
	{
		#print "### $1\n" ;
		$cnt = $1 ;
	}
	for ( $i = 0 ; $i < $cnt ; ++$i )
	{
		my $id = $cols[$i*4+1]." ".$cols[$i*4+2]." ".$cols[$i*4+3]." ".$cols[$i*4+4] ;
		push @blocksId, $id ;
	}

	for ( $i = 0 ; $i < $cnt - 1 ; ++$i )
	{
		my $l = <FP1> ;
		chomp $l ;

		push @connections, $l ;
	}

	my @keepConnection ; # whether we want to keep connection between block $i and $i+1.
	for ( $i = 0 ; $i < $cnt - 1 ; ++$i )
	{
		push @keepConnection, 1 ;
	}
	push @keepConnection, 0 ; # the last contig has not extension.

	for ( $i = 0 ; $i < $cnt - 1 ; ++$i )
	{
		my @cols1 = split /\s+/,  $connections[$i] ;
		#print STDERR $connections[$i], " ", $i, " ",  @cols1[2], "\n" if ( 1 ) ;
		next if ( !defined( $blocksValidation{ $cols1[2] } ) ) ;
		my $validation = $blocksValidation{ $cols1[2] } ;
		my $len = $blocksLength{ $cols1[2] } ;	

		for ( $j = 0 ; $j < $len - 1 ; ++$j )
		{
			if ( $validation != 1 )
			{
				$keepConnection[$i + $j] = 0 ;
			}
		}
		$i = $i + $j ;
	}

	for ( $i = 0 ; $i < $cnt ; ++$i ) 
	{
		for ( $j = $i ; $j < $cnt ; ++$j )
		{
			last if ( $keepConnection[$j] == 0 ) ;
		}

		next if ( $j == $i ) ; # a single contig after break the connection, no need to output it now.

		print $j - $i + 1, ": " ;
		my $k ;
		for ( $k = $i ; $k <= $j ; ++$k )
		{
			print $blocksId[$k]." " ;
		}
		print "\n" ;

		for ( $k = $i ; $k < $j ; ++$k )
		{
			print $connections[$k], "\n" ;
		}

		for ( $i = $j + 1 ; $i < $cnt ; ++$i )
		{
			last if ( $keepConnection[$i] != 0 ) ;
		}
		--$i ;
	}
}
close FP1 ;

