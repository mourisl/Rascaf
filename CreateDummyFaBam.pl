#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl rascaf.out prefix\n" if ( @ARGV == 0 ) ;

my $prefix = $ARGV[1] ;
open FP1, $ARGV[0] ;
open FP2, ">$prefix.fa" ;

my @chrId ;
my @chrLength ;
while ( <FP1> )
{
	last if ( /^command/ ) ;
	next if ( !/^[0-9]/ ) ;

	my @cols = split ;
	if ( scalar( @chrId ) == 0 || $chrId[-1] ne $cols[1] )
	{
		push @chrId, $cols[1] ;
		push @chrLength, $cols[3] ;

		print FP2 ">", $chrId[-1], "\n" ;
		print FP2 "A"x($cols[3] - $cols[2] + 1), "\n" ;
	}
	else
	{
		print FP2 "N"x($cols[2] - $chrLength[-1] - 1), "\n" ;
		print FP2 "A"x($cols[3] - $cols[2] + 1), "\n" ;
		$chrLength[-1] = $cols[3] ;
	}
}
#print scalar( @chrId ), " ", scalar( @chrLength ), "\n" ;
close FP1 ;
close FP2 ;

open FP1, ">$prefix.sam" ;
#@HD	VN:1.0	SO:coordinate
#@SQ	SN:chr1	LN:249250621
print FP1 "\@HD\tVN:1.0\tSO:coordinate\n" ;
my $i ;
for ( $i = 0 ; $i < scalar( @chrId ) ; ++$i )
{
	print FP1 "\@SQ\tSN:", $chrId[$i], "\tLN:", $chrLength[$i], "\n" ;
}
close FP1 ;
