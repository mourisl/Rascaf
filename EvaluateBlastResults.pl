#!/bin/perl

# evaluate the blast output for the connection sequences

use strict ;

die "usage: a.pl blast.out rascaf_cs.out [-diffScaf -minIdt 0]\n" if ( @ARGV == 0 ) ;

my %split ;
my %csIds ;

my $onlyUseDifferentScaffoldConnection = 0 ;
my $minIdentity = 0 ;
for ( my $i = 2 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "-diffScaf" )
	{
		$onlyUseDifferentScaffoldConnection = 1 ;
	}
	elsif ( $ARGV[$i] eq "-minIdt" )
	{
		$minIdentity = $ARGV[$i + 1] ;
		++$i ;
	}
	else
	{
		die "wrong option.\n" ;
	}
}

open FP1, $ARGV[1] ;
while ( <FP1> )
{
	next if (!/^>/ ) ;
	chomp ;
	#my $line = $_ ;
	my @cols = split ;
	my $csId = substr( $cols[0], 1 ) ;

	next if ( $onlyUseDifferentScaffoldConnection == 1 && $cols[2] == 0 ) ;
	my $sum = 1 ; # 1-based index?
	$csIds{ $csId } = 1 ;
	for ( my $i = 4 ; $i < scalar( @cols ) ; ++$i )
	{
		push @{ $split{ $csId } }, $sum ;
		$sum += $cols[$i] ;
	} 
	push @{ $split{ $csId } }, $sum ;

}
close FP1 ;

open FP1, $ARGV[0] ;
my @prevCols ;

my @blastLines ;
my %csIdToBlastId ;
my %blastIdToCsId ;
my %blastIdResults ; # 0:uncertain, 1:successful, -1:fail due to orientation, -2:fail due to wrong order
my %csIdResults ; # summary of results from blastIdResults for each connected sequence
my %blastIdSplits ; # record what are the splits associated with the blast hit.

# Obtain the list of splits (gene blocks) overlapped with the range
sub GetSplitList
{
	my $csId = $_[0] ;
	my $start = $_[1] ;
	my $end = $_[2] ;
	my @ret ;

	my @splitList = @{ $split{ $csId } } ;
	for ( my $i = 0 ; $i < scalar( @splitList ) - 1 ; ++$i )
	{
		my $s = $splitList[$i] ;
		my $e = $splitList[$i + 1] ;

		next if ( $s >= $end - 20 || $e <= $start + 20 ) ; # must have 20bp overlap
		push @ret, $i ;
	}
	return @ret ;
}

sub GetOverlapSize
{
	my $s1 = $_[0] ;
	my $e1 = $_[1] ;
	my $s2 = $_[2] ;
	my $e2 = $_[3] ;

	my $ret = 0 ;
}

sub ProcessBlastHit
{
	my @blastLines = @{ $_[0] } ;
	@blastLines = sort { (split /\s+/, $a )[6] <=> ( split /\s+/, $b)[6] } @blastLines ;
	#print @blastLines, "\n" ;

	my @plusCount ;
	my @minusCount ;
	my @hitOrder ;

	my $csId = (split /\s+/, $blastLines[0] )[0] ;
	if ( !(defined $csIds{ $csId } ) )
	{
		return -5 ;
	}
	my $blastId = $csId."_".(split /\s+/, $blastLines[0] )[1] ;
	
	my $splitCnt = scalar( @{ $split{ $csId } } ) - 1 ;

	$blastIdToCsId{ $blastId } = $csId ;
	push @{ $csIdToBlastId{ $csId } }, $blastId ;

	for ( my $i = 0 ; $i < $splitCnt ; ++$i )
	{
		$plusCount[$i] = $minusCount[$i] = 0 ;
	}

	my $flagOrient = 0 ;
	my $flagOrder = 0 ;
	my $prevOrient = 0 ;
	my @prevSplitList ;

	# compute the length of hit on each orientation
	for ( my $i = 0 ; $i < @blastLines ; ++$i )
	{	
		my $localOrient = 0 ;
		my @cols = split /\s+/, $blastLines[$i] ;
		if ( ( $cols[6]<=>$cols[7] ) == ( $cols[8]<=>$cols[9] ) )
		{
			$localOrient = 1 ;
		}
		else
		{
			$localOrient = -1 ;
		}

		my @splitList = GetSplitList( $csId, $cols[6], $cols[7] ) ;
		if ( $prevOrient != 0 && $prevOrient != $localOrient 
			&& $prevSplitList[-1] != $splitList[0] )
		{
			$flagOrient = 1 ;
		}
		$prevOrient = $localOrient ;
		@prevSplitList = @splitList ;
	}

	if ( $flagOrient == 1 )
	{
		return -1 ;
	}

	# Then test the positions for each hit
	for ( my $i = 0 ; $i < scalar( @blastLines ) - 1 ; ++$i )
	{
		my @cols1 = split /\s+/, $blastLines[$i] ;
		my @cols2 = split /\s+/, $blastLines[$i + 1] ;
		my $tmpFlag = 0 ;
		
		my $prevOrient = -1 ; 
		if ( ( $cols1[6]<=>$cols1[7] ) == ( $cols1[8]<=>$cols1[9] ) )
		{
			$prevOrient = 1 ;
		}
		my $nextOrient = -1 ;
		if ( ( $cols2[6]<=>$cols2[7] ) == ( $cols2[8]<=>$cols2[9] ) )
		{
			$nextOrient = 1 ;
		}

		if ( $cols1[7] >= $cols2[7] )
		{
			# the next blast hit is a subset of current hit.
			next ;
		}

		# if the orientation changes, we ignore their relative order.
		next if ( $prevOrient != $nextOrient ) ;

		if ( $nextOrient == 1 && $cols1[8] >= $cols2[8] )
		{
			$tmpFlag = 1 ;
		}
		elsif ( $nextOrient == -1 && $cols1[8] <= $cols2[8] )
		{
			$tmpFlag = 1 ;
		}

		# if this happens within a contig, then we ignore it
		if ( $tmpFlag )
		{
			if ( scalar( GetSplitList( $csId, $cols1[6], $cols2[7] ) ) == 1 )
			{
				$tmpFlag = 0 ;
			}
		}

		if ( $tmpFlag == 1 )
		{
			$flagOrder = 1 ;
		}
	}

	if ( $flagOrder == 1 )
	{
		return -2 ;
	}

	# Obtain what are the splits (gene blocks) get involved
	my @contained ;
	for ( my $i = 0 ; $i < $splitCnt ; ++$i )
	{
		$contained[ $i ] = 0 ;
	}
	for ( my $i = 0 ; $i < scalar( @blastLines ) ; ++$i )
	{
		my @cols = split /\s+/, $blastLines[$i] ;
		my @splitList = GetSplitList( $csId, $cols[6], $cols[7] ) ;    	
		for ( my $j = 0 ; $j < scalar( @splitList ) ; ++$j )
		{
			#print "$i $j\n" ;
			$contained[ $splitList[$j] ] = 1 ;
		}
	}
	my $cnt = 0 ;
	for ( my $i = 0 ; $i < $splitCnt ; ++$i )
	{
		++$cnt if ( $contained[$i] == 1 ) ;
	}
	#print "$blastId $cnt $splitCnt\n" ;
	if ( $cnt == $splitCnt )
	{
		return 1 ;
	}
	
	#@{ $blastIdSplits{ $blastId } } = GetSplitList( $csId, $cols[6], $cols[7] ) ;
	for ( my $i = 0 ; $i < $splitCnt ; ++$i )
	{
		if ( $contained[$i] == 1 )
		{
			push @{ $blastIdSplits{ $blastId } }, $i ; 
		}
	} 

	return 0 ;
}

while ( <FP1> )
{
	chomp ;
#rascaf_CS_0	gi|694415396|ref|XM_009337601.1|	96.10	2255	24	4	1	2255	29	2219	0.0	 3710
	my @cols = split ;
	my $line = $_ ;
	next if ( $cols[2] < $minIdentity ) ;
	if ( scalar( @blastLines ) > 0 && ( $cols[1] ne (split /\s+/, $blastLines[0] )[1] ) )
	{
		#print $line, "\n" ;
		# sort the lines according to the start position of queries
		my @cols = split /\s+/, $blastLines[0] ;
		my $blastId = $cols[0]."_".$cols[1] ;
		my $result = ProcessBlastHit( \@blastLines ) ;
		#print "$blastId $result\n" ;
		$blastIdResults{ $blastId } = $result ;
		undef @blastLines ;
	}
	push @blastLines, $line ;
}

if ( scalar( @blastLines ) > 0 )
{
	my @cols = split /\s+/, $blastLines[0] ;
	my $blastId = $cols[0]."_".$cols[1] ;
	my $result = ProcessBlastHit( \@blastLines ) ;
	$blastIdResults{ $blastId } = $result ;
}
close FP1 ;

# summarize the results
foreach my $key (keys %blastIdResults ) 
{
	my $csId = $blastIdToCsId{ $key } ;
	my $result = $blastIdResults{$key} ;
	#print "$csId $key $result\n" ;
	if ( !( defined $csIdResults{$csId} ) )
	{	
		$csIdResults{ $csId } = $result ;
	}
	my $recordResult = $csIdResults{ $csId } ;
	if ( $result == 1 )
	{
		$csIdResults{ $csId } = 1 ;
	}
	elsif ( $result < 0 && $recordResult == 0 )
	{
		$csIdResults{ $csId } = $result ;
	}
}

foreach my $key (keys %csIds )
{
	if ( !( defined $csIdResults{$key} ) )
	{	
		$csIdResults{ $key } = -4 ;
	}
}

# Make sure the class-0 connections are not from two separate pieces
foreach my $key ( keys %csIdResults )
{
	next if ( $csIdResults{ $key } != 0 ) ;
	my @blastIdList = @{ $csIdToBlastId{ $key } } ;
	my $splitCnt = scalar( @{ $split{ $key } } ) - 1 ;
	my @contained ;
	my $i ;
	for ( $i = 0 ; $i < $splitCnt ; ++$i )
	{
		$contained[$i] = 0 ;
	}
	for ( $i = 0 ; $i < scalar( @blastIdList ) ; ++$i )
	{
		#print "$key ", $blastIdList[$i], " ", $blastIdSplits{ $blastIdList[$i] }, " ", $blastIdResults{ $blastIdList[$i]}, " ", $csIdResults{$key}, "\n" ;
		if ( defined $blastIdSplits{ $blastIdList[$i] } )
		{
			# if might not be defined if the blast hit is too short and overlap is less 20bp on both segment
			my @splitList = @{ $blastIdSplits{ $blastIdList[$i] } } ; 
			for ( my $j = 0 ; $j < scalar( @splitList ) ; ++$j )
			{
				$contained[ $splitList[$j] ] = 1 ;
			}
		}
	}

	for ( $i = 0 ; $i < $splitCnt ; ++$i )
	{
		last if ( $contained[$i] == 0 ) ;
	}

	if ( $i >= $splitCnt )
	{
		$csIdResults{ $key } = -3 ;
	}
}

# Output the results
foreach my $key ( sort { ( split /_+/, $a )[-1]<=>(split /_+/, $b)[-1] } ( keys %csIds ) )
{
	print "$key ", $csIdResults{$key}, "\n" ;
}
