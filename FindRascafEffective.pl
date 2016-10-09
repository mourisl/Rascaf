# Find the scaffolds in rascaf_scaffold.info where the rascaf joined two different scaffolds
#!/bin/perl

use strict ;

die "usage: a.pl rascaf_scaffold.info" if ( @ARGV == 0 ) ;

open FP1, $ARGV[0] ;

while ( <FP1> )
{
	my $line = $_ ;
	my @cols = split ;
	#>scaffold_1 (scaffold2 7 +) (scaffold2 8 +)
	my $prevScafId = "--" ;
	for ( my $i = 1 ; $i < @cols ; $i += 3 )
	{
		my $scafId = substr( $cols[$i], 1 ) ;
		if ( $i != 1 && $scafId ne $prevScafId )
		{
			print $line ;
			last ;
		}

		$prevScafId = $scafId ;
	}
}

