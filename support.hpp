// The class that deals the number of alignments supporting blocks or edges
#ifndef _LSONG_RSCAF_SUPPORT_HEADER
#define _LSONG_RSCAF_SUPPORT_HEADER

extern int minimumSupport ; 

class Support
{
private:
	int uniqSupport ; 
	int multiSupport ;


	int clipSupport ;
	int nmSum ;
	double multiSupportCoefficient ;

	int64_t leftPos, rightPos ;
	int coordCnt ; // Record how many coordinates showed up.
	int64_t prevCoord ;
		
	bool IsSignificantDifferent( double a, double b )
	{
		if ( a - 6 * sqrt((double)a) <= b  
				&& b <= a + 6 * sqrt( (double)a ) )
		{
			return false ;
		}

		if ( b - 6 * sqrt( (double)b ) <= a && a <= b + 6 * sqrt( (double)b ) )
		{
			return false ;
		}

		return true ;
	}
public:
	double plusSupport ;
	double minusSupport ;
	Support()
	{
		uniqSupport = multiSupport = 0 ;
		plusSupport = minusSupport = 0 ;
		clipSupport = 0 ;
		leftPos = rightPos = -1 ;
		multiSupportCoefficient = 1.0 ;
		nmSum = 0 ;
		coordCnt = 0 ;
		prevCoord = -1 ;
	}

	~Support()
	{
	} 

	void Add( Alignments &align, bool ignoreCoord = false )
	{
		if ( align.IsUnique() )
			++uniqSupport ;
		else
			++multiSupport ;

		int strand = align.GetStrand() ;
		if ( strand == 1 )
			plusSupport += align.GetStrandWeight() ;
		else if ( strand == -1 )
			minusSupport += align.GetStrandWeight() ;
		int nm = align.GetFieldI( "NM" ) ;
		if ( nm >= 0 )
			nmSum += nm ;

		if ( !ignoreCoord )
		{
			if ( leftPos == -1 || align.segments[0].a < leftPos )
				leftPos = align.segments[0].a ;
			if ( align.segments[ align.segCnt - 1 ].b > rightPos )
				rightPos = align.segments[ align.segCnt - 1 ].b ;

			if ( align.segments[0].a != prevCoord ) // This makes sense when the coordinates are sorted
				++coordCnt ;
			else if ( align.GetFieldZ( "SA" ) != NULL )
				++coordCnt ;
			prevCoord = align.segments[0].a ;
		}
	}

	void Add( Support &in )
	{
		uniqSupport += in.uniqSupport ;
		multiSupport += in.multiSupport ;

		plusSupport += in.plusSupport ;
		minusSupport += in.minusSupport ;

		nmSum += in.nmSum ;

		clipSupport += in.clipSupport ;

		if ( in.leftPos != -1 && ( leftPos == -1 || in.leftPos < leftPos ) )
			leftPos = in.leftPos ;
		if ( in.rightPos > rightPos )
			rightPos = in.rightPos ;

		prevCoord = prevCoord > in.prevCoord ? prevCoord : in.prevCoord ;
		coordCnt += in.coordCnt ;
	}

	bool IsGood()
	{
		if ( uniqSupport < 0.1 * ( uniqSupport + multiSupport ) || uniqSupport == 0 )
			return false ;
		if ( nmSum / ( uniqSupport + multiSupport ) >= 2.5 ) //&& !IsUnique() )
			return false ;
		//else if ( nmSum / ( uniqSupport + multiSupport ) >= 3.5  )
		//	return false ;
		return true ;
	}

	bool IsUnique()
	{
		if ( uniqSupport < 0.95 * ( uniqSupport + multiSupport ) )
			return false ;
		return true ;
	}

	int GetCount()
	{
		return (int)( uniqSupport + multiSupportCoefficient * multiSupport + 1e-6) ;
	}

	int GetUniqCount()
	{
		return uniqSupport ;
	}

	int GetStrand()
	{
		if ( !IsSignificantDifferent( plusSupport, minusSupport ) && plusSupport > 1 && minusSupport > 1 )
			return 0 ;
		else if ( plusSupport <= 1 && minusSupport <= 1 )
			return 0 ;
		else if ( plusSupport == minusSupport )
			return 0 ;
		else if ( plusSupport > minusSupport )
			return 1 ;
		else 
			return -1 ;
	}

	bool HasStrandSupport()
	{
		if ( plusSupport >= minimumSupport || minusSupport >= minimumSupport )
			return true ;
		return false ;
	}

	int GetLeftMostPos()
	{
		return leftPos ;
	}

	int GetRightMostPos()
	{
		return rightPos ;
	}

	int GetCoordCnt()
	{
		return coordCnt ;
	}

	void operator=( const Support &in )
	{
		uniqSupport = in.uniqSupport ;
		multiSupport = in.multiSupport ;

		plusSupport = in.plusSupport ;
		minusSupport = in.minusSupport ;
		
		clipSupport += in.clipSupport ;
		
		nmSum = in.nmSum ;
		
		leftPos = in.leftPos ;
		rightPos = in.rightPos ;

		multiSupportCoefficient = in.multiSupportCoefficient ;
		
		prevCoord = in.prevCoord ;
		coordCnt = in.coordCnt ;
	}

	void Clear()
	{
		uniqSupport = multiSupport = 0 ;
		plusSupport = minusSupport = 0 ;
		clipSupport = 0 ;
		nmSum = 0 ;
		leftPos = rightPos = -1 ;
		multiSupportCoefficient = 1.0 ;

		coordCnt = 0 ;
		prevCoord = -1 ;
	}
} ;
#endif
