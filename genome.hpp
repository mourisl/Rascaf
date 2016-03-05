// The class handles the genome file in fasta file

#ifndef _LSONG_RSCAF_GENOME_HEADER
#define _LSONG_RSCAF_GENOME_HEADER

#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <stdio.h>

#include <inttypes.h>

#include "alignments.hpp"
#include "KmerCode.hpp" 
#include "defs.h"

extern char nucToNum[26] ;
extern char numToNuc[26] ;
extern bool VERBOSE ;
extern FILE *fpOut ;

struct _contig
{
	int64_t start, end ;
	int chrId ; // scaffold id
	int id ;
} ;

// Use bit to store the sequence
// The non-ACGT character's bits are non-defined.
class BitSequence 
{
private:
	int len ;
	//std::vector<uint32_t> sequence ;
	uint32_t *sequence ;

public:
	BitSequence() { len = 0 ; sequence = NULL ;} 
	BitSequence( int l )
	{
		len = 0 ;
		sequence = new uint32_t[ l / 16 + 1 ] ;
	}
	
	~BitSequence() 
	{
		//if ( sequence != NULL )
		//	delete[] sequence ;
	}

	int GetLength()
	{
		return len ;
	}

	void Append( char c )
	{
		if ( ( len & 15 ) == 0 )
		{
			sequence[ len / 16 ] = 0 ;
		}
		++len ;
		//printf( "%d %c\n", len, c ) ;
		Set( c, len - 1 ) ;
	}
	
	// pos is 0-based coordinate
	// notice that the order within one 32 bit butcket is reversed
	void Set( char c, int pos ) 
	{
		if ( pos >= len )		
			return ;

		if ( c >= 'a' && c <= 'z' )
		{
			c = c - 'a' + 'A' ;
		}
		if ( c == 'N' )
			c = 'A' ;

		int ind = pos >> 4 ;
		int offset = pos & 15 ;
		int mask = ( (int)( nucToNum[c - 'A'] & 3 ) ) << ( 2 * offset ) ;
		sequence[ind] = sequence[ind] | mask ;
		//if ( c != 'A' )
		//	printf( "%d: %c %c %d %d : %d\n", pos, c, Get(pos), ind, offset, mask ) ;
		//Print() ;
	}

	char Get( int pos )
	{
		if ( pos >= len )
			return 'N' ;

		int ind = pos >> 4 ;
		int offset = pos & 15 ;
		//printf( "%d: %d\n", pos, sequence[ind] ) ;
		return numToNuc[ ( ( sequence[ind] >> ( 2 * offset ) ) & 3 ) ] ;
	}

	void Print()
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
			std::cout<<Get( i ) ;
		std::cout<<"\n" ;
	}

	void Print( FILE *fp, int start, int end, bool rc )
	{	
		if ( !rc )
		{
			for ( int i = start ; i <= end ; ++i )
				fprintf( fp, "%c", Get( i ) ) ;
		}
		else
		{
			for ( int i = end ; i >= start ; --i )
			{
				char c = Get( i ) ;
				if ( c == 'A' )
					c = 'T' ;
				else if ( c == 'C' )
					c = 'G' ; 
				else if ( c == 'G' )
					c = 'C' ; 
				else //if ( c == 'T' )
					c = 'A' ; 
				fprintf( fp, "%c", c ) ;
			} 
		}
	}
} ;


class Genome
{
private:
	std::vector<BitSequence> genomes ;
	bool isOpen ;
	std::vector<struct _contig> contigs ;
	std::vector<struct _pair> contigRanges ;

public:
	Genome() { isOpen = false ;}
	~Genome() {}

	void Open( Alignments &alignments, char *fa )
	{
		std::ifstream fp ;	
		fp.open( fa ) ;
		std::string line ;
		int chrId = -1 ;
		std::string s( "" ) ;

		struct _contig tmpContig ;
		struct _pair tmpContigRange ;
		int64_t offset = 0 ;

		tmpContig.start = -1 ;
		tmpContigRange.a = 0 ;
		while ( getline( fp, line ) )
		{
			/*std::cout<<line<<"\n" ;
			
			fflush( stdout ) ;*/
			if ( line[0] == '>' )
			{
				//printf( "%d\n", alignments.GetChromIdFromName( "chr20_1" ) ) ;
				if ( tmpContig.start != -1 )
				{
					tmpContig.end = offset - 1 ;
					tmpContig.id = contigs.size() ;
					contigs.push_back( tmpContig ) ;
					
					tmpContig.start = -1 ;
				}

				if ( chrId != -1 )
				{
					tmpContigRange.b = contigs.size() - 1 ;
					//printf( "%d %s (%d %d)\n", chrId, alignments.GetChromName( chrId ), (int)tmpContigRange.a, (int)tmpContigRange.b ) ;

					if ( (int)contigRanges.size() <= chrId )
					{
						while ( (int)contigRanges.size() <= chrId )
							contigRanges.push_back( tmpContigRange ) ;
					}
					else
					{
						contigRanges[ chrId ] = tmpContigRange ;
					}
					tmpContigRange.a = tmpContigRange.b + 1 ;
				}
				
				char *s = strdup( line.c_str() + 1 ) ;
				//std::cout<<line<<"\n" ;	
				for ( int i = 0 ; s[i] ; ++i )
					if ( s[i] == ' ' || s[i] == '\t' )
					{
						s[i] = '\0' ;
						break ;
					}
				chrId = alignments.GetChromIdFromName( s ) ;
				
				if ( (int)genomes.size() <= chrId )
				{
					while ( (int)genomes.size() < chrId )
					{
						BitSequence bs( 0 ) ;
						genomes.push_back( bs ) ;
					}
					//printf( "%d %s %d\n", chrId, alignments.GetChromName( chrId ), alignments.GetChromLength( chrId ) ) ;

					BitSequence bs( alignments.GetChromLength( chrId ) );
					//printf( "%s %d %d\n", s, chrId, alignments.GetChromLength( chrId ) ) ;
					genomes.push_back( bs ) ;
				}
				else
				{
					BitSequence bs( alignments.GetChromLength( chrId ) );
					//printf( "%s %d %d\n", s, chrId, alignments.GetChromLength( chrId ) ) ;
					genomes[ chrId ] = bs ;
				}
				free( s ) ;

				tmpContig.chrId = chrId ;
				tmpContig.start = -1 ;
				offset = 0 ;
			}
			else
			{
				int len = line.length() ;
				int i ;
				BitSequence &bs = genomes[chrId] ;
				for ( i = 0 ; i < len ; ++i, ++offset )
				{
					if ( ( line[i] >= 'A' && line[i] <= 'Z' ) ||
						( line[i] >= 'a' && line[i] <= 'z' ) )
					{
						bs.Append( line[i] ) ;
						//if ( line[i] != 'N' )
						//	printf( "%c %c\n", line[i], bs.Get(i) ) ;

						if ( line[i] == 'n' || line[i] == 'N' )
						{
							if ( tmpContig.start != -1 )
							{
								tmpContig.end = offset - 1 ;
								tmpContig.id = contigs.size() ;
								contigs.push_back( tmpContig ) ;
								//printf( "(%d %d)\n", (int)tmpContig.start, (int)tmpContig.end ) ;

								tmpContig.start = -1 ;
							}
						}
						else
						{
							if ( tmpContig.start == -1 )
							{
								tmpContig.start = offset ;
								//printf( "==(%lld %lld) %lld %d %c\n", tmpContig.start, tmpContig.end, offset, i, line[i - 1] ) ;
							}
						}
					}
				}
				//bs.Print() ;
			}
		}

		if ( tmpContig.start != -1 )
		{
			tmpContig.end = offset - 1 ;
			tmpContig.id = contigs.size() ;
			contigs.push_back( tmpContig ) ;

			tmpContig.start = -1 ;
		}

		if ( chrId != -1 )
		{
			tmpContigRange.b = contigs.size() - 1 ;
			//printf( "%d (%d %d)\n", chrId, (int)tmpContigRange.a, (int)tmpContigRange.b ) ;

			if ( (int)contigRanges.size() <= chrId )
			{
				while ( (int)contigRanges.size() <= chrId )
					contigRanges.push_back( tmpContigRange ) ;
			}
			else
			{
				contigRanges[ chrId ] = tmpContigRange ;
			}
			tmpContigRange.a = tmpContigRange.b + 1 ;
		}

		fp.close() ;
		isOpen = true ;
		
		if ( fpOut != NULL )
		{
			fprintf( fpOut, "Contigs\n" ) ;
			for ( int i = 0 ; i < contigs.size() ; ++i )
				fprintf( fpOut, "%d: %s %"PRId64" %"PRId64"\n", i, alignments.GetChromName( contigs[i].chrId ), contigs[i].start + 1, contigs[i].end + 1 ) ;
		}
	}

	bool IsOpen()
	{
		return isOpen ;
	}

	void SetIsOpen( bool in )
	{
		isOpen = in ;
	}


	int ScaffoldToContigCoord( int chrId, int64_t pos, int &contigId, int64_t &cpos )
	{
		if ( !isOpen )
		{
			contigId = chrId ;
			cpos = pos ;
			return chrId ;
		}

		//TODO: some caching
		int l, r, m ;
		l = contigRanges[ chrId ].a ;
		r = contigRanges[ chrId ].b ;

		while ( l <= r )
		{
			m = ( l + r ) / 2 ;
			if ( pos < contigs[m].start )
				r = m - 1 ;
			else if ( pos > contigs[m].end )
				l = m + 1 ;
			else
			{
				contigId = m ;
				cpos = pos - contigs[m].start ;
				//printf( "%s %d\n", __func__, contigId ) ;
				return contigId ;
			}
		}

		//contigId = cpos = -1 ;
		contigId = l ;
		cpos = contigs[l].start ;

		return contigId ;
	}

	int GetChrIdFromContigId( int contigId )
	{
		if ( !isOpen )
			return contigId ;
		//printf( "%s: %d %d\n", __func__, contigId, contigs[ contigId ].chrId ) ;
		return contigs[ contigId ].chrId ;
	}

	int GetContigId( int chrId, int64_t pos )
	{
		if ( !isOpen )
		{
			return chrId ;
		}

		int contigId = -1 ;
		int64_t cpos = -1 ;

		ScaffoldToContigCoord( chrId, pos, contigId, cpos ) ;
		//printf( "%d %d\n", chrId, contigId ) ;
		return contigId ;
	}

	int GetContigCount()
	{
		return contigs.size() ;
	}

	int GetChrCount()
	{
		return genomes.size() ;
	}

	int GetChrLength( int chrId )
	{
		return genomes[chrId].GetLength() ;
	}
	
	int GetChrContigRange( int chrId, int &from, int &to )
	{
		if ( !isOpen )
		{
			from = chrId ;
			to = chrId ;
			return 1 ;
		}

		from = contigRanges[ chrId ].a ;
		to = contigRanges[ chrId ].b ;
		return to - from + 1 ;
	}

	// 0-based pos
	char GetNucleotide( int chrId, int pos )
	{
		if ( !isOpen )	
			return '\0' ;
		//printf( "%c\n", genomes[chrId].Get(pos) ) ;
		return genomes[chrId].Get( pos ) ;	
	}

	struct _contig GetContigInfo( int contigId )
	{
		if ( isOpen )
			return contigs[ contigId ] ;
		else
		{
			struct _contig ret ;
			ret.chrId = contigId ;
			ret.id = contigId ;
			ret.start = 0 ;
			ret.end = 0 ;//genomes[contigId].GetLength() - 1 ;

			return ret ;
		}
	}

	void PrintContig( FILE *fp, int contigId, bool reverseComplement )
	{
		if ( isOpen )
			genomes[ contigs[ contigId ].chrId ].Print( fp, contigs[contigId ].start, contigs[ contigId ].end, reverseComplement ) ;
		else
			genomes[ contigId ].Print( fp, 0, genomes[ contigId].GetLength() - 1, reverseComplement ) ;
	}

	// The function handle kmers========================================================
	void AddKmer( int chrId, int from, int to, int kl, std::map<uint64_t, int> &kmers ) 
	{
		KmerCode code( kl ) ;
		int i ;
		BitSequence &s = genomes[chrId] ;
		for ( i = from ; i < from + kl - 1 ; ++i )
		{
			//std::cout<<s.Get(i)<<"\n" ;
			code.Append( s.Get( i ) ) ;
		}

		for ( ; i <= to ; ++i )
		{
			if ( code.IsValid() )
			{
				uint64_t key = code.GetCanonicalKmerCode() ;
				if ( kmers.count( key ) > 0 )
					++kmers[key] ;
				else
					kmers[key] = 1 ;
			}
			code.Append( s.Get( i ) ) ;
		}
		
		if ( code.IsValid() )
		{
			uint64_t key = code.GetCanonicalKmerCode() ;
			if ( kmers.count( key ) > 0 )
				++kmers[key] ;
			else
				kmers[key] = 1 ;
		}
	}

	int CountStoredKmer( int chrId, int from, int to, int kl, std::map<uint64_t, int> &kmers, bool test = false )
	{
		KmerCode code( kl ) ;
		int i ;
		int ret = 0 ;
		BitSequence &s = genomes[chrId] ;
		for ( i = from ; i < from + kl - 1 ; ++i )
			code.Append( s.Get( i ) ) ;

		for ( ; i <= to ; ++i )
		{
			if ( code.IsValid() )
			{
				std::map<uint64_t, int>::iterator it = kmers.find( code.GetCanonicalKmerCode() ) ;
				if ( it != kmers.end() )
				{
					if ( test )
						printf( "%d %d %d\n", i, kmers[ code.GetCanonicalKmerCode()], it->second ) ;
					ret += it->second ;
				}
			}
			code.Append( s.Get( i ) ) ;
		}
		
		if ( code.IsValid() )
		{
			std::map<uint64_t, int>::iterator it = kmers.find( code.GetCanonicalKmerCode() ) ;
			if ( it != kmers.end() )
			{
				//if ( test )
				//	printf( "%d %d\n", i, kmers[ code.GetCanonicalKmerCode()] ) ;
				if ( test )
					printf( "%d %d %d\n", i, kmers[ code.GetCanonicalKmerCode()], it->second ) ;
				ret += it->second ;
			}
		}
		return ret ;
	}

	int CompareKmerSets( std::map<uint64_t, int> &kmersA, std::map<uint64_t, int> &kmersB )
	{
		int ret = 0 ;
		std::map<uint64_t, int>::iterator itA = kmersA.begin() ; 
		std::map<uint64_t, int>::iterator itB ;
		for ( ; itA != kmersA.end() ; ++itA )
		{
			int a = itA->second ;
			itB = kmersB.find( itA->first ) ;
			if ( itB == kmersB.end() )
				continue ;
			int b = itB->second ;
			ret += ( a < b ? b : a ) ; // use the larger one
		}
		return ret ;
	}
} ;

#endif
