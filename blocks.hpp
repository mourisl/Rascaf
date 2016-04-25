// The class manage the blocks
// Li Song

#ifndef _LSONG_RSCAF_BLOCKS_HEADER
#define _LSONG_RSCAF_BLOCKS_HEADER


//#define DEBUG
#define DEBUG_U 2196
#define DEBUG_V 2201

#include <stdlib.h> 
#include <vector>
#include <map>
#include <assert.h>
#include <math.h>
#include <set>
#include <inttypes.h>


#include "defs.h"
#include "support.hpp" 
#include "genome.hpp"

extern int minimumSupport ;
extern int minimumEffectiveLength ;
extern int kmerSize ;
extern bool VERBOSE ;
extern FILE *fpOut ;
extern bool aggressiveMode ;

struct _block
{
	int chrId ;
	int contigId ;
	int64_t start, end ;
	int64_t leftSplice, rightSplice ; // record the leftmost and rightmost coordinates of a splice site within this block.
	Support support ;
} ;

struct _geneBlock
{
	std::vector<int> exonBlockIds ;
	int64_t start, end ;
	int chrId ;
	int contigId ;
	Support support ;
} ;

// a plain edge
struct _edge
{
	int u, v ; // node u and v are connected
	Support support ;
} ;

// bi-directional edge
struct _mateEdge 
{
	int u, v ;  
	
	// left-bit for u, right-bit for v.
	// 1-need flip (alignment is reversed).
	Support support[4] ; 

	int supportUse ; // Indicate which support number we should use.

	bool valid ;
	bool semiValid ; // The valid information from one direction. 
			 // Only when two semiValid is false, the whole connection fails.
} ;

struct _geneBlockBubble
{
	int u ;
	int j1, j2 ; 
} ;

class Blocks
{
	private:
		std::vector<struct _block> exonBlocks ;
		std::map<int, int> exonBlocksChrIdOffset ;
		std::map<int, struct _pair> geneBlocksChrIdOffset ;

		int GetFather( int tag, int *father )
		{
			if ( father[tag] == tag )
				return tag ;
			return father[tag] = GetFather( father[tag], father ) ;
		}

		int *repeatFather ;

		// Get the exonic length of a gene block
		int GetEffectiveLength( int tag )
		{
			int cnt = geneBlocks[tag].exonBlockIds.size() ;
			int i ;
			int ret = 0 ;
			for ( i = 0 ; i < cnt ; ++i )
				ret += exonBlocks[ geneBlocks[tag].exonBlockIds[i] ].end - exonBlocks[ geneBlocks[tag].exonBlockIds[i] ].start + 1 ;
			return ret ;
		}

		void Split( const char *s, char delimit, std::vector<std::string> &fields )
		{
			int i ;
			fields.clear() ;
			if ( s == NULL )
				return ;

			std::string f ;
			for ( i = 0 ; s[i] ; ++i )
			{
				if ( s[i] == delimit || s[i] == '\n' )	
				{
					fields.push_back( f ) ;
					f.clear() ;
				}
				else
					f.append( 1, s[i] ) ;
			}
			fields.push_back( f ) ;
			f.clear() ;
		}
		
		// Test whether the j1th and j2th connection from gene block u form a simple bubble
		// simple bubble:
		// u - j1 -j2
		//  \______/
		// If so, return which one is closer
		int IsSimpleBubble( int u, int j1, int j2 )
		{
			int v1 = geneBlockGraph[u][j1].v ;
			int v2 = geneBlockGraph[u][j2].v ;
			
			int size ;
			int j ;

			size = geneBlockGraph[v1].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( geneBlockGraph[v1][j].v == v2 )
				{
					// Think about how the supportUse corresponding to the
					// dummy node.
					int su = geneBlockGraph[u][j1].supportUse & 2 ;
					su |= ( geneBlockGraph[v1][j].supportUse & 1 ) ;
				
					if ( geneBlockGraph[u][j2].supportUse == su ) 
						return j1 ;
					else
						return -1 ;
				}
			}

			size = geneBlockGraph[v2].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( geneBlockGraph[v2][j].v == v1 )
				{
					// Think about how the supportUse corresponding to the
					// dummy node.
					int su = geneBlockGraph[u][j2].supportUse & 2 ;
					su |= ( geneBlockGraph[v2][j].supportUse & 1 ) ;
				
					if ( geneBlockGraph[u][j1].supportUse == su ) 
						return j2 ;
					else
						return -1 ;
				}
			}
			return -1 ;
		}
	public:
		std::vector<struct _geneBlock> geneBlocks ;
		std::vector<struct _mateEdge> *geneBlockGraph ;

		int readLength ;
		int fragLength ;
		int fragStd ;

		Blocks() { geneBlockGraph = NULL ; } 	
		~Blocks() 
		{
			if ( geneBlockGraph != NULL )
				delete[] geneBlockGraph ;
		}

		int BuildExonBlocks( Alignments &alignments, Genome &genome )
		{
			unsigned int tag = 0 ;
			while ( alignments.Next() )
			{
				int i, j, k ;
				int segCnt = alignments.segCnt ;
				struct _pair *segments = alignments.segments ;

				while ( tag < exonBlocks.size() && ( exonBlocks[tag].end < segments[0].a - 1 
							|| exonBlocks[tag].chrId != alignments.GetChromId() ) )  
				{
					++tag ;
				}

				for ( i = 0 ; i < segCnt ; ++i )
				{
					//if ( i == 0 )
					//	printf( "hi %d %d %d\n", i, segments[i].a, segments[i].b ) ;
					for ( j = tag ; j < (int)exonBlocks.size() ; ++j )
					{
						if ( exonBlocks[j].end >= segments[i].a - 1 )
							break ;
					}

					if ( j >= (int)exonBlocks.size() )
					{
						// Append a new block
						struct _block newSeg ;
						newSeg.chrId = alignments.GetChromId() ;
						newSeg.start = segments[i].a ;
						newSeg.end = segments[i].b ;
						newSeg.support.Add( alignments ) ;

						newSeg.leftSplice = -1 ;
						newSeg.rightSplice = -1 ;
						if ( i > 0 )
							newSeg.leftSplice = segments[i].a ; 
						if ( i < segCnt - 1 )
							newSeg.rightSplice = segments[i].b ;

						exonBlocks.push_back( newSeg ) ;
					}
					else if ( exonBlocks[j].end < segments[i].b || 
							( exonBlocks[j].start > segments[i].a && exonBlocks[j].start <= segments[i].b + 1 ) ) 
					{
						// If overlaps with a current exon block, so we extend it
						if ( exonBlocks[j].end < segments[i].b )
						{
							// extends toward right 
							exonBlocks[j].end = segments[i].b ;
							exonBlocks[j].support.Add( alignments ) ;
							if ( i > 0 && ( exonBlocks[j].leftSplice == -1 || segments[i].a < exonBlocks[j].leftSplice ) )
								exonBlocks[j].leftSplice = segments[i].a ;
							if ( i < segCnt - 1 && segments[i].b > exonBlocks[j].rightSplice )
								exonBlocks[j].rightSplice = segments[i].b ;

							// Merge with next few exon blocks
							for ( k = j + 1 ; k < (int)exonBlocks.size() ; ++k )
							{
								if ( exonBlocks[k].start <= exonBlocks[j].end + 1 )
								{
									if ( exonBlocks[k].end > exonBlocks[j].end )
										exonBlocks[j].end = exonBlocks[k].end ;

									if ( exonBlocks[k].leftSplice != -1 && ( exonBlocks[j].leftSplice == -1 || exonBlocks[k].leftSplice < exonBlocks[j].leftSplice ) )
										exonBlocks[j].leftSplice = exonBlocks[k].leftSplice ;
									if ( exonBlocks[k].rightSplice != -1 && exonBlocks[k].rightSplice > exonBlocks[j].rightSplice )
										exonBlocks[j].rightSplice = exonBlocks[k].rightSplice ;

									exonBlocks[j].support.Add( exonBlocks[k].support ) ;
								}
								else
									break ;
							}

							if ( k > j + 1 )
							{
								// Remove the merged blocks
								int a, b ;
								for ( a = j + 1, b = k ; b < (int)exonBlocks.size() ; ++a, ++b )
									exonBlocks[a] = exonBlocks[b] ;
								for ( a = 0 ; a < k - ( j + 1 ) ; ++a )
									exonBlocks.pop_back() ;
							}
						}
						else if ( exonBlocks[j].start > segments[i].a && exonBlocks[j].start <= segments[i].b + 1 ) 
						{
							// extends toward left
							exonBlocks[j].start = segments[i].a ;
							exonBlocks[j].support.Add( alignments ) ;
							if ( i > 0 && ( exonBlocks[j].leftSplice == -1 || segments[i].a < exonBlocks[j].leftSplice ) )
								exonBlocks[j].leftSplice = segments[i].a ;
							if ( i < segCnt - 1 && segments[i].b > exonBlocks[j].rightSplice )
								exonBlocks[j].rightSplice = segments[i].b ;

							// Merge with few previous exon blocks
							for ( k = j - 1 ; k >= 0 ; --k )
							{
								if ( exonBlocks[k].end >= exonBlocks[k + 1].start - 1 )
								{
									if ( exonBlocks[k + 1].start < exonBlocks[k].start )
									{
										exonBlocks[k].start = exonBlocks[k + 1].start ;
									}

									if ( exonBlocks[k].leftSplice != -1 && ( exonBlocks[j].leftSplice == -1 || exonBlocks[k].leftSplice < exonBlocks[j].leftSplice ) )
										exonBlocks[j].leftSplice = exonBlocks[k].leftSplice ;
									if ( exonBlocks[k].rightSplice != -1 && exonBlocks[k].rightSplice > exonBlocks[j].rightSplice )
										exonBlocks[j].rightSplice = exonBlocks[k].rightSplice ;

									exonBlocks[k].support.Add( exonBlocks[k + 1].support ) ;
								}
								else
									break ;
							}

							if ( k < j - 1 )
							{
								int a, b ;
								for ( a = k + 2, b = j + 1 ; b < (int)exonBlocks.size() ; ++a, ++b )
									exonBlocks[a] = exonBlocks[b] ;
								for ( a = 0 ; a < ( j - 1 ) - k ; ++a )
									exonBlocks.pop_back() ;

							}
						}
					}
					else if ( exonBlocks[j].start > segments[i].b + 1 )
					{
						int size = exonBlocks.size() ;
						int a ;
						// No overlap, insert a new block
						struct _block newSeg ;
						newSeg.chrId = alignments.GetChromId() ;
						newSeg.start = segments[i].a ;
						newSeg.end = segments[i].b ;
						newSeg.support.Add( alignments ) ;

						newSeg.leftSplice = -1 ;
						newSeg.rightSplice = -1 ;
						if ( i > 0 )
							newSeg.leftSplice = segments[i].a ; 
						if ( i < segCnt - 1 )
							newSeg.rightSplice = segments[i].b ;

						// Insert at position j
						exonBlocks.push_back( newSeg ) ;	
						for ( a = size ; a > j ; --a )
							exonBlocks[a] = exonBlocks[a - 1] ;
						exonBlocks[a] = newSeg ;
					}
					else
					{
						// The segment is contained in j
						exonBlocks[j].support.Add( alignments ) ;
					}
				}
			}

			/*for ( int i = 0 ; i < (int)exonBlocks.size() ; ++i )
			  {
			  printf( "%d %d %d\n", exonBlocks[i].start, exonBlocks[i].end, exonBlocks[i].support.GetCount() ) ;
			  }*/	

			if ( exonBlocks.size() > 0 )
			{
				// Build the map for the offsets of chr id in the exonBlock list.
				exonBlocksChrIdOffset[ exonBlocks[0].chrId] = 0 ;
				int cnt = exonBlocks.size() ;
				for ( int i = 1 ; i < cnt ; ++i )
				{
					if ( exonBlocks[i].chrId != exonBlocks[i - 1].chrId )
						exonBlocksChrIdOffset[ exonBlocks[i].chrId ] = i ;
				}

				// Put the contig id.
				if ( genome.IsOpen() )
				{
					for ( int i = 0 ; i < cnt ; ++i )
					{
						int id1 = genome.GetContigId( exonBlocks[i].chrId, exonBlocks[i].start ) ;
						int id2 = genome.GetContigId( exonBlocks[i].chrId, exonBlocks[i].end ) ;

						if ( id1 == -1 && id2 == -1 )
							exonBlocks[i].contigId = -1 ; // TODO: should not happen, remove it if necessary
						else if ( id1 != -1 )
							exonBlocks[i].contigId = id1 ;
						else
							exonBlocks[i].contigId = id2 ;
					}
				}
				else
				{
					for ( int i = 0 ; i < cnt ; ++i )
						exonBlocks[i].contigId = exonBlocks[i].chrId ;
				}
			}
			return exonBlocks.size() ;
		}

		// Extend the exon blocks with other blocks.
		// Notice that how to control the number of support
		int ExtendExonBlocks( Blocks &otherBlocks )
		{
			std::vector<struct _block> origExonBlocks( exonBlocks ) ; 
			std::vector<struct _block> &otherExonBlocks = otherBlocks.exonBlocks ;


			exonBlocks.clear() ;
			int i = 0, j = 0, k = 0 ;
			int ret = -1 ;
			int origSize = origExonBlocks.size() ;
			int otherSize = otherExonBlocks.size() ;
			bool containOrig = false ;
			while ( 1 )
			{
				if ( i >= origSize && j >= otherSize )
					break ;

				if ( i < origSize && ( j >= otherSize || origExonBlocks[i].chrId < otherExonBlocks[j].chrId ||
						( origExonBlocks[i].chrId == otherExonBlocks[j].chrId && origExonBlocks[i].start <= otherExonBlocks[j].start ) ) )
				{
					if ( ret >= 0 && exonBlocks[ret].chrId == origExonBlocks[i].chrId && origExonBlocks[i].start <= exonBlocks[ret].end )	
					{
						// merge into the current block
						if ( !containOrig )
							exonBlocks[ret].support = origExonBlocks[i].support ;
						else
						{
							exonBlocks[ret].support.Add( origExonBlocks[i].support ) ;
						}
						if ( origExonBlocks[i].end > exonBlocks[ret].end )
							exonBlocks[ret].end = origExonBlocks[i].end ;
						containOrig = true ;
					}
					else
					{
						exonBlocks.push_back( origExonBlocks[i] ) ;
						containOrig = true ;
						++ret ;
					}
					++i ;
				}
				else
				{
					if ( ret >= 0 && exonBlocks[ret].chrId == otherExonBlocks[j].chrId && otherExonBlocks[j].start <= exonBlocks[ret].end )	
					{
						// merge into the current block
						if ( !containOrig )
							exonBlocks[ret].support.Add( otherExonBlocks[j].support ) ;
						if ( otherExonBlocks[j].end > exonBlocks[ret].end )
							exonBlocks[ret].end = otherExonBlocks[j].end ;
						containOrig = false ;
					}
					else
					{
						exonBlocks.push_back( otherExonBlocks[j] ) ;
						//exonBlocks[ret + 1].support.Clear() ;
						containOrig = false ;
						++ret ;
					}
					++j ;
				}
			}

			// Rebuild the map for the offsets of chr id in the exonBlock list.
			if ( ret >= 0 )
			{
				exonBlocksChrIdOffset.clear() ;

				exonBlocksChrIdOffset[ exonBlocks[0].chrId] = 0 ;
				int cnt = exonBlocks.size() ;
				for ( int i = 1 ; i < cnt ; ++i )
				{
					if ( exonBlocks[i].chrId != exonBlocks[i - 1].chrId )
						exonBlocksChrIdOffset[ exonBlocks[i].chrId ] = i ;
				}
			}

			return ret + 1 ;
		}

		void GetAlignmentsInfo( Alignments &alignments )
		{
			int tag = 0 ;
			int64_t totalReadLength = 0 ;
			int readCnt = 0 ;
			int64_t sum = 0 ;
			int64_t sqSum = 0 ;
			int exonBlockCnt = exonBlocks.size() ;

			while ( alignments.Next() && readCnt < 100000 )	
			{
				int i, j, k ;
				int segCnt = alignments.segCnt ;
				struct _pair *segments = alignments.segments ;

				if ( tag < exonBlockCnt && alignments.GetChromId() != exonBlocks[tag].chrId )
				{
					std::map<int, int>::iterator it ;

					it = exonBlocksChrIdOffset.find( alignments.GetChromId() ) ;

					if ( it != exonBlocksChrIdOffset.end() )
					{
						if ( tag < it->second )
							tag = it->second ;
						else
							continue ;
					}
					else
						continue ; // skip this read
				}

				while ( tag < exonBlockCnt && exonBlocks[tag].end < segments[0].a && 
						exonBlocks[tag].chrId == alignments.GetChromId() )
				{
					++tag ;
				}

				int64_t mPos ;
				int mChrId ;

				k = -1 ;
				for ( j = tag ; j < exonBlockCnt && exonBlocks[j].end < segments[segCnt - 1].b 
						&& exonBlocks[j].chrId == alignments.GetChromId() ; ++j )	
					;
				if ( j < exonBlockCnt && exonBlocks[j].chrId == alignments.GetChromId() && exonBlocks[j].start <= segments[ segCnt - 1 ].a )
				{
					k = j ;
				}
				if ( k == -1 )
					continue ;
				alignments.GetMatePosition( mChrId, mPos ) ;

				if ( mChrId != alignments.GetChromId() || mPos < exonBlocks[j].start || mPos > exonBlocks[j].end || mPos < segments[0].a )
					continue ;
				int rl = 0 ;
				for ( i = 0 ; i < segCnt ; ++i )
				{
					rl += segments[i].b - segments[i].a + 1 ;
				}
				totalReadLength += rl ;
				int tmp = 2 * rl + mPos - segments[ segCnt - 1 ].a - 1 ;
				sum += tmp ;
				sqSum += tmp * tmp ;
				++readCnt ;
			}

			assert( readCnt > 30 ) ;

			readLength = totalReadLength / readCnt ;
			fragLength = sum / readCnt ;
			fragStd = sqrt( (double)sqSum / readCnt - fragLength * fragLength ) ;
			//fprintf( stderr, "Fragment length: %d std: %d\n", (int)fragLength, (int)fragStd ) ;
			/*readLength = 100 ;
			  fragLength = 200 ;
			  fragStd = 50 ;*/
		}


		int BuildGeneBlocks( Alignments &alignments )
		{
			int i, j, k ;
			int exonBlockCnt = exonBlocks.size() ;
			int tag = 0 ;


			//struct _adjGraph *adj = (struct _adjGraph *)malloc( exonBlockCnt * sizeof( *adj ) )  ;
			std::vector<struct _edge> *adj = new std::vector<struct _edge>[exonBlockCnt] ;
			int *father = new int[ exonBlockCnt] ;

			for ( i = 0 ; i < exonBlockCnt ; ++i )
			{
				father[i] = i ;
			}

			while ( alignments.Next() )
			{
				int i, j, k ;
				int segCnt = alignments.segCnt ;
				struct _pair *segments = alignments.segments ;

				if ( tag < exonBlockCnt && alignments.GetChromId() != exonBlocks[tag].chrId )
				{
					std::map<int, int>::iterator it ;

					it = exonBlocksChrIdOffset.find( alignments.GetChromId() ) ;

					if ( it != exonBlocksChrIdOffset.end() )
					{
						if ( tag < it->second )
							tag = it->second ;
						else
							continue ;
					}
					else
						continue ; // skip this read
				}

				while ( tag < exonBlockCnt && exonBlocks[tag].end < segments[0].a && 
						exonBlocks[tag].chrId == alignments.GetChromId() )
				{
					++tag ;
				}

				//printf( "%s %d\n", __func__, tag ) ;
				// Find which blocks this read is compatible with
				int segmentBlocks[MAX_SEG_COUNT] ;
				//printf( "%d\n", segCnt ) ;
				k = tag ;
				for ( i = 0 ; i < segCnt ; ++i )
				{
					segmentBlocks[i] = -1 ;
					for ( j = k ; j < exonBlockCnt && exonBlocks[j].end < segments[i].b 
							&& exonBlocks[j].chrId == alignments.GetChromId() ; ++j )	
						;
					if ( j < exonBlockCnt && exonBlocks[j].chrId == alignments.GetChromId() )
					{
						if ( exonBlocks[j].start <= segments[i].a && exonBlocks[j].end >= segments[i].b )
						{
							segmentBlocks[i] = j ;
						}
					}
					k = j ;
				}
				
				/*if ( !strcmp( alignments.GetReadId(),  "Id_5023627" ) )
				{
					for ( i = 0 ; i < segCnt ; ++i )
					{
						printf( "%d\n", segmentBlocks[i] ) ;
					}
				}*/
				// Connect those exon blocks
				// Notice that the segmentBlocks are already sorted.
				//continue ;
				for ( i = 0 ; i < segCnt - 1 ; ++i )
				{
					if ( segmentBlocks[i] == -1 || segmentBlocks[i + 1] == -1 )
						continue ;

					int size = adj[ segmentBlocks[i] ].size() ;
					for ( j = size - 1 ; j >= 0 ; --j )
					{
						if ( adj[ segmentBlocks[i] ][j].v == segmentBlocks[i + 1] )
						{
							adj[segmentBlocks[i]][j].support.Add( alignments ) ; 
							break ;
						}
					}

					if ( j < 0 )
					{
						struct _edge newP ;
						newP.u = segmentBlocks[i] ;
						newP.v = segmentBlocks[i + 1] ;
						newP.support.Add( alignments ) ;

						adj[ segmentBlocks[i] ].push_back( newP ) ;
					}
				}

				// Also connect the exon blocks by the mate-pair information.
				for ( i = 0 ; i < segCnt ; ++i )
				{
					if ( segmentBlocks[i] != -1 )
						break ;
				}
				
				if ( segmentBlocks[i] != -1 )
				{
					int mChrId ;
					int64_t mPos ;

					alignments.GetMatePosition( mChrId, mPos ) ;
					
					if ( mChrId == alignments.GetChromId() 
						&& alignments.IsReverse() != alignments.IsMateReverse() && mPos < segments[0].a )
					{
						int k = segmentBlocks[i] ;
						for ( ; k >= 0 ; --k )
						{
							if ( exonBlocks[k].start <= mPos && mPos <= exonBlocks[k].end )
							{
								int size = adj[ segmentBlocks[i] ].size() ;
								for ( j = size - 1 ; j >= 0 ; --j )
								{
									if ( adj[ segmentBlocks[i] ][j].v == k )
									{
										adj[segmentBlocks[i]][j].support.Add( alignments ) ; 
										break ;
									}
								}

								if ( j < 0 )
								{
									struct _edge newP ;
									newP.u = segmentBlocks[i] ;
									newP.v = k ;
									newP.support.Add( alignments ) ;
									adj[ segmentBlocks[i] ].push_back( newP ) ;
								}
								break ;
							}
							if ( exonBlocks[k].chrId != mChrId || mPos > exonBlocks[k].end )
								break ;
						}
					}
				}
			}
			// Clustering the exon blocks
			for ( i = 0 ; i < exonBlockCnt ; ++i )
			{
				int size = adj[i].size() ;
				for ( j = 0 ; j < size ; ++j ) 
				{
					//Filter out the connections 
					if ( exonBlocks[i].chrId != exonBlocks[j].chrId && 
						( !adj[i][j].support.IsGood() /*|| adj[i][j].support.GetCount() < minimumSupport */ ) )
						continue ;
					
					int a = GetFather( i, father ) ;
					int b = GetFather( adj[i][j].v, father ) ;
					if ( a <= b )
						father[b] = a ;
					else
						father[a] = b ;
				}
			}
			//printf( "father: %d %d\n", GetFather( 19564, father) , GetFather( 19565, father ) ) ;
			//filter out gene blocks.
			// We do this in this stage because we will break the gene blocks
			// when we actually build them due to interleaved genes.
			//Support *geneBlockSupports = new Support[ exonBlockCnt ] ;
			struct _pair *geneBlockGoodCnt = new struct _pair[ exonBlockCnt ] ;
			bool *validGeneBlock = new bool[ exonBlockCnt ] ;

			memset( geneBlockGoodCnt, 0, sizeof( geneBlockGoodCnt[0] ) * exonBlockCnt ) ;

			for ( i = 0 ; i < exonBlockCnt ; ++i )
			{
				int f = GetFather( i, father ) ;
				// TODO: we can also look at each exon block.
				//geneBlockSupports[f].Add( exonBlocks[i].support ) ;
				if ( exonBlocks[i].support.IsGood() )
					++geneBlockGoodCnt[f].a ;
				++geneBlockGoodCnt[f].b ;

				//if ( f == 5 )
				//	printf( "%d %d %d\n",i, geneBlockGoodCnt[f].a, geneBlockGoodCnt[f].b ) ;
			}

			for ( i = 0 ; i < exonBlockCnt ; ++i )
			{
				if ( father[i] == i )
				{
					if ( geneBlockGoodCnt[i].a >= 1 
							&&  geneBlockGoodCnt[i].a >= 0.2 * geneBlockGoodCnt[i].b )
						validGeneBlock[i] = true ;
					else
						validGeneBlock[i] = false ;

					// Possible pseudo-gene
					/*if ( geneBlockGoodCnt[i].b == 1 && !exonBlocks[i].support.IsUnique() )
					  {
					  validGeneBlock[i] = false ;
					  }*/
				}
				//if ( i == 5 )
				//	printf( "%d %d %d: %d\n",i, geneBlockGoodCnt[i].a, geneBlockGoodCnt[i].b, validGeneBlock[i] ) ;
			}
			delete []geneBlockGoodCnt ;

			int geneBlockCnt = 0 ;
			int currentFather = -1 ;
			for ( i = 0 ; i < exonBlockCnt ; ++i )
			{
				// Ignore the invalid gene blocks.
				if ( !validGeneBlock[ GetFather(i, father ) ] )
					continue ;
				// We split the cluster if the blocks are interleaved.
				if ( GetFather( i, father ) != currentFather || 
						( i > 0 && exonBlocks[i].contigId != exonBlocks[i - 1].contigId ) )
				{
					// Create a new gene block
					struct _geneBlock newBlock ;

					newBlock.chrId = exonBlocks[i].chrId ;
					newBlock.contigId = exonBlocks[i].contigId ;
					newBlock.start = exonBlocks[i].start ;
					newBlock.end = exonBlocks[i].end ;
					newBlock.exonBlockIds.push_back( i ) ;
					geneBlocks.push_back( newBlock ) ;

					++geneBlockCnt ;
					currentFather = GetFather( i, father ) ; 
				}
				else
				{
					// Update the gene block
					/*int f = GetFather( i, father ) ;
					  for ( j = geneBlockCnt - 1 ; j >= 0 ; --j )
					  {
					  if ( geneBlocks[j].exonBlockIds[0] == f )
					  break ;
					  }
					 */
					j = geneBlockCnt - 1 ;
					assert( j >= 0 ) ;
					geneBlocks[j].exonBlockIds.push_back( i ) ;
					if ( exonBlocks[i].start < geneBlocks[j].start )
						geneBlocks[j].start = exonBlocks[i].start ;
					if ( exonBlocks[i].end > geneBlocks[j].end )
						geneBlocks[j].end = exonBlocks[i].end ;
				}
			}
			delete[] adj ;
			delete[] father ;
			delete[] validGeneBlock ;

			// Determine the strand of the gene block
			for ( i = 0 ; i < geneBlockCnt ; ++i )
			{
				int cnt = geneBlocks[i].exonBlockIds.size() ;
				Support support ;
				for ( j = 0 ; j < cnt ; ++j )
					support.Add( exonBlocks[ geneBlocks[i].exonBlockIds[j] ].support ) ;
				geneBlocks[i].support = support ;
			}

			geneBlocksChrIdOffset[ geneBlocks[0].chrId ].a = 0 ;
			for (  i = 1 ; i < geneBlockCnt ; ++i )
			{
				if (  geneBlocks[i].chrId != geneBlocks[i - 1].chrId )
				{
					geneBlocksChrIdOffset[ geneBlocks[i].chrId ].a = i ;
					geneBlocksChrIdOffset[ geneBlocks[i - 1].chrId].b = i - 1 ;
				}
			}
			geneBlocksChrIdOffset[ geneBlocks[i - 1].chrId ].b = i - 1 ;

			if ( VERBOSE )
			{
				fprintf( fpOut, "Gene blocks:\n" ) ;
				for ( i = 0 ; i< geneBlocks.size() ; ++i )
				{
					fprintf( fpOut, "geneblock %d: %s %"PRId64" %"PRId64"\n", i, alignments.GetChromName( geneBlocks[i].chrId ), geneBlocks[i].start + 1, geneBlocks[i].end + 1 ) ; 
					int size = geneBlocks[i].exonBlockIds.size() ;
					for ( j = 0 ; j < size ; ++j )
					{
						int id = geneBlocks[i].exonBlockIds[j] ;
						fprintf( fpOut, "\t%d: %"PRId64" %"PRId64"\n", j, exonBlocks[id].start + 1, exonBlocks[id].end + 1 ) ;
					}
				}
			}
			return geneBlocks.size() ;
		}

		int GetGeneBlockEffectiveLength( int ind )
		{
			return GetEffectiveLength( ind ) ;
		}

		int GetGeneBlockEffectiveSequence( int ind, int flip, char *buffer, Genome &genome )
		{
			std::vector<int> &exonBlockIds = geneBlocks[ind].exonBlockIds ;
			int ecnt = exonBlockIds.size() ;
			int i ;
			int ret = 0 ;

			if ( flip == 0 )
			{
				for ( i = 0 ; i < ecnt ; ++i )
				{
					struct _block &eblock = exonBlocks[ exonBlockIds[i] ] ;
					int chrId = eblock.chrId ;

					for ( int l = eblock.start ; l <= eblock.end ; ++l, ++ret )
					{
						buffer[ret] = genome.GetNucleotide( chrId, l ) ;
					}
					//printf( "++ %d %d %d\n", ret, (int)eblock.start, (int)eblock.end ) ;
				}
			}
			else
			{
				for ( i = ecnt - 1 ; i >= 0 ; --i )
				{
					struct _block &eblock = exonBlocks[ exonBlockIds[i] ] ;
					int chrId = eblock.chrId ;

					for ( int l = eblock.end ; l >= eblock.start ; --l, ++ret )
					{
						//if ( chrId == 74619 || chrId == 91834 )
						//	printf( "%c", genome.GetNucleotide( chrId, l ) ) ;
						buffer[ret] = genome.GetNucleotide( chrId, l ) ;
						if ( nucToNum[ buffer[ret] - 'A' ] != -1 )
							buffer[ret] = numToNuc[ 3 - nucToNum[ buffer[ret] - 'A' ] ] ; // The complementary nucleotide
						else
							buffer[ret] = 'N' ;
					}
					//printf( "-- %d %d %d\n", ret, (int)eblock.start, (int)eblock.end ) ;
				}

			}
			buffer[ret] = '\0' ;
			//if ( exonBlocks[ exonBlockIds[0] ].chrId == 74619 || exonBlocks[ exonBlockIds[0] ].chrId == 91834 )
			//	printf( "(%d %lld)\n%s\n", ret, buffer, buffer ) ;
			return ret ;
		}

		int GetGeneBlockCount()
		{
			return geneBlocks.size() ;
		}

		void BuildGeneBlockGraph( Alignments &alignments )	
		{
			int i, j, k ;
			int tag = 0 ;
			int geneBlockCnt = geneBlocks.size() ;
			
			std::vector<struct _edge> *repeatGraph = new std::vector<struct _edge>[ geneBlockCnt ] ;
			geneBlockGraph = new std::vector<struct _mateEdge>[ geneBlockCnt ] ;
			repeatFather = new int[geneBlockCnt ] ;

			for ( i = 0 ; i < geneBlockCnt ; ++i )
			{
				repeatFather[i] = i ;	
			}

			while ( alignments.Next() )
			{
				struct _pair *segments = alignments.segments ;
				if ( tag < geneBlockCnt && alignments.GetChromId() != geneBlocks[tag].chrId )
				{
					std::map<int, struct _pair>::iterator it ;

					it = geneBlocksChrIdOffset.find( alignments.GetChromId() ) ;

					if ( it != geneBlocksChrIdOffset.end() )
					{
						if ( tag < it->second.a )
							tag = it->second.a ;
						else
							continue ;
					}
					else
						continue ; // skip this read
				}

				while ( tag < geneBlockCnt && geneBlocks[tag].end < segments[0].a && 
						geneBlocks[tag].chrId == alignments.GetChromId() )
				{
					++tag ;
				}

				/*if ( !strcmp( alignments.GetReadId(), "Id_10718637" ) )
				  {
				  printf( "%d %lld: %d %lld %lld %d\n", alignments.GetChromId(), segments[0].a, tag, geneBlocks[tag].start, geneBlocks[tag].end, geneBlocks[tag].chrId ) ;
				  }*/
				if ( tag >= geneBlockCnt || geneBlocks[tag].chrId != alignments.GetChromId() )
					continue ;
				int tagE = GetExonBlockInGeneBlock( tag, geneBlocks[tag].chrId, segments[0].a ) ;
				if ( tagE == -1 )
					continue ;
				// Test the mates
				// Notice that, we add the information twice.
				int mChrId ;
				int64_t mPos ;
				int kE ;
				alignments.GetMatePosition( mChrId, mPos ) ;
				k = FindGeneBlock( mChrId, mPos ) ;

				/*if ( !strcmp( alignments.GetReadId(), "Id_30314187" ) )
				  {
				  printf( "m: %d %lld: %d %lld %lld %d\n", mChrId, mPos, k, geneBlocks[k].start, geneBlocks[k].end, geneBlocks[k].chrId ) ;
				  }*/
				// Skip the read if it does not compatible, or the two mates are in the
				// same gene block.
				if ( k == -1 || k == tag ) 
					continue ;
				kE = GetExonBlockInGeneBlock( k, mChrId, mPos ) ;
				if ( kE == -1 )
					continue ;
				//if ( tag == 281 && k == 279 )
				//{
				//		printf( "%d %d %s\n", tag, k, alignments.GetReadId() ) ;		
				//}

				// Add the edge
				int cnt = geneBlockGraph[tag].size() ;
				int directionTag = 0 ;

				// TODO: may need to change to hash table for edge representation
				directionTag = alignments.IsReverse() ? 2 : 0 ;
				directionTag |= ( alignments.IsMateReverse() ? 1 : 0 ) ;

				// Test the strand. 
				// Defered to clean block graph stage.
				/*if ( geneBlocks[tag].support.GetStrand() != 0 && geneBlocks[k].support.GetStrand() != 0 )	
				  {
				  int s1 = geneBlocks[tag].support.GetStrand() ;
				  int s2 = geneBlocks[k].support.GetStrand() ;

				  if ( alignments.IsReverse() == alignments.IsMateReverse() )
				  {
				  if ( s1 == s2 )
				  continue ;
				  }
				  else
				  {
				  if ( s1 != s2 )
				  continue ;
				  }
				  }*/

				// Test the insert size
				int insert1 = GetGeneBlockResidual( tag, tagE, segments[0].a, alignments.IsReverse() ? -1 : 1 ) ;
				int insert2 = GetGeneBlockResidual( k, kE, mPos, alignments.IsMateReverse() ? -1 : 1 ) ;
				/*if ( !strcmp( alignments.GetReadId(), "Id_10718637" ) )
				{
					printf( "%d %d\n", insert1 + insert2, fragLength + 2 * fragStd ) ;
				}*/
				if ( insert1 + insert2 <= fragLength + 2 * fragStd ) // 400 here is to take short alternative splicing events into account.
				{
					// Add or update the edge
					/*if ( tag == 1084 && k == 392 )
					  {
					  printf( "%s\n", alignments.GetReadId() ) ;
					  }*/
					for ( i = 0 ; i < cnt ; ++i )
					{
						if ( geneBlockGraph[tag][i].v == k )
						{
							geneBlockGraph[tag][i].support[ directionTag ].Add( alignments ) ;
							break ;
						}
					}

					if ( i >= cnt )
					{
						struct _mateEdge newE ;
						newE.u = tag ;
						newE.v = k ;
						//memset( newE.support, 0, sizeof( newE.support ) ) ;

						newE.support[ directionTag ].Add( alignments ) ;
						newE.valid = true ;

						geneBlockGraph[tag].push_back( newE ) ;
					}
				}

				// The part for the repeat graph, using CC and CP field
				int64_t cPos ;
				int cChrId ;
				if ( alignments.GetRepeatPosition( cChrId, cPos ) )
				{
					k = FindGeneBlock( cChrId, cPos ) ;

					// Skip the read if it does not compatible, or the two mates are in the
					// same gene block.
					if ( k != -1 && GetExonBlockInGeneBlock( k, mChrId, mPos ) != -1 ) 
					{
						// Add or update the edge
						cnt = repeatGraph[tag].size() ;
						for ( i = 0 ; i < cnt ; ++i )
						{
							if ( repeatGraph[tag][i].v == k )
							{
								repeatGraph[tag][i].support.Add( alignments ) ;
								break ;
							}
						}

						if ( i >= cnt )
						{
							struct _edge newE ;
							newE.u = tag ;
							newE.v = k ;
							//memset( newE.support, 0, sizeof( newE.support ) ) ;

							newE.support.Add( alignments ) ;
							repeatGraph[tag].push_back( newE ) ;
						}

					}
				}
			}

			// Get the sets of repeated gene blocks
			for ( i = 0 ; i < geneBlockCnt ; ++i )
			{
				int cnt = repeatGraph[i].size() ;
				for ( j = 0 ; j < cnt ; ++j )
				{
					//TODO: check the threshold
					if ( repeatGraph[i][j].support.GetCount() >= minimumSupport ) 
					{
						int fi = GetFather( i, repeatFather ) ;
						int fj = GetFather( repeatGraph[i][j].v, repeatFather ) ;
						repeatFather[fj] = fi ;
					}
				}
			}

			delete[] repeatGraph ;
		}

		// Use the clipped alignment information to add the weight 
		// of the gene block graph
		void AddGeneBlockGraphByClippedAlignments( Alignments alignments )
		{
			int i, k ;
			int tag = 0 ;
			int geneBlockCnt = geneBlocks.size() ;
			while ( alignments.Next() )
			{
				struct _pair *segments = alignments.segments ;
				int start = segments[0].a ;
				bool isReverse = alignments.IsReverse() ;
				bool reverseMateRole = false ; // triggered when letting the mate connect to the supplementary alignment

				if ( tag < geneBlockCnt && alignments.GetChromId() != geneBlocks[tag].chrId )
				{
					std::map<int, struct _pair>::iterator it ;

					it = geneBlocksChrIdOffset.find( alignments.GetChromId() ) ;

					if ( it != geneBlocksChrIdOffset.end() )
					{
						if ( tag < it->second.a )
							tag = it->second.a ;
						else
							continue ;
					}
					else
						continue ; // skip this read
				}

				while ( tag < geneBlockCnt && geneBlocks[tag].end < start && 
						geneBlocks[tag].chrId == alignments.GetChromId() )
				{
					++tag ;
				}

				if ( tag >= geneBlockCnt || geneBlocks[tag].chrId != alignments.GetChromId() )
					continue ;
				int tagG = tag ;
				int tagE = GetExonBlockInGeneBlock( tagG, geneBlocks[tag].chrId, start ) ;
				if ( tagE == -1 )
					continue ;

				// Test the mates
				// Notice that, we add the information twice.
				int mChrId ;
				int64_t mPos ;
				int kE ;
				alignments.GetMatePosition( mChrId, mPos ) ;
				k = FindGeneBlock( mChrId, mPos ) ;

				/*if ( !strcmp( alignments.GetReadId(), "Id_30063560" ) )
				  {
				  printf( "m: %d %lld: %d %lld %lld %d\n", mChrId, mPos, k, geneBlocks[k].start, geneBlocks[k].end, geneBlocks[k].chrId ) ;
				  }*/
				// Skip the read if it does not compatible, or the two mates are in the
				// same gene block.
				//if ( k == -1 ) 
				//	continue ;
				bool skippable = false ;
				if ( k == -1 || k == tagG )
					skippable = true ;
				if ( skippable || geneBlocks[tagG].contigId == geneBlocks[k].contigId )
				{
					if ( skippable && alignments.IsSupplementary() )
						continue ;
					// Look at other aligned part, if it aligned to
					// a different gene block, then use that as the mate.
					char *SA = alignments.GetFieldZ( "SA" ) ;
					if ( SA != NULL )
					{
						std::vector< std::string > fragments ;
						std::vector< std::string > hit ;

						Split( SA, ';', fragments ) ;
						int fcnt = fragments.size() ;
						int f ;
						for ( f = 0 ; f < fcnt ; ++f )
						{
							//SA:Z:chr20_343,89495,+,36M64S,60,0;
							if ( fragments[f].c_str()[0] == '\0' )
								continue ;
							Split( fragments[f].c_str(), ',', hit ) ;
							int hChrId = alignments.GetChromIdFromName( hit[0].c_str() ) ;
							int hPos = atoi( hit[1].c_str() ) ;
							int gb = FindGeneBlock( hChrId, hPos ) ;
							
							if ( gb == -1 )
								continue ;
							if ( ( skippable && gb != tagG ) || 
								( !skippable && geneBlocks[gb].contigId != geneBlocks[k].contigId ) )
							{
								start = hPos ;
								tagG = gb ;
								tagE = GetExonBlockInGeneBlock( tagG, geneBlocks[tagG].chrId, start ) ;
								isReverse = ( hit[2].c_str() )[0] == '+' ? false : true ;
								if ( tagE == -1 )
									continue ;
								reverseMateRole = true ;
								break ;
							}
						}
						if ( f >= fcnt )
							continue ;
					}
					else
						continue ;
				}
				kE = GetExonBlockInGeneBlock( k, mChrId, mPos ) ;
				if ( kE == -1 )
					continue ;
				//if ( tagG == 281 && k == 279 )
				//{
				//		printf( "%d %d %s\n", tagG, k, alignments.GetReadId() ) ;		
				//}

				// Add the edge
				int directionTag = 0 ;

				// TODO: may need to change to hash table for edge representation
				directionTag = isReverse ? 2 : 0 ;
				directionTag |= ( alignments.IsMateReverse() ? 1 : 0 ) ;

				// Test the insert size
				int insert1 = GetGeneBlockResidual( tagG, tagE, start, isReverse ? -1 : 1 ) ;
				int insert2 = GetGeneBlockResidual( k, kE, mPos, alignments.IsMateReverse() ? -1 : 1 ) ;
				
				if ( reverseMateRole )
				{
					int tmp ;
					tmp = tagG ;
					tagG = k ;
					k = tmp ;
					directionTag ^= 3 ;
				}
				int cnt = geneBlockGraph[tagG].size() ;
				if ( insert1 + insert2 <= fragLength + 2 * fragStd )
				{
					// Add or update the edge
					/*if ( tagG == 1162 && k == 1165 )
					  {
					  printf( "%s\n", alignments.GetReadId() ) ;
					  }*/
					for ( i = 0 ; i < cnt ; ++i )
					{
						if ( geneBlockGraph[tagG][i].v == k )
						{
							// TODO: should I use a different counter for clipped reads?
							geneBlockGraph[tagG][i].support[ directionTag ].Add( alignments, reverseMateRole ) ;
							break ;
						}
					}

					if ( i >= cnt )
					{
						struct _mateEdge newE ;
						newE.u = tagG ;
						newE.v = k ;
						//memset( newE.support, 0, sizeof( newE.support ) ) ;

						newE.support[ directionTag ].Add( alignments, reverseMateRole ) ;
						newE.valid = true ;

						geneBlockGraph[tagG].push_back( newE ) ;
					}
				}
			}
		}

		// -1: if we can not find one
		int GetExonBlockInGeneBlock( int geneBlockInd, int32_t chrId, int64_t pos )
		{
			int i ;
			int exonCnt = geneBlocks[ geneBlockInd ].exonBlockIds.size() ;

			if ( chrId != geneBlocks[geneBlockInd].chrId )
				return -1 ;

			for ( i = 0 ; i < exonCnt ; ++i )
			{
				int k = geneBlocks[ geneBlockInd ].exonBlockIds[i] ;
				if ( exonBlocks[k].start > pos ) 
					return -1 ;
				else if ( pos >= exonBlocks[k].start && pos <= exonBlocks[k].end )	
					return i ;
			}
			return -1 ;
		}

		int GetGeneBlockResidual( int gid, int eid, int pos, int direction )
		{
			// TODO: cache the previous result
			int i ;
			int cnt = geneBlocks[gid].exonBlockIds.size() ;
			int ret = 0 ;
			for ( i = eid + direction ; i >= 0 && i < cnt ; i += direction )
			{
				struct _block &e = exonBlocks[ geneBlocks[gid].exonBlockIds[i] ] ; 
				int l = e.leftSplice ;
				int r = e.rightSplice ;
				if ( l == -1 )
					l = e.start ;
				if ( r == -1 )
					r = e.end ;
				ret += r - l + 1 ; 
			}


			return ret ;
		}


		int FindGeneBlock( int chrId, int64_t pos ) 
		{
			int l, r, m ;

			// TODO: cache the previous result.
			std::map<int, struct _pair>::iterator it = geneBlocksChrIdOffset.find( chrId ) ;

			if ( it == geneBlocksChrIdOffset.end() )
				return -1 ;

			l = it->second.a ;
			r = it->second.b ;

			while ( l <= r )
			{
				m = ( l + r ) / 2 ;
				if ( geneBlocks[m].start <= pos && geneBlocks[m].end >= pos )
					return m ;

				if ( geneBlocks[m].start > pos )
					r = m - 1 ;
				else
					l = m + 1 ;
			}

			return -1 ;
		}

		// Take the exons into account, compute the length of exons betwee
		int GetGeneBlockEffectiveCoverage( int geneBlockId, int start, int end )
		{
			std::vector<int> &exonIds = geneBlocks[ geneBlockId ].exonBlockIds ;
			int size = exonIds.size() ;
			int i ;
			int ret = 0 ;
			for ( i = 0 ; i < size ; ++i )
			{
				int s = exonBlocks[ exonIds[i] ].start ;
				int e = exonBlocks[ exonIds[i] ].end ;
				if ( s <= start && e >= start )
				{
					ret += e - start + 1 ;
				}
				else if ( start <= s && end >= e )
				{
					ret += e - s + 1 ;
				}
				if ( s <= end && e >= end  )
				{
					ret += end - s + 1 ;
					break ;
				}
			}
			return ret ;
		}

		// Clean up the graph.
		void CleanGeneBlockGraph( Alignments &alignments, Genome &genome )
		{
			int blockCnt = geneBlocks.size() ;
			int i, j, k ;
			bool *validGeneBlock = new bool[blockCnt] ;
			std::vector< struct _geneBlockBubble > bubbles ;

			for ( i = 0 ; i < blockCnt ; ++i )
			{
				int cnt = geneBlockGraph[i].size() ;
				//printf( "%d\n", cnt ) ;

				for ( j = 0 ; j < cnt ; ++j )
				{
					int max = -1 ;
					int maxtag = 0 ;
#ifdef DEBUG
					if ( ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V ) )
					{
						for ( k = 0 ; k < 4 ; ++k )
							printf( "%d(%d) ", geneBlockGraph[i][j].support[k].GetCount(), geneBlockGraph[i][j].support[k].IsGood() ) ;
						printf( "\n" ) ;
					}
#endif

					for ( k = 0 ; k < 4 ; ++k )
					{
						if ( geneBlockGraph[i][j].support[k].IsGood() && geneBlockGraph[i][j].support[k].GetCount() > max )
						{
							max = geneBlockGraph[i][j].support[k].GetCount() ;
							maxtag = k ;	
						}
					}


					if ( max == -1 )
					{
						maxtag = -1 ;
						geneBlockGraph[i][j].supportUse = maxtag ;
						continue ;
					}

					// The best choice should be much better than the other choice
					int max2 = 0 ;
					for ( k = 0 ; k < 4 ; ++k )
					{
						if ( k == maxtag )
							continue ;
						if ( geneBlockGraph[i][j].support[k].IsGood() && geneBlockGraph[i][j].support[k].GetCount() > max2 )
							max2 = geneBlockGraph[i][j].support[k].GetCount() ;
					}

					if ( max2 > 0 && !IsSignificantDifferent( max, fragLength - 2 * readLength, max2, fragLength - 2 * readLength ) )
					{
						maxtag = -1 ;
					}	

					//if ( i == 1084 && geneBlockGraph[i][j].v == 392 )

					geneBlockGraph[i][j].supportUse = maxtag ;
				}
			}

			if ( VERBOSE )
			{
				fprintf( fpOut, "Raw gene block graph\n" ) ;
				OutputGeneBlockGraph( alignments ) ;
			}

			// Remove unqualified edges
			struct _pair *geneBlockInfo = new struct _pair[blockCnt] ;
			for ( i = 0 ; i < blockCnt ; ++i )
			{
				int cnt = geneBlocks[i].exonBlockIds.size() ;
				int len = 0 ;
				int count = 0 ;
				for ( j = 0 ; j < cnt ; ++j )
				{
					k = geneBlocks[i].exonBlockIds[j] ;
					len += exonBlocks[k].end - exonBlocks[k].start + 1 ;
					count += exonBlocks[k].support.GetCount() ;
				}
				assert( j > 0 ) ;
				geneBlockInfo[i].a = count ;
				geneBlockInfo[i].b = len ;
				//printf( "%lf: %d %d\n", avgDepth[i], count, len ) ;

			}

			for ( i = 0 ; i < blockCnt ; ++i )
			{
				int cnt = geneBlockGraph[i].size() ;
				int leni = 0 ;
				if ( cnt == 0 )
					continue ;

				if ( genome.IsOpen() )
				{
					k = geneBlocks[i].exonBlockIds.size() ;
					for ( j = 0 ; j < k ; ++j )
					{
						int ii = geneBlocks[i].exonBlockIds[j] ;
						leni += exonBlocks[ii].end - exonBlocks[ii].start + 1 ;
					}
				}
				
				// Test ambiguous extension
				for ( int test = 0 ; test <= 2 ; test += 2 )
				{
					int max = 0, max2 = 0 ;
					bool valid = true ;
					int maxtag = -1, max2tag ;
					int maxAll = 0, maxAllTag ;
					// Get the max from all the connection	
					for ( j = 0 ; j < cnt ; ++j )
					{
						int su = geneBlockGraph[i][j].supportUse ;
						if ( su == -1 )
							continue ;
						if ( ( su & 2 ) != test )
							continue ;
						
						if ( geneBlockGraph[i][j].support[ su ].GetCount() > maxAll )
						{
							maxAllTag = j ;
							maxAll = geneBlockGraph[i][j].support[ su ].GetCount() ;
						}
					}
					// Get the max from the connection to different scaffold.
					for ( j = 0 ; j < cnt ; ++j )
					{
						int su = geneBlockGraph[i][j].supportUse ;
						if ( su == -1 )
							continue ;
						if ( ( su & 2 ) != test )
							continue ;
						// ignore the connection within the same scaffold
						if ( geneBlocks[ geneBlockGraph[i][j].v ].chrId == geneBlocks[ i ].chrId )
							continue ;

						if ( geneBlockGraph[i][j].support[ su ].GetCount() > max )
						{
							max2 = max ;
							max2tag = maxtag ;
							max = geneBlockGraph[i][j].support[ su ].GetCount() ;
							maxtag = j ;
						}
						else if ( geneBlockGraph[i][j].support[ su ].GetCount() > max2 ) 
						{
							max2 = geneBlockGraph[i][j].support[ su ].GetCount() ;
							max2tag = j ;
						}
					}
					
					// If the connection within the same scaffold is signifcant better,
					// then we should consider it.
					if ( maxAllTag != maxtag )
					{
						if ( max <= 0 || IsSignificantDifferent( max, fragLength - 2 * readLength, maxAll, fragLength - 2 * readLength ) )
						{
							max = maxAll ;
							maxtag = maxAllTag ;
						}
					}
					if ( max > 0 && max2 > 0 && ( !IsSignificantDifferent( max, fragLength - 2 * readLength, max2, fragLength - 2 * readLength )
								) )//|| ( max > 100 && max2 > 100 ) ) )
					{
						// Test whether it is bubble-like
						int bubbleRet = IsSimpleBubble( i, maxtag, max2tag ) ;
						if ( bubbleRet != -1 )
						{
							struct _geneBlockBubble bb ;
							bb.u = i ;
							if ( bubbleRet == max2tag )
							{
								max = max2 ;
								maxtag = max2tag ;

								bb.j1 = max2tag ;
								bb.j2 = maxtag ;
							}
							else
							{
								bb.j1 = maxtag ;
								bb.j2 = max2tag ;
							}
							// Maybe no rescue is better.
							//bubbles.push_back( bb ) ;
								
						}
						else if ( geneBlockGraph[i][maxtag].supportUse != geneBlockGraph[i][max2tag].supportUse 
							|| geneBlocks[ geneBlockGraph[i][maxtag].v ].chrId !=  geneBlocks[ geneBlockGraph[i][max2tag].v ].chrId )
							valid = false ;
					}
					//if ( max == 0 )
					//	valid = false ;
#ifdef DEBUG
					if ( i == DEBUG_U ) //geneBlockGraph[i][j].v == 583 )
						printf( "ambiguous test: %d %d: %d\n", max, max2, valid ) ;
#endif 

					for ( j = 0 ; j < cnt ; ++j )
					{
						int su = geneBlockGraph[i][j].supportUse ;
						if ( su == -1 )
							continue ;
						if ( ( su & 2 ) != test )
							continue ;
						//if ( i == 19160 ) //geneBlockGraph[i][j].v == 583 )
						//	printf( "%d: %d %d %d: %d\n", __LINE__, j, geneBlockGraph[i][j].v, max, geneBlockGraph[i][j].support[ su ].GetCount() ) ;
						//if ( geneBlockGraph[i][j].support[ su ].GetCount() < max )
						if ( max > 0 && j != maxtag )
							geneBlockGraph[i][j].valid = false ;

						//if ( i == 19160 ) //geneBlockGraph[i][j].v == 583 )
						//	printf( "%d\n", geneBlockGraph[i][0].valid ) ;

						if ( valid == false )
							geneBlockGraph[i][j].valid = false ;
						else if ( maxtag != -1 && j != maxtag && su == geneBlockGraph[i][maxtag].supportUse && 
							geneBlocks[ geneBlockGraph[i][j].v ].chrId == geneBlocks[ geneBlockGraph[i][maxtag].v ].chrId)
							geneBlockGraph[i][maxtag].support[ geneBlockGraph[i][maxtag].supportUse ].Add( geneBlockGraph[i][j].support[ su ] ) ;
							
					}
					//if ( i == 19160 ) //geneBlockGraph[i][j].v == 583 )
					//	printf( "%d\n", geneBlockGraph[i][0].valid ) ;

					// If this gene block is too short, and contains an ambiguous extension
					// Then remove it totally
					//if ( i == 10314 )
					//	printf( "%d %d %d\n", valid, max, max2 ) ;
					if ( valid == false && leni < 2 * readLength )
					{
						for ( j = 0 ;j < cnt ; ++j )
						{
							geneBlockGraph[i][j].valid = false ;
						}

					}
				}
			}

			// Rescue the connection from bubbles
			int bubbleCnt = bubbles.size() ;
			for ( i = 0 ; i < bubbleCnt ; ++i )
			{
				int bu = bubbles[i].u ;
				int bj1 = bubbles[i].j1 ;
				int bj2 = bubbles[i].j2 ;
			
				// First, determine whether the connections from the bubbles are clean
				if ( geneBlockGraph[bu].size() >= 5 
					|| geneBlockGraph[ geneBlockGraph[ bu][ bj1 ].v ].size() >= 5 
					|| geneBlockGraph[ geneBlockGraph[ bu][ bj2 ].v ].size() >= 5 ) 
					continue ;	

				int u, j1 ;
				for ( k = 0 ; k < 4 ; ++k )
				{
					int size ;			
					if ( k == 0 )
					{
						// u -> j1 ;
						u = bu ;
						j1 = bj1 ;
					}
					else if ( k == 1 ) 
					{
						// j1 -> u
						u = geneBlockGraph[bu][ bj1 ].v ;
						size = geneBlockGraph[u].size() ;
						for ( j = 0 ; j < size ; ++j )
							if ( geneBlockGraph[u][j].v == bu )
								break ;
						if ( j >= size )
							j1 = -1 ;
						else
							j1 = j ;
					}
					else if ( k == 2 )
					{
						// j1->j2 
						u = geneBlockGraph[bu][bj1].v ;
						size = geneBlockGraph[u].size() ;
						for ( j = 0 ; j < size ; ++j )
							if ( geneBlockGraph[u][j].v == bj2 )
								break ;
						if ( j >= size )
							j1 = -1 ;
						else
							j1 = j ;
					}
					else if ( k == 3 )
					{
						// j2->j1
						u = geneBlockGraph[bu][bj2].v ;
						size = geneBlockGraph[u].size() ;
						for ( j = 0 ; j < size ; ++j )
							if ( geneBlockGraph[u][j].v == bj1 )
								break ;
						if ( j >= size )
							j1 = -1 ;
						else
							j1 = j ;
					}

					if ( j1 == -1 )
						continue ;
					int su = geneBlockGraph[u][j1].supportUse ;
					size = geneBlockGraph[u].size() ;
					for ( j = 0 ; j < size ; ++j )
					{
						if ( geneBlockGraph[u][j].valid == true 
								&& ( geneBlockGraph[u][j].supportUse & su & 2 ) != 0 )
							break ;
					}
					if ( j >= size )
						geneBlockGraph[u][j1].valid = true ;
				}
			}

			for ( i = 0 ; i < blockCnt ; ++i )
			{
				int cnt = geneBlockGraph[i].size() ;
				int leni = 0 ;
				if ( cnt == 0 )
					continue ;

				std::map<uint64_t, int> kmers ;
				if ( genome.IsOpen() )
				{
					k = geneBlocks[i].exonBlockIds.size() ;
					for ( j = 0 ; j < k ; ++j )
					{
						int ii = geneBlocks[i].exonBlockIds[j] ;
						genome.AddKmer( geneBlocks[i].chrId, exonBlocks[ii].start, exonBlocks[ii].end, kmerSize, kmers ) ;
						leni += exonBlocks[ii].end - exonBlocks[ii].start + 1 ;
					}
				}
				
				for ( j = 0 ; j < cnt ; ++j )
				{
					geneBlockGraph[i][j].semiValid = true ;

					//if ( geneBlockGraph[i][j].valid == false )
					//	continue ;
					bool valid = geneBlockGraph[i][j].valid ;

					// In aggressive mode, we allow the ambiguous extension
					if ( aggressiveMode == true )
						valid = true ;			

					//if ( i == 1084 && geneBlockGraph[i][j].v == 392 )
					/*if ( i == 19160 )
					  {
					  printf( "%d %d\n", geneBlockGraph[i][j].v, geneBlockGraph[i][j].supportUse ) ;
					  }*/
#ifdef DEBUG
					if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
						printf( "1: %d\n", valid ) ;
#endif
					// The types of support is ambiguous
					int v = geneBlockGraph[i][j].v ;
					if ( geneBlockGraph[i][j].supportUse == -1 )
					{
						valid = false ;
					}

#ifdef DEBUG
					if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
						printf( "2: %d\n", valid ) ;
#endif

					//printf( "%d %d\n", i, geneBlockGraph[i][j].v ) ;
					if ( valid == true && geneBlockGraph[i][j].support[ geneBlockGraph[i][j].supportUse ].GetCount() < minimumSupport )
						valid = false ;

					//geneBlockGraph[i][j].valid = valid ;
					//continue ;
#ifdef DEBUG
					if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
						printf( "3: %d\n", valid ) ;
#endif



					// Test whether the two side has too much different expression level
					if ( 0 )
					{
						if ( valid == true && IsSignificantDifferent_SimpleTest( geneBlockInfo[i].a, geneBlockInfo[i].b, 
									geneBlockInfo[v].a, geneBlockInfo[v].b ) )
							valid = false ;
					}


#ifdef DEBUG
					if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
						printf( "4: %d\n", valid ) ;
#endif
					// Test the strand
					if ( valid == true && geneBlocks[i].support.GetStrand() != 0 && geneBlocks[v].support.GetStrand() != 0 )
					{
						int su = geneBlocks[i].support.GetStrand() ;
						int sv = geneBlocks[v].support.GetStrand() ;
						int supportUse = geneBlockGraph[i][j].supportUse ;

						if ( su == sv && ( supportUse == 0 || supportUse == 3 ) )
							valid = false ;
						else if ( su != sv && ( supportUse == 1 || supportUse == 2 ) )
							valid = false ;
					}

#ifdef DEBUG
					if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
						printf( "5: %d\n", valid ) ;
#endif

					// TODO: a better way to decide this threshold
					// Test whether the connection is too few comparing with the expression level
					if ( 0 ) //valid == true && fragLength > 2 * readLength )
					{
						int su = geneBlockGraph[i][j].supportUse ;
						double c = geneBlockGraph[i][j].support[ su ].GetCount() / (double)( fragLength - readLength )  ;
						int cnt = geneBlockGraph[i][j].support[ su ].GetCount() ;
						/*if ( ( i == 11500 && geneBlockGraph[i][j].v == 37818 )
						  || ( i == 37818 && geneBlockGraph[i][j].v == 11500  ) )
						  {
						  printf( "%d %d: %lf\n", i, geneBlockGraph[i][j].v, c ) ;
						  }*/
						if ( IsSignificantDifferent_SimpleTest( geneBlockInfo[i].a, geneBlockInfo[i].b, cnt, fragLength - 2 * readLength ) ||
								IsSignificantDifferent_SimpleTest( geneBlockInfo[v].a, geneBlockInfo[v].b, cnt, fragLength - 2 * readLength ) )
						{
							//printf( "%lf\n", c) ;
							//printf( "hi\n" ) ;
							if ( !( cnt / ( fragLength - 2 * readLength ) > 10 && !geneBlockGraph[i][j].support[ su ].IsUnique() ) )
								valid = false ;
						}
					}
#ifdef DEBUG
					if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
						printf( "6: %d: %lf %lf: %d %d\n", valid,  
							(double)geneBlockInfo[i].a / geneBlockInfo[i].b,
							(double)geneBlockInfo[v].a / geneBlockInfo[v].b, 
							geneBlockGraph[i][j].support[ geneBlockGraph[i][j].supportUse ].GetCount(), ( fragLength - 2 * readLength ) ) ;
#endif

					// Gene family - use read to fileter
					if ( valid == true )
					{
						int fi = GetFather( i, repeatFather ) ;
						int fj = GetFather( geneBlockGraph[i][j].v, repeatFather ) ;
						if ( fi == fj )
							valid = false ;

						/*if ( ( i == 11500 && geneBlockGraph[i][j].v == 37818 )
						  || ( i == 37818 && geneBlockGraph[i][j].v == 11500  ) )
						  {
						  printf( "%d %d: %d\n", i, geneBlockGraph[i][j].v, valid ) ;
						  }*/
					}
					// Gene family - use kmer to filter
					if ( genome.IsOpen() && valid == true )
					{
						int v = geneBlockGraph[i][j].v ;
						int m, n ;
						n = geneBlocks[v].exonBlockIds.size() ;
						int lenj = 0 ;
						int kmerCoverage = 0 ;

						for ( m = 0 ; m < n ; ++m )
						{
							int ii = geneBlocks[v].exonBlockIds[m] ;
							//genome.AddKmer( geneBlocks[v].chrId, exonBlocks[ii].start, exonBlocks[ii].end, kmerSize, kmersV ) ;
							kmerCoverage += genome.GetKmerCoverage( geneBlocks[v].chrId, exonBlocks[ii].start, exonBlocks[ii].end, kmerSize, kmers ) ; 
							lenj += exonBlocks[ii].end - exonBlocks[ii].start + 1 ;
						}

						//cnt = genome.CompareKmerSets( kmers, kmersV ) ;
						/*printf( "%d %d: (%s: %d-%d) (%s: %d-%d): %d %d %d\n", i, v, alignments.GetChromName( geneBlocks[i].chrId ), geneBlocks[i].start, geneBlocks[i].end,
						  alignments.GetChromName( geneBlocks[v].chrId ), geneBlocks[v].start, geneBlocks[v].end,
						  cnt, leni, lenj ) ;*/

						//if ( cnt > 10 || ( cnt > 1 && ( cnt > 0.1 * leni || cnt > 0.1 * lenj ) ) )
						int su = geneBlockGraph[i][j].supportUse ;
						int span = GetGeneBlockEffectiveCoverage(i, geneBlockGraph[i][j].support[su].GetLeftMostPos() 
								,geneBlockGraph[i][j].support[su].GetRightMostPos() ) + readLength - kmerSize + 1 ;
						if ( ( kmerCoverage > int( 1.5 * kmerSize ) && ( !geneBlockGraph[i][j].support[ su ].IsUnique() ) ) || 
							( kmerCoverage >= kmerSize + 2 && ( kmerCoverage > 0.1 * leni || kmerCoverage > 0.1 * lenj || kmerCoverage > 30 * span / readLength ) ) )
						{
							valid = false ;
						}
#ifdef DEBUG
						if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
							printf( "7: %d: %d %d\n", valid, kmerCoverage, span ) ;
#endif
					}

					if (  valid == true )
					{
						// Test whether the connection is from the middle of two gene blocks. 
						// We expect to see the read close to the boundary.
						// Actually, the test of insert size should take care of this, but the test here makes it more stringent
						int su = geneBlockGraph[i][j].supportUse ;

						if ( su & 2 )
						{
							// should close to left boundary
							struct _block &eblock = exonBlocks[ geneBlocks[i].exonBlockIds[0] ] ;
							int boundary = eblock.leftSplice ;
							if ( boundary == -1 )
								boundary = eblock.start ;
							int pos = geneBlockGraph[i][j].support[su].GetLeftMostPos() ;

							if ( pos > boundary + readLength )
								geneBlockGraph[i][j].semiValid = false ;
						}
						else
						{
							// should close to right boundary
							struct _block &eblock = exonBlocks[ geneBlocks[i].exonBlockIds[ geneBlocks[i].exonBlockIds.size() - 1 ] ] ;
							int boundary = eblock.rightSplice ;
							if ( boundary == -1 )
								boundary = eblock.end ;
							int pos = geneBlockGraph[i][j].support[su].GetRightMostPos() ;

							if ( pos < boundary - readLength )
								geneBlockGraph[i][j].semiValid = false ;

						}

						// If all the support mates are from the same position
						if ( geneBlockGraph[i][j].support[su].GetCoordCnt() == 1 )
							geneBlockGraph[i][j].semiValid = false ;
					}

					geneBlockGraph[i][j].valid = valid ;
#ifdef DEBUG
					if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
						printf( "8: %d\n", valid ) ;
#endif
				}
			}

			// Remove the one-side edge that makes the graph non-symmetric
			for ( i = 0 ; i < blockCnt ; ++i )
			{
				int cnt = geneBlockGraph[i].size() ;
				for ( j = 0 ; j < cnt ; ++j )
				{
					if ( geneBlockGraph[i][j].valid == true )
					{
						int v = geneBlockGraph[i][j].v ;
						int cnt = geneBlockGraph[v].size() ;
						for ( k = 0 ; k < cnt ; ++k )
						{
							if ( geneBlockGraph[v][k].v == i )
							{
								int su = geneBlockGraph[v][k].supportUse ;
								int tmp = ( su >> 1 ) | ( ( su & 1 ) << 1 ) ;
								if ( tmp == geneBlockGraph[i][j].supportUse )
									break ;
								else
									k = cnt ;
							}
						}

						if ( k >= cnt )
							geneBlockGraph[i][j].valid = false ;
#ifdef DEBUG
						if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
							printf( "9: %d\n", geneBlockGraph[i][j].valid ) ;
#endif
					}
				}
			}

			//printf( "valid 8: %d\n", geneBlockGraph[85185][0].valid ) ;

			// Mark the invalid edges caused by semiValid flag
			for ( i = 0 ; i < blockCnt ; ++i )
			{
				int cnt = geneBlockGraph[i].size() ;
				for ( j = 0 ; j < cnt ; ++j )
				{
					if ( geneBlockGraph[i][j].valid == true && geneBlockGraph[i][j].semiValid == false )
					{
						int v = geneBlockGraph[i][j].v ;
						int cnt = geneBlockGraph[v].size() ;
						bool otherSemiValid = true ;
						for ( k = 0 ; k < cnt ; ++k )
						{
							if ( geneBlockGraph[v][k].v == i )
							{
								otherSemiValid = geneBlockGraph[v][k].semiValid ;
							}
						}

						if ( geneBlockGraph[i][j].semiValid == false && otherSemiValid == false )
							geneBlockGraph[i][j].valid = false ;
					}
#ifdef DEBUG
					if ( i == DEBUG_U && geneBlockGraph[i][j].v == DEBUG_V )
						printf( "10: %d\n", geneBlockGraph[i][j].valid ) ;
#endif
				}
			}
			//printf( "valid 9: %d\n", geneBlockGraph[85185][0].valid ) ;

			// Make sure the graph is symmetric
			for ( i = 0 ; i < blockCnt ; ++i )
			{	
				int cnt = geneBlockGraph[i].size() ;
				for ( j = 0 ; j < cnt ; ++j )
				{
					if ( geneBlockGraph[i][j].valid == false )
					{
						int v = geneBlockGraph[i][j].v ;
						int cnt = geneBlockGraph[v].size() ;
						for ( k = 0 ; k < cnt ; ++k )
						{
							if ( geneBlockGraph[v][k].v == i )
							{
								geneBlockGraph[v][k].valid = false ;
							}
						}
					}
				}
			}
			//printf( "valid 10: %d\n", geneBlockGraph[85185][0].valid ) ;

			// Remove unqualified single-exons or short gene blocks that is not in the middle
			for ( i = 0 ; i < blockCnt ; ++i )
			{
				validGeneBlock[i] = true ;

				int len = GetEffectiveLength( i ) ; ;
				if ( ( geneBlocks[i].exonBlockIds.size() == 1 && !geneBlocks[i].support.IsUnique() ) ||
						( len < minimumEffectiveLength && !geneBlocks[i].support.IsUnique() ) ) //|| len <= 100 )
				{
					k = geneBlockGraph[i].size() ;
					int cnt[2] = {0, 0} ;
					for ( j = 0 ; j < k ; ++j )
					{
						if ( geneBlockGraph[i][j].valid == false )
							continue ;

						if ( geneBlockGraph[i][j].supportUse == -1 )
							continue ;
						if ( geneBlockGraph[i][j].supportUse & 2 )
							++cnt[1] ;
						else
							++cnt[0] ;
					}
					//if ( i == 10746 || i == 10745 )
					//	printf( "%d %d %d\n", len, cnt[0], cnt[1] ) ;
					if ( cnt[0] == 0 || cnt[1] == 0 )
					{
						geneBlockGraph[i].clear() ;
						validGeneBlock[i] = false ;
					}
				}

				// If the two blocks are short and they are the only connection
				/*if ( len < minimumEffectiveLength )
				  {
				  k = geneBlockGraph[i].size() ;
				  if ( k == 1 )
				  {
				//int u = geneBlockGraph[i][0].u ;
				int v = geneBlockGraph[i][0].v ;

				if ( geneBlockGraph[v].size() == 1 )
				{
				int lenv = GetEffectiveLength( v ) ;
				if ( lenv < minimumEffectiveLength )
				validGeneBlock[i] = false ; 
				}
				}
				}*/

				/*if ( len < minimumEffectiveLength )
				  {
				  geneBlockGraph[i].clear() ;
				  validGeneBlock[i] = false ;
				  }*/
			}

#ifdef DEBUG
			printf( "11: validGeneBlock[i]=%d\n", validGeneBlock[ DEBUG_U ] ) ;
#endif

			// Clean up the graph
			for ( i = 0 ; i < blockCnt ; ++i )
			{
				if ( !validGeneBlock[i] )
				{
					geneBlockGraph[i].clear() ;
					continue ;
				}
				int cnt = geneBlockGraph[i].size() ;
				for ( j = 0, k = 0 ; j < cnt ; ++j )
				{
					if ( geneBlockGraph[i][j].valid && validGeneBlock[ geneBlockGraph[i][j].v ] ) 
					{
						geneBlockGraph[i][k] = geneBlockGraph[i][j] ;
						++k ;
					}
				}
				for ( j = cnt - 1 ; j >= k ; --j )
				{
					geneBlockGraph[i].pop_back() ;
				}
			}
			delete[] geneBlockInfo ;
			delete[] validGeneBlock ;

			if ( VERBOSE )
			{
				fprintf( fpOut, "Gene block graph:\n" ) ;
				OutputGeneBlockGraph( alignments ) ;
			}
		}

		bool IsSignificantDifferent( double a, double b )
		{
			if ( a - 6 * sqrt(a) <= b  
					&& b <= a + 6 * sqrt( a ) )
			{
				return false ;
			}

			if ( b - 6 * sqrt( b ) <= a && a <= b + 6 * sqrt( b ) )
			{
				return false ;
			}

			return true ;
		}
		bool IsSignificantDifferent_SimpleTest( int cntA, int lenA, int cntB, int lenB ) 
		{
			double rA = (double)cntA/ lenA ;
			double rB = (double)cntB/ lenB ;

			return IsSignificantDifferent( rA, rB ) ;
		}

		bool IsSignificantDifferent( int cntA, int lenA, int cntB, int lenB ) 
		{
			//return IsSignificantDifferent( (double)cntA / lenA, (double)cntB / lenB ) ;

			// The code below is wrong. fixed by log transformation
			double rA = (double)cntA/ lenA ;
			double rB = (double)cntB/ lenB ;

			rA = sqrt( rA ) ;
			rB = sqrt( rB ) ;

			//double rNULL = (double)(rA * lenA + rB * lenB ) / ( lenA + lenB ) ;
			double rNULL = sqrt( (double)( cntA + cntB ) / ( lenA + lenB ) ) ;

			/*double cntNA = rNULL * lenA ;
			  double cntNB = rNULL * lenB ;


			  double chi = ( cntA - cntNA ) * ( cntA - cntNA ) / cntNA + 
			  ( cntB - cntNB ) * ( cntB - cntNB ) / cntNB ;*/

			double cA = rA * lenA ;
			double cB = rB * lenB ;
			double eA = rNULL * lenA ;
			double eB = rNULL * lenB ;
			double chi = ( cA - eA ) * ( cA - eA ) / eA + ( cB - eB ) * ( cB - eB ) / eB ;
			//double chi = ( rA - rNULL ) * ( rA - rNULL ) / rNULL + ( rB - rNULL ) * ( rB - rNULL ) / rNULL ;
			//bool tmp = IsSignificantDifferent( (double)cntA / lenA, (double)cntB / lenB ) ;
			bool ret = false ;
			if ( chi > 3.841 )
				ret = true ;
			//printf( "%d %d %d %d: %d\n", cntA, lenA, cntB, lenB, ret ) ;
			return ret ;
		}

		void OutputGeneBlockGraph( Alignments &alignments )
		{
			int i, j ;
			int bcnt = geneBlocks.size() ;

			for ( i = 0 ; i < bcnt ; ++i )
			{
				fprintf( fpOut, "%d %s %d (%"PRId64" %"PRId64"): ", i, alignments.GetChromName( geneBlocks[i].chrId ), geneBlocks[i].contigId, geneBlocks[i].start + 1, geneBlocks[i].end + 1 ) ;
				int cnt = geneBlockGraph[i].size() ;
				for ( j = 0 ; j < cnt ; ++j )
				{
					int v = geneBlockGraph[i][j].v ;
					int s ;
					if ( geneBlockGraph[i][j].supportUse == -1 )
						s = -1 ;
					else
						s = geneBlockGraph[i][j].support[ geneBlockGraph[i][j].supportUse ].GetCount() ;
					fprintf( fpOut, "(%d %s %d (%"PRId64" %"PRId64")):(%d %d) ", v, alignments.GetChromName( geneBlocks[v].chrId ), geneBlocks[v].contigId,
							geneBlocks[v].start + 1, geneBlocks[v].end + 1,
							geneBlockGraph[i][j].supportUse, s ) ;
				}
				fprintf( fpOut, "\n" ) ;
			}
		}
} ;

#endif
