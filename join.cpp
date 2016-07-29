// Use the connections to create scaffold
// Li Song
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <vector>

#include "defs.h"
#include "genome.hpp" 
#include "alignments.hpp"
#include "ContigGraph.hpp"

const char *usage = "usage: rascaf_join [OPTIONS]\n"
		"OPTIONS:\n"
		"\t-r STRING: the path to the rascaf output. Can use multiple of -r. (required)\n"
		"\t-o STRING: the prefix of the output file. (default: rascaf_scaffold)\n"
		"\t-ms INT: minimum support alignments for the connection (default: 2)\n"
		"\t-ignoreGap: ignore the gap size, which do not consider the number of Ns between contigs (default: not used)\n";

bool VERBOSE = false ;
int minSupport = 2 ;
bool ignoreGap = false ; // if ignore gaps, then we will ignore the gap size, which will be aggressive
FILE *fpOut = NULL ;
int MAX_NEIGHBOR ;
int breakN ;

struct _part
{
	int chrId ;
	int contigId ;
	char strand ;
} ;

struct _contigList
{
	int contigId ;
	int strand ;

	int prevContigId ;
	int nextContigId ;

	int contigInScaf ;
} ;

std::vector< std::vector<struct _part> > connects ;

int CompScaffold( const void *p1, const void *p2 )
{
	struct _pair *a = (struct _pair*)p1 ;
	struct _pair *b = (struct _pair*)p2 ;
	if ( a->b < b->b )
		return 1 ;
	else if ( a->b > b->b )
		return -1 ;
	else
		return 0 ;
}

int AddConnection( char *s, Alignments &alignments, std::vector<struct _part> &parts ) 
{
	int n ;
	int i ; 
	int cnt = 0 ;
	for ( i = 0 ; s[i] != ':'; ++i )
	{
		n = n * 10 + s[i] - '0' ;
	}
	
	//printf( "%s\n", s ) ;
	while ( 1 )
	{
		if ( s[i] == '(')
		{
			struct _part np ;
			int tag ;
			int stage = 0 ; 
			tag = i + 1 ;
			for ( i = i + 1 ; s[i] != ')' ; ++i )
			{
				if ( s[i] == ' ' )
				{
					if ( stage == 0 )
					{
						s[i] = '\0' ;
						np.chrId = alignments.GetChromIdFromName( &s[tag] ) ;
						s[i] = ' ' ;
						tag = i + 1 ;
					}
					else if ( stage == 1 )
					{
						tag = i + 1 ;
					}
					else if ( stage == 2 )
					{
						s[i] = '\0' ;
						np.contigId = atoi( &s[tag] ) ;
						//printf( "%d\n", np.contigId ) ;
						s[i] = ' ' ;
					}
					++stage ;
				}
			}
			np.strand = s[i - 1] ;

			parts.push_back( np ) ;	
		}
		else if ( s[i] == '\0' )
			break ;
		++i ;
	}
	return n ;
}

void ForwardSearch( int u, int inDummy, int time, int *visitTime, int *counter, int *visitDummy, ContigGraph &contigGraph )
{
	if ( visitTime[u] == 2 * time )
		return ;
	visitTime[u] = 2 * time ;
	counter[u] = 1 ;
	visitDummy[u] = inDummy ;
	struct _pair *buffer = new struct _pair[ MAX_NEIGHBOR ] ;
	int ncnt ;
	int i ;
	
	ncnt = contigGraph.GetNeighbors( u, 1 - inDummy, buffer, MAX_NEIGHBOR ) ;
	for ( i = 0 ; i < ncnt ; ++i )
	{
		/*if ( time == 639 )
		{
			printf( "forwardsearch: (%d %d)=>(%d %d)\n", u, 1 - inDummy, buffer[i].a, buffer[i].b ) ;
		}*/
		ForwardSearch( buffer[i].a, buffer[i].b, time, visitTime, counter, visitDummy, contigGraph ) ;	
	}
	delete[] buffer ;
}
void BackwardSearch( int u, int inDummy, int time, int *visitTime, int *counter, ContigGraph &contigGraph, int chosenNodes[], int &chosenCnt )
{
	if ( visitTime[u] == 2 * time + 1 )
		return ;
	if ( visitTime[u] != 2 * time )
		counter[u] = 0 ;
	visitTime[u] = 2 * time + 1 ;
	++counter[u] ;
	if ( counter[u] == 2 )
	{
		chosenNodes[ chosenCnt ] = u ;
		++chosenCnt ;
	}
	struct _pair *buffer = new struct _pair[ MAX_NEIGHBOR ];
	int ncnt ;
	int i ;
	ncnt = contigGraph.GetNeighbors( u, 1 - inDummy, buffer, MAX_NEIGHBOR ) ;
	for ( i = 0 ; i < ncnt ; ++i )
		BackwardSearch( buffer[i].a, buffer[i].b, time, visitTime, counter, contigGraph, chosenNodes, chosenCnt ) ;
	delete[] buffer ;
}

void BackwardSearchForTriangularCycle( int u, int inDummy, int time, int *visitTime, int *counter, int * visitDummy, ContigGraph &contigGraph, int chosenNodes[], int &chosenCnt )
{
	if ( visitTime[u] == 2 * time + 1 )
		return ;
	if ( visitTime[u] == 2 * time && visitDummy[u] == inDummy )
	{
		visitTime[u] = 2 * time + 1 ;
		chosenNodes[ chosenCnt ] = u ;
		++chosenCnt ;
		return ;
	}
	visitTime[u] = 2 * time + 1 ;
	
	struct _pair *buffer = new struct _pair[ MAX_NEIGHBOR ];
	int ncnt ;
	int i ;
	
	ncnt = contigGraph.GetNeighbors( u, 1 - inDummy, buffer, MAX_NEIGHBOR ) ;
	for ( i = 0 ; i < ncnt ; ++i )
		BackwardSearchForTriangularCycle( buffer[i].a, buffer[i].b, time, visitTime, counter, visitDummy, contigGraph, chosenNodes, chosenCnt ) ;	
	delete[] buffer ;
}

// Search one dangling paths
void SearchDangling( int u, int inDummy, bool *used, int time, int *visitTime, ContigGraph &contigGraph, bool add, int chosenNodes[], int chosenDummyNodes[], int &chosenCnt, Genome &genome )
{
	//if ( u == 10 )
	//	printf( "hi\n" ) ;
	if ( visitTime[u] == time )
		return ;
	visitTime[u] = time ;
	if ( add ) // Only one path will be added.
	{
		chosenNodes[ chosenCnt ] = u ;
		chosenDummyNodes[ chosenCnt ] = inDummy ;
		++chosenCnt ;
	}
	//if ( u == 10246 )
	//	printf( "hi2 %d\n", chosenCnt ) ;
	struct _pair *neighbors = new struct _pair[ MAX_NEIGHBOR ] ;
	int ncnt ;
	int i ;
	ncnt = contigGraph.GetNeighbors( u, 1 - inDummy, neighbors, MAX_NEIGHBOR ) ; 
	bool newAdd = true ;
	//if ( u == 10 )
	//	printf( "hi3 %d %d\n", ncnt, 1 - inDummy ) ;
	
	// Firstly, bias towards a connection onto another scaffold
	for ( i = 0 ; i < ncnt ; ++i )		
	{
		if ( used[ neighbors[i].a ] )
			continue ;
		if ( genome.GetChrIdFromContigId( u ) != genome.GetChrIdFromContigId( neighbors[i].a ) ) 
		{
			SearchDangling( neighbors[i].a, neighbors[i].b, used, time, visitTime, 
				contigGraph, newAdd, chosenNodes, chosenDummyNodes, chosenCnt, genome ) ;
			newAdd = false ;
		}
	}

	// Lastly, bias towards the closest contig in the raw assembly, since we know the connections are on the same scaffold
	int direction ;
	if ( inDummy == 1 )
		direction = -1 ;
	else
		direction = 1 ;

	int min = genome.GetContigCount() + 1 ;
	int mintag = -1 ;

	for ( i = 0 ; i < ncnt ; ++i )
	{
		//if ( u == 10 )
		//	printf( "%d\n", i ) ;
		if ( used[ neighbors[i].a ] )
			continue ;
		//if ( u == 34674 && neighbors[i].a == 144159 )
		//	printf( "hi\n" ) ;
		//SearchDangling( neighbors[i].a, neighbors[i].b, used, time, visitTime, contigGraph, newAdd, chosenNodes, chosenDummyNodes, chosenCnt, genome ) ;
		//newAdd = false ;
		int tmp = direction * ( neighbors[i].a - u ) ;
		if ( tmp < min )
		{
			min = tmp ;
			mintag = i ;
		}
	}

	if ( mintag != -1 )
	{
		SearchDangling( neighbors[ mintag ].a, neighbors[ mintag ].b, used, time, visitTime, 
			contigGraph, newAdd, chosenNodes, chosenDummyNodes, chosenCnt, genome ) ;
		newAdd = false ;
	}
	delete[] neighbors ;
}

int main( int argc, char *argv[] )
{
	Alignments alignments ;
	Genome genome ;
	std::vector<int> rascafFileId ; 
	
	char line[2048] ;
	char prefix[512] = "rascaf_scaffold" ;
	int rawAssemblyInd = 1 ;
	FILE *rascafFile ;
	bool contigLevel = false ;
	int i ;
	FILE *outputFile ;
	FILE *infoFile ;

	breakN = 1 ;

	if ( argc < 2 )
	{
		fprintf( stderr, "%s", usage ) ;
		exit( 1 ) ;
	}
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-o", argv[i] ) )
		{
			strcpy( prefix, argv[i + 1 ] ) ;
			++i ;
		}
		else if ( !strcmp( "-ms", argv[i] ) )
		{
			minSupport = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-ignoreGap", argv[i] ) )
		{
			ignoreGap = true ;
		}
		else if ( !strcmp( "-r", argv[i] ) )
		{
			rascafFileId.push_back( i + 1 ) ; 
			++i ;
		}
		else
		{
			fprintf( stderr, "Unknown option: %s\n", argv[i] ) ;
			exit( EXIT_FAILURE ) ;
		}

	}

	if ( rascafFileId.size() == 0 )
	{
		fprintf( stderr, "Must use -r to specify rascaf output file.\n" ) ;
		exit( EXIT_FAILURE ) ;
	}
	
	MAX_NEIGHBOR = 1 + rascafFileId.size() ;
	
	// Get the bam file.
	rascafFile = fopen( argv[ rascafFileId[0] ], "r" ) ;
	while ( fgets( line, sizeof( line ), rascafFile ) != NULL ) 
	{
		if ( strstr( line, "command line:" ) )
		{
			char *p ;
			char buffer[512] ;

			p = strstr( line, "-breakN" ) ;
			if ( p != NULL )
			{
				p += 7 ;
				while ( *p == ' ' )
					++p ;
				for ( i = 0 ; *p && *p != ' ' ; ++p, ++i )
					buffer[i] = *p ;
				buffer[i] = '\0' ;
				breakN = atoi( buffer ) ;
			}

			p = strstr( line, "-b" ) ;
			if ( p == NULL )
			{
				fprintf( stderr, "Could not find the bam file specified by -b in Rascaf.\n" ) ;
				exit( 1 ) ;
			}
			p += 2 ;
			while ( *p == ' ' )
				++p ;
			for ( i = 0 ; *p && *p != ' ' ; ++p, ++i )
				buffer[i] = *p ;
			buffer[i] = '\0' ;
			alignments.Open( buffer ) ;

			p = strstr( line, "-f") ;
			if ( p == NULL )
			{
				fprintf( stderr, "Could not find the raw assembly file specified by -f in Rascaf.\n" ) ;
				exit( 1 ) ;
			}
			p += 2 ;
			while ( *p == ' ' )
				++p ;
			for ( i = 0 ; *p && *p != ' ' && *p != '\n' ; ++p, ++i )
				buffer[i] = *p ;
			buffer[i] = '\0' ;
			fprintf( stderr, "Found raw assembly file: %s\n", buffer ) ;
			genome.Open( alignments, buffer ) ;
			
			break ;
		}
	}
	fclose( rascafFile ) ;
	// Parse the input.
	for ( unsigned int fid = 0 ; fid < rascafFileId.size() ; ++fid )
	{
		rascafFile = fopen( argv[ rascafFileId[fid] ], "r" ) ;
		bool start = false ;
		int tag ;
		while ( fgets( line, sizeof( line ), rascafFile ) != NULL ) 
		{
			if ( strstr( line, "command line:" ) )
			{
				start = true ;
				if ( strstr( line, "-f" ) )
				{
					contigLevel = true ;
				}

				continue ;
			}

			if ( !start )
				continue ;

			if ( !strcmp( line, "WARNINGS:\n" ) )
				break ;

			std::vector<struct _part> nparts  ;
			if ( line[0] >= '0' && line[0] <= '9' )
			{
				AddConnection( line, alignments, nparts ) ;
				connects.push_back( nparts ) ;
				tag = 0 ;
			}
			else if ( line[0] == '\t' || line[0] == ' ' )
			{
				// Break the nparts if the support is too low.
				int num = 0 ;
				for ( i = 0 ; line[i] < '0' || line[i] > '9' ; ++i )
					;
				for ( ; line[i] >= '0' && line[i] <= '9' ; ++i )
					num = num * 10 + line[i] - '0' ;
				++tag ;
				if ( num < minSupport )
				{
					nparts = connects.back() ;
					connects.pop_back() ;
					int size = nparts.size() ;

					std::vector<struct _part> newNParts ;
					for ( i = 0 ; i < tag ; ++i )
						newNParts.push_back( nparts[i] ) ;
					if ( newNParts.size() > 1 )
						connects.push_back( newNParts ) ;

					newNParts.clear() ;
					for ( ; i < size ; ++i  ) 
						newNParts.push_back( nparts[i] ) ;
					if ( newNParts.size() > 1 )
						connects.push_back( newNParts ) ;

					tag = 0 ;
				}
			}
		}
		fclose( rascafFile ) ;
	}
	if ( contigLevel == false )
	{
		genome.SetIsOpen( contigLevel ) ;
	}

	// Build the graph
	int contigCnt = genome.GetContigCount() ;
	int edgeCnt = 0 ;
	int csize = connects.size() ;

	for ( i = 0 ; i < csize ; ++i )
		edgeCnt += connects[i].size() ;

	ContigGraph contigGraph( contigCnt, contigCnt + edgeCnt ) ;
	for ( i = 0 ; i < contigCnt - 1 ; ++i )
	{
		if ( genome.GetChrIdFromContigId( i ) == genome.GetChrIdFromContigId( i + 1 ) )
		{
			contigGraph.AddEdge( i, 1, i + 1, 0 ) ;
		}
	}
	for ( i = 0 ; i < csize ; ++i )	
	{
		std::vector<struct _part> &parts = connects[i] ;
		int size = parts.size() ;
		for ( int j = 0 ; j < size - 1 ; ++j )
		{
			struct _part &a = parts[j] ;
			struct _part &b = parts[j + 1] ;
			
			// Two dummy nodes for each contig. Left is 0, right is 1
			int dummyU = 0 ;
			int dummyV = 0 ;
			if ( a.strand == '+' )
				dummyU = 1 ;
			if ( b.strand == '-' )
				dummyV = 1 ;
			contigGraph.AddEdge( a.contigId, dummyU, b.contigId, dummyV, true ) ;
		}
	}

	// Check the cycles in the contig graph. This may introduces when combining different rascaf outputs.
	int *visitTime = new int[contigCnt] ;
	struct _pair *neighbors = new struct _pair[ MAX_NEIGHBOR ] ;

	bool *isInCycle = new bool[contigCnt] ;
	std::vector<int> cycleNodes ;
	memset( visitTime, -1, sizeof( int ) * contigCnt ) ;
	memset( isInCycle, false, sizeof( bool ) * contigCnt ) ;
	for ( i = 0 ; i < contigCnt ; ++i )
	{
		if ( isInCycle[i] )
			continue ;
		if ( contigGraph.IsInCycle( i, cycleNodes, visitTime ) )	
		{
			int cnt = cycleNodes.size() ;
			//printf( "===\n") ;
			for ( int j = 0 ; j < cnt ; ++j )
			{
				//printf( "In cycle %d\n", cycleNodes[j] ) ;
				isInCycle[ cycleNodes[j] ] = true ;
			}
		}
	}
	//exit( 1 ) ; 
	// Remove the connected edges involving the nodes in the cycle
	for ( i = 0 ; i < contigCnt ; ++i )
	{
		if ( isInCycle[i] )
		{	
			for ( int dummy = 0 ; dummy <= 1 ; ++dummy )
			{
				int ncnt = contigGraph.GetNeighbors( i, dummy, neighbors, MAX_NEIGHBOR ) ;
				for ( int j = 0 ; j < ncnt ; ++j )
				{
					if ( neighbors[j].a == i + 2 * dummy - 1 && neighbors[j].b != dummy 
						&& genome.GetChrIdFromContigId( i ) == genome.GetChrIdFromContigId( neighbors[j].a ) )
						continue ; // the connection created by the raw assembly
					else
						contigGraph.RemoveEdge( i, dummy, neighbors[j].a, neighbors[j].b ) ;
				}
			}
		}
	}
	//delete[] isInCycle ;
	//printf( "hi: %d %d\n", __LINE__, contigCnt ) ;
	//printf( "%d %d\n", contigGraph.GetNeighbors( 0, 0, neighbors, MAX_NEIGHBOR ), contigGraph.GetNeighbors( 0, 1, neighbors, MAX_NEIGHBOR ) ) ;
	// Sort the scaffolds from fasta file, so that longer scaffold come first
	int scafCnt = genome.GetChrCount() ;
	struct _pair *scafInfo = new struct _pair[scafCnt] ;
	memset( scafInfo, -1, sizeof( struct _pair) * scafCnt ) ;
	for ( i = 0 ; i < contigCnt ; ++i )	
	{
		int chrId = genome.GetChrIdFromContigId( i ) ;
		if ( scafInfo[chrId].a == -1 )
		{
			scafInfo[ chrId ].a = i ;
			scafInfo[ chrId ].b = genome.GetChrLength( chrId ) ;
		}
	}
	qsort( scafInfo, scafCnt, sizeof( struct _pair ), CompScaffold ) ;

	// Merge the branches and build the scaffold
	ContigGraph scaffold( contigCnt, 2 * contigCnt ) ;
	
	// Use a method similar to topological sort
	bool *used = new bool[contigCnt] ;
	int *degree = new int[2 *contigCnt] ;
	int *danglingVisitTime = new int[contigCnt] ;
	int *counter = new int[contigCnt] ;
	int *visitDummy = new int[ contigCnt ] ;
	int *buffer = new int[contigCnt] ;
	int *buffer2 = new int[contigCnt] ;
	bool *isInQueue = new bool[ contigCnt ] ;
	int *chosen = new int[contigCnt] ;
	int chosenCnt ;

	memset( isInCycle, false, sizeof( bool ) * contigCnt ) ;
	memset( visitTime, -1, sizeof( int ) * contigCnt ) ;	
	memset( visitDummy, -1, sizeof( int ) * contigCnt ) ;	
	memset( counter, -1, sizeof( int ) * contigCnt ) ;

	// Use those memory to remove triangular cycles
	for ( i = 0 ; i < scafCnt ; ++i )
	{
		int from, to ;
		if ( scafInfo[i].a == -1 )
			continue ;
		genome.GetChrContigRange( genome.GetChrIdFromContigId( scafInfo[i].a ), from, to ) ;
		ForwardSearch( from, 0, i, visitTime, counter, visitDummy, contigGraph ) ;
		chosenCnt = 0 ;
		BackwardSearchForTriangularCycle( to, 1, i, visitTime, counter, visitDummy, contigGraph, chosen, chosenCnt ) ;

		for ( int j = 0 ; j < chosenCnt ; ++j )
		{
			//printf( "%d\n", chosen[j] ) ;
			isInCycle[ chosen[j] ] = true ;
		}
	}

	for ( i = 0 ; i < contigCnt ; ++i )
	{
		if ( isInCycle[i] )
		{	
			for ( int dummy = 0 ; dummy <= 1 ; ++dummy )
			{
				int ncnt = contigGraph.GetNeighbors( i, dummy, neighbors, MAX_NEIGHBOR ) ;
				for ( int j = 0 ; j < ncnt ; ++j )
				{
					if ( neighbors[j].a == i + 2 * dummy - 1 && neighbors[j].b != dummy 
						&& genome.GetChrIdFromContigId( i ) == genome.GetChrIdFromContigId( neighbors[j].a ) )
						continue ; // the connection created by the raw assembly
					else
						contigGraph.RemoveEdge( i, dummy, neighbors[j].a, neighbors[j].b ) ;
				}
			}
		}
	}


	memset( used, false, sizeof( bool ) * contigCnt ) ;
	memset( visitTime, -1, sizeof( int ) * contigCnt ) ;	
	memset( visitDummy, -1, sizeof( int ) * contigCnt ) ;	
	memset( danglingVisitTime, -1, sizeof( int ) * contigCnt ) ;
	memset( counter, -1, sizeof( int ) * contigCnt ) ;
	memset( isInQueue, false, sizeof( bool ) * contigCnt ) ;
	ContigGraph newGraph( contigCnt, edgeCnt ) ;		
	
	// Compute the gap size
	int *gapSize = new int[contigCnt] ;
	for ( i = 0 ; i < contigCnt - 1 ; ++i )
	{	
		if ( genome.GetChrIdFromContigId( i ) == genome.GetChrIdFromContigId( i + 1 ) )
		{
			struct _contig c1 = genome.GetContigInfo( i ) ;
			struct _contig c2 = genome.GetContigInfo( i + 1 ) ;
			gapSize[i] = c2.start - c1.end - 1 ;
		}
		else
			gapSize[i] = -1 ;
	}

	// Start search
	int ncnt ;
	struct _pair *queue = new struct _pair[ contigCnt ] ;
	int head = 0, tail ;
	int danglingTime = 0 ;
		
	// Pre-allocate the subgraph.
	ContigGraph subgraph( contigCnt, 3 * contigCnt ) ; 
	for ( i = 0 ; i < scafCnt ; ++i )
	{
		//if ( used[144281] == true )
		//	printf( "changed %d %d\n", i, scafInfo[i - 1].a ) ;
		if ( scafInfo[i].a == -1 )
			continue ;
		int from, to ;
		genome.GetChrContigRange( genome.GetChrIdFromContigId( scafInfo[i].a ), from, to ) ;
		//printf( "%d: %d %d %d\n", i, scafInfo[i].b, from, to ) ;
		ForwardSearch( from, 0, i, visitTime, counter, visitDummy, contigGraph ) ;
		chosenCnt = 0 ;
		BackwardSearch( to, 1, i, visitTime, counter, contigGraph, chosen, chosenCnt ) ;

		/*printf( "%s %d (%d %d) %d\n", alignments.GetChromName( genome.GetChrIdFromContigId( scafInfo[i].a ) ), i, from, to, chosenCnt ) ;
		if ( chosenCnt > 1 )
		{
			printf( "=== " ) ;
			for ( int j = 0 ; j < chosenCnt ; ++j )
				printf( "%d ", chosen[j] ) ;
			printf( "\n" ) ;
		}*/

		for ( int j = 0 ; j < chosenCnt ; ++j )
		{
			ncnt = contigGraph.GetNeighbors( chosen[j], 0, neighbors, MAX_NEIGHBOR ) ;
			//printf( "%d %d %d: %d %d %d\n", j, chosen[j], ncnt, neighbors[0].a, visitTime[ neighbors[0].a ], 
			//	counter[neighbors[0].a ] ) ;
			for ( int k = 0 ; k < ncnt ; ++k )
			{
				//if ( i == 639 )
				//	printf( "Neighbor from 0 %d: %d %d\n", k, neighbors[k].a, neighbors[k].b ) ;
				if ( visitTime[ neighbors[k].a ] == 2 * i + 1 && counter[neighbors[k].a ] == 2 )
				{
					subgraph.AddEdge( chosen[j], 0, neighbors[k].a, neighbors[k].b, true ) ;
					//printf( "subgraph: (%d %d)=>(%d %d)\n", chosen[j], 0, neighbors[k].a, neighbors[k].b ) ;
				}
			}

			ncnt = contigGraph.GetNeighbors( chosen[j], 1, neighbors, MAX_NEIGHBOR ) ;
			for ( int k = 0 ; k < ncnt ; ++k )
			{
				//if ( i == 639 )
				//	printf( "Neighbor from 1 %d: %d %d\n", k, neighbors[k].a, neighbors[k].b ) ;
				if ( visitTime[ neighbors[k].a ] == 2 * i + 1 && counter[neighbors[k].a ] == 2 )
				{
					subgraph.AddEdge( chosen[j], 1, neighbors[k].a, neighbors[k].b, true ) ;
					//printf( "subgraph: (%d %d)=>(%d %d)\n", chosen[j], 1, neighbors[k].a, neighbors[k].b ) ;
				}
			}
		}

		// Initialize the degree counter
		for ( int j = 0 ; j < chosenCnt ; ++j )
		{
			for ( int l = 0 ; l < 2 ; ++l )
			{
				/*if ( i == 6145 )
				{
					std::vector<struct _pair> neighbors ;
					ncnt = subgraph.GetNeighbors( chosen[j], l, neighbors ) ;
					printf( "%d ncnt=%d\n", l, ncnt ) ;
				}*/

				ncnt = subgraph.GetNeighbors( chosen[j], l, neighbors, MAX_NEIGHBOR ) ;
				degree[ 2 * chosen[j] + l ] = ncnt ;
			} 
		}

		// "topological" sort
		head = 0 ;
		isInQueue[from] = true ;
		queue[0].a = from ;
		queue[0].b = 0 ;
		tail = 1 ;
		int prevTag = -1 ;
		int *prevAdd = buffer ; // reuse counter to save some memory.
		int *nextAdd = buffer2 ;
		int firstAdd = -1 ;

		while ( head < tail )
		{
			int tailTag = tail ;
			for ( int j = head ; j < tailTag ; ++j )
			{
				nextAdd[j] = -1 ;
				if ( !used[ queue[j].a ] )	
				{
					used[ queue[j].a ] = true ;
					if ( prevTag != -1 )
					{
						scaffold.AddEdge( queue[ prevTag].a, 1 - queue[prevTag].b, queue[j].a, queue[j].b ) ;
						nextAdd[ prevTag ] = j ;
					
						/*if ( i == 639 )
							printf( "(%lld %lld)=>(%lld %lld)\n", queue[ prevTag].a, 1 - queue[prevTag].b, queue[j].a, queue[j].b ) ;*/
					}
					else
						firstAdd = j ;
				
					prevTag = j ;
				}
				prevAdd[j] = prevTag ; // the most recent(<=) queue id when added to scaffold. 

				ncnt = subgraph.GetNeighbors( queue[j].a, 1 - queue[j].b, neighbors, MAX_NEIGHBOR ) ;
				for ( int k = 0 ; k < ncnt ; ++k )
				{
					--degree[ 2 * neighbors[k].a + neighbors[k].b ] ;
					if ( degree[ 2 * neighbors[k].a + neighbors[k].b ] == 0 && !isInQueue[neighbors[k].a] )
					{
						isInQueue[ neighbors[k].a ] = true ;
						queue[ tail ] = neighbors[k] ; // Interesting assignment, I think.
						++tail ;
						/*if ( i == 639 )
							printf( "pushed in queue: %d\n", neighbors[k].a ) ;*/
						// Put the consecutive contigs together.
						struct _pair testNeighbors[ MAX_NEIGHBOR ] ;
						struct _pair tag ;
						tag = neighbors[k] ;
						while ( 1 )  
						{
							if ( contigGraph.GetNeighbors( tag.a, 1 - tag.b, testNeighbors, MAX_NEIGHBOR ) != 1 )
								break ;
							int n = subgraph.GetNeighbors( tag.a, 1 - tag.b, testNeighbors, MAX_NEIGHBOR ) ;
							if ( n != 1 )
								break ;
							//printf( "%d %d\n", n, testNeighbors[0].a ) ;

							struct _pair backNeighbors[ MAX_NEIGHBOR ] ;
							if ( contigGraph.GetNeighbors( testNeighbors[0].a, testNeighbors[0].b, backNeighbors, MAX_NEIGHBOR ) != 1 )
								break ;
							n = subgraph.GetNeighbors( testNeighbors[0].a, testNeighbors[0].b, backNeighbors, MAX_NEIGHBOR ) ;
							if ( n != 1 )
								break ;
							isInQueue[ testNeighbors[0].a ] = true ;
							queue[tail] = testNeighbors[0] ;
							++tail ;
							/*if ( i == 639 )
								printf( "pushed in queue: %d\n", testNeighbors[0].a ) ;*/
							tag = testNeighbors[0] ;
						}
					}
				}
			}

			head = tailTag ;
		}
		// Remove the effect on the subgraph. 
		/*if ( tail != chosenCnt )
		{
			printf( "WARNING: not matched\n" ) ;
			exit( 1 ) ;
		}*/
		for ( int j = 0 ; j < tail ; ++j )
		{
			visitDummy[ queue[j].a ] = -1 ;
			counter[ queue[j].a ] = -1 ;
			subgraph.RemoveAdjacentEdges( queue[j].a ) ;
			isInQueue[ queue[j].a ] = false ;
		}
		subgraph.ResetEdgeUsed() ;

		// no point is picked
		if ( prevTag == -1 )
		{
			continue ;
		}

		// Update the gap size
		prevTag = -1 ;
		for ( int j = 0 ; j < tail - 1 ; ++j )
		{
			if ( genome.GetChrIdFromContigId( queue[j].a ) == genome.GetChrIdFromContigId( from ) )
				prevTag = queue[j].a ;
			else if ( prevTag != -1 )
			{
				struct _contig c = genome.GetContigInfo( queue[j].a ) ;
				gapSize[prevTag] -= ( c.end - c.start + 1) ;
			}
		}
		// Add the dangling contigs. Use the fact that the queue holding the contigs in the same order as in the scaffold.
		// 5'->3' dangling
		int *chosenDummy = degree ; 
		for ( int j = tail - 1 ; j >= 0 ; --j )
		{
			//if ( j < tail - 1 ) 
			//	continue ;
			chosenCnt = 0 ;
			//if ( queue[j].a == 0 )
			//	printf( "Dummy: %d %d %d\n", j, queue[j].b, 1 - queue[j].b ) ;
			SearchDangling( queue[j].a, queue[j].b, used, danglingTime, danglingVisitTime, contigGraph, false, chosen, chosenDummy, chosenCnt, genome ) ;
			++danglingTime ;
			int prevTag = prevAdd[j] ;
			/*if ( queue[j].a == 0 )
			{
				struct _pair neighbors[5] ;
				int ncnt = contigGraph.GetNeighbors( queue[j].a, 1 - queue[j].b, neighbors, 5 ) ;
				printf( "%d %d %d %d: %d %d\n", queue[j].b, chosenCnt, prevTag, ncnt, neighbors[0].a, used[ neighbors[0].a ] ) ;
			}*/	
			if ( prevTag == -1 )
				break ;

			// Trim the dangling list
			int k = chosenCnt - 1 ;
			if ( j > 0 && j < tail - 1 )
			{
				for ( k = chosenCnt - 1 ; k >= 1 ; --k )
					if ( genome.GetChrIdFromContigId( chosen[k] ) != genome.GetChrIdFromContigId( chosen[k - 1] ) )
						break ;
			}
			
			// Test the gap size
			int len = 0 ;
			for ( int l = 0 ; l <= k ; ++l )
			{
				struct _contig c = genome.GetContigInfo( chosen[k] ) ;
				len += c.end - c.start + 1 ;
			}

			if ( j < tail - 1 )
			{
				int l ;
				for ( l = j ; l >= 0 ; --l )
					if ( genome.GetChrIdFromContigId( queue[l].a ) == genome.GetChrIdFromContigId( from ) )
						break ;

				if ( !ignoreGap && len >= gapSize[ queue[l].a ] + 100 )
					continue ;
				else
					gapSize[ queue[l].a ] -= len ;
			}


			for ( ; k >= 0 ; --k )
			{
				used[ chosen[k] ] = true ;
				//printf( "Dangling 1: %d=>%d\n", queue[prevTag].a, chosen[k] ) ;
				scaffold.InsertNode( queue[ prevTag ].a, 1 - queue[ prevTag ].b, chosen[k], chosenDummy[k] ) ;
			}
		}

		// 3'->5' dangling
		for ( int j = 0 ; j < tail ; ++j )
		{
			//if ( j > 0 )
			//	continue ;
			chosenCnt = 0 ;
			SearchDangling( queue[j].a, 1 - queue[j].b, used, danglingTime, danglingVisitTime, contigGraph, false, chosen, chosenDummy, chosenCnt, genome ) ;
			++danglingTime ;

			int prevTag = prevAdd[j] ;
			int nextTag ;
			if ( prevTag == -1 || j <= firstAdd )
				nextTag = firstAdd ;
			else if ( j == prevTag )
				nextTag = j ;
			else 
				nextTag = nextAdd[ prevTag ] ;
			if ( nextTag == -1 ) 
				break ;
			/*if ( queue[j].a == 37549 )
			{
				struct _pair neighbors[5] ;
				int ncnt = contigGraph.GetNeighbors( queue[j].a, queue[j].b, neighbors, 5 ) ;
				fprintf( stderr, "%d %d %d: %d %d %d: %d %d %d\n", j, queue[j].a, queue[j].b, chosenCnt, nextTag, ncnt, chosen[0], chosenDummy[0], used[ chosen[0] ] ) ;
			}*/
			
			// trim the danling list
			int k = chosenCnt - 1 ;
			if ( j < tail - 1 && j > 0 )
			{
				for ( k = chosenCnt - 1 ; k >= 1 ; --k )
					if ( genome.GetChrIdFromContigId( chosen[k] ) != genome.GetChrIdFromContigId( chosen[k - 1] ) )
						break ;
			}

			// Test the gap size
			int len = 0 ;
			for ( int l = 0 ; l <= k ; ++l )
			{
				struct _contig c = genome.GetContigInfo( chosen[k] ) ;
				len += c.end - c.start + 1 ;
			}

			if ( j > 0 )
			{
				int l ;
				for ( l = j - 1 ; l >= 0 ; --l ) // Notice the j-1 here, because we want the gap strictly before current contig
					if ( genome.GetChrIdFromContigId( queue[l].a ) == genome.GetChrIdFromContigId( from ) )
						break ;

				if ( !ignoreGap && len >= gapSize[ queue[l].a ] + 100 )
					continue ;
				else
					gapSize[ queue[l].a ] -= len ;
			}

			for ( ; k >= 0 ; --k )
			{
				used[ chosen[k] ] = true ;
				scaffold.InsertNode( queue[nextTag].a, queue[nextTag].b, chosen[k], chosenDummy[k] ) ;
				//printf( "Dangling 2: %d<=%d\n", queue[nextTag].a, chosen[k] ) ;
				//if ( chosen[k] == 10246 )
				//	printf( "hi %d %d %d %d\n", j, queue[j].a, k, chosen[k] ) ;
			}
		}
	}
	//return 0 ;
	
	// Output the scaffold
	int id = 0 ;
	char infoFileName[512] ;
	char outputFileName[512] ;
	sprintf( infoFileName, "%s.info", prefix ) ;
	sprintf( outputFileName, "%s.fa", prefix ) ;

	outputFile = fopen( outputFileName, "w" ) ;
	infoFile = fopen( infoFileName, "w") ;

	memset( used, false, sizeof( bool ) * contigCnt ) ;
	for ( i = 0 ; i < contigCnt ; ++i )
	{
		//printf( "%d (%s)\n", i, alignments.GetChromName( genome.GetChrIdFromContigId( i ) )  ) ; fflush( stdout ) ;
		/*if ( i == 10246 )
		{
			std::vector<struct _pair> neighbors ;
			scaffold.GetNeighbors( i, 0, neighbors ) ;
			printf( "%u\n", neighbors.size() ) ;
		}*/
		if ( used[i] )
			continue ;
		int ncnt1 = scaffold.GetNeighbors( i, 0, neighbors, MAX_NEIGHBOR ) ;
		int ncnt2 = scaffold.GetNeighbors( i, 1, neighbors, MAX_NEIGHBOR ) ;
		if ( ncnt1 == 0 || ncnt2 == 0 ) // The end of a scaffold
		{
			fprintf( outputFile, ">scaffold_%d\n", id) ;
			fprintf( infoFile, ">scaffold_%d", id ) ;
			++id ;
			int p = i ;
			int dummyP = 1 ;
			if ( ncnt1 == 0 )
				dummyP = 0 ;
			
			used[i] = true ;
			genome.PrintContig( outputFile, i, dummyP ) ;
			fprintf( infoFile, " (%s %d %c)", alignments.GetChromName( genome.GetChrIdFromContigId( p ) ), p, dummyP == 0 ? '+' : '-' ) ;
			while ( 1 )
			{
				ncnt = scaffold.GetNeighbors( p, 1 - dummyP, neighbors, MAX_NEIGHBOR ) ;		
				if ( ncnt == 0 )
					break ;
				// ncnt must be 1
				int insertN = 17 ;
				
				if ( genome.GetChrIdFromContigId( p ) == genome.GetChrIdFromContigId( neighbors[0].a ) )
				{
					struct _contig cp, cna ;
					cp = genome.GetContigInfo( p ) ;
					cna = genome.GetContigInfo( neighbors[0].a ) ;
					if ( p < neighbors[0].a )
						insertN = cna.start - cp.end - 1 ;
					else if ( p > neighbors[0].a )
						insertN = cp.start - cna.end - 1 ;

				}

				p = neighbors[0].a ;
				dummyP = neighbors[0].b ;
				for ( int j = 0 ; j < insertN ; ++j )	
					fprintf( outputFile, "N" ) ;
				used[p] = true ;
				genome.PrintContig( outputFile, p, dummyP ) ;

				fprintf( infoFile, " (%s %d %c)", alignments.GetChromName( genome.GetChrIdFromContigId( p ) ), p, dummyP == 0 ? '+' : '-' ) ;
			}
			fprintf( outputFile, "\n" ) ;
			fprintf( infoFile, "\n" ) ;
		}
	}

	for ( i = 0 ; i < contigCnt ; ++i ) 
		if ( !used[i] )
		{
			fprintf( stderr, "Unreported contig %d.\n", i ) ;
		}
	fclose( outputFile ) ;
	fclose( infoFile ) ;

	delete[] buffer ;
	delete[] buffer2 ;
	delete[] chosen ;
	delete[] queue ;
	delete[] counter ;
	delete[] visitTime ;
	delete[] used ;
	delete[] scafInfo ;
	delete[] isInQueue ;
	delete[] gapSize ;

	//fclose( rascafFile ) ;
	return 0 ;
}
