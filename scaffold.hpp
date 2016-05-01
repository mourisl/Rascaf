// The class predicting the connections
// Li Song

#ifndef _LSONG_RSCAF_SCAFFOLD_HEADER
#define _LSONG_RSCAF_SCAFFOLD_HEADER

#include <vector>
#include <algorithm>

#include <inttypes.h>

#include "defs.h"
#include "blocks.hpp"
#include "alignments.hpp"
#include "ContigGraph.hpp"

extern char *prefix ;
extern bool outputConnectionSequence ;

struct _scaffold
{
	int contigId ;
	int flip ;
} ;

struct _blockScaffoldWrapper
{
	int support ;
	std::vector<struct _mateEdge> *s ; // Use point so the swapping when sorting can be efficiently done
} ;

struct _misassembledInfo
{
	int u, v ;

	// For type 3, we change the location of the blocks and flip the orientation of v
	int type ; // 0-deletion, 1-transversion u, 2-transversion v, 3-translocation, 4-misassembly in scaffold 
} ;

// Conceptually, each chromosome has two dummy nodes
// left node(0) leads to the 3'->5' neighbour(means need flip), right node(1) says the 5'->3' neighbor
// neighborType tells the information of the neighbor. 
// each node corresponds to a contig
struct _scaffoldNode
{
	//int contigId ;
	int neighbor[2] ;  // The id of the neighbor is mapped chr id.
	int neighborType[2] ;

	struct _mateEdge mateEdges[2] ; // The information of the edges we use
} ;

bool CompEdge( struct _mateEdge e1, struct _mateEdge e2 )
{
	return e1.support[ e1.supportUse ].GetCount() > e2.support[ e2.supportUse ].GetCount() ; 
}

bool CompScaffold( struct _blockScaffoldWrapper s1, struct _blockScaffoldWrapper s2 )
{
	return s1.support > s2.support ;	
}

class Scaffold
{
private:
	Blocks &blocks ;
	Genome &genome ;

	std::vector< std::vector<int> > components ;
	std::vector<bool> blockUsed ;

	void SearchComponent( int tag, bool *visited, std::vector<int> &list )
	{
		if ( visited[tag] )
			return ;
		visited[tag] = true ;
		list.push_back( tag ) ;

		std::vector< struct _mateEdge > &edges = blocks.geneBlockGraph[tag] ;
		std::vector<struct _mateEdge>::iterator it ;
		for ( it = edges.begin() ; it != edges.end() ; ++it )
		{
			if ( !visited[ it->v ] )
				SearchComponent( it->v, visited, list ) ;
		}
	}

	std::vector< struct _blockScaffoldWrapper > blockScaffolds ; // a-id of the block, b-direction
	
	int contigCnt ;

	struct _scaffoldNode *scaffoldNodes ;

	std::vector< struct _mateEdge > problematicMates ; 

	struct _misassembledInfo *misassembledInfo ; // record which two gene blocks resulting in reporting the misassembly
	std::vector< std::vector<int> > misassembledCycles ;
	std::vector<int> misassembledContigs ;

	int GetFather( int tag, int *father )
	{
		if ( father[tag] == tag )
			return tag ;
		return father[tag] = GetFather( father[tag], father ) ;
	}
	
	// Find a circle 
	bool FindContigGraphCircle( int startGeneBlock, int currGeneBlock, int timeStamp, int *visitTime, 
		std::vector<int> *geneBlocksInContig, std::vector<int> *geneBlockScaffoldGraph )
	{
		int j ;
		visitTime[currGeneBlock] = timeStamp ;
		int currContigId = blocks.geneBlocks[ currGeneBlock ].contigId ;
		int cnt = geneBlockScaffoldGraph[currGeneBlock].size() ;

		for ( j = 0 ; j < cnt ; ++j )
		{
			if ( visitTime[ geneBlockScaffoldGraph[ currGeneBlock ][j] ] == timeStamp )
				continue ;

			if ( currContigId != blocks.geneBlocks[ startGeneBlock ].contigId && 
				geneBlockScaffoldGraph[ currGeneBlock ][j] != startGeneBlock &&
				blocks.geneBlocks[ geneBlockScaffoldGraph[ currGeneBlock ][j] ].contigId == blocks.geneBlocks[ startGeneBlock ].contigId )
			{
				misassembledInfo[ blocks.geneBlocks[ startGeneBlock ].contigId ].u = startGeneBlock ;
				misassembledInfo[ blocks.geneBlocks[ startGeneBlock ].contigId ].v = geneBlockScaffoldGraph[ currGeneBlock ][j] ;
				misassembledInfo[ blocks.geneBlocks[ startGeneBlock ].contigId ].type = 0 ;
				return true ;
			}
			if ( FindContigGraphCircle( startGeneBlock, geneBlockScaffoldGraph[ currGeneBlock ][j], timeStamp, visitTime, 
					geneBlocksInContig, geneBlockScaffoldGraph ) )
				return true ;
		}

		if ( currContigId != blocks.geneBlocks[ startGeneBlock ].contigId )
		{
			cnt = geneBlocksInContig[ currContigId ].size() ;
			for ( j = 0 ; j < cnt ; ++j )
			{
				int tag = geneBlocksInContig[ currContigId ][j] ;
				if ( tag == currGeneBlock )
					continue ;

				if ( visitTime[tag] == timeStamp )
					continue ;

				if ( FindContigGraphCircle( startGeneBlock, tag, timeStamp, visitTime, 
							geneBlocksInContig, geneBlockScaffoldGraph ) )
					return true ;
			}
		}
		return false ;
	}

	void FindMisassemblies()
	{
		std::vector<int> *geneBlockScaffoldGraph = new std::vector<int>[ blocks.geneBlocks.size() ] ; // This graph is undirected
									// the index are for the gene block 
		std::vector<int> *geneBlocksInContig = new std::vector<int>[ contigCnt ] ; // The other gene blocks in the same contig

		int bscafCnt = blockScaffolds.size() ;
		int *visitTime = new int[ blocks.geneBlocks.size() ] ; 
		misassembledInfo = new struct _misassembledInfo[ contigCnt ] ;	
		int i, j, k ;
		int edgeCnt = 0 ;

		memset( misassembledInfo, -1, sizeof( *misassembledInfo ) * contigCnt ) ;
		memset( visitTime, -1, sizeof( int ) * blocks.geneBlocks.size() ) ;
		
		// First, try to find the misassemblies that involves other contigs
		for ( i = 0 ; i < bscafCnt ; ++i )
		{
			int cnt = blockScaffolds[i].s->size() ;
			std::vector<struct _mateEdge> &s = *( blockScaffolds[i].s ) ;
			edgeCnt += cnt ;

			for ( j = 0 ; j < cnt ; ++j )
			{
				// Since the gene block scaffold are disjoint paths,
				// there is no need to worry about repeat edges
				geneBlockScaffoldGraph[ s[j].u ].push_back( s[j].v ) ;
				geneBlockScaffoldGraph[ s[j].v ].push_back( s[j].u ) ;
			}
		}

		k = blocks.geneBlocks.size() ;
		for ( i = 0 ; i < k ; ++i )
			geneBlocksInContig[ blocks.geneBlocks[i].contigId ].push_back( i ) ;


		for ( i = 0 ; i < k ; ++i )
		{
			FindContigGraphCircle( i, i, i, visitTime, geneBlocksInContig, geneBlockScaffoldGraph ) ;
			/*{
				//printf( "Misassembly from gene block: %lld-%lld", blocks.geneBlocks[i].start, blocks.geneBlocks[i].end ) ;
				//isMisassembled[ blocks.geneBlocks[i].contigId ] = true ;
			}*/
		}

		// Then, find the misassemblies that can directly determined by the geneblocks within the contig
		for ( i = 0 ; i < bscafCnt ; ++i )	
		{
			int cnt = blockScaffolds[i].s->size() ;
			std::vector<struct _mateEdge> &s = *( blockScaffolds[i].s ) ;

			for ( j = 1 ; j < cnt ; ++j )
			{
				if ( blocks.geneBlocks[ s[j].u ].contigId != blocks.geneBlocks[ s[j].v ].contigId )
					continue ;
				// Test the direction
				if ( s[j].supportUse == 0 )
				{
					misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].u = s[j].u ;
					misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].v = s[j].v ;

					if ( blocks.geneBlocks[ s[j].u ].start <= blocks.geneBlocks[ s[j].v ].start )
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].type = 2 ;
					else
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].type = 1 ;
						
				}
				else if ( s[j].supportUse == 3 )
				{
					misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].u = s[j].u ;
					misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].v = s[j].v ;

					if ( blocks.geneBlocks[ s[j].u ].start <= blocks.geneBlocks[ s[j].v ].start )
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].type = 2 ;
					else
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].type = 1 ;
						
				}
				else if ( s[j].supportUse == 1 )
				{
					// from left to right
					if ( blocks.geneBlocks[ s[j].u ].start > blocks.geneBlocks[ s[j].v ].start )
					{
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].u = s[j].u ;
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].v = s[j].v ;
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].type = 3 ;
					}
				}
				else // supportUsed == 2
				{
					// from right to left
					if ( blocks.geneBlocks[ s[j].u ].start < blocks.geneBlocks[ s[j].v ].start )
					{
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].u = s[j].u ;
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].v = s[j].v ;
						misassembledInfo[ blocks.geneBlocks[ s[j].u ].contigId ].type = 3 ;
					}
				}
			}
		}

		delete[] visitTime ;
		delete[] geneBlockScaffoldGraph ;
		delete[] geneBlocksInContig ;
	}

public:
	Scaffold( Blocks &in, Genome &inGenome ):blocks(in), genome( inGenome )
	{ 
		for ( std::vector< struct _geneBlock >::iterator it = blocks.geneBlocks.begin() ; it != blocks.geneBlocks.end() ; ++it )
		{
			blockUsed.push_back( false ) ;
		}
		misassembledInfo = NULL ;
	} 
	//Scaffold( ) {} 
	~Scaffold() 
	{
		int i ;
		int cnt = blockScaffolds.size() ;
		for ( i = 0 ; i < cnt ; ++i )
			delete blockScaffolds[i].s ;

		if ( misassembledInfo != NULL )
			delete[] misassembledInfo ;
	}
	
	int BuildComponent()	
	{
		int i, j ;	
		int geneBlockCnt = blocks.geneBlocks.size() ;
		bool *visited = new bool[geneBlockCnt] ;

		for ( i = 0 ; i < geneBlockCnt ; ++i )
			visited[i] = false ;

		for ( i = 0 ; i < geneBlockCnt ; ++i )
		{
			if ( visited[i] )
				continue ;
			std::vector<int> newComponent ;

			SearchComponent( i, visited, newComponent ) ;
			if ( newComponent.size() <= 1 )
				continue ;
			components.push_back( newComponent ) ;
		}
		delete[] visited ;

		return components.size() ;
	}

	// scaffold one component of gene block graph.
	void ScaffoldComponent( int componentId )
	{
		std::vector<int> &component = components[ componentId ] ;

		int nodeCnt = component.size() ;
		int edgeCnt ;
		int i, j, k ;
		//printf( "%s %d\n", __func__, nodeCnt ) ;
		// A greedy heuristic method
		// TODO: enumeration method for simple component
		std::vector< struct _mateEdge > edges ;

		for ( i = 0 ; i < nodeCnt ; ++i )
		{
			int cnt = blocks.geneBlockGraph[ component[i] ].size() ;
			for ( j = 0 ; j < cnt ; ++j )
			{
				edges.push_back( blocks.geneBlockGraph[ component[i] ][j] ) ;
			}
		}

		std::sort( edges.begin(), edges.end(), CompEdge ) ;
		edgeCnt = edges.size() ;

		while ( 1 )
		{
			// Find the anchor
			int tmp ;
			for ( i = 0 ; i < edgeCnt ; ++i )		
			{
				if ( blockUsed[ edges[i].u ] || blockUsed[ edges[i].v ] )
					continue ;
				break ;
			}
			if ( i >= edgeCnt )
				break ;

			blockUsed[ edges[i].u ] = true ;
			blockUsed[ edges[i].v ] = true ;

			// Extend the anchor until the support drops a lot or connect with other scaffoled part
			std::vector<struct _mateEdge> chain[2] ;
			chain[1].push_back( edges[i] ) ; // chain2 are used to search towards right.
			
			tmp = blocks.geneBlockGraph[ edges[i].v ].size() ;
			for ( j = 0 ; j < tmp ; ++j )
			{
				if ( blocks.geneBlockGraph[ edges[i].v ][j].v == edges[i].u )
				{
					chain[0].push_back( blocks.geneBlockGraph[ edges[i].v ][j] ) ;
					break ;
				}
			}
			//printf( "%d %d\n", edges[i].u, edges[i].v ) ;
			assert( chain[0].size() > 0 ) ;
			//printf( "%d %d\n", edges[i].supportUse, chain[0][0].supportUse ) ;
			// Extend.
			for ( k = 0 ; k < 2 ; ++k )
			{
				while ( 1 )
				{
					int csize = chain[k].size() ;
					int u ;

					u = chain[k][csize - 1].v ;

					int cnt = blocks.geneBlockGraph[u].size() ;
					int max = 0, maxtag ;
					int max2 = 0 ;
					for ( j = 0 ; j < cnt ; ++j )
					{
						struct _mateEdge &edge = blocks.geneBlockGraph[u][j] ;
						// The first part makes sure we are exiting the current block.
						if ( ( edge.supportUse / 2 != ( chain[k][ csize - 1 ].supportUse & 1 ) ) )
						{ 
							int c = edge.support[ edge.supportUse ].GetCount() ;
							if ( c >= max )
							{
								max2 = max ;
								max = edge.support[ edge.supportUse ].GetCount() ;
								maxtag = j ;
							}
							else if ( c > max2 )
								max2 = c ;
						}
					}

					if ( max == 0 )
						break ;
					// We already make sure the extension is unambiguous when cleaning the graph.
					//if ( max2 != 0 && !blocks.IsSignificantDifferent( max, 100, max2, 100 ) ) // The extension is ambiguous
					//	break ;

					struct _mateEdge &edge = blocks.geneBlockGraph[u][maxtag] ;
					if ( blockUsed[ edge.v ]  )
						break ;
					if ( edge.support[ edge.supportUse ].GetCount() < chain[k][csize - 1].support[ chain[k][csize - 1].supportUse ].GetCount() / 50 ) // stop if the support drops too much
					//if ( blocks.IsSignificantDifferent( edge.support[ edge.supportUse ].GetCount(), 100, chain[k][csize - 1].support[ chain[k][csize - 1].supportUse ].GetCount(), 100 ) )
						break ;
					blockUsed[ blocks.geneBlockGraph[u][maxtag].v ] = true ;
					chain[k].push_back( edge ) ;
				}
			}
			// Concatenate chain1 and chain2.
			// Add the whole chain to the scaffolding list.
			struct _blockScaffoldWrapper scaffoldW ;
			scaffoldW.s = new std::vector< struct _mateEdge > ;
			std::vector< struct _mateEdge > *scaffold = scaffoldW.s ;
			int cnt = chain[0].size() ;
			for ( i = cnt - 1 ; i >= 0 ; --i )
			{
				scaffold->push_back( chain[0][i] ) ;
			}
			cnt = chain[1].size() ;
			for ( i = 1 ; i < cnt ; ++i )
			{
				scaffold->push_back( chain[1][i] ) ;
			}
			
			// If this just connect two single-exon gene blocks
			// TODO: make it less stringent
			/*if ( scaffold->size() == 1 )
			{
				int u = edges[i].u ;
				int v = edges[i].v ;
				if ( blocks.geneBlocks[u].exonBlockIds.size() == 1  &&
					blocks.geneBlocks[v].exonBlockIds.size() == 1 )
				{
					delete scaffoldW.s ;
					continue ;
				}
			}*/
			
			// Compute the support for the scaffold
			// Possible alternative, the sum of supports including those spanning several 
			// blocks
			scaffoldW.support = 100000000 ;
			for ( k = 0 ; k < 2 ; ++k )
			{
				cnt = chain[k].size() ;
	
				for ( i = 0 ; i < cnt ; ++i )
				{
					if ( chain[k][i].support[ chain[k][i].supportUse ].GetCount() < scaffoldW.support )
						scaffoldW.support = chain[k][i].support[ chain[k][i].supportUse ].GetCount() ;
				}
			}
			blockScaffolds.push_back( scaffoldW ) ;
		}
	}

	void ScaffoldGenome()
	{
		int i, j ;
		int bscafCnt ;
		contigCnt = 0 ;

		std::sort( blockScaffolds.begin(), blockScaffolds.end(), CompScaffold ) ;					
		// The sort and mapping can make scaffolds with higher support comes first
		bscafCnt = blockScaffolds.size() ;

		if ( genome.IsOpen() )
		{
			contigCnt = genome.GetContigCount()  ;
		}
		else
		{
			int cnt = blocks.geneBlocks.size() ;
			int max = -1 ;
			for ( i = 0 ; i < cnt ; ++i )
				if ( blocks.geneBlocks[i].contigId > max )
					max = blocks.geneBlocks[i].contigId ;
			contigCnt = max + 1 ;
		}

		/*for ( i = 0 ; i < bscafCnt ; ++i )
		{
			std::vector<struct _mateEdge> &s = *( blockScaffolds[i].s ) ;
			int cnt = s.size() ;
			for ( j = 0 ; j < cnt ; ++j )
				printf( "(%d=>%d) ", s[j].u, s[j].v ) ;
			printf( "\n" ) ;
		}*/

		FindMisassemblies() ;	


		scaffoldNodes = new struct _scaffoldNode[contigCnt] ;

		for ( i = 0 ; i < contigCnt ; ++i )
		{
			memset( &scaffoldNodes[i], -1, sizeof( scaffoldNodes[i] ) ) ;
		}

		for ( i = 0 ; i < bscafCnt ; ++i )
		{
			// TODO: stop early if the block scaffold is not well
			std::vector<struct _mateEdge> &s = *( blockScaffolds[i].s ) ;

			int cnt = s.size() ;
			for ( j = 0 ; j < cnt ; ++j )
			{
				int u = s[j].u ;
				int v = s[j].v ;

				int contigIdu = blocks.geneBlocks[u].contigId ;
				int contigIdv = blocks.geneBlocks[v].contigId ;
				int dummyNodeu = 1 ;
				int dummyNodev = 1 ;
				if ( s[j].supportUse & 2 )	
					dummyNodeu = 0 ;
				if ( s[j].supportUse & 1 )
					dummyNodev = 0 ;
				
				if ( misassembledInfo[ contigIdu ].u != -1 || misassembledInfo[ contigIdv ].v != -1 || contigIdu == contigIdv )
					continue ;

				if ( scaffoldNodes[ contigIdu ].neighbor[ dummyNodeu ] == -1 &&
					scaffoldNodes[ contigIdv ].neighbor[ dummyNodev ] == -1 )
				{
					//printf( "connect: %d %d %d\n", contigIdu, contigIdv, s[j].support[ s[j].supportUse ].GetCount() ) ;
					scaffoldNodes[ contigIdu ].neighbor[ dummyNodeu ] = contigIdv ;
					scaffoldNodes[ contigIdu ].neighborType[ dummyNodeu ] = dummyNodev ;
					scaffoldNodes[ contigIdu ].mateEdges[ dummyNodeu ] = s[j] ;
			
					scaffoldNodes[ contigIdv ].neighbor[ dummyNodev ] = contigIdu ;
					scaffoldNodes[ contigIdv ].neighborType[ dummyNodev ] = dummyNodeu ;
					scaffoldNodes[ contigIdv ].mateEdges[ dummyNodev ] = s[j] ;
				}
				else
				{
					//printf( "error: %d %d %d\n", contigIdu, contigIdv, s[j].support[ s[j].supportUse ].GetCount() ) ;
					// Break the originally assignment if the support is about the same.
					bool reported = false ;
					if ( scaffoldNodes[ contigIdu ].neighbor[ dummyNodeu ] >= 0 )
					{
						struct _mateEdge &me = scaffoldNodes[ contigIdu ].mateEdges[ dummyNodeu ] ;
						if ( !blocks.IsSignificantDifferent( me.support[ me.supportUse ].GetCount(), 100, s[j].support[ s[j].supportUse ].GetCount(), 100 ) &&
							( me.support[ me.supportUse ].GetCount() < 20 || !me.support[ me.supportUse ].IsUnique() ) )
						{
							reported = true ;
							scaffoldNodes[ contigIdu ].neighbor[ dummyNodeu ] = -2 ;
							problematicMates.push_back( me ) ;
						}
					}
					if ( scaffoldNodes[ contigIdv ].neighbor[ dummyNodev ] >= 0 )
					{
						struct _mateEdge &me = scaffoldNodes[ contigIdv ].mateEdges[ dummyNodev ] ;
						if ( !blocks.IsSignificantDifferent( me.support[ me.supportUse ].GetCount(), 100, s[j].support[ s[j].supportUse ].GetCount(), 100 ) &
							( me.support[ me.supportUse ].GetCount() < 20 || !me.support[ me.supportUse ].IsUnique() ) )
						{
							scaffoldNodes[ contigIdv ].neighbor[ dummyNodev ] = -2 ;
							if ( !reported )
								problematicMates.push_back( me ) ;
						}
					}

					problematicMates.push_back( s[j] ) ;
				}
			}
		}

		for ( i = 0 ; i < contigCnt ; ++i )
		{
			if ( scaffoldNodes[i].neighbor[0] < -1 )
				scaffoldNodes[i].neighbor[0] = -1 ;
			if ( scaffoldNodes[i].neighbor[1] < -1 )
				scaffoldNodes[i].neighbor[1] = -1 ;
		}

		// Remove cycle.
		// which records for the misassemblies in the scaffolding.

		// Build graph
		int edgeCnt = 0 ;
		for ( i = 0 ; i < contigCnt ; ++i )
		{
			if ( scaffoldNodes[i].neighbor[0] != -1 )
				++edgeCnt ;
			if ( scaffoldNodes[i + 1].neighbor[1] != -1 )
				++edgeCnt ;
		}
		ContigGraph contigGraph( contigCnt, contigCnt + edgeCnt ) ;
		for ( i = 0 ; i < contigCnt - 1 ; ++i )
		{
			if ( genome.GetChrIdFromContigId( i ) == genome.GetChrIdFromContigId( i + 1 ) )
				contigGraph.AddEdge( i, 1, i + 1, 0 ) ;
		}

		for ( i = 0 ; i < contigCnt ; ++i )	
		{
			if ( scaffoldNodes[i].neighbor[0] != -1 )	
				contigGraph.AddEdge( i, 0, scaffoldNodes[i].neighbor[0], scaffoldNodes[i].neighborType[0] ) ;	
			if ( scaffoldNodes[i].neighbor[1] != -1 )	
				contigGraph.AddEdge( i, 1, scaffoldNodes[i].neighbor[1], scaffoldNodes[i].neighborType[1] ) ;	
		}

		// Record cycles
		int *visitTime = new int[2 * contigCnt + 2] ; 
		memset( visitTime, -1, sizeof( int ) * ( 2 * contigCnt + 2 ) ) ;
		
		std::vector<int> cycleNodes ;
		for ( i = 0 ; i < contigCnt ; ++i )	
		{
			if ( misassembledInfo[i].u != -1 || misassembledInfo[i].v != -1 )
				continue ;
			//printf( "%d\n", i ) ;
			cycleNodes.clear() ;

			if ( contigGraph.IsInCycle( i, cycleNodes, visitTime ) )
			{
				int ccnt = cycleNodes.size() ;
				for ( j = 0 ; j < ccnt ; ++j )
				{
					//printf( "hi %d %d: %d %d\n", i, cycleNodes[j], misassembledInfo[cycleNodes[j]].u, misassembledInfo[ cycleNodes[j] ].v ) ;
					if ( misassembledInfo[ cycleNodes[j] ].u == -1 && misassembledInfo[ cycleNodes[j] ].v == -1 )
					{
						misassembledInfo[ cycleNodes[j] ].u = cycleNodes[j] ;
						misassembledInfo[ cycleNodes[j] ].v = misassembledCycles.size() ;
						misassembledInfo[ cycleNodes[j] ].type = 4 ;
						//printf( "hi2 %d %d: %d %d\n", i, cycleNodes[j], misassembledInfo[cycleNodes[j]].u, misassembledInfo[ cycleNodes[j] ].v ) ;
						// Remove the effects on the nieghbors
						int v, dummyV ;
						for ( int l = 0 ; l < 2 ; ++l )
						{
							v = scaffoldNodes[ cycleNodes[j] ].neighbor[l] ;
							if ( v != -1 )
							{
								dummyV = scaffoldNodes[ cycleNodes[j] ].neighborType[l] ;
								scaffoldNodes[v].neighbor[ dummyV ] = -1 ;
								scaffoldNodes[v].neighborType[ dummyV ] = -1 ;
							}
						}
						scaffoldNodes[ cycleNodes[j] ].neighbor[0] = scaffoldNodes[ cycleNodes[j] ].neighbor[1] = -1 ;
						scaffoldNodes[ cycleNodes[j] ].neighborType[0] = scaffoldNodes[ cycleNodes[j] ].neighborType[1] = -1 ;
					}
				}
				misassembledCycles.push_back( cycleNodes ) ;
			}
		}

		delete[] visitTime ;
		
		for ( i = 0 ; i < contigCnt ; ++i )
		{
			if ( misassembledInfo[i].u != -1 )
			{
				//struct _contig c = genome.GetContigInfo( i ) ;
				misassembledContigs.push_back( i ) ;
			}
		}
	}

	void Output( FILE *fpOut, Alignments &alignments ) 
	{
		int i, j ;
		std::vector<struct _scaffold> scaffold ;		
		std::vector<int> visit ;
		bool *used = new bool[ contigCnt ] ;
		memset( used, false, sizeof( bool ) * contigCnt ) ;

		FILE *fpCS ; // the genomic sequence involved in connection 
		int headerCnt = 0 ;
		if ( outputConnectionSequence )
		{
			char buffer[256] ;
			sprintf( buffer, "%s_cs.fa", prefix ) ;
			fpCS = fopen( buffer, "w" ) ;
		}

		// Output misassembled contig information
		int mcCnt = misassembledContigs.size() ;
		for ( i = 0 ; i < mcCnt ; ++i )
		{
			struct _contig c = genome.GetContigInfo( misassembledContigs[i] ) ;			
			if ( misassembledInfo[c.id].type <= 3 )
			{
				fprintf( fpOut, "Contig %d is misassembled. (%s %" PRId64 "-%" PRId64 "): (%"PRId64"-%"PRId64") (%"PRId64"-%"PRId64"): ", 
						c.id, alignments.GetChromName( c.chrId ), c.start + 1, c.end + 1,
						blocks.geneBlocks[ misassembledInfo[ c.id ].u ].start + 1, blocks.geneBlocks[ misassembledInfo[ c.id ].u ].end + 1, 
						blocks.geneBlocks[ misassembledInfo[ c.id ].v ].start + 1, blocks.geneBlocks[ misassembledInfo[ c.id ].v ].end + 1 ) ;

				if ( misassembledInfo[ c.id ].type == 0 )
				{
					fprintf( fpOut, "Found deletion between.\n") ;	
				}
				else if ( misassembledInfo[ c.id ].type == 1 )
				{
					fprintf( fpOut, "Reverse first block.\n" ) ;
				}
				else if ( misassembledInfo[ c.id ].type == 2 )
				{
					fprintf( fpOut, "Reverse second block.\n" ) ;
				}
				else if ( misassembledInfo[ c.id ].type == 3 )
				{
					fprintf( fpOut, "Swap the two blocks.\n" ) ;
				}
			}
			else if ( misassembledInfo[ c.id ].type == 4 )
			{
				if ( !used[ misassembledInfo[ c.id ].v ] )
				{
					used[ misassembledInfo[ c.id ].v ] = true ;
					std::vector<int> &cycleNodes = misassembledCycles[ misassembledInfo[ c.id ].v ] ;
					int cnt = cycleNodes.size() ;
					fprintf( fpOut, "Misassembled scaffolds found among %d contigs: ", cnt ) ;

					for ( j = 0 ; j < cnt - 1 ; ++j )
					{
						struct _contig c= genome.GetContigInfo( cycleNodes[j] ) ;
						fprintf( fpOut, "(%d %s %"PRId64"-%"PRId64") ",  c.id, alignments.GetChromName( c.chrId ), c.start + 1, c.end + 1 ) ;
					}
					struct _contig c= genome.GetContigInfo( cycleNodes[j] ) ;
					fprintf( fpOut, "(%d %s %"PRId64"-%"PRId64")\n",  c.id, alignments.GetChromName( c.chrId ), c.start + 1, c.end + 1 ) ;
				}
				
			}
		}
		memset( used, false, sizeof( bool ) * contigCnt ) ;

		//printf( "%d\n", chrIdCnt ) ;	
		for ( i = 0 ; i < contigCnt ; ++i )
		{
			if ( used[i] )
				continue ;
			if ( scaffoldNodes[i].neighbor[0] != -1 && scaffoldNodes[i].neighbor[1] != -1 )
				continue ;
			if ( scaffoldNodes[i].neighbor[0] == -1 && scaffoldNodes[i].neighbor[1] == -1 ) // single-contig scaffold
				continue ;

			int p = i ;
			int direction = 0 ;
			if ( scaffoldNodes[i].neighbor[0] == -1 )
				direction = 1 ;
			else
				direction = -1 ;

			scaffold.clear() ;
			visit.clear() ;
			while ( p != -1 )
			{
				//if ( used[p] ) // avoid a cycle
				//	break ;
				used[p] = true ;
				struct _scaffold newS ;
				newS.contigId = p ;//scaffoldNodes[p].contigId ;

				//printf( "%d %d\n", p, scaffoldNodes[p].chrId ) ;
				if ( direction == -1 )
					newS.flip = 1 ;
				else
					newS.flip = 0 ;
				scaffold.push_back( newS ) ;
				visit.push_back( p ) ;

				int tag ;
				if ( direction == -1 )
					tag = 0 ;
				else
					tag = 1 ;
					
				if ( scaffoldNodes[p].neighborType[tag] == 1 )
					direction = -1 ;
				else
					direction = 1 ;
				p = scaffoldNodes[p].neighbor[tag] ;
			}

			// Output the connections
			int cnt = scaffold.size() ;
			fprintf( fpOut, "%d: ", cnt ) ;
			for ( j = 0 ; j < cnt ; ++j )
			{
				char orientation = ( scaffold[j].flip ? '-' : '+' ) ;
				struct _contig contigInfo = genome.GetContigInfo( scaffold[j].contigId ) ;
				fprintf( fpOut, "(%s %"PRId64" %d %c) ", alignments.GetChromName( genome.GetChrIdFromContigId( scaffold[j].contigId ) ), 
						//alignments.GetChromLength( genome.GetChrIdFromContigId( scaffold[j].contigId ) ), 
						contigInfo.end - contigInfo.start + 1, 
						scaffold[j].contigId, orientation ) ;
			}
			fprintf( fpOut, "\n" ) ;
			
			for ( j = 0 ; j < cnt - 1 ; ++j )
			{
				int tag ;
				if ( scaffold[j].flip )
					tag = 0 ;
				else
					tag = 1 ;

				struct _mateEdge &edge = scaffoldNodes[ visit[j] ].mateEdges[tag] ;
				int64_t a, b, c, d ;
				if ( blocks.geneBlocks[ edge.u ].contigId == visit[j] )
				{
					a = blocks.geneBlocks[ edge.u ].start ;
					b = blocks.geneBlocks[ edge.u ].end ;
					c = blocks.geneBlocks[ edge.v ].start ;
					d = blocks.geneBlocks[ edge.v ].end ;
				}
				else
				{
					a = blocks.geneBlocks[ edge.v ].start ;
					b = blocks.geneBlocks[ edge.v ].end ;
					c = blocks.geneBlocks[ edge.u ].start ;
					d = blocks.geneBlocks[ edge.u ].end ;
				}
				//fprintf( fpOut, "%d: (%s %" PRI64 "-%" PRI64 ") (%s %" PRI64 "-%" PRI64 ")\n", 
				fprintf( fpOut, "\t%d: (%s:%"PRId64"-%"PRId64") (%s:%"PRId64"-%"PRId64")\n", 
						edge.support[ edge.supportUse ].GetCount(), alignments.GetChromName( genome.GetChrIdFromContigId( scaffold[j].contigId ) ), a + 1, b + 1,
						alignments.GetChromName( genome.GetChrIdFromContigId( scaffold[j + 1].contigId ) ), c + 1, d + 1 ) ;
			}

			// Output the genomic sequence associated with the connection
			if ( outputConnectionSequence )
			{
				std::vector<struct _pair> geneBlockIds ;
				for ( j = 0 ; j < cnt - 1 ; ++j )		
				{
					int tag ;
					if ( scaffold[j].flip )
						tag = 0 ;
					else
						tag = 1 ; 

					struct _mateEdge &edge = scaffoldNodes[ visit[j] ].mateEdges[tag] ;
				
					if ( blocks.geneBlocks[ edge.u ].contigId == visit[j] )
					{
						struct _pair id ;
						id.a = edge.u ;
						id.b = scaffold[j].flip ? 1 : 0  ;
						geneBlockIds.push_back( id ) ;

						id.a = edge.v ;
						id.b = scaffold[j + 1].flip ? 1 : 0 ;
						geneBlockIds.push_back( id ) ;
					}
					else
					{
						struct _pair id ;
						id.a = edge.v ;
						id.b = scaffold[j].flip ? 1 : 0  ;
						geneBlockIds.push_back( id ) ;

						id.a = edge.u ;
						id.b = scaffold[j + 1].flip ? 1 : 0 ;
						geneBlockIds.push_back( id ) ;
					}
				}

				int idCnt = geneBlockIds.size() ;
				int tag = 0 ;
				char header[1023] ; 
				//printf( "%d\n", idCnt ) ;
				for ( j = 2 ; j <= idCnt ; j += 2 )
				{
					if ( j < idCnt && geneBlockIds[j].a == geneBlockIds[j - 1].a )
					{
						// The connection spanning multiple gene blocks
						continue ;
					}
					
					/*if ( geneBlockIds[j].a == 693 )
					{
						printf( "== %d: %d %d %d %d\n", j, tag, idCnt, geneBlockIds[j - 1].a, geneBlockIds[j].a ) ;
					}*/

					header[0] = '\0' ;
					int len = 0, headerLen = 0 ;
					struct _geneBlock &gb = blocks.geneBlocks[ geneBlockIds[tag].a ] ; 
					int differentScaffoldConnection = 0 ;
					// Test whether this is a connection from the contigs from the same scaffold
					int k ;
					for ( k = tag + 1 ; k < j ; ++k )
					{
						if ( blocks.geneBlocks[ geneBlockIds[k].a ].chrId != blocks.geneBlocks[geneBlockIds[k - 1].a ].chrId )
							break ;
					} 
					if ( k < j )
						differentScaffoldConnection = 1 ;

					sprintf( header, ">rascaf_CS_%d (%s:%"PRId64"-%"PRId64") %d %d", headerCnt, 
						alignments.GetChromName( gb.chrId ), gb.start + 1, gb.end + 1, 
						differentScaffoldConnection, ( j - tag ) / 2 + 1 ) ;	
					headerLen = strlen( header ) ;
					int tmp ;
					for ( int k = tag ; k < j ; k += 2 )
					{
						tmp = blocks.GetGeneBlockEffectiveLength( geneBlockIds[k].a ) ; 
						len += tmp ;
						sprintf( header + headerLen, " %d", tmp ) ;
						for ( ; header[ headerLen ] ; ++ headerLen )
							;
					}
					tmp = blocks.GetGeneBlockEffectiveLength( geneBlockIds[j - 1].a ) ;
					len += tmp ;
					sprintf( header + headerLen, " %d", tmp ) ;
					//printf( "%s\n", header ) ;
					char *seq = ( char * )malloc( sizeof( char ) * ( len + 3 ) ) ;
					seq[0] = '\0' ;
					int seqLen = 0 ;
					// Create the header
					for ( int k = tag ; k < j ; k += 2 )
					{
						int l = blocks.GetGeneBlockEffectiveSequence( geneBlockIds[k].a, geneBlockIds[k].b, seq + seqLen, genome ) ;
						seqLen += l ;
					}
					blocks.GetGeneBlockEffectiveSequence( geneBlockIds[j - 1].a, geneBlockIds[j - 1].b, seq + seqLen, genome ) ;
					
					++headerCnt ;
					//if ( seq[0] == '\0' )
					//	printf( "%s | %d %d %d %d %d\n", header, tag, j, seqLen, len, geneBlockIds[0].a ) ;
									
					fprintf( fpCS, "%s\n%s\n", header, seq ) ;

					free( seq ) ;
					tag = j ;
				}
			}
		}
		delete[] used ;

		if ( outputConnectionSequence )
			fclose( fpCS ) ;

		// Print out the problematic mates.
		fprintf( fpOut, "WARNINGS:\n" ) ;
		int pcnt = problematicMates.size() ;
		for ( i = 0 ; i < pcnt ; ++i )
		{
			char du, dv ;
			struct _mateEdge &edge = problematicMates[i] ;
			if ( edge.supportUse & 2 )
				du = '-' ;
			else
				du = '+' ;

			if ( edge.supportUse & 1 )
				dv = '+' ;
			else
				dv = '-' ;
			
			struct _geneBlock &u = blocks.geneBlocks[ edge.u ] ;
			struct _geneBlock &v = blocks.geneBlocks[ edge.v ] ;
			fprintf( fpOut, "%d: (%s:%"PRId64"-%"PRId64" %d %c) (%s:%"PRId64"-%"PRId64" %d %c)\n", 
				edge.support[ edge.supportUse ].GetCount(), 
				alignments.GetChromName( u.chrId ), u.start + 1, u.end + 1, u.contigId, du, 
				alignments.GetChromName( v.chrId ), v.start + 1, v.end + 1, v.contigId, dv ) ;
		}
	}
} ;

#endif
