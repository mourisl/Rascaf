#ifndef _MOURISL_CONTIG_GRAPH_HEADER
#define _MOURISL_CONTIG_GRAPH_HEADER

#include <stdio.h>
#include <vector>
#include <inttypes.h>

#include "defs.h"

struct _contigGraphEdge 
{
	int u, v ; 
	int dummyU, dummyV ; // use left dummy node (0) or the right dummy node (1) ;
	int next ;
} ;

class ContigGraph
{
private:
	int n ; // node count
	int m ; // max edge count
	int edgeUsed ;
	int *contigScafIds ;

	bool SearchCycle( int node, int inDummy, int startNode, std::vector<int> &cycleNodes, int time, int visitTime[] )
	{
		if ( visitTime[node] == time )
			return false ;
		//if ( startNode == 21330 )
		//	printf( "hi: %d %d %d\n", node, inDummy ) ;
		visitTime[node] = time ;

		int p = contigGraph[node].next ;
		while ( p != -1 )
		{
			//printf( "%d %d %d\n", p, contigGraph[p].u, contigGraph[p].next ) ;
			if ( contigGraph[p].u != node || contigGraph[p].dummyU == inDummy )
			{
				p = contigGraph[p].next ;
				continue ;
			}
			
			if ( contigGraph[p].v == startNode )
			{
				cycleNodes.push_back( node ) ;
				return true ;
			}
			
			if ( SearchCycle( contigGraph[p].v, contigGraph[p].dummyV, startNode, cycleNodes, time, visitTime ) )
			{
				cycleNodes.push_back( node ) ;
				return true ;
			}
			p = contigGraph[p].next ;
		}
		return false ;
	}
	struct _contigGraphEdge *contigGraph ; // Use adjacent list to represent the graph
public:

	ContigGraph( int nodeCnt, int edgeCnt ) 
	{
		int i ;
		n = nodeCnt ;
		m = edgeCnt ;
		contigGraph = new struct _contigGraphEdge[nodeCnt + 2 * edgeCnt+4] ;
		contigScafIds = new int[n] ;
		for ( i = 0 ; i < n ; ++i )
		{
			contigGraph[i].next = -1 ;
			contigScafIds[i] = -1 ;
		}
		edgeUsed = i ;

	}

	~ContigGraph() 
	{
		delete[] contigGraph ;
		delete[] contigScafIds ;
	}

	void AddEdge( int u, int dummyU, int v, int dummyV, bool checkDuplicate = false ) 
	{
		if ( checkDuplicate )
		{
			std::vector<struct _pair> neighbors ;
			int ncnt = GetNeighbors( u, dummyU, neighbors ) ;
			for ( int i = 0 ; i < ncnt ; ++i )
			{
				if ( neighbors[i].a == v && neighbors[i].b == dummyV )
					return ;
			}
		}
		contigGraph[ edgeUsed ].u = u ;
		contigGraph[ edgeUsed ].v = v ;
		contigGraph[ edgeUsed ].dummyU = dummyU ;
		contigGraph[ edgeUsed ].dummyV = dummyV ;
		contigGraph[ edgeUsed ].next = contigGraph[u].next ;
		contigGraph[u].next = edgeUsed ;
		++edgeUsed ;
		
		contigGraph[edgeUsed].u = v ;
		contigGraph[edgeUsed].v = u ;
		contigGraph[edgeUsed].dummyU = dummyV ;
		contigGraph[edgeUsed].dummyV = dummyU ;
		contigGraph[ edgeUsed ].next = contigGraph[v].next ;
		contigGraph[v].next = edgeUsed ;
		++edgeUsed ;

		if ( edgeUsed >= n + 2 * m + 4 )
		{
			fprintf( stderr, "Too many edges.\n" ) ;
			exit( 1 ) ;
		}
	}

	void RemoveEdge( int u, int dummyU, int v, int dummyV )
	{
		int prev ;
		int p  ;
	
		prev = u ;
		p = contigGraph[u].next ;
		while ( p != -1 )
		{
			if ( contigGraph[p].dummyU == dummyU && contigGraph[p].v == v 
				&& contigGraph[p].dummyV == dummyV )
			{
				contigGraph[prev].next = contigGraph[p].next ;
			}
			
			prev = p ;
			p = contigGraph[p].next ;
		}
			
		prev = v ;
		p = contigGraph[v].next ;
		while ( p != -1 )
		{
			if ( contigGraph[p].dummyU == dummyV && contigGraph[p].v == u 
				&& contigGraph[p].dummyV == dummyU )
			{
				contigGraph[prev].next = contigGraph[p].next ;
			}
			
			prev = p ;
			p = contigGraph[p].next ;
		}
	}

	bool IsInCycle( int node, std::vector<int> &cycleNodes, int visitTime[] )
	{
		cycleNodes.clear() ;
		if ( SearchCycle( node, 0, node, cycleNodes, 2 * node, visitTime ) )
			return true ;
		cycleNodes.clear() ;
		if ( SearchCycle( node, 1, node, cycleNodes, 2 * node + 1, visitTime ) )
			return true ;
		
		return false ;
	}
	
	int GetNeighbors( int u, int dummyU, std::vector<struct _pair> &neighbors )  
	{
		int p = contigGraph[u].next ;
		int cnt = 0 ;
		while ( p != -1 )
		{
			if ( contigGraph[p].dummyU == dummyU )
			{
				struct _pair np ;
				np.a = contigGraph[p].v ;
				np.b = contigGraph[p].dummyV ;
				neighbors.push_back( np ) ;
				++cnt ;
			}
			p = contigGraph[p].next ;
		}
		return cnt ;
	}

	int GetNeighbors( int u, int dummyU, struct _pair neighbors[], int maxN )  
	{
		int p = contigGraph[u].next ;
		int cnt = 0 ;
		while ( p != -1 )
		{
			if ( contigGraph[p].dummyU == dummyU )
			{
				if ( cnt >= maxN )
				{
					std::vector<struct _pair> neighbors ;
					int ncnt = GetNeighbors( u, dummyU, neighbors ) ;
					fprintf( stderr, "(%d %d) has too many neighbors. %d>=%d: ", u, dummyU, ncnt, maxN ) ;
					for ( int i = 0 ; i < ncnt ; ++i )
					{
						fprintf( stderr, " (%"PRId64" %"PRId64")", neighbors[i].a, neighbors[i].b ) ;
					}
					fprintf( stderr, "\n") ;
					exit( 1 ) ;
				}

				struct _pair np ;
				np.a = contigGraph[p].v ;
				np.b = contigGraph[p].dummyV ;
				neighbors[cnt] = np ;
				++cnt ;
			}
			p = contigGraph[p].next ;
		}
		return cnt ;
	}

	void AdjustNeighbor( int u, int dummyU, int oldV, int oldDummyV, int newV, int newDummyV )
	{
		int p = contigGraph[u].next ;
		while ( p != -1 )
		{
			if ( contigGraph[p].dummyU == dummyU && contigGraph[p].v == oldV && contigGraph[p].dummyV == oldDummyV )
			{
				contigGraph[p].v = newV ; 
				contigGraph[p].dummyV = newDummyV ;
			}
			p = contigGraph[p].next ;
		}
		return ;
	}

	// connect dummyU and dummyV, and move the original edges connected
	// to dummyU to 1-dummyV
	void InsertNode( int u, int dummyU, int v, int dummyV )
	{
		int p = contigGraph[u].next ;
		int cnt = 0 ;
		int tagU ;
		int tagV ;
		contigGraph[u].next = -1 ;
		while ( p != -1 )
		{
			int nextp = contigGraph[p].next ;
			if ( contigGraph[p].dummyU == dummyU )
			{
				contigGraph[p].u = v ;
				contigGraph[p].dummyU = 1 - dummyV ;
				contigGraph[p].next = contigGraph[v].next ; 
				contigGraph[v].next = p ;

				AdjustNeighbor( contigGraph[p].v, contigGraph[p].dummyV, u, dummyU, v, 1- dummyV ) ;
			}
			else
			{
				contigGraph[p].next = contigGraph[u].next ;
				contigGraph[u].next = p ;
			}
			p = nextp ;
		}
		AddEdge( u, dummyU, v, dummyV ) ;
	}

	void SetContigScafId( int cid, int sid )
	{
		contigScafIds[cid] = sid ;
	}

	int GetContigScafId( int cid )
	{
		return contigScafIds[cid] ;
	}

	void RemoveAdjacentEdges( int u ) 
	{
		contigGraph[u].next = -1 ;
	}

	void ResetEdgeUsed() 
	{
		edgeUsed = n ; 
	}
} ;

#endif
