// The program for scaffolding with RNA-seq alignments
// Li Song
// email: lsong10@jhu.edu

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alignments.hpp"
#include "blocks.hpp"
#include "scaffold.hpp"
#include "genome.hpp"

char usage[] = "usage: a.out [options]\n"
	       "options:\n"
	       "\t-b STRING (required): the path to the alignment BAM file\n"
	       "\t-f STRING (recommended): the paths to the raw assembly fasta file(default: not used)\n"
	       "\t-o STRING : prefix of the output file (default: rascaf)\n"
	       "\t-bc STRING: the path to the alignment BAM file allowing clipping from non-spliced aligner (default: not used)\n"
	       "\t-ms INT: minimum support for connecting two contigs(default: 2)\n"
	       "\t-ml INT: minimum exonic length(default: 200)\n"
	       "\t-k INT: the size of a kmer(<=32. default: 21)\n"
	       "\t-cs : output the genomic sequence involved in connections in file $prefix_cs.fa (default: not used)\n"
	       //"\t-aggressive: make connection decisions more aggressively, may introduce much more misassemblies. (default: not used)\n"
	       "\t-v : verbose mode (default: false)\n" ;

int minimumSupport ;
int minimumEffectiveLength ;
int kmerSize ;
bool outputConnectionSequence ;
bool aggressiveMode ;
char *prefix ;
bool VERBOSE ;
FILE *fpOut ;

int main( int argc, char *argv[] )
{
	int i ;
	int ret ;

	Alignments alignments ;
	Alignments clippedAlignments ;
	Blocks blocks ;
	Genome genome ;
	char *genomeFile = NULL ;
	
	if ( argc < 2 )
	{
		printf( "%s", usage ) ;
		exit( 0 ) ;
	}

	minimumSupport = 2 ;
	minimumEffectiveLength = 200 ;
	kmerSize = 21 ;
	prefix = NULL ;
	VERBOSE = false ;
	outputConnectionSequence = false ;
	aggressiveMode = false ;

	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-b", argv[i] ) )
		{
			alignments.Open( argv[i + 1]) ;
			++i ;
		}
		else if ( !strcmp( "-o", argv[i] ) )
		{
			prefix = argv[i + 1] ;
			++i ;
		}
		else if ( !strcmp( "-f", argv[i] ) )
		{
			genomeFile = argv[i + 1] ;
			++i ;
		}
		else if ( !strcmp( "-ms", argv[i] ) )
		{
			minimumSupport = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-ml", argv[i] ) )
		{
			minimumEffectiveLength = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-k", argv[i] ) )
		{
			kmerSize = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-v", argv[i] ) )
		{
			VERBOSE = true ;
		}
		else if ( !strcmp( "-cs", argv[i] ) )
		{
			outputConnectionSequence = true ;
		}
		/*else if ( !strcmp( "-aggressive", argv[i] ) )
		{
			aggressiveMode = true ;
		}*/
		else if ( !strcmp( "-bc", argv[i] ) )
		{
			// So far, assume the input is from BWA mem
			clippedAlignments.Open( argv[i + 1] ) ;
			clippedAlignments.SetAllowSupplementary( true ) ;
			++i ;
		}
		else
		{
			fprintf( stderr, "Unknown parameter: %s\n", argv[i] ) ;
			exit( 1 ) ;
		}
	}

	if ( !alignments.IsOpened() )
	{
		printf( "Must use -b to specify the bam file." ) ;
		return 0 ;
	}

	if ( prefix != NULL )
	{
		char buffer[255] ;
		sprintf( buffer, "%s.out", prefix ) ;
		fpOut = fopen( buffer, "w" ) ;
	}
	else
	{
		char buffer[255] ;
		prefix = strdup( "rascaf" ) ;
		sprintf( buffer, "%s.out", prefix ) ;
		fpOut = fopen( buffer, "w" ) ;
	}
	
	if ( genomeFile != NULL )
	{
		genome.Open( alignments, genomeFile ) ;
		alignments.Rewind() ;
	}

	if ( outputConnectionSequence == true && genomeFile == NULL )
	{
		fprintf( stderr, "Must use -f to specify assembly file when using -cs\n" ) ;	
		exit( EXIT_FAILURE ) ;
	}
	// 74619
	//printf( "%c\n", genome.GetNucleotide( 74619, 4 ) ) ;
	//exit(0) ;
	// Build the graph
	ret = blocks.BuildExonBlocks( alignments, genome ) ;
	alignments.Rewind() ;
	fprintf( stderr, "Found %d exon blocks.\n", ret ) ;
	if ( clippedAlignments.IsOpened() )
	{
		fprintf( stderr, "Extend exon blocks with clipped alignments.\n" ) ;
		Blocks extendBlocks ;
		extendBlocks.BuildExonBlocks( clippedAlignments, genome ) ;
		clippedAlignments.Rewind() ;

		ret = blocks.ExtendExonBlocks( extendBlocks ) ;
		fprintf( stderr, "Found %d exon blocks after extension.\n", ret ) ;
	}

	blocks.GetAlignmentsInfo( alignments ) ;
	alignments.Rewind() ;

	ret = blocks.BuildGeneBlocks( alignments ) ;
	alignments.Rewind() ;
	fprintf( stderr, "Found %d gene blocks.\n", ret ) ;
	
	blocks.BuildGeneBlockGraph( alignments ) ;
	blocks.AddGeneBlockGraphByClippedAlignments( clippedAlignments ) ; 
	
	// Cleaning
	blocks.CleanGeneBlockGraph( alignments, genome ) ;

	// Scaffolding
	Scaffold scaffold( blocks, genome ) ;
	//scaffold.Init( blocks ) ;
	int componentCnt = scaffold.BuildComponent() ;
	fprintf( stderr, "Found %d non-trivial gene block components.\n", componentCnt ) ;
	// Possible for parallelization
	for ( i = 0 ; i < componentCnt ; ++i )
	{
		scaffold.ScaffoldComponent( i ) ;
	}
	
	scaffold.ScaffoldGenome() ;
	
	// Output the command line
	fprintf( fpOut, "command line: " ) ;
	char *fullpath = (char *)malloc( sizeof( char ) * 4096 ) ;
	for ( i = 0 ; i < argc ; ++i )
	{
		char c = ' ' ;
		if ( i == argc - 1 )
			c = '\n' ;
		if ( i > 0 && !strcmp( argv[i - 1], "-b" ) )
		{
			if ( realpath( argv[i], fullpath ) == NULL )
			{
				fprintf( stderr, "Failed to resolve the path of file %s.\n", argv[i] ) ;
				exit( 1 ) ;
			}
			fprintf( fpOut, "%s%c", fullpath, c ) ;
		}
		else if ( i > 0 && !strcmp( argv[i - 1], "-f" ) )
		{
			if ( realpath( argv[i], fullpath ) == NULL )
			{
				fprintf( stderr, "Failed to resolve the path of file %s.\n", argv[i] ) ;
				exit( 1 ) ;
			}
			fprintf( fpOut, "%s%c", fullpath, c ) ;
		}
		else
			fprintf( fpOut, "%s%c", argv[i], c ) ;
	}
	free( fullpath ) ;
	scaffold.Output( fpOut, alignments ) ;
	return 0 ;
}
