Rascaf
=======

Described in: 
 
Song, L., Shankar, D. and Florea, L. "Rascaf: Improving Genome Assembly with RNA-seq Data", Plant and Animal Genome XXIV (PAG-2015) Conference; manuscript in preparation.

Copyright (C) 2015-2016, and GNU GPL, by Li Song, Liliana Florea

Rascaf includes source code from:

SAMtools - Copyright (C) 2008-2009, Genome Research Ltd, Heng Li 

### What is Rascaf?

Rascaf (RnA-seq SCAFfolder) uses continuity and order information from paired-end RNA-seq reads to improve a draft assembly, particularly in the gene regions. It takes as input an assembly and one or several RNA-seq data sets aligned to the genome, and recruits additional contigs into the assembly, potentially adjusting some scaffolds to better fit the data and to create longer gene models. Rascaf works in three stages. It first computes a set of candidate contig connections from the raw (original) assembly that are supported by the RNA-seq data. Then, in an optional step, the user can choose to validate and filter the connections by searching the merged gene sequences against public sequence databases. Finally, Rascaf uses these connections to select and/or re-arrange additional contigs within scaffolds and chromosomes. When Rascaf is run with multiple RNA-seq data sets, it first generates a set of connections for each set independently. Rascaf then reconciles all connections during a 'join' step that detects and resolves any conflicts. 

### Install

1. Clone the [GitHub repo](https://github.com/mourisl/rascaf), e.g. with `git clone https://github.com/mourisl/rascaf.git`
2. Run `make` in the repo directory

### Usage
Rascaf is comprised of two executable files, "rascaf" and "rascaf-join". "rascaf" identifies the connections from a single RNA-seq data set. "rascaf-join" uses the connections found by "rascaf" to build the scaffolds and, if applicable, to combine different data sets.

For "rascaf":    

	Usage: ./rascaf [OPTIONS]
	OPTIONS:
	Required:
		-b STRING: path to the BAM file for the alignment	
	Recommended:
		-f STRING: path to the raw assembly fasta file
	Others:
		-o STRING : prefix of the output file (default: rascaf)
		-ms INT: minimum support for connecting two contigs(default: 2)
		-ml INT: minimum exonic length if no intron (default: 200)
		-k INT: the size of a kmer(<=32. default: 21)
		-cs : output the genomic sequence involved in connections in file $prefix_cs.fa (default: not used)
		-v : verbose mode (default: false)


For "rascaf-join":

	Usage: ./rascaf-join [OPTIONS]
	OPTIONS:
	Required:
		-r STRING: path to the rascaf connection file. Can use multiple -r to specify multiple connection files 
	Others:
		-o STRING: prefix of the output file (default: rascaf_scaffold)
		-ms INT: minimum support alignments for the connection (default: 2)
		-ignoreGap: ignore the gap size, which do not consider the number of Ns between contigs (default: not used)		

### Output

1. For "rascaf":

Rascaf will output a list of contig connections determined from the RNA-seq data into a file '$prefix.out' (default: "rascaf.out"). This file is recognized and handled automatically by the downstream analysis software ("rascaf-join"); details below are provided for debugging purposes:

Each contig connection is recorded into a connection block:

	N: (scaff_id1 contig_len1 contig_id1 contig_ori1) ... (scaff_idN contig_lenN contig_idN contig_oriN)
	    M: (scaff_id1:start1-end1) (scaff_id2:start2-end2)        # pairwise connections 1 .. N-1
	    ....
	where:
	   scaff_id, contig_len, contig_id, contig_ori - scaffold ID, etc.
	   N  -  number of contigs being connected
	   M  -  number of read alignments supporting the connection
	   start, end - coordinates of the gene block within the scaffold

Example of connection between two scaffolds:  

	2: (scaffold00015 963557 219 +) (scaffold01887 956 115731 -)
	    23: (scaffold00015:332817-333132) (scaffold01887:76-854)

The list of connections may be preceded in the file by a number of messages regarding likely errors detected in the raw assembly, and may be followed by one or several warnings on the newly identified contig connections.

NOTE: If desired, the user can manually edit the connections in "rascaf.out" to filter weakly supported or unwanted connections before running "rascaf-join".

2. For "rascaf-join":

Rascaf_join will output the scaffolds in the file '$prefix.fa' (default: "rascaf_scaffold.fa"), along with reports on how the new scaffolds are built from the contigs in the original draft assembly (file '$prefix.info' file).

More specifically, the scaffolds info file will have one line for each output scaffold: 

	\>out_scaff_id (contig_id1 orig_scaff_id1 contig_ori1) ... (contig_idN orig_scaff_idN contig_oriN)
	...
	where:
	    out_scaff_id  - ID of output scaffold
	    contig_id1  - ID of first contig in the output scaffold
	    orig_scaff_1 - ID of original scaffold containing contig_id1
	    contig_ori1 - orientation of contid_id1 in the output scaffold 

### Example

Suppose we have two sorted BAM files from two RNA-seq data sets, A.bam and B.bam, and that the raw assembly is in the file "assembly.fa". To build the scaffold, use the commands:

	>./rascaf -b A.bam -f assembly.fa -o A
	>./rascaf -b B.bam -f assembly.fa -o B
	>./rascaf-join -r A.out -r B.out -o assembly_scaffold

The new scaffolds can be found in the "assembly_scaffold.fa" file.


### Terms of use

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received (LICENSE.txt) a copy of the GNU General
Public License along with this program; if not, you can obtain one from
http://www.gnu.org/licenses/gpl.txt or by writing to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
### Support

Create a [GitHub issue](https://github.com/mourisl/rascaf/issues).
