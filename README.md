
# WHAT IS KCLUST?

kClust is a fast and sensitive clustering method for the clustering of protein sequences. It is able to cluster large protein databases down to 20-30% sequence identity.
kClust generates a clustering where each cluster is represented by its longest sequence (representative sequence).

Since the development of kClust we have developed a new software from scratch based on similar ideas, MMseqs[2]. We recommend you to use MMseqs2 instead of kClust.
[Download MMseqs2 here](https://github.com/soedinglab/MMseqs2)


## HOW TO COMPILE KCLUST?

To compile, change into the directory where you have kClust source files and type 

    make all
    
It may be necessary to assign different values to CXX (Compiler) and CXXFLAGS (Compiler options) variables in the Makefile that are appropriate for your system.

## GETTING STARTED

To get started with the clustering, type

    kClust -i /home/maria/your_db_in_fasta_format.fas -d /path/to/output/directory

This call clusters the database down to 30% sequence identity. To adjust the minimum sequence identity in the clusters, use the -s option (score per column).

## KCLUST OPTIONS LIST

Typing

    kClust
    
provides a full list of kClust options and the correspondance of the minimum sequence identity in the cluster and score per column (see -s option).

## OUTPUT FILES

kClust produces following files as output:

    db_sorted.fas       : The input database sorted by sequence length.
    headers.dmp         : The mapping sequence index -> sequence header.
    representatives.fas : All representative sequences.
    clusters.dmp        : Mapping sequence index -> index of the representative sequence of its cluster. Each representative sequence in a cluster is therefore mapped to itself.

If the --write-time-benchmark option is set, kClust will additionally write the files tb_prefiltering.dat, tb_kDP.dat and tb_all.dat which contain the sequences that need the most computation time in the prefiltering step, kDP step and overall, respectively.

## GENERATING CLUSTER ALIGNMENTS

For generating one multiple sequence alignment file for each cluster, please use kClust_mkAln. Type

    kClust_mkAln

to get the full list of options. 
To get started, type

    kClust_mkAln -c '/usr/bin/kalign -q -i $infile -o $outfile' -d /path/to/output/directory

This will generate an 'alignments' directory with 100 subdirectories containing the cluster files with the multiple sequence alignment files. The cluster names have no meaning, their purpose is an equal distribution of the cluster files over the subdirectories.

I warmly recommend to use following options:

    --no-pseudo-headers               : keeps original headers in the alignments files, instead of replacing them with the sequence index
    --write-add-header                : writes an additional header describing the cluster
    --merge-[ncbi|uniprot]-headers    : the additional header (see --write-add-header option) is compressed, so it is clearly visible what the cluster is about, and which organisms the sequences stem from. The merging is only implemented for the NCBI and Uniprot header formats.

!WARNING! The -p option for the parallel calculation of multiple sequence alignments on the computer cluster works at the time only on our computer cluster and will not work anywhere else!

## ITERATIVE KCLUST

kClust can take the results from its previous run as input and further merge clusters. This merging is based on the profile-consensus sequence comparison of the clusters and is therefore more sensitive than simple sequence-sequence comparison. As basis for the merging, you have to calculate the multiple sequence alignments of the clusters (with kClust_mkAln). Afterwards, you have to convert the MSAs into a3m format (for example, with hhconsensus program) and compute HMMs in hhm format (for example, with hhmake). Then, you can type

    kClust -i /path/to/output/directory -d /path/to/output/directory_second_round -P
  
and kClust will merge clusters from the first round and write the results into /path/to/output/directory_second_round directory in the usual kClust format.

## References

    [1] Maria Hauser, Christian E Mayer and Johannes Söding (2013).
        kClust: fast and sensitive clustering of large protein sequence databases.
        BMC Bioinformatics, 2013, 14:248. DOI: 10.1186/1471-2105-14-248
    [2] Hauser M, SteineggerM, Soeding J. (2016).
        MMseqs software suite for fast and deep clustering and searching of large protein sequence sets
        Bioinformatics, 2016, doi: 10.1093/bioinformatics/btw006

