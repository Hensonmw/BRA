#BRA
    __________ __________    _____   
    \______   \\______   \  /  _  \  
     |    |  _/ |       _/ /  /_\  \
     |    |   \ |    |   \/    |    \
     |______  / |____|_  /\____|__  /
            \/         \/         \/
...Bacterial Read Assembly....

========

BRA is a program to assemble Sanger sequence generated .fasta and .abi files. BRA takes input from the user of multiple sequences, creates a reverse complemnt, assembles the overlapping reads using [PEAR](http://sco.h-its.org/exelixis/web/software/pear/), and then BLAST individual assembled contigs.  

INSTALLATION
------------  

BRA was written in Python 3. For ease, the author recommends that you use the
[conda](https://www.continuum.io/downloads) package manager for python. Anaconda installation instruction can be found here https://docs.continuum.io/anaconda/install.

BRA also utilizes the following modules developed for Python:

[BioPython](http://biopython.org/wiki/Documentation):

    conda install biopython

BRA requires PEAR to be download and located in your path. PEAR does have dependencies
that need to be installed if they are not already on your computer:

    apt-get install build-essential autoconf automake libtool
    git clone https://github.com/xflouris/PEAR.git
    cd PEAR
    ./autogen.sh
    ./configure
    make
    sudo make install

    Zhang, J., Kobert, K., Flouri, T. and Stamatakis, A., 2014. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics, 30(5), pp.614-620.


Alternatively, MACOSX users can install it using [homebrew](http://brew.sh/):

    brew install pear

You will also need to have blastn installed in your $PATH. MACOSX user can install it using:

    brew install blastn

Alternatively, you can install it by downloading the latest version from
NCBI at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ . Then following their commands
for installing on your local computer.

Once installed, you will need to create a database. This can be done by following these commands:

    makeblastdb -in databasename.fna -dbtype nucl -out outputname.fna

You will need to give this PATH of these files with the --db_loc flag.


BRA can be download from the github repository:

    git clone https://github.com/Hensonmw/BRA.git

USAGE
-----

BRA requires the following user input:

--D the path to the directory where you would like the output placed

--F_fa the path to the forward sequences in Fasta file format

--F_qual the path to the forward sequences quality score in .abi format

--R_fa the path to the reverse sequences in Fasta file format

--R_qual the path to the reverse sequences quality score in .abi format

--o the output file name for the assebmled contigs.

--db the name of the database file

--db_loc the path to the database file

The short version is:

	BRA.py --D path/for/directory --F_fa path/to/Forward.fasta --F_qual path/to/Forward.abi --R_fa path/to/Reverse.fasta --R_qual path/to/Reverse.abi --o outputname --db databasename.fna --db_loc path/to/database

To help users with running BRA, sample data is provided (F_fa, F_qual, R_fa, R_qual, LSUCC_allBORisolate.fna).

    BRA.py --D ~/ --F_fa test_files/F_fa --F_qual test_files/F_qual --R_fa test_files/R_fa --R_qual test_files/R_qual --o output --db LSUCC_allBORisolate.fna --db_loc ~/test_files

**Output will be provided in the directory called Output. Within Output, there where
be:**

output.assembled.fastq          -->Contains the assembled contigs

output.assembled.fastq_corr     -->This file will exist if there were any byte characters
                                were found in the output.assembled.fastq

output.discarded.fastq          -->Unassembled reads

output.unassembled.forward.fastq-->Unassembled forward reads

output.unassembled.reverse.fastq-->Unassembled reverse reads

F.fastq                         -->Forward read in Fastq format

R.fastq                         -->Reverse read in Fastq format

RC.fastq                        -->Reverse complement read in Fastq format

query_id.fasta                  -->Fasta file of the contig

query_id_blastresult            -->Blast results in text format

**Files containing the parsed BLAST results will be returned with the query ID in text format**

Example BLAST OUTPUT will look like:

    > LSUCC267
    Length=1333

     Score = 957 bits (518),  Expect = 0.0
     Identities = 592/627 (94%), Gaps = 8/627 (1%)
     Strand=Plus/Minus

    Query  168   TTCACCGCGGCAAAGTCG-TCCACGTATAGCT-GCGATTCCACGATGATGCCCTCCACTG  225
                 ||||||||||| | |  | ||| || ||  || |||||||| |  | |||||||| | |
    Sbjct  1258  TTCACCGCGGC-ATGCTGATCCGCG-ATTACTAGCGATTCCGCCTTCATGCCCTCGAGTT  1201

    Query  226   GCAGAGGACAATCCGAAC-GACTCCTACTTTTGGAGA-CAGTGTCACCCTTGCGGGGTCG  283
                 |||||||||||||||||| ||   | |||||||||||  ||  |||||||||||||||||
    Sbjct  1200  GCAGAGGACAATCCGAACTGA-GACGACTTTTGGAGATTAG-CTCACCCTTGCGGGGTCG  1143

    Query  284   CTGCTCACTGTCATCGCCATTGTAGCACGTGTGTAGCCCAGCCTGTAAGGGCCATGAGGA  343
                 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    Sbjct  1142  CTGCTCACTGTCATCGCCATTGTAGCACGTGTGTAGCCCAGCCTGTAAGGGCCATGAGGA  1083



FAQ
---

### Does BRA assemble contigs with mismatches?
At this time, BRA will only assembly contigs that have no mismatches. Future versions will hopefully
use the quality scores.
