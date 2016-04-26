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

BRA requires PEAR to be download and located in your path. PEAR does have dependacies
that need to be installed if they are not already on your computer:

    apt-get install build-essential autoconf automake libtool
    git clone https://github.com/xflouris/PEAR.git
    cd PEAR
    ./autogen.sh
    ./configure
    make
    sudo make install


Alternatively, MACOSX users can install it using [homebrew](http://brew.sh/):

    brew install pear


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

The short version is:

	BRA.py --D path/for/directory --F_fa path/to/Forward.fasta --F_qual path/to/Forward.abi --R_fa path/to/Reverse.fasta --R_qual path/to/Reverse.abi --o outputname

To help users with running BRA, sample data is provided (F_fa, F_qual, R_fa, R_qual).
    BRA.py --D ~/Documents --F_fa test_files/F_fa --F_qual test_files/F_qual --R_fa test_files/R_fa --R_qual test_files/R_qual --o output

Output will be provided in the directory called Output. Within Output, there where
be:

output.assembled.fastq          -->Contains the assembled contigs

output.assembled.fastq_corr     -->This file will exist if there were any byte characters
                                were found in the output.assembled.fastq

output.discarded.fastq          -->Unassembled reads

output.unassembled.forward.fastq-->Unassembled forward reads

output.unassembled.reverse.fastq-->Unassembled reverse reads

F.fastq                         -->Forward read in Fastq format

R.fastq                         -->Reverse read in Fastq format

RC.fastq                        -->Reverse complement read in Fastq format

blast_results.xml               -->Unparsed BLAST results

query_id                        -->Parsed BLAST results

**Files containing the parsed BLAST results will be returned with the query ID in text format**

example BLAST OUTPUT will look like:

    gi|347440079|gb|JN119202.1| Uncultured marine bacterium clone GD-C5 16S ribosomal RNA gene, partial sequence
    1146.0

    gi|268308535|gb|FJ820391.1| Uncultured bacterium clone S0008 16S ribosomal RNA gene, partial sequence
    1146.0

    106.0

    gi|268308529|gb|FJ820385.1| Uncultured bacterium clone H0046 16S ribosomal RNA gene, partial sequence
    1146.0

    106.0






FAQ
---

### Does BRA assemble contigs with mismatches?
At this time, BRA will only assembly contigs that have no mismatches. Future versions will hopefully
use the quality scores.
