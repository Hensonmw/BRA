#!/usr/bin/env python
# encoding: utf-8


"""
Copyright 2016 Michael W Henson. All rights reserved.

Michael W. Henson || thethrashlab.com
"""


import argparse
import os
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
import subprocess
import shlex
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time
import pdb


def askingforfiles():
    parser = argparse.ArgumentParser(
        description="Directory for your output folder")
    parser.add_argument(
        "--D",
        required=True,
        help="Provide the desired directory path please ",
        type=str
    )
    parser.add_argument(
        "--F_fa",
        required=True,
        help="Provide the forward .fasta file",
        type=str
    )
    parser.add_argument(
        "--R_fa",
        required=True,
        help="Provide the reverse .fasta file",
        type=str
    )
    parser.add_argument(
        "--F_qual",
        required=True,
        help="Provide the forward .abi file",
        type=str
    )
    parser.add_argument(
        "--R_qual",
        required=True,
        help="Provide the reverse .abi file",
        type=str
    )
    parser.add_argument(
        "--o",
        required=True,
        help="Provide the outfile name for your assembled contigs",
        type=str
    )
    parser.add_argument(
        "--db",
        required=True,
        help="Provide the database file name",
        type=str
        )
    parser.add_argument(
        "--db_loc",
        required=True,
        help="Provide the database file path",
        type=str
        )
    return parser.parse_args()


def makedirectory(directory_file):
    # Making a new directory for the output from this program
    newpath = os.path.join(directory_file, 'Output')
    print("Check here for your outputs")
    print(newpath)
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    return newpath


def R_qual(R_read, R_qual):
    # Creating a list to store our sequences in
    sequences = []
    # Creating a fastq file from the fasta and quality scores
    records = PairedFastaQualIterator(open(R_read), open(R_qual))
    for read in records:
        if "n" in read.seq:
            raise IOError("Sequence {} has ambigous bp".format(read.id))
        else:
            sequences.append(read)
    return sequences


def reverse(R_read, F_read):
    ''' Because we are working with DNA, we need to take the reverse read from our
    sanger sequences and create the reverse complement (RC). This will allow us
    to align our sequences.
    '''
    sequences = []
    # Creating a list file of the reads
    for read in R_read:
        x = read.reverse_complement(read.seq)
        sequences.append(x)
    return sequences


def F_qual(F_read, F_qual):
    sequences = []
    records = PairedFastaQualIterator(open(F_read), open(F_qual))
    for read in records:
        if "n" in read.seq:
            raise IOError("Sequence {} has ambigous bp".format(read.id))
        else:
            sequences.append(read)
    return sequences


def contigmaker(r):
    f_reads = "F.fastq"
    rc_reads = "RC.fastq"
    output = r
    text_command = 'pear -f {} -r {} -o {} -s 1'.format(
        f_reads,
        rc_reads,
        output
        )
    cmd = shlex.split(text_command)
    my_proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
        )
    if my_proc.returncode == 0:
        pass
    else:
        raise IOError("An Error has been raised\n{}".format(my_proc.stderr))


def correcting(output):
    file_name = output+".assembled.fastq"
    new = open(file_name+"_corr", "w")
    '''
     Because of issues with a byte character being read in from pear
     We need to remove them from the file and replace them.
    '''
    with open(file_name, encoding="utf-8", errors='backslashreplace') as assembled:
        for line in assembled:
            try:
                x = line.replace("\\xa9", "I").replace("\\xa0", "I").replace("\\xa3","I")
                new.write(x)
            except:
                new.write(line)
    new.close()


def fnafiles(fastq):
    seqid = []
    for record in SeqIO.parse(fastq+".assembled.fastq_corr", 'fastq'):
        seqid.append(record.id)
        with open(record.id+".fasta", "w") as fna:
            SeqIO.write(record, fna, "fasta")
    return seqid


def BestBlastHit(seqid, db, output, loc, D):
    os.chdir(loc)
    for seq in seqid:
        text_command = 'blastn -db {} -query {} -out {} -max_hsps 2'.format(
            db,
            D+"/"+seq+".fasta",
            D+"/"+seq+"_blastresult"
            )
        cmd = shlex.split(text_command)
        my_proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
            )
        if my_proc.returncode == 0:
            pass
        else:
            raise IOError("An Error has been raised\n{}".format(my_proc.stderr))


def main():
    time.sleep(1)
    path = askingforfiles()
    X = makedirectory(path.D)
    R_fq = R_qual(path.R_fa, path.R_qual)
    F_fq = F_qual(path.F_fa, path.F_qual)
    RC = reverse(R_fq, F_fq)
    os.chdir(X)
    R_record = open("R.fastq", "w")
    SeqIO.write(R_fq, R_record, "fastq")
    F_record = open("F.fastq", "w")
    SeqIO.write(F_fq, F_record, "fastq")
    RC_output = open("RC.fastq", "w")
    SeqIO.write(RC, RC_output, "fastq")
    RC_output.close()
    F_record.close()
    contigmaker(path.o)
    correcting(path.o)
    seqid = fnafiles(path.o)
    BestBlastHit(seqid, path.db, path.o, path.db_loc, X)


if __name__ == '__main__':
    main()
