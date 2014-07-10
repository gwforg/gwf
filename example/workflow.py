
# This is an example workflow for read mapping using bwa and samtools.

from gwf import *
from gwf.bwa import bwa_index, bwa_map
from gwf.samtools import samtools_sort


target('IndexGenome') << bwa_index(refGenome='ponAbe2')
target('MapReads')    << bwa_map(refGenome='ponAbe2', R1='Masala_R1.fastq.gz', R2='Masala_R2.fastq.gz', bamfile='Masala.unsorted.bam')
target('SortBAM')     << samtools_sort(name='Masala')

