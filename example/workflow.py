
# This is an example workflow for read mapping using bwa and samtools.

from gwf import *
from gwf.bwa import bwa_index, bwa_map
from gwf.samtools import samtools_sort

target('UnZipGenome', input='ponAbe2.fa.gz', output='ponAbe2.fa', walltime="00:00:05") << '''
zcat ponAbe2.fa.gz > ponAbe2.fa
'''
target('IndexGenome', walltime="00:59:00") << bwa_index(refGenome='ponAbe2')
target('MapReads', walltime="00:59:00")    << bwa_map(refGenome='ponAbe2', R1='Masala_R1.fastq.gz', R2='Masala_R2.fastq.gz', bamfile='Masala.unsorted.bam')
target('SortBAM', walltime="00:59:00")     << samtools_sort(name='Masala')

