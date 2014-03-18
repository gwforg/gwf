
# This is an example workflow for read mapping using bwa and samtools.

from gwf import *
from gwf.bwa import *
from gwf.samtools import *

target('IndexGenome') << bwa_index(refGenome='ponAbe2')

samples = ['Buschi', 'Suma', 'Moni', 'Masala']

for sample in samples:
    R1 = '{name}_R1.fastq.gz'.format(name=sample)
    R2 = '{name}_R2.fastq.gz'.format(name=sample)
    bamfile = '{name}.bam'.format(name=sample)
    target('Map'+sample)  << bwa_map(refGenome='ponAbe2', R1=R1, R2=R2, bamfile=bamfile)
    target('Sort'+sample) << samtools_sort(name=sample)


from gwf_workflow.workflow import build_workflow
build_workflow()
