
# This is an example workflow for read mapping using bwa and samtools.

from gwf import *

bwa_index = template(input='{refGenome}.fa.gz', 
                     output=['{refGenome}.amb', '{refGenome}.ann', '{refGenome}.pac']) \
    << 'bwa index -p {refGenome} -a bwtsw {refGenome}.fa.gz'
                     
bwa_map = template(input=['{name}_R1.fastq.gz {name}_R2.fastq.gz', 
                          '{refGenome}.amb {refGenome}.ann {refGenome}.pac'],
                   output='{name}.unsorted.bam',
                   pbs=['-l nodes=1:ppn=16', '-l walltime=1:0:0']) << '''

bwa mem -t 16 {refGenome} {name}_R1.fastq.gz {name}_R2.fastq.gz | \
    samtools view -Shb - > {name}.unsorted.bam

'''

samtools_sort = template(input='{name}.unsorted.bam', output='{name}.sorted.rmdup.bam',
                         pbs='-l walltime=1:0:0') << '''
                         
samtools sort -o {name}.unsorted.bam /scratch/$PBS_JOBID/{name} | \
	 samtools rmdup -s - {name}.sorted.rmdup.bam

'''

target('IndexGenome') << bwa_index(refGenome='ponAbe2')
target('MapReads')    << bwa_map(refGenome='ponAbe2', name='Masala')
target('SortBAM')     << samtools_sort(name='Masala')