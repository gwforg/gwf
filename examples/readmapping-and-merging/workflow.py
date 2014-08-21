
# This is an example workflow for read mapping using bwa and samtools.

from gwf import *


# Templates
unzip = template(input='{refGenome}.fa.gz', output='{refGenome}.fa') << '''
zcat {refGenome}.fa.gz > {refGenome}.fa
'''

bwa_index = template(input='{refGenome}.fa', 
                     output=['{refGenome}.amb', 
                             '{refGenome}.ann', 
                             '{refGenome}.pac']) << '''

bwa index -p {refGenome} -a bwtsw {refGenome}.fa
'''

bwa_map = template(input=['{R1}', '{R2}', 
                          '{refGenome}.amb', 
                          '{refGenome}.ann', 
                          '{refGenome}.pac'],
                   output='{bamfile}',
                   cores=16) << '''

bwa mem -t 16 {refGenome} {R1} {R2} | \
    samtools view -Shb - > /scratch/$GWF_JOBID/unsorted.bam

samtools sort -o /scratch/$GWF_JOBID/unsorted.bam /scratch/$GWF_JOBID/sort | \
    samtools rmdup -s - {bamfile}

'''


# Workflow
R1files = ['Masala_{}_R1.fastq.gz'.format(i) for i in range(1,3)]
R2files = ['Masala_{}_R2.fastq.gz'.format(i) for i in range(1,3)]
bamfiles = ['Masala_{}.sorted.rmdup.bam'.format(i) for i in range(1,3)]

target('UnZipGenome', walltime="00:00:05") << unzip(refGenome='ponAbe2')
target('IndexGenome', walltime="00:59:00") << bwa_index(refGenome='ponAbe2')

for i in range(1,3):
	target('MapReads_{}'.format(i), walltime="00:59:00") << \
		bwa_map(refGenome='ponAbe2', 
			R1='Masala_{}_R1.fastq.gz'.format(i), 
			R2='Masala_{}_R2.fastq.gz'.format(i),
		 	bamfile='Masala_{}.unsorted.bam'.format(i))

