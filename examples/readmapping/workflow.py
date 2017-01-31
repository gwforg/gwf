
# This is an example workflow for read mapping using bwa and samtools.

from gwf import template, Workflow

gwf = Workflow()

# Templates
bwa_index = template(inputs=['{refGenome}.fa'],
                     outputs=['{refGenome}.amb', '{refGenome}.ann', '{refGenome}.pac']) \
    << 'bwa index -p {refGenome} -a bwtsw {refGenome}.fa'

bwa_map = template(inputs=['{R1}', '{R2}', 
	                       '{refGenome}.amb', '{refGenome}.ann', 
	                       '{refGenome}.pac'],
                   outputs=['{bamfile}'], cores=16) << '''
bwa mem -t 16 {refGenome} {R1} {R2} | \
    samtools sort | \
    samtools rmdup -s - {bamfile}
'''


gwf.target('UnZipGenome', 
	       inputs=['ponAbe2.fa.gz'], 
	       outputs=['ponAbe2.fa'], 
	       walltime="00:00:05") << '''
gzcat ponAbe2.fa.gz > ponAbe2.fa
'''
gwf.target('IndexGenome', walltime="00:59:00") << \
	bwa_index(refGenome='ponAbe2')
gwf.target('MapReads', walltime="00:59:00") << \
	bwa_map(refGenome='ponAbe2', 
		    R1='Masala_R1.fastq.gz', 
		    R2='Masala_R2.fastq.gz', 
		    bamfile='Masala.bam')


