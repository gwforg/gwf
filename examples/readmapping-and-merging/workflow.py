
# This is an example workflow for read mapping using bwa and samtools.

from gwf import Workflow, template

gwf = Workflow()


# Templates
unzip = template(inputs='{refGenome}.fa.gz', outputs='{refGenome}.fa') << '''
zcat {refGenome}.fa.gz > {refGenome}.fa
'''

bwa_index = template(inputs='{refGenome}.fa',
                     outputs=['{refGenome}.amb',
                             '{refGenome}.ann', 
                             '{refGenome}.pac']) << '''

bwa index -p {refGenome} -a bwtsw {refGenome}.fa
'''

bwa_map = template(inputs=['{R1}', '{R2}',
                          '{refGenome}.amb', 
                          '{refGenome}.ann', 
                          '{refGenome}.pac'],
                   outputs='{bamfile}',
                   cores=16) << '''

bwa mem -t 16 {refGenome} {R1} {R2} | \
    samtools view -Shb - > /scratch/$GWF_JOBID/unsorted.bam

samtools sort -o /scratch/$GWF_JOBID/unsorted.bam /scratch/$GWF_JOBID/sort | \
    samtools rmdup -s - {bamfile}

'''

def merge(individual):
	inputfiles = ['{}_{}.sorted.rmdup.bam'.format(individual, i) for i in range(1,3)]
	outputfile = '{}.bam'.format(individual)
	shell_spec = '''

	samtools merge - {inputbams} | samtools rmdup -s - {name}.bam

	'''.format(inputbams = ' '.join(inputfiles), name=individual)
	options = {'input': inputfiles, 'output': outputfile}

	return (options, shell_spec)


# Workflow
R1files = ['Masala_{}_R1.fastq.gz'.format(i) for i in range(1,3)]
R2files = ['Masala_{}_R2.fastq.gz'.format(i) for i in range(1,3)]
bamfiles = ['Masala_{}.sorted.rmdup.bam'.format(i) for i in range(1,3)]

gwf.target('UnZipGenome', walltime="01:00:00") << unzip(refGenome='ponAbe2')
gwf.target('IndexGenome', walltime="01:00:00") << bwa_index(refGenome='ponAbe2')

for i in range(1,3):
	gwf.target('MapReads_{}'.format(i), walltime="01:00:00") << \
		bwa_map(refGenome='ponAbe2', 
				R1='Masala_{}_R1.fastq.gz'.format(i), 
				R2='Masala_{}_R2.fastq.gz'.format(i),
		 		bamfile='Masala_{}.sorted.rmdup.bam'.format(i))

gwf.target('Merge', walltime="01:00:00") << merge(individual='Masala')
