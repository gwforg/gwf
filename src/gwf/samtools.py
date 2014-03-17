from gwfp import *

samtools_sort = template(input='{name}.unsorted.bam', output='{name}.sorted.rmdup.bam',
                         pbs='-l walltime=1:0:0') << '''
                         
samtools sort -o {name}.unsorted.bam /scratch/$PBS_JOBID/{name} | \
	 samtools rmdup -s - {name}.sorted.rmdup.bam

'''