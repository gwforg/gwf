from gwf import *

samtools_sort = template(input='{name}.unsorted.bam', output='{name}.sorted.rmdup.bam') << '''
                         
samtools sort -o {name}.unsorted.bam /scratch/$PBS_JOBID/{name} | \
	 samtools rmdup -s - {name}.sorted.rmdup.bam

'''