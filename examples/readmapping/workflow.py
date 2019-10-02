"""This is an example workflow for read-mapping using bwa and samtools."""

from gwf import Workflow, AnonymousTarget

gwf = Workflow()


def unzip(inputfile, outputfile):
    """A template for unzipping files."""
    inputs = [inputfile]
    outputs = [outputfile]
    options = {
        'cores': 1,
        'memory': '2g',
    }

    spec = '''
    gzcat {} > {}
    '''.format(inputfile, outputfile)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def bwa_index(ref_genome):
    """Template for indexing a genome with `bwa index`."""
    inputs = ['{}.fa'.format(ref_genome)]
    outputs = ['{}.amb'.format(ref_genome),
               '{}.ann'.format(ref_genome),
               '{}.pac'.format(ref_genome),
               '{}.bwt'.format(ref_genome),
               '{}.sa'.format(ref_genome),
               ]
    options = {
        'cores': 16,
        'memory': '1g',
    }

    spec = """
    bwa index -p {ref_genome} -a bwtsw {ref_genome}.fa
    """.format(ref_genome=ref_genome)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def bwa_map(ref_genome, r1, r2, bamfile):
    """Template for mapping reads to a reference genome with `bwa` and `samtools`."""
    inputs = [r1, r2,
              '{}.amb'.format(ref_genome),
              '{}.ann'.format(ref_genome),
              '{}.pac'.format(ref_genome),
             ]
    outputs = [bamfile]
    options = {
        'cores': 16,
        'memory': '1g',
    }

    spec = '''
    bwa mem -t 16 {ref_genome} {r1} {r2} | \
    samtools sort | \
    samtools rmdup -s - {bamfile}
    '''.format(ref_genome=ref_genome, r1=r1, r2=r2, bamfile=bamfile)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


gwf.target_from_template(
    name='UnzipGenome',
    template=unzip(
        inputfile='ponAbe2.fa.gz',
        outputfile='ponAbe2.fa'
    )
)

gwf.target_from_template(
    name='IndexGenome',
    template=bwa_index(
        ref_genome='ponAbe2'
    )
)

gwf.target_from_template(
    name='MapReads',
    template=bwa_map(
        ref_genome='ponAbe2',
        r1='Masala_R1.fastq.gz',
        r2='Masala_R2.fastq.gz',
        bamfile='Masala.bam'
    )
)
