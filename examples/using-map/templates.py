import os.path

from gwf import AnonymousTarget


def filter_chromosomes(sample_file, output_dir):
    name = os.path.basename(sample_file)
    inputs = {'sample_file': sample_file}
    outputs = {
        'sample_file': os.path.join(
            output_dir, '{sample_file}.filtered'.format(sample_file=name)
        )
    }
    options = {}
    spec = """
    grep -E 'chrom1|chrom3' {inputs[sample_file]} > {outputs[sample_file]}
    """.format(inputs=inputs, outputs=outputs)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def compute_stats(sample_file, output_dir):
    inputs = {'sample_file': sample_file}

    name = os.path.basename(sample_file)
    outputs = {
        't_count_file': os.path.join(output_dir, '{name}.tcounts'.format(name=name)),
        'seq_len_file': os.path.join(output_dir, '{name}.seqlengths'.format(name=name)),
    }
    options = {}
    spec = """
    python scripts/compute_stats.py \
        {inputs[sample_file]} \
        {outputs[t_count_file]} \
        {outputs[seq_len_file]}
    """.format(inputs=inputs, outputs=outputs)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def sum_ts(t_count_files, output_dir):
    inputs = {'t_count_files': t_count_files}
    outputs = {'sum_file': os.path.join(output_dir, 'sum.txt')}
    options = {}
    spec = """
    cat {t_count_files} | awk -f scripts/sum_ts.awk > {outputs[sum_file]}
    """.format(t_count_files=' '.join(t_count_files), outputs=outputs)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def compute_avg_len(seq_len_file, output_dir):
    name = os.path.basename(seq_len_file)
    inputs = {'seq_len_file': seq_len_file}
    outputs = {'avg_file': os.path.join(output_dir, name + '.avg')}
    options = {}
    spec = """
    cat {inputs[seq_len_file]} | awk -f scripts/compute_avg.awk > {outputs[avg_file]}
    """.format(inputs=inputs, outputs=outputs)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)