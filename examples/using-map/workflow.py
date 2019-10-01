import os

from gwf import Workflow
from gwf.workflow import collect, select

from templates import (
    filter_chromosomes,
    compute_stats,
    sum_ts,
    compute_avg_len,
)

FOO_FILES = ['data/sample{}'.format(i) for i in range(3)]

os.makedirs('results', exist_ok=True)


gwf = Workflow()

A = gwf.map(filter_chromosomes, FOO_FILES, extra={'output_dir': 'results'})
B = gwf.map(compute_stats, A.outputs, extra={'output_dir': 'results'})

C = gwf.target_from_template(
    'sum_counts',
    sum_ts(**collect(B.outputs, ['t_count_file']), output_dir='results')
)

D = gwf.map(compute_avg_len, select(B.outputs, ['seq_len_file']), extra={'output_dir': 'results'})
