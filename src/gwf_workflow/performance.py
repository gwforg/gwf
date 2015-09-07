import re
import sys

def job_cpu_performance(job):
	'''

		The CPU performance of a job is calculated as how big a percentage the
		used cpu time is of the spent walltime. Note that the total walltime
		equals the number of cores times the walltime reported by SLURM.

	'''

	cores = int(job['cores'])

	used_cpu_time = time_in_seconds(job['used_cpu_time'])
	used_walltime = time_in_seconds(job['used_walltime'])

	# The jobs finished within the blink of an eye. Makes no sense
	# to calculate performance.
	if used_walltime == 0 or used_cpu_time == 0:
		return 'NA'

	# Return the performance as a well-formatted string.
	else:
		total_walltime = cores * used_walltime
		return '{:5.1f}%'.format(100*float(used_cpu_time)/total_walltime)

def job_memory_performace(job):
	'''

		The memory performance of a job is calcuated as how much the maximum memory used
		is compared to the memory reserved. Note that the formula

			max_memory_used / memory_reserved

		is true no matter how many nodes are used: In the case of a single node, it reflects
		memory used per node. In case of multiple nodes, it reflects the memory used per core.
		In both cases, this is exactly what we want.

	'''


	max_memory_used = memory_in_bytes(job['max_memory_used'])
	memory_reserved = memory_in_bytes(job['memory_reserved'])
	
	# Trivial jobs. Makes no sense to calculate performance.
	if max_memory_used == 0.0:
		return 'NA'

	# Return the performance as a well-formatted string.	
	else:
		return '{:5.1f}%'.format(100*max_memory_used/memory_reserved)

def time_in_seconds(time):
	'''	Converts a time string DAYS-HOURS:MINUTES:SECONDS into number of seconds '''

	time_parts_re = re.compile(r'(((?P<days>\d+)-)?(?P<hours>\d\d):)?(?P<minutes>\d\d):(?P<seconds>\d\d(\.\d+)?)')
	time_parts = time_parts_re.match(time)

	if time_parts == None:
		return 0

	seconds = int(time_parts.group('seconds'))
	minutes	= int(time_parts.group('minutes'))
	hours	= int(time_parts.group('hours') or '0')
	days 	= int(time_parts.group('days') or '0')

	return seconds + 60 * minutes + 60 * 60 * hours + 24 * 60 * 60 * days

def memory_in_bytes(memory):
	''' Converts compact representation of memory into bytes '''

	if memory[-1].isdigit():
		return float(memory)
	elif memory[-1] == 'K':
		return 1024 * float(memory[:-1])
	elif memory[-1] == 'M':
		return 1024 * 1024 * float(memory[:-1])
	elif memory[-1] == 'G':
		return 1024 * 1024 * 1024 * float(memory[:-1])
	elif memory[-1] == 'T':
		return 1024 * 1024 * 1024 * 1024 * float(memory[:-1])
	else:
		print 'INVALID MEMORY FORMAT: {}'.format(memory)
		sys.exit(1)