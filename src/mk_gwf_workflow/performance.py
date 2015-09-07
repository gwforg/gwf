import re
import sys

def job_cpu_performance(job):

	cores = int(job['cores'])

	used_cpu_time = time_in_seconds(job['used_cpu_time'])
	used_walltime = time_in_seconds(job['used_walltime'])

	if used_walltime == 0:
		return '100%'

	else:
		total_walltime = cores * used_walltime
		return '{:5.1f}%'.format(100*float(used_cpu_time)/total_walltime)

def job_memory_performace(job):
	
	max_memory_used = job['max_memory_used']
	memory_reserved = job['memory_reserved']
	print max_memory_used, memory_reserved, memory_in_bytes(max_memory_used), memory_in_bytes(memory_reserved)

def time_in_seconds(time):

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

	if memory[-1].isdigit():
		return float(memory)
	elif memory[-1] == 'M':
		return 1024 * 1024 * float(memory[:-1])
	elif memory[-1] == 'G':
		return 1024 * 1024 * 1024 * float(memory[:-1])
	else:
		print 'INVALID MEMORY FORMAT: {}'.format(memory)
		sys.exit(1)