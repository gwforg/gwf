import re

def job_cpu_performance(job):

	used_cpu_time = job['used_cpu_time']
	used_walltime = job['used_walltime']

	print used_cpu_time, time_in_seconds(used_cpu_time)
	print used_walltime, time_in_seconds(used_walltime)

	if used_cpu_time == '--':
		return 'NA'

	else:
		print job['id'], used_cpu_time, used_walltime

def job_memory_performace(job):

	pass

def time_in_seconds(time):

	time_parts_re = re.compile(r'(((?P<days>\d+)-)?(?P<hours>\d\d):)?(?P<minutes>\d\d):(?P<seconds>\d\d(\.\d+)?)')
	time_parts = time_parts_re.match(time)

	if time_parts == None:
		return 0.0

	seconds = float(time_parts.group('seconds'))
	minutes	= float(time_parts.group('minutes'))
	hours	= float(time_parts.group('hours') or '0.0')
	days 	= float(time_parts.group('days') or '0.0')

	return seconds + 60 * minutes + 60 * 60 * hours + 24 * 60 * 60 * days