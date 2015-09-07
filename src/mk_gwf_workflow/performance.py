import re

def job_cpu_performance(job):

	cores = int(job['cores'])

	used_cpu_time = time_in_seconds(job['used_cpu_time'])
	used_walltime = time_in_seconds(job['used_walltime'])

	

def job_memory_performace(job):

	pass

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