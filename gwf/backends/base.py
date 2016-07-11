class Backend(object):

    def get_state_of_jobs(self):
        raise NotImplementedError()

    def submit_command(self, target, script_name, dependent_ids):
        raise NotImplementedError()

    def build_cancel_command(self, job_ids):
        raise NotImplementedError()
