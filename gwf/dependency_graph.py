'''Graph for dependency relationships of targets.'''


class Node:
    '''A node in the dependencies DAG.'''
    
    def __init__(self, name, task, dependencies):
        self.name = name
        self.task = task
        self.dependencies = dependencies
        
        # The node needs to be run if what it contains needs to be run
        # or any of the upstream nodes need to be run...
        self.should_run = self.task.should_run or \
                any(dep.should_run for _,dep in dependencies)
        

class DependencyGraph:
    '''A complete dependency graph, with code for scheduling a workflow.'''

    def __init__(self, workflow):
        self.workflow = workflow
        self.nodes = dict()
        for name, target in workflow.targets.items():
            if name not in self.nodes:
                self.nodes[name] = self.build_DAG(target)
    
        
    def add_node(self, task, dependencies):
        '''Create a graph node for a workflow task, assuming that its
        dependencies are already wrapped in nodes. The type of node
        depends on the type of task.
        
        Here we call the objects from the workflow for "tasks" and the
        nodes in the dependency graph for "nodes" and we represent each
        "task" with one "node", keeping the graph structure captured by
        nodes in the dependency DAG and the execution logic in the
        objects they refer to.'''
        
        node = Node(task.name, task, dependencies)
        self.nodes[task.name] = node
        return node
        
    
    def build_DAG(self, task):
        '''Run through all the dependencies for "target" and build the
        nodes that are needed into the graph, returning a new node with
        dependencies.
        
        Here we call the objects from the workflow for "tasks" and the
        nodes in the dependency graph for "nodes" and we represent each
        "task" with one "node", keeping the graph structure captured by
        nodes in the dependency DAG and the execution logic in the
        objects they refer to.'''
        
        def dfs(task):
            if task.name in self.nodes:
                return self.get_node(task.name)
                
            else:
                deps = [(fname, dfs(dep)) for fname, dep in task.dependencies]
                return self.add_node(task, deps)
        
        return dfs(task)

    def has_node(self, name):
        return name in self.nodes

    def get_node(self, name):
        return self.nodes[name]
   	
    	
    	

    def schedule(self, target_name):
        '''Linearize the targets to be run.
        
        Returns a list of tasks to be run (in the order they should run or
        be submitted to the cluster to make sure dependences are handled
        correctly) and a set of the names of tasks that will be scheduled
        (to make sure dependency flags are set in the qsub command).
        
        '''
        
        root = self.nodes[target_name]
        
        processed = set()
        scheduled = set()
        schedule = []
        
        def dfs(node):
            if node in processed:
                # we have already processed the node, and
                # if we should run the target name is scheduled
                # otherwise it isn't.
                return node.name in scheduled

            # schedule all dependencies before we schedule this task
            for _,dep in node.dependencies:
                dfs(dep)
            
            # If this task needs to run, then schedule it
            if node.should_run or node.task.job_in_queue:
                schedule.append(node)
                scheduled.add(node.name)
                
            processed.add(node)

        dfs(root)
            
        return schedule, scheduled
