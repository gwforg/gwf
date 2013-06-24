'''Graph for dependency relationships of targets.'''

class Node:
    '''A node in the dependencies DAG.'''

    def __init__(self, target, dependencies):
        self.target = target
        self.dependencies = dependencies
        
        self.should_run = target.should_run() or \
        					any(dep.should_run for dep in dependencies)

        

class DependencyGraph:
    '''A complete dependency graph, with code for scheduling a workflow.'''

    def __init__(self, workflow):
        self.nodes = dict()

        for name, target in workflow.targets.items():
            if name not in self.nodes:
                self.nodes[name] = self.build_DAG(target)
    
    def add_node(self, target, dependencies):
        node = Node(target, dependencies)
        self.nodes[target.name] = node
        return node
    
    def build_DAG(self, target):
        '''Run through all the dependencies for "target" and build the
        nodes that are needed into the graph, returning a new node with
        dependencies.'''
        
        def dfs(targ):
            if targ.name in self.nodes:
                return self.get_node(targ.name)
            else:
                deps = [(fname, dfs(dep)) 
                        for fname, dep in targ.dependencies]
                return self.add_node(targ, deps)
        
        return dfs(target)

    def has_node(self, name):
        return name in self.nodes

    def get_node(self, name):
        return self.nodes[name]

    def print_workflow_graph(self, out):
    	'''Print the workflow to file object out in graphviz format.'''
    	
    	print >> out, 'digraph workflow {'
    	
    	# Handle nodes
    	for node in self.nodes.values():
    		print >> out, node.target.name, ';' # FIXME annotate with run info
    		
    	# FIXME: handle system files ... can't right here yet
    	
    	for src in self.nodes.values():
    	    for fname,dst in src.dependencies:
    	        print >> out, src.target.name, '->', dst.target.name,
    	        print >> out, '[label="%s"]' % fname,
    	        print >> out, ';'
    	
    	print >> out, '}'
    	
    	
    	
    	

    def schedule(self, target_name):
        '''Linearize the targets to be run.
        
        Returns a list of tasks to be run (in the order they should run or
        be submitted to the cluster to make sure dependences are handled
        correctly) and a set of the names of tasks that will be scheduled
        (to make sure dependency flags are set in the qsub command).
        
        '''
        
        # FIXME: After refactoring, this code can be simplified
        
        root = self.nodes[target_name]
        
        processed = set()
        scheduled = set()
        schedule = []
        
        def dfs(node):
            if node in processed:
                # we have already processed the node, and
                # if we should run the target name is scheduled
                # otherwise it isn't.
                return node.target.name in scheduled
            
            # Process dependencies and find out if any upstream tasks
            # needs to run.
            upstream_runs = False
            for _,dep in node.dependencies:
                upstream_runs |= dfs(dep)
            
            # If this task needs to run, then schedule it

            if upstream_runs or node.target.should_run():
                schedule.append(node)
                scheduled.add(node.target.name)
                to_run = True
            else:
                to_run = False
                
            processed.add(node)
            return to_run

        dfs(root)
            
        return schedule, scheduled
