'''Graph for dependency relationships of targets.'''

from exceptions import NotImplementedError

import workflow

class Node:
    '''A node in the dependencies DAG.'''
    
    def __init__(self, name, dependencies):
        self.name = name
        self.dependencies = dependencies
        
        # The node needs to be run if what it contains needs to be run
        # or any of the upstream nodes need to be run...
        self.should_run = self.node_should_run() or \
                any(dep.should_run for _,dep in dependencies)
    
    def node_should_run(self):
        raise NotImplementedError() # Must be defined in concrete classes
        
    
    def print_graphviz(self, out):
        pass
    

class TargetNode(Node):
    def __init__(self, target, dependencies):
        # unconventional order, but I need what I refer to when calling
        # the constructor for Node...
        self.target = target
        Node.__init__(self, target.name, dependencies)
        
    def node_should_run(self):
        return self.target.should_run()

    def print_graphviz(self, out):
    
        if self.target.should_run():
            col = 'color = red, style=bold'
        elif self.should_run:
            col = 'color = red'
        else:
            col = 'color = darkgreen'

        print >> out, '"%s"'%self.name, 
        print >> out, '[',
        print >> out, ','.join(['shape = doubleoctagon', col]),
        print >> out, ']',
        print >> out, ';'


        
class FileNode(Node):
    def __init__(self, sysfile):
        # unconventional order, but I need what I refer to when calling
        # the constructor for Node...
        self.sysfile = sysfile
        Node.__init__(self, sysfile.name, [])

    def node_should_run(self):
        return not self.sysfile.file_exists()

    def print_graphviz(self, out):
    
        if self.sysfile.file_exists():
            col = 'color = darkgreen'
        else:
            col = 'color = red, style=bold'
    
        print >> out, '"%s"'%self.name, 
        print >> out, '[',
        print >> out, ','.join(['shape = note', col]),
        print >> out, ']',
        print >> out, ';'
        

class DependencyGraph:
    '''A complete dependency graph, with code for scheduling a workflow.'''

    def __init__(self, workflow):
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
        
        # FIXME: find a better, OO, way of dispatching based on type,
        # 'cause this is definitely ugly and is not going to work once I start
        # making new types of tasks!
        
        if isinstance(task, workflow.Target):
            node = TargetNode(task, dependencies)
            self.nodes[task.name] = node
            return node
        elif isinstance(task, workflow.SystemFile):
            node = FileNode(task)
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

    def print_workflow_graph(self, out):
    	'''Print the workflow to file object out in graphviz format.'''
    	
    	print >> out, 'digraph workflow {'
    	
    	# Handle nodes
    	for node in self.nodes.values():
    	    node.print_graphviz(out)
    		
    	# FIXME: handle system files ... can't right here yet
    	
    	for src in self.nodes.values():
    	    for fname,dst in src.dependencies:
    	        print >> out, '"%s"'%src.name, '->', '"%s"'%dst.name,
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
