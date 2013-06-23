'''Graph for dependency relationships of targets.'''

class Node:
    '''A node in the dependencies DAG.'''

    def __init__(self, target, dependencies):
        self.target = target
        self.dependencies = dependencies

        

class DependencyGraph:
    '''A complete dependency graph, with code for scheduling a workflow.'''

    def __init__(self):
        self.nodes = dict()
        self.root = None

    def add_node(self, name, target, dependencies):
        node = Node(target, dependencies)
        self.nodes[name] = node
        return node

    def has_node(self, name):
        return name in self.nodes

    def get_node(self, name):
        return self.nodes[name]

    def set_root(self, node):
        '''Make "node" the root of the graph. The root is used when computing
        the scheduled scripts.'''
        self.root = node

    def print_dependency_graph(self):
        '''Prints the graph to stdout.  A very simple function mostly useful
        for debugging and not really for user output.'''
        
        assert self.root is not None
        printed = set()
        def dfs(node, indent=''):
            print indent, node.target.name, node.target.should_run(),
            if node in printed:
                print '[...]'
                return
            else:
                print # add newline if we recurse...
                
            printed.add(node)
            for dep in node.dependencies:
                dfs(dep, indent+'\t')

        dfs(self.root)
                

    def schedule(self):
        '''Linearize the targets to be run.
        
        Returns a list of tasks to be run (in the order they should run or
        be submitted to the cluster to make sure dependences are handled
        correctly) and a set of the names of tasks that will be scheduled
        (to make sure dependency flags are set in the qsub command).
        
        '''
        
        assert self.root is not None
        
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
            for dep in node.dependencies:
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

        dfs(self.root)
            
        return schedule, scheduled
