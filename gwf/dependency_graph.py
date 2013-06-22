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
        self.root = node

    def print_dependency_graph(self):
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
        '''Linearize the targets to be run.'''
        
        assert self.root is not None
        scheduled = set()
        schedule = []
        
        def dfs(node):
            if node in scheduled:
                return
            if node.target.should_run():
                for dep in node.dependencies:
                    dfs(dep)
                schedule.append(node)
                scheduled.add(node)

        if self.root.target.should_run():
            dfs(self.root)
        return schedule
