
from workflow import SystemFile, Target

def task_shape(task):
    '''Dispatching types of tasks to graphviz shapes...
    
    There's a tight coupling here that I'm not totally happy about,
    but the alternative is to have tasks contain code for giving the shapes.
    '''
    if isinstance(task, SystemFile):
        return 'note'
    
    if isinstance(task, Target):
        if len(task.input) == 0 and len(task.output) > 0:
            return 'invhouse'
        if len(task.output) == 0 and len(task.input) > 0:
            return 'house'
        return 'octagon'
        
    return 'box'

def print_node(node, out):
    '''Print the graphviz description of this node to "out".'''
    
    shape = 'shape = %s' % task_shape(node.task)
    
    if node.task.job_in_queue:
        col = 'color = orange, style=bold'
    elif node.task.should_run:
        col = 'color = red, style=bold'
    elif node.should_run:
        col = 'color = red'
    else:
        col = 'color = darkgreen'

    print >> out, '"%s"'%node.name, 
    print >> out, '[',
    print >> out, ','.join([shape, col]),
    print >> out, ']',
    print >> out, ';'

def print_graphviz(graph, out):
    '''Print the dependency graph, graph, to output stream out.'''
    
    print >> out, 'digraph workflow {'
    	
    # Handle nodes
    for node in graph.nodes.values():
        print_node(node, out)
    		    	
    # Handle edges
    for src in graph.nodes.values():
        for fname,dst in src.dependencies:
            print >> out, '"%s"'%dst.name, '->', '"%s"'%src.name,
            print >> out, '[label="%s"]' % fname,
            print >> out, ';'
    	        

    print >> out, '}'

def print_graphviz_dependencies(graph, target, out):
    '''Print the dependency graph for a given target.'''
    
    # collect the dependencies
    dependencies = set()
    def dfs(node):
        if node in dependencies:
            return # already processed
        
        dependencies.add(node)
        for _,dep in node.dependencies:
            dfs(dep)
    dfs(target)
       
    print dependencies
    print >> out, 'digraph workflow {'
    	
    # Handle nodes
    for src in dependencies:
        print_node(src, out)
        for fname,dst in src.dependencies:
            print >> out, '"%s"'%dst.name, '->', '"%s"'%src.name,
            print >> out, '[label="%s"]' % fname,
            print >> out, ';'

    print >> out, '}'
    