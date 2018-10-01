"""
Create a graph object in order to find the most propriate path.
"""
from Constant import g

class Graph():
    def __init__(self, graph_dict = None):
        if graph_dict == None: graph_dict = {}
        self.__graph_dict = graph_dict
    def vertices(self):
        return(list(self.__graph_dict.keys()))
    def find_all_paths(self, start_vertex, end_vertex, path=[]):
        graph = self.__graph_dict
        path  = path + [start_vertex]
        if start_vertex == end_vertex: return([path])
        if start_vertex not in graph: return([])
        paths = []

        if len(paths) > 0: print(paths)
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_paths = self.find_all_paths(vertex, end_vertex, path)
                for p in extended_paths:
                    #print(len(p),p)
                    paths.append(p)
        return(paths)

graph_test = Graph(g)

for key in g:
    print(key, g[key])
    print("inside grapy.py")

print(graph_test.find_all_paths("VioC_VO.pdb", "5LSQ.pdb"))



#print(graph.vertices()) # To get the vertices/Node of the graph.
