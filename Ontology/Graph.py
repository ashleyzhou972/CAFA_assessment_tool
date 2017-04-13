# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
Module containing abstract representation of graph.
"""

from functools import total_ordering

class DiGraph(object):
    """
    Base class for directed graph representation.

    Nodes' labels can be any hashable objects.

    """
    
    _REACHABLE = "reachable"
    _IS_VISITED = "is_visited"
    
    def __init__(self, edges=None):
        """
        Initialize graph with edges.
    
        Parameters
        ----------
        data - list of edges in a graph.
        """
        
        self.cycles = []
        self.attrs = {}
        self.nodes = {}
        if edges != None:
            for (u, v) in edges:
                self.add_edge(u, v)
    
    def add_edge(self, u, v, data = None):
        """
        Adds an edge u->v to the graph

        Parameters
        ----------
        u,v - nodes connected by the edge
        """
        if u not in self.nodes:
            self.add_node(u)
        if v not in self.nodes:
            self.add_node(v)
        
        u_node = self.nodes[u]
        v_node = self.nodes[v]
        v_node.pred.add(DiEdge(u_node, data))
        u_node.succ.add(DiEdge(v_node, data))
        
    def add_node(self, u, data = None):
        """
        Adds node to the graph

        Parameters
        ----------
        u - node to add
        data - node data
        """
        
        self.nodes[u] = DiNode(u, data)
        
    def node_exists(self, u):
        return u in self.nodes

    def update_node(self, u, data):
        """
        Updates node data. If node does not exists it is added.

        Parameters
        ----------
        u - node to update
        data - node data
        """
        if self.node_exists(u):
            self.nodes[u].data = data
        else:
            self.add_node(u, data)


    def get_node(self, u):
        """
        Gets node from the graph

        Parameters
        ----------
        u - node id
        """
        return self.nodes[u]
    
    def get_induced_subgraph(self, nodes_ids):
        """
        Gets graph induced by given nodes labels.
        
        Parameters
        ----------
        nodes_ids - list of nodes labels
        
        >>> g = DiGraph([(1,2), (2,3), (3,4), (3,5), (5,2), (5,6), (6,8), (6,7), (2,9), (9,2)])
        >>> print(g)
        DiGraph(nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9], edges = [1->2, 2->9, 2->3, 3->4, 3->5, 5->2, 5->6, 6->8, 6->7, 9->2])
        >>> ig = g.get_induced_subgraph([2, 3, 4, 5, 8, 9])
        >>> print(ig)
        DiGraph(nodes = [2, 3, 4, 5, 8, 9], edges = [2->9, 2->3, 3->4, 3->5, 5->2, 9->2])
        """
        
        igraph = DiGraph()
        id_set = set(nodes_ids)
        
        for nid in nodes_ids:
            # copying wanted nodes
            n = self.nodes[nid]
            igraph.update_node(nid, n.data)
            for edge in n.succ:
                if edge.to_node.label in id_set:
                    # copying edges to wanted nodes
                    igraph.add_edge(nid, edge.to_node.label, edge.data)
        return igraph
    
    def _get_reachable(self, node):
        """
        Gets ids of all nodes reachable from given node. Finds cycles along the way.

        Parameters
        ----------
        node - node which descendants we want to obtain.
        """
        if DiGraph._REACHABLE in node.attr:
            return ([], node.attr[DiGraph._REACHABLE])
        else:
            my_set = set()
            
            if DiGraph._IS_VISITED in node.attr:
                in_cycle = [node.label]
            else:
                node.attr[DiGraph._IS_VISITED] = None
                
                in_cycle = []
                for edge in node.succ:
                    up_cycle, up_set = self._get_reachable(edge.to_node)
                    if len(up_cycle) > 0:
                        if up_cycle[0] == node.label:
                            # we closed the cycle
                            self.cycles.append(up_cycle)
                        else:
                            up_cycle.append(node.label)
                            in_cycle = up_cycle
                        up_set |= my_set
                        my_set = up_set # we are binding sets of reachable nodes of every node in cycle
                    else:
                        my_set |= up_set
                    my_set.add(edge.to_node.label)
                node.attr[DiGraph._REACHABLE] = my_set
                node.attr.pop(DiGraph._IS_VISITED)
            return (in_cycle, my_set)
    
    def __repr__(self):
        first = True
        result = "DiGraph(nodes = " + str(list(self.nodes.keys())) + ", edges = ["
        for n in self.nodes.values():
            for e in n.succ:
                if not first:
                    result += ", "
                else:
                    first = False
                result += str(n.label) + str(e)
        return result + "])"
    
class DiEdge(object):
    """
    Class representing an edge in the graph.
    """
    
    def __init__(self, to_node, data):
        self.to_node = to_node
        self.data = data

    def __eq__(self, other):
        return self.to_node == other.to_node

    def __hash__(self):
        return hash(self.to_node)

    def __str__(self):
        res = ""
        if self.data != None:
            res += "-"
            res += str(self.data)
        return res + "->" + str(self.to_node.label)

    def __repr__(self):
        return str(self.to_node) + "(" + str(self.data) + ")"

class DiNode(object):
    """
    Class containing information about graph structure. Only used
    internally in DiGraph.

    Nodes with the same label are not distinguishable.
    
    >>> a = DiNode(1)
    >>> b = DiNode(1)
    >>> a == b
    True
    >>> a
    DiNode(label = 1)
    """

    def __init__(self, label, data = None):
        """
        Initialize node with data.

        Parameters
        ----------
        label - node label
        data - internal node data
        """
        self.label = label
        self.data = data
        self.pred = set()
        self.succ = set()
        self.attr = {}

    def __eq__(self, other):
        if other is None:
            return False
        else:
            return self.label == other.label

    def __hash__(self):
        return hash(self.label)

    def __str__(self):
        return "Node: " + str(self.label)

    def __repr__(self):
        return "DiNode(label = " + repr(self.label)+ ")"

if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

