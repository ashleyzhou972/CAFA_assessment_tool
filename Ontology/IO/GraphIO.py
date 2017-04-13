# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
I/O and visualization operations on graphs.
"""
from __future__ import print_function
from Bio._py3k import basestring

import sys
from .Interfaces import OntoWriter

_INDENT = "  "
        
class GmlWriter(OntoWriter):
    """
    Writes graph to gml format.
    
    Simple gml writer. Writes gml files given DiGraph instance. It stores nodes
    and edges with their respective data attribute. Note that data attribute
    is stored only if it is instance of dict class.

    >>> from Bio.Ontology.Graph import DiGraph
    >>> from Bio._py3k import StringIO
    >>> out = StringIO()
    
    Constructing the graph:
    >>> g = DiGraph([(1, 2), (2,3), (3,2)])
    >>> g.update_node(2, {"a" : {"b" : "string"}, "c" : 3})

    Writing the graph:
    >>> writer = GmlWriter(out)
    >>> writer.write(g)
    
    >>> print(out.getvalue())
    graph [
      directed 1
      node [
        id 0
        label 1
      ]
      node [
        id 1
        label 2
        a [
          b "string"
        ]
        c 3
      ]
      node [
        id 2
        label 3
      ]
      edge [
        source 0
        target 1
      ]
      edge [
        source 1
        target 2
      ]
      edge [
        source 2
        target 1
      ]
    ]
    """
    
    def __init__(self, file_handle):
        self.handle = file_handle
                
    def data_to_gml(self, label, data, indent):
        lines = []
        if isinstance(data, dict):
            lines.append(_INDENT * indent + str(label) + " [")
            for k, v in data.items():
                lines += self.data_to_gml(k, v, indent + 1)
            lines.append(_INDENT * indent + "]")
        elif data != None:
            if isinstance(data, basestring) or isinstance(data, list):
                rdata = '"' + str(data) + '"'
            else:
                rdata = str(data)
            lines.append(_INDENT * indent + label + " " + rdata)
        return lines
    
    def get_lines(self, graph):
        i = 0
        node_ids = {}
        edges = []
        
        lines = ["graph [", _INDENT + "directed 1"]
        
        for k, v in graph.attrs.items():
            lines += self.data_to_gml(k, v, 1)
            
        for label, node in graph.nodes.items():
            lines.append(_INDENT + "node [")
            lines.append(_INDENT * 2 + "id " + str(i))
            lines.append(_INDENT * 2 + "label \"" + str(label) + '"')
            if isinstance(node.data, dict):
                for k, v in node.data.items():
                    lines += self.data_to_gml(k, v, 2)
            lines.append(_INDENT + "]")
            
            for edge in node.succ:
                edges.append((i, edge.to_node.label, edge.data))
            node_ids[label] = i
            i += 1
    
        for src, target_label, data in edges:
            lines.append(_INDENT + "edge [")
            lines.append(_INDENT * 2  + "source " + str(src))
            lines.append(_INDENT * 2  + "target " + str(node_ids[target_label]))
            if isinstance(data, dict):
                for k, v in data.items():
                    lines += self.data_to_gml(k, v, 2)
            else:
                # by default it will be label
                lines += self.data_to_gml("label", data, 2)
            lines.append(_INDENT + "]")
        lines.append("]")
        return lines
    
    def write(self, graph):
        self.handle.write("\n".join(self.get_lines(graph)))

class GraphVisualizer(OntoWriter):
    """
    Stores graph in png format using graphviz library.
    """
    
    def __init__(self, file_handle, dpi = 96, color = "#00bd28"):
        self.handle = file_handle
        self.dpi = dpi
        self.color = color
    
    def term_to_label(self, term):
        return "{0}\n{1}".format(term.id, term.name)

    def to_printable_graph(self, graph):
        try:
            import pygraphviz
        except ImportError:
            print("Error while printing. To use this functionality you need to have pygraphviz installed.", file=sys.stderr)
        else:
            viz_graph = pygraphviz.AGraph()
            viz_graph.graph_attr.update(dpi = str(self.dpi))
            viz_graph.node_attr.update(shape="box", style="rounded,filled")
            viz_graph.edge_attr.update(shape="normal", color="black", dir="back")
            
            entry_labels = {}
            
            for node in graph.nodes.values():
                new_label = self.term_to_label(node.data)
                viz_graph.add_node(new_label, fillcolor = self.color)
                entry_labels[node.label] = new_label
                
            for label, node in graph.nodes.items():
                for edge in node.succ:
                    viz_graph.add_edge(entry_labels[edge.to_node.label], entry_labels[label],  label=edge.data)
                    
            return viz_graph
        
    def write(self, graph):
        vg = self.to_printable_graph(graph)
        vg.draw(self.handle, prog="dot")


if __name__ == '__main__':
    from Bio._utils import run_doctest
    run_doctest()

