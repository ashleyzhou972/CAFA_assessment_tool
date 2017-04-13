# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

from __future__ import print_function
from Bio._py3k import range

from Ontology.IO.GraphIO import GmlWriter
from Ontology.Graph import DiGraph
from Ontology.Stats import corrections_labels
import sys

def rgb_to_triple(rgb):
    """
    Returns triple of ints from rgb color in format #xxxxxx.
    
    >>> rgb_to_triple("#00bd28")
    (0, 189, 40)
    
    """
    if rgb[0] == '#':
        return (int(rgb[1:3], 16), int(rgb[3:5], 16), int(rgb[5:8], 16))
    else:
        raise ValueError("Not an rgb value.")

def triple_to_rgb(triple):
    """
    Returns rgb color in format #xxxxxx from triple of ints
    
    >>> triple_to_rgb((0, 189, 40))
    '#0bd28'
    
    """
    r, g, b = triple
    return "#{0}{1}{2}".format(hex(r)[2:], hex(g)[2:], hex(b)[2:])

def get_gradient(color_a, color_b, k):
    """
    Returns gradient of colors from a to b with k steps
    
    >>> get_gradient("#aabbcc", "#001122", 10)
    ['#aabbcc', '#99aabb', '#8899aa', '#778899', '#667788', '#556677', '#445566', '#334455', '#223344', '#112233']
    
    """
    r1, g1, b1 = rgb_to_triple(color_a)
    r2, g2, b2 = rgb_to_triple(color_b)
    d = float(k)
    r2 = (r2 - r1) / d
    g2 = (g2 - g1) / d
    b2 = (b2 - b1) / d
    grad = []
    for _ in range(k):
        grad.append(triple_to_rgb((int(r1), int(g1), int(b1))))
        r1 += r2
        g1 += g2
        b1 += b2
    return grad

def get_gradient_index(val, min_val, max_val, steps):
    """
    Returns the index in gradient given p-value, minimum value, maximum value
    and number of steps in gradient.
    
    >>> get_gradient_index(0.002, 0.001, 0.0029, 10)
    5
    >>> get_gradient_index(0.0011, 0.001, 0.0029, 10)
    0
    >>> get_gradient_index(0.002899999, 0.001, 0.0029, 10)
    9
    >>> get_gradient_index(0.0029, 0.001, 0.0029, 10)
    9
    >>> get_gradient_index(0.001, 0.001, 0.001, 10)
    0
    """
    if max_val == min_val:
        return 0
    else:
        d = (max_val - min_val) / steps
        r = int((val - min_val) / d)
        return r if r < steps else r - 1

def print_enrichment_chart(file_handle, vals, title):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Error while printing. To use this functionality you need to have matplotlib installed.", file=sys.stderr)
    else:
        fig, ax1 = plt.subplots()
        
        xs = list(range(len(vals)))
        ys =  vals
        
        ax1.plot(xs, ys)
        
        bar_ys = [int(ys[0] > 0)]
        for i in range(1, len(ys)):
            bar_ys.append(int(ys[i] > ys[i - 1]))
        bar_ys = [bar_ys]
        
        pos = ax1.axes.get_position()
        
        ax0 = fig.add_axes([pos.x0, pos.y1, pos.width, 0.1])
        
        ax0.imshow(bar_ys, cmap=plt.cm.Blues, interpolation='nearest')
        ax0.axes.get_yaxis().set_visible(False)
        ax0.axes.get_xaxis().set_visible(False)
        ax0.set_title(title)
        
        plt.savefig(file_handle, bbox_inches=0)
        plt.close()

class PrettyPrinter(object):
    """
    Base class for printers.
    """
    
    def __init__(self, file_handle):
        self.handle = file_handle
        
    def pretty_print(self, enrichment, graph):
        raise NotImplementedError("pretty_print not implmented yet")


class GmlPrinter(PrettyPrinter):
    """
    Stores found enrichments as a graph in gml format.
    """
    
    def __init__(self, file_handle, gradient_step = 10, color_a = "#7eff00",
                 color_b = "#f8ff8d", color_none = "#c3c3c3"):
        
        self.handle = file_handle
        self.gradient_step = gradient_step
        self.color_none = color_none
        self.gradient = get_gradient(color_a, color_b, self.gradient_step)
        self.min_p = 1.0
        self.max_p = 0.0
        
    def to_printable_data(self, e_entry):
        return {"name" : e_entry.name,
                "pvalue" : e_entry.p_value,
                "graphics" : { "fill" : self.gradient[get_gradient_index(e_entry.p_value, self.min_p, self.max_p, self.gradient_step)]
                              }
                }
    def term_to_printable(self, term):
        return {"name" : term.name,
                "graphics" : { "fill" : self.color_none
                              }
                }
    def entry_to_label(self, entry):
        return str(entry.id)
    
    def term_to_label(self, term):
        return str(term.id)
    
    def to_printable_graph(self, enrichment, graph):
        viz_graph = DiGraph()
        viz_graph.attrs["defaultnodesize"] = "labelsize"
        viz_graph.attrs["label"] = str(enrichment)
        
        entry_labels = {}
        
        for entry in enrichment.entries:
            if entry.p_value < self.min_p:
                self.min_p = entry.p_value
            if entry.p_value > self.max_p:
                self.max_p = entry.p_value
        
        for entry in enrichment.entries:
            new_label = self.entry_to_label(entry)
            viz_graph.add_node(new_label, self.to_printable_data(entry))
            entry_labels[entry.id] = new_label
        
        for label, node in graph.nodes.items():
            if label not in entry_labels:
                new_label = self.term_to_label(node.data)
                viz_graph.add_node(new_label, self.term_to_printable(node.data))
                entry_labels[label] = new_label
        
        for label, u in graph.nodes.items():
            for edge in u.succ:
                viz_graph.add_edge(entry_labels[label], entry_labels[edge.to_node.label])
        
        return viz_graph
    
    def pretty_print(self, enrichment, graph):
        nodes_ids = set()
        for x in enrichment.entries:
            nodes_ids.add(x.id)
            nodes_ids = nodes_ids.union(graph.get_ancestors(x.id))
        
        g = graph.get_induced_subgraph(nodes_ids)
        
        vg = self.to_printable_graph(enrichment, g)
        GmlWriter(self.handle).write(vg)

class GraphVizPrinter(PrettyPrinter):
    """
    Stores found enrichments as visualization in png format using graphviz library.
    """
    
    def __init__(self, file_handle, gradient_step = 10, color_a = "#7eff00",
                 color_b = "#f8ff8d", color_none = "#c3c3c3", dpi = 96):
        self.handle = file_handle
        self.color_none = color_none
        self.gradient_step = gradient_step
        self.gradient = get_gradient(color_a, color_b, self.gradient_step)
        self.dpi = str(dpi)
        self.min_p = 1.0
        self.max_p = 0.0
        
    def entry_to_label(self, entry):
        return "{0}\n{1}\np:{2}".format(entry.id, entry.name, entry.p_value)

    def term_to_label(self, term):
        return "{0}\n{1}".format(term.id, term.name)
    
    def to_printable_graph(self, enrichment, graph):
        try:
            import pygraphviz
        except ImportError:
            print("Error while printing. To use this functionality you need to have pygraphviz installed.", file=sys.stderr)
        else:
            viz_graph = pygraphviz.AGraph()
            viz_graph.graph_attr.update(dpi = self.dpi)
            viz_graph.node_attr.update(shape="box", style="rounded,filled")
            viz_graph.edge_attr.update(shape="normal", color="black", dir="back")
            
            entry_labels = {}
            
            for entry in enrichment.entries:
                if entry.p_value < self.min_p:
                    self.min_p = entry.p_value
                if entry.p_value > self.max_p:
                    self.max_p = entry.p_value
            
            for entry in enrichment.entries:
                new_label = self.entry_to_label(entry)
                col = self.gradient[get_gradient_index(entry.p_value, self.min_p, self.max_p, self.gradient_step)]
                viz_graph.add_node(new_label, fillcolor = col)
                entry_labels[entry.id] = new_label
            
            for label, node in graph.nodes.items():
                if label not in entry_labels:
                    new_label = self.term_to_label(node.data)
                    viz_graph.add_node(new_label, fillcolor = self.color_none)
                    entry_labels[label] = new_label
            
            for label, u in graph.nodes.items():
                for edge in u.succ:
                    viz_graph.add_edge(entry_labels[edge.to_node.label], entry_labels[label],  label=edge.data)
                    
            return viz_graph
            
            
        
    def pretty_print(self, enrichment, graph):
        nodes_ids = set()
        for x in enrichment.entries:
            nodes_ids.add(x.id)
            nodes_ids = nodes_ids.union(graph.get_ancestors(x.id))
        g = graph.get_induced_subgraph(nodes_ids)
        
        vg = self.to_printable_graph(enrichment, g)
        vg.draw(self.handle, prog="dot")

class TxtPrinter(PrettyPrinter):
    """
    Prints found enrichments to txt file.
    """
    
    def __init__(self, file_handle):
        self.handle = file_handle
        
    def pretty_print(self, enrichment, graph):
        self.handle.write("Enrichments found using {0} method.\n\nEnrichments:\n\n"
                          .format(enrichment.method))
        sorted_entries = sorted(enrichment.entries, key = lambda x: x.p_value)
        for x in sorted_entries:
            self.handle.write(str(x))
            self.handle.write("\n\n")
        if (len(enrichment.warnings) > 0):
            self.handle.write("Warnings:\n")
            for x in enrichment.warnings:
                self.handle.write(str(x))

_DEFAULT_STYLE = """<style type="text/css">
.warning
{
    color:#D13F31;
}
body
{
    margin:45px;
}
h1
{
    color:#669;
}
table
{
    font-family: Sans-Serif;
    font-size: 14px;
    background: #fff;
    margin: 45px 0px;
    border-collapse: collapse;
    text-align: left;
}
th
{
    font-size: 16px;
    font-weight: normal;
    color: #009;
    padding: 12px 10px;
    border-bottom: 2px solid #6678b0;
}
td
{
    color: #669;
    padding: 8px 10px;
    border-bottom: 1px solid #ccc;
}
tbody tr:hover td
{
    color: #009;
}
</style>
"""

class HtmlPrinter(PrettyPrinter):
    """
    Prints found enrichments to html file.
    """
    
    def __init__(self, file_handle, go_to_url = "#",
                 style = _DEFAULT_STYLE):
        self.handle = file_handle
        self.go_to_url = go_to_url
        self.style = style
    
    def write_tag(self, tag, text, attrs = None):
        self.open_tag(tag, attrs)
        self.handle.write(text)
        self.close_tag(tag)
    
    def open_tag(self, tag, attrs = None):
        ot = "<" + tag
        if attrs != None:
            for k, v in attrs.items():
                ot += ' {0}="{1}"'.format(k,v)
        ot += ">\n"
        self.handle.write(ot)
        
    def close_tag(self, tag):
        self.handle.write("</" + tag + ">\n")
    
    def pretty_print(self, enrichment, graph):
        self.handle.write("<!DOCTYPE html>")
        
        self.open_tag("html")
        self.open_tag("head")
        self.handle.write(self.style)
        self.close_tag("head")
        self.open_tag("body")
        self.write_tag("h1", "Enrichments found using {0} method."
                          .format(enrichment.method))
        
        sorted_entries = sorted(enrichment.entries, key = lambda x: x.p_value)
        
        self.open_tag("table")
        
        self.open_tag("tr")
        headers = ["ID", "name", "p-value"]
        for x in enrichment.corrections:
            if x in corrections_labels:
                headers.append(corrections_labels[x])
            else:
                headers.append(x)
        for header in headers:
            self.write_tag("th", header)
        self.close_tag("tr")
        
        for x in sorted_entries:
            self.open_tag("tr")
            self.open_tag("td")
            self.write_tag("a", str(x.id), {"href" : self.go_to_url + str(x.id)})
            self.close_tag("td")
            self.write_tag("td", str(x.name))
            self.write_tag("td", str(x.p_value))
            for cr in enrichment.corrections:
                self.write_tag("td", str(x.corrections[cr]))
            self.close_tag("tr")
            
        self.close_tag("table")
        
        if (len(enrichment.warnings) > 0):
            self.write_tag("h1", "Warnings:", {"class" : "warning"})
            self.open_tag("ul", {"class" : "warning"})
            for x in enrichment.warnings:
                self.write_tag("li", str(x))
            self.close_tag("ul")
        
        self.close_tag("body")
        self.close_tag("html")

if __name__ == '__main__':
    from Bio._utils import run_doctest
    run_doctest()
