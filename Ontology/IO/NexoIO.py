# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
I/O operations for NeXO.
"""

import xml.sax
import collections
import ast
from Ontology.Data import OntologyGraph, OntologyTerm, TermAssociation, GeneAnnotation
from .Interfaces import OntoReader

_SKIP = 0
_NODE = 1
_EDGE = 2

class NexoContentHandler(xml.sax.ContentHandler):
    def __init__(self, get_all_attrs, annotation_source):
        xml.sax.ContentHandler.__init__(self)
        
        self.get_all_attrs = get_all_attrs
        self.annotation_source = annotation_source
        
        self.state = _SKIP
        self.old_state = _SKIP
        
        self.nodes = {}
        self.annotations = collections.defaultdict(list)
        self.edges = []
        
        self.current_term = None
        self.current_edge = None
    
    def _split_list(self, val):
        if val.startswith('['):
            return ast.literal_eval(val)
        else:
            return [val]
        
    def startElement(self, name, attrs):
        if name == "node":
            self.state = _NODE
            term_id = attrs["label"]
            self.current_term = OntologyTerm(term_id, term_id, {})
            self.nodes[attrs["id"]] = self.current_term

        elif name == "edge":
            self.state = _EDGE
            self.current_edge = [attrs["source"], attrs["target"]]
            self.edges.append(self.current_edge)
            
        elif name == "graphics":
            self.old_state = self.state
            self.state = _SKIP
            
        elif name == "att":
            if self.state == _NODE:
                if attrs["name"] == "Term":
                    val = attrs.get("value")
                    if val != None:
                        self.current_term.name = val
                elif (attrs["name"] == "Assigned Genes" and self.annotation_source == "genes")\
                 or (attrs["name"] == "Assigned Orfs" and self.annotation_source == "orfs"):
                    val = attrs.get("value")
                    if val != None:
                        for gene in self._split_list(val):
                            self.annotations[gene].append(self.current_term.id)
                elif self.get_all_attrs:
                    val = attrs.get("value")
                    if val != None:
                        self.current_term.attrs[attrs["name"]] = val
            elif self.state == _EDGE:
                if attrs["name"] == "NeXO relation type":
                    self.current_edge.append(attrs.get("value"))
            
    def endElement(self, name):
        if name == "node" or name == "edge":
            self.state = _SKIP
        elif name == "graphics":
            self.state = self.old_state
 
    def characters(self, content):
        pass
    
class NexoReader(OntoReader):
    """
    Class for reading Nexo xgmml network.
    """


    def __init__(self, file_handle,  get_all_attrs = False, annotation_source = "genes"):
        self.handle = file_handle
        self.get_all_attrs = get_all_attrs
        self.annotation_source = annotation_source
    
    def read(self):
        """
        Returns gene annotation list and ontology graph read from nexo file.
        """
        
        content_handler = NexoContentHandler(self.get_all_attrs, self.annotation_source)
        xml.sax.parse(self.handle, content_handler)
        
        annotations = []
        for obj, assocs in content_handler.annotations.items():
            annotations.append(GeneAnnotation(obj,
                        associations = [TermAssociation(x) for x in assocs]))
        
        graph = OntologyGraph()

        for _, node in content_handler.nodes.items():
            graph.add_node(node.id, node)
        
        edge_types = set()
        for edge in content_handler.edges:
            source = content_handler.nodes[edge[0]].id
            target = content_handler.nodes[edge[1]].id
            graph.add_edge(target, source, edge[2]) # in our representation it is inverted
            edge_types.add(edge[2])
            
        for edge_type in edge_types:
            graph.typedefs[edge_type] = {"id" : edge_type}
            
        return (annotations, graph)
        
