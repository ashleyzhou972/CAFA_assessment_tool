# Copyright 2013 by Kamil Koziara. All rights reserved
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Additions by Iddo Friedberg, 2016
# 2/2016: added support for edge types in get_parents and 
#         _get_reachable

"""
Module with classes representing ontologies and annotation data.
"""
from __future__ import print_function

#from Bio.Ontology.Graph import DiGraph
from Ontology.Graph import DiGraph
import copy
import sys
from collections import defaultdict

class OntologyGraph(DiGraph):
    """
    Represents ontology graph.
    """

    def __init__(self):
        DiGraph.__init__(self)
        self.typedefs = {}
        self.synonyms = {}
#        self.alt_ids = defaultdict(set)
        self.alt_ids = {}
        self.namespace = {}
        
    def get_node(self, u):
        x = self.nodes.get(u)
        if x is None:
            return self.nodes.get(self.alt_ids[u])
#            return self.synonyms.get(u)
        else:
            return x
    
    def node_exists(self, u):
        return u in self.nodes or u in self.alt_ids
    
    def get_term(self, oid):
        return self.get_node(oid).data
    
    def get_ancestors(self, oid):
        node = self.get_node(oid)
        _, res = self._get_reachable(node)
        return copy.copy(res) # return copy so user won't disturb cached values

    def get_ids(self):
        return list(self.nodes) + list(self.alt_ids)

    def get_namespace(self,oid):
        return self.namespace[oid]

#    def get_parents(self, oid):
        # Old get_parents, takes any edge type
#        node = self.get_node(oid)
#        res = set()
#        for edge in node.succ:
#            res.add(edge.to_node.label)
#        return res
    
    def get_parents(self, oid, accepted_edges=set({"is_a","part_of"})):
        """
        Gets parent node, along accepted edges.
        If the accepted_edges set is empty, all edge types are acceptable.
        """
        node = self.get_node(oid)
        res = set()
        for edge in node.succ:
            if (not accepted_edges) or (accepted_edges and (edge.data in accepted_edges)):
                res.add(edge.to_node.label)
        return res

    def get_relationship_types(self):
        return list(self.typedefs.keys())
        
    def trim(self, kept_edges):
        """
        Returns graph with only given edges left.
        
        Params:
        kept_edges - iterable specyfing edges we want to preserve
        """
        
        filter_set = set(kept_edges)
        
        fgraph = OntologyGraph()
        
        for label, node in self.nodes.items():
            # copying all the nodes
            fgraph.update_node(label, node.data)
            for edge in node.succ:
                # copying only the wanted edges
                if edge.data in filter_set:
                    fgraph.add_edge(label, edge.to_node.label, edge.data)
        
        fgraph.synonyms = dict(self.synonyms)
        
        # copying wanted relationships definitions
        for rel, typedef in self.typedefs.items():
            if rel in filter_set:
                fgraph.typedefs[rel] = typedef
            
        return fgraph
    
    def get_induced_subgraph(self, nodes_ids):
        """
        Returns graph with only given nodes left
        """
        
        idg = super(OntologyGraph, self).get_induced_subgraph(nodes_ids)
        
        igraph = OntologyGraph()
        igraph.nodes = idg.nodes
        igraph.synonyms = copy.copy(self.synonyms)
        igraph.typedefs = copy.copy(self.typedefs)
        
        return igraph
    
    def to_networkx(self, annotations = None):
        """
        Exports OntologyGraph to networkx DiGraph object
        
        Parameters:
        - annotations - (optional) annotations from term to list of genes in
            form of a dictionary {"term_id" : ['gene_a', ...]}
        """
        
        try:
            import networkx as nx
        except ImportError:
            print("Error while exporting. To use this functionality you need to have networkx installed.", file=sys.stderr)
        else:
            # Copy the graph structure
            nxgraph = nx.classes.DiGraph()
            for v in self.nodes.values():
                attrs = dict(v.data.attrs)
                attrs['name'] = v.data.name
                nxgraph.add_node(v.label, **attrs)
                for edge in v.succ:
                    nxgraph.add_edge(v.label, edge.to_node.label, relation = edge.data)
            # Add annotations
            if annotations != None:
                for k, v in annotations.items():
                    if k in nxgraph.node:
                        nxgraph.node[k]['annotated_genes'] = v
            # Add rest of the data
            nxgraph.graph['typedefs'] = self.typedefs
            nxgraph.graph['synonyms'] = self.synonyms
            
            return nxgraph
            
    def _get_reachable(self, node, accepted_edges=set({"part_of","is_a"})):
        """
        Gets ids of all nodes reachable from given node. Finds cycles along the way.
        Uses only accepted edges
        If the accepted_edges set is empty or None, all edges are accepted.

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
                    if accepted_edges and (edge.data not in accepted_edges):
                        continue
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

class OntologyTerm(object):
    """
    Represents ontology term.
    """
    
    def __init__(self, term_id, term_name, attrs = {}):
        self.id = term_id
        self.name = term_name
        self.attrs = attrs
        
    def __str__(self):
        s = "[Term]\n"
        s += "id: " + self.id + "\n"
        s += "name: " + self.name + "\n"
        for k, v in self.attrs.items():
            for vi in v:
                s+= k + ": " + vi + "\n"
        return s
            
    def __repr__(self):
        return "OntologyTerm(id = " + self.id + ", name = " + self.name + ")" 

class GeneAnnotation(object):
    """
    Represents one generic gene annotation object
    """
    
    def __init__(self, gene_id, associations = [], attrs = {}):
        self.id = gene_id
        self.associations = associations
        self.attrs = attrs
            
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __repr__(self):
        return "GeneAnnotation(db_object_id = {0})".format(self.id)
    
    def __str__(self):
        s = "DB Object ID: " + self.id + "\n"
        for k, v in self.attrs.items():
            s += k + ": " + str(v)+ "\n"       
        if len(self.associations) > 0:
            s += "\nAssociations:\n"
            for a in self.associations:
                s += str(a) + "\n"
        return s
    
class TermAssociation(object):
    """
    Represents one gene to ontology term association
    """
    
    def __init__(self, term_id, attrs = {}):
        self.term_id = term_id
        self.attrs = attrs
        
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __repr__(self):
        return "TermAssociation(id = {0})".format(self.term_id)
    
    def __str__(self):
        s = "ID: " + self.term_id + "\n"
        for k, v in self.attrs.items():
            s += k + ": " + str(v) + "\n"
        return s
