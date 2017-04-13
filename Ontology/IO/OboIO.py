# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
I/O operations for ontologies.
"""

import collections
import re, shlex
from Ontology.Data import OntologyTerm, OntologyGraph
from .Interfaces import OntoIterator, OntoReader, OntoWriter

_START = 0
_STANZA = 1
_READ_STANZA = 2
_EOF = 3

_IS_A_TYPE = {"id" : "is_a",
              "name" : "is_a",
              "range": "OBO:TERM_OR_TYPE",
              "domain" : "OBO:TERM_OR_TYPE",
              "def" : "The basic subclassing relationship [OBO:defs]" }

class OboWriter(OntoWriter):
    """
    Writes obo files.
    
    Writes OntologyTerms to obo files.
    """
    
    def __init__(self, file_handle, version = "1.2"):
        self.handle = file_handle
        self.version = version
    
    def write(self, terms_list):
        # now only terms are valid for writing
        self.handle.write("format-version:" + self.version + "\n")
        for term in terms_list:
            if isinstance(term, OntologyTerm):
                self.handle.write("\n[Term]\n")
                self.handle.write("id: " + term.id + "\n")
                self.handle.write("name: " + term.name + "\n")
                for k, v in term.attrs.items():
                    for vi in v:
                        self.handle.write(k + ": " + vi + "\n")
    
class OboIterator(OntoIterator):
    """
    Parses obo files.
    """
    
    def __init__(self, file_handle):
        self._handle = file_handle
        self._reg = re.compile("^\[([\w]+)\]$")
        self._state = _START
        self._dict = None
        self._found_stanza_type = None
    
    def _strip_comment(self, line):
        pos = line.find('!')
        if pos > 0:
            return line[0:pos]
        else:
            return line
    
    def _read_line(self):
        iterate = True
        result = ''
        while iterate:
            line = self._handle.readline()
            if line:
                line = self._strip_comment(line)
                pos = line.find("\\\n")
                if pos > 0:
                    line = line[:pos]
                else:
                    iterate = False
                result += line
            else:
                return None
        return result
    
    def _split_tag(self, line):
        pos = line.find(':')
        if pos < 0:
            raise ValueError(("Invalid obo file: Incorrect tag: ':'"
                              " expected in line '{0}'.").format(line))
        return (line[0:pos].strip(), line[(pos+1):].strip())
    
    def _read_stanza(self):
        while self._state != _STANZA and self._state != _EOF:
            line = self._read_line()
            if line is None:
                self._state = _EOF
            else:
                line = line.strip()
                match = self._reg.match(line)
                if match != None:
                    self._state = _STANZA
                    self._found_stanza_type, = match.groups()
                elif len(line) > 0 and self._state == _READ_STANZA:
                    k, v = self._split_tag(line)
                    self._dict[k].append(v)


    def __next__(self):
        return self.next()

    def next(self):
        self._read_stanza()
        if self._state == _EOF:
            raise StopIteration
        self._dict = collections.defaultdict(list)
        self._state = _READ_STANZA
        stanza_type = self._found_stanza_type
        self._read_stanza()
        return (stanza_type, self._dict)

def terms_to_graph(terms):
    """
    Crates OntologyGraph from terms list obtained from OboIterator
    """
    
    g = OntologyGraph()
    defined_relations = set()
    found_relations = set()
    
    for (term_type, data) in terms:
        if term_type == "Term": # Add only terms and typedefs for now
            nid = data.pop("id")[0]
            name = data.pop("name")[0]
            term = OntologyTerm(nid, name, data)
            if g.node_exists(nid):
                g.update_node(nid, term)
            else:
                g.add_node(nid, term)
            if "is_a" in data:
                for edge in data["is_a"]:
                    g.add_edge(nid, edge, "is_a")
            if "synonym" in data:
                node = g.get_node(nid)
                for synonym in data["synonym"]:
                    g.synonyms[shlex.split(synonym)[0]] = node
            if "relationship" in data:
                for edge in data["relationship"]:
                    p = edge.split()
                    if len(p) == 2:
                        g.add_edge(nid, p[1], p[0])
                        found_relations.add(p[0])
                    else:
                        raise ValueError("Incorrect relationship: " + edge)
            if "alt_id" in data:
                for alt_id in data["alt_id"]:
                    g.alt_ids[alt_id] = nid
            if "namespace" in data:
                g.namespace[nid] = data.pop("namespace")[0]

        elif term_type == "Typedef":
            rid = data["id"][0]
            g.typedefs[rid] = data
            defined_relations.add(rid)
    
    # validate whether all relationships were defined
    not_defined = found_relations.difference(defined_relations)
    if len(not_defined) > 0:
        raise ValueError("Undefined relationships found: " + str(not_defined))
    
    g.typedefs["is_a"] = _IS_A_TYPE
    for alt_id in g.alt_ids:
        g.namespace[alt_id] = g.namespace[g.alt_ids[alt_id]]
    return g

class OboReader(OntoReader):
    """
    Reads obo file to OntologyGraph.
    """

    def __init__(self, file_handle):
        self._handle = file_handle
    
    
    def read(self):
        """
        Returns obo file representation as OntologyGraph.
        """
        return terms_to_graph(OboIterator(self._handle))
