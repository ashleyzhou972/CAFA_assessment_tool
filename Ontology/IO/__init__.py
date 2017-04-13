# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
Module giving acces to I/O modules through generic methods.
"""


from Bio.File import as_handle
from Bio._py3k import basestring

from . import OboIO
from . import GoaIO
from . import GraphIO
from . import PrettyIO
from . import NexoIO
from . import EnrichmentIO

_FormatToIterator = { "obo" : OboIO.OboIterator,
                      "tsv" : GoaIO.TsvIterator}

_FormatToReader = { "nexo" : NexoIO.NexoReader,
                    "obo"  : OboIO.OboReader,
                    "etsv"  : EnrichmentIO.EnrichmentReader,
                    "gaf" : GoaIO.GafReader }

_FormatToWriter = { "png" : GraphIO.GraphVisualizer,
                    "etsv" : EnrichmentIO.EnrichmentWriter}

_FormatToPrinter = {"gml" : PrettyIO.GmlPrinter,
                    "png" : PrettyIO.GraphVizPrinter,
                    "txt" : PrettyIO.TxtPrinter,
                    "html": PrettyIO.HtmlPrinter}

def write(data, handle, file_format, **params):
    """
    Writes given data to file.

    Parameters:
     - data - data to write to a file,
     - handle - File handle object to write to, or filename as string
                   (note older versions of Biopython only took a handle),
     - file_format - lower case string describing the file format to write,
         Formats:
             - png - writes picture of graph to png format (this feature needs
               pygraphviz to be installed)
             - etsv
     - params - additional parameters
     
    You should close the handle after calling this function.

    """
    
    if not isinstance(file_format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not file_format:
        raise ValueError("Format required (lower case string)")

    with as_handle(handle, 'w') as fp:
        #Map the file format to a writer class
        if file_format in _FormatToWriter:
            writer_class = _FormatToWriter[file_format]
            writer_class(fp, **params).write(data)
        else:
            raise ValueError("Unknown format '%s'" % file_format)

def read(handle, file_format, **params):
    """
    Read file in given format.
    
    Parameters:
     - handle - File handle object to read from, or filename as a string,
     - file_format - lower case string describing the file format to write,
         Formats:
             - nexo
             - obo
             - etsv
             - gaf
     - params - additional parameters

    You should close the handle after calling this function.
    """

    if not isinstance(file_format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not file_format:
        raise ValueError("Format required (lower case string)")          
    if file_format != file_format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    with as_handle(handle, 'rU') as fp:
        if file_format in _FormatToReader:
            reader_generator = _FormatToReader[file_format]
            return reader_generator(fp, **params).read()
        else:
            raise ValueError("Unknown format '%s'" % file_format)

def parse(handle, file_format):
    """
    Iterate over a gene ontology file.
    
    Parameters:
     - handle - File handle object to read from, or filename as a string,
     - file_format - lower case string describing the file format to write,
         Formats:
             - obo
             - tsv
             
    You should close the handle after calling this function.
    """

    if not isinstance(file_format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not file_format:
        raise ValueError("Format required (lower case string)")          
    if file_format != file_format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    with as_handle(handle, 'rU') as fp:
        if file_format in _FormatToIterator:
            iterator_generator = _FormatToIterator[file_format]
            it = iterator_generator(fp)

            for el in it:
                yield el
        else:
            raise ValueError("Unknown format '%s'" % file_format)

def pretty_print(enrichment, graph, handle, file_format, **params):
    """
    Print results returned by enrichment finder in a specified format.
    
     Parameters:
     - enrichment - result from EnrichmentFinder
     - graph - OntologyGraph with containing enriched nodes
     - handle - File handle object to read from, or filename as a string,
     - file_format - lower case string describing the file format to write,
         Formats:
             - gml
             - png
             - txt
             - html
     - params - additional parameters
     
    You should close the handle after calling this function.
    """
    
    if not isinstance(file_format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not file_format:
        raise ValueError("Format required (lower case string)")

    with as_handle(handle, 'w') as fp:
        #Map the file format to a writer class
        if file_format in _FormatToPrinter:
            writer_class = _FormatToPrinter[file_format]
            writer = writer_class(fp, **params)
            writer.pretty_print(enrichment, graph)
        else:
            raise ValueError("Unknown format '%s'" % file_format)
