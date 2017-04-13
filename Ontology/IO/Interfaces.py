# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
Contains base classes for IO modules that can be used through Bio.Ontology.IO
module.
"""

class OntoWriter(object):
    """
    Base class for writers.
    """
    
    def __init__(self, file_handle):
        self.handle = file_handle
        
    def write(self, data):
        raise NotImplementedError("write not implmented yet")
    
class OntoReader(object):
    """
    Base class for readers.
    """
    
    def __init__(self, file_handle):
        self.handle = file_handle
        
    def read(self):
        raise NotImplementedError("read not implmented yet")
    
class OntoIterator(object):
    """
    Base class for iterators.
    """
    
    def __init__(self, file_handle):
        self.handle = file_handle
        
    def __iter__(self):
        return self
    
    def __next__(self):
        raise NotImplementedError("next not implmented yet")

    def next(self):
        raise NotImplementedError("next not implmented yet")
