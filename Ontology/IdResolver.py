# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
Module containing resolvers for ambigious entities ids.
"""

import collections

class Resolver(object):
    """
    Resolver which simply returns key given for disambiguation.
    """
    
    def __init__(self, annotations):
        pass
    
    def resolve(self, oid):
        return oid

class SetPickerResolver(Resolver):
    """
    Resolver which picks the key which belongs to association from given
    mapping of the keys to synonyms.
    """
    def __init__(self, synonyms, annotations):
        self.base_keys = set()
        for obj in annotations:
            self.base_keys.add(obj.id)
        self.synonyms = synonyms
        
    def resolve(self, oid):
        if oid in self.base_keys or oid not in self.synonyms:
            return oid
        else:
            for x in self.synonyms[oid]:
                if x in self.base_keys:
                    return x
            return oid
        
class FirstOneResolver(Resolver):
    """
    Resolver which picks first possible key for ambiguous entry.
    """
    
    def __init__(self, annotations):
        self.base_keys = set()
        alter = collections.defaultdict(set)
        for obj in annotations:
            self.base_keys.add(obj.id)
            if 'Synonym' in obj.attrs:
                for aid in obj.attrs['Synonym']:
                    alter[aid].add(obj.id)
        self.alter_keys = dict([(k,list(v)) for (k, v) in alter.items()])
        
    
    def resolve(self, oid):
        if oid in self.base_keys or oid not in self.alter_keys:
            return oid
        else:
            return self.alter_keys[oid][0]
