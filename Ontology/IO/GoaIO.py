# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
I/O operations for gene annotation files.
"""
from __future__ import print_function

import sys
import csv
import collections
from Ontology.Data import GeneAnnotation, TermAssociation
from .Interfaces import OntoIterator, OntoReader

class TsvIterator(OntoIterator):
    """
    Parses TSV files
    """

    def __init__(self, file_handle):
        self._reader = csv.reader(file_handle, delimiter='\t')

    def __iter__(self):
        return self._reader

    def __next__(self):
        return next(self._reader)
    
    def next(self):
        return next(self._reader)

# GAF version 2.0
GAF20FIELDS = ['DB',
        'DB_Object_ID',
        'DB_Object_Symbol',
        'Qualifier',
        'GO_ID',
        'DB:Reference',
        'Evidence',
        'With',
        'Aspect',
        'DB_Object_Name',
        'Synonym',
        'DB_Object_Type',
        'Taxon_ID',
        'Date',
        'Assigned_By',
        'Annotation_Extension',
        'Gene_Product_Form_ID']

# GAF version 1.0
GAF10FIELDS = ['DB',
        'DB_Object_ID',
        'DB_Object_Symbol',
        'Qualifier',
        'GO_ID',
        'DB:Reference',
        'Evidence',
        'With',
        'Aspect',
        'DB_Object_Name',
        'Synonym',
        'DB_Object_Type',
        'Taxon_ID',
        'Date',
        'Assigned_By'] 

GAF_VERSION = { "1.0" : GAF10FIELDS,
                "2.0" : GAF20FIELDS}

def _split_multi(value):
    if len(value) > 0:
        return value.split('|')
    else:
        return []

def _to_goa(obj_rows, version):
    row = obj_rows[0]
    
    obj_id = row[1]
    obj_attrs = {GAF20FIELDS[0] : row[0],
                 GAF20FIELDS[2] : row[2],
                 GAF20FIELDS[9] : row[9],
                 GAF20FIELDS[10] : _split_multi(row[10]),
                 GAF20FIELDS[11] : row[11],
                 GAF20FIELDS[12]: _split_multi(row[12])}
    
    if version == "1.0":
        row_len = 15
    else:
        row_len = 17
        obj_attrs[GAF20FIELDS[15]] = _split_multi(row[15])
        obj_attrs[GAF20FIELDS[16]] = row[16]
        
    assocs = []
    for row in obj_rows:
        if len(row) == row_len:
            assocs.append(TermAssociation(row[4],
                                       {GAF20FIELDS[3] : _split_multi(row[3]),
                                        GAF20FIELDS[5] : _split_multi(row[5]),
                                        GAF20FIELDS[6] : row[6],
                                        GAF20FIELDS[7] :_split_multi(row[7]),
                                        GAF20FIELDS[8] : row[8],
                                        GAF20FIELDS[13] : row[13],
                                        GAF20FIELDS[14] : row[14]}
                                          ))
        else:
            raise ValueError("Invalid gaf file: Incorrect row length.")
    
    return GeneAnnotation(obj_id, assocs, obj_attrs)

class GafReader(OntoReader):
    """
    Reads GAF files into list of GeneAnnotation.
    
    GAF file is list of tab separated values in the following order:
        'DB', 'DB Object ID', 'DB Object Symbol', 'Qualifier', 'GO ID',
        'DB:Reference', 'Evidence Code', 'With (or) From', 'Aspect',
        'DB Object Name', 'DB Object Synonym', 'DB Object Type',
        'Taxon', 'Date', 'Assigned By', 'Annotation Extension',
        'Gene Product Form ID'
    """

    _ID_IDX = 1
    
    def __init__(self, file_handle, assoc_format = "dict"):
        """
        Parameters:
        ----------
        - assoc_format - states format of returned association:
          o "dict" - as a dictionary (faster)
          o "in_mem_sql" - as dict-like object with underlying in-memory database
                         (more memory efficient)
         
        """
        self.handle = file_handle
        self.assoc_format = assoc_format

    
    def read(self):
        first = self.handle.readline()
        if first and first.startswith('!gaf-version:'):
            version = first[(first.find(':') + 1):].strip()
        else:
            raise ValueError("Invalid gaf file: No version specified.")
        if version not in GAF_VERSION:
            raise ValueError("Incorrect version.")
        
        tsv_iter = TsvIterator(self.handle)
        if self.assoc_format == "dict":
            raw_records = collections.defaultdict(list)
            for row in tsv_iter:
                first = row[0]
                if not first.startswith('!'):
                    raw_records[row[self._ID_IDX]].append(row)
            return dict([(k, _to_goa(v, version)) for k, v in raw_records.items()]) # Possible py2 slow down
        elif self.assoc_format == "in_mem_sql":
            try:
                sqla = InSqlAssoc(GAF_VERSION[version], [1,4], lambda x:  _to_goa(x, version))
            except ImportError:
                print("Error: To use in_mem_sql association you need to have sqlite3 bindings installed.", file=sys.stderr)
            else:
                for row in tsv_iter:
                    if not row[0].startswith('!'):
                        sqla.add_row(row)
                return sqla
        else:
            raise ValueError("Incorrect assoc_format parameter.")
            
    
class InSqlAssoc(object):
    """
    Immutable dictionary-like structure for storing annotations.
    
    It provides slower access but is more memory efficient thus more suitable
    for big annotations files.
    """
    
    def __init__(self, fields, index, selection_to_obj_fun, db_path = ":memory:"):
        """
        Parameters:
        ----------
        - fields - name of the columns in db representation
        - index - pair of fields indexing associations: (gene_id, ontology_term_id)
        - selection_to_obj_fun - function transforming list of rows into
          GeneAssociation
        - db_path - path to database file, special value ":memory:" creates
          database in memory
         
        """
        import sqlite3
        
        self.fields = fields
        self.fun = selection_to_obj_fun
        self.con = sqlite3.connect(db_path)
        self.index = index
        
        cur = self.con.cursor()
        
        query = 'CREATE TABLE assocs ("' + self.fields[0] + '" '
        for field in fields[1:]:
            query += ', "' + field + '" '
        query += ');'
        cur.execute(query)
        cur.execute('CREATE INDEX obj_idx ON assocs ({0});'.format(self.fields[index[0]]))
        self.con.commit()
        
    def add_row(self, row):
        if len(row) != len(self.fields):
            raise TypeError("Incorrect number of fields in a row.")
        else:
            cur = self.con.cursor()
            cur.execute("INSERT INTO assocs VALUES (?" + (",?" * (len(self.fields) - 1)) + ");", row)
            self.con.commit()
          
    def __len__(self):
        cur = self.con.cursor()
        cur.execute('SELECT COUNT(DISTINCT "' + self.fields[self.index[0]] + '") FROM assocs;')
        return cur.fetchone()[0]
    
    def __contains__(self, key):
        cur = self.con.cursor()
        cur.execute('SELECT * FROM assocs WHERE "' + self.fields[self.index[0]]\
                     + '" = ?;', [key])
        return len(list(cur)) > 0 #TODO sth prettier
    
    def __getitem__(self, key):
        cur = self.con.cursor()
        cur.execute('SELECT * FROM assocs WHERE "' + self.fields[self.index[0]]\
                     + '" = ?;', [key])
        return self.fun(list(cur))
    
    def __iter__(self):
        cur = self.con.cursor()
        cur.execute('SELECT * FROM assocs ORDER BY "{0}"'.format(self.fields[self.index[0]]))
        cur_id = None
        row_list = []
        for row in cur:
            if cur_id and cur_id != row[self.index[0]]:
                obj = self.fun(row_list)
                row_list = [row]
                cur_id = row[self.index[0]]
                yield (cur_id, obj)
            else:
                cur_id = row[self.index[0]]
                row_list.append(row)
        yield (cur_id, self.fun(row_list))

    def itervalues(self):
        for _, v in self:
            yield v
            
    def iterkeys(self):
        for k, _ in self:
            yield k
            
    def keys(self):
        return self.keys()
    
    def values(self):
        return self.values()
