# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:15:35 2016
This file contains functions to preprocess CAFA related files
including: 
-split prediction file by ontology
-split obo file by ontology
-spllit ancestor file by ontology
-split benchmark file by species


@author: Ashley Zhou
"""

import os
#os.chdir('/home/nzhou/git')
from Ontology.IO import OboIO
from precrec.GOPred import GOPred

def prediction_ontology_split_write(pred_path, obo_path):
    """
    Separate the prediction file into the different ontologies
    pred_path should be a handle!!!!!!!!!!
    This should work
    Original code by Dr. Friedberg
    Last edited by Ashley 05/11/2016
    """
    
    all_pred = GOPred()
    #pred_path = open('/home/nzhou/git/CAFAAssess/precrec/M1HS.74.Homo_sapiens.txt')
    all_pred.read(pred_path)
    #obo_path = '/home/nzhou/git/Ontology/go-basic.obo'
    go_graph = OboIO.OboReader(open(obo_path)).read()
    mfo_out = open("%s_MFO.txt" % os.path.splitext(pred_path.name)[0],"w")
    bpo_out = open("%s_BPO.txt" % os.path.splitext(pred_path.name)[0],"w")
    cco_out = open("%s_CCO.txt" % os.path.splitext(pred_path.name)[0],"w")
    for protein in all_pred.data.items():
        for u in protein[1]:
            if go_graph.get_namespace(u['term']) == 'molecular_function': 
                mfo_out.write("%s\t%s\t%.2f\n" % (protein[0], u['term'], u['confidence']))
            elif go_graph.get_namespace(u['term']) == 'biological_process': 
                bpo_out.write("%s\t%s\t%.2f\n" % (protein[0], u['term'], u['confidence']))
            elif go_graph.get_namespace(u['term']) == 'cellular_component': 
                cco_out.write("%s\t%s\t%.2f\n" % (protein[0], u['term'], u['confidence']))
            else:
                raise ValueError ("Term %s not found in any ontology" % u['term'])
    mfo_out.close()
    bpo_out.close()
    cco_out.close()
    
    
def go_ontology_ancestors_split_write(obo_path):
    """
    Input: an OBO file
    Output: 3 files with ancestors
    by Dr. Friedberg
    """
    obo_mfo_out = open("%s_ancestors_mfo.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_bpo_out = open("%s_ancestors_bpo.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_cco_out = open("%s_ancestors_cco.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_parser = OboIO.OboReader(open(obo_path))
    go = obo_parser.read()
    mfo_terms, bpo_terms, cco_terms = go_ontology_split(go)
    for term in mfo_terms:
        ancestors = go.get_ancestors(term)
        obo_mfo_out.write("%s\t%s\n" % (term,",".join(ancestors)))
    for term in bpo_terms:
        ancestors = go.get_ancestors(term)
        obo_bpo_out.write("%s\t%s\n" % (term,",".join(ancestors)))
    for term in cco_terms:
        ancestors = go.get_ancestors(term)
        obo_cco_out.write("%s\t%s\n" % (term,",".join(ancestors)))

    obo_mfo_out.close()
    obo_bpo_out.close()
    obo_cco_out.close()

def go_ontology_split_write(obo_path):
    """
    Split a GO obo file into three files with different namespaces
    by Dr.Friedberg
    """
    obo_mfo_out = open("%s_mfo.obo" % os.path.splitext(obo_path)[0],"w")
    obo_bpo_out = open("%s_bpo.obo" % os.path.splitext(obo_path)[0],"w")
    obo_cco_out = open("%s_cco.obo" % os.path.splitext(obo_path)[0],"w")
    obo_parser = OboIO.OboReader(open(obo_path))
    go = obo_parser.read()
    mfo_terms, bpo_terms, cco_terms = go_ontology_split(go)
    for term in mfo_terms:
        obo_mfo_out.write("%s\n" % term)
    for term in bpo_terms:
        obo_bpo_out.write("%s\n" % term)
    for term in cco_terms:
        obo_cco_out.write("%s\n" % term)

    obo_mfo_out.close()
    obo_bpo_out.close()
    obo_cco_out.close()

def go_ontology_split(ontology):
    """
    Split an GO obo file into three ontologies
    by Dr. Friedberg
    """
    mfo_terms = set({})
    bpo_terms = set({})
    cco_terms = set({})
    for node in ontology.get_ids(): # loop over node IDs and alt_id's
        if ontology.namespace[node] == "molecular_function":
            mfo_terms.add(node)
        elif ontology.namespace[node] == "biological_process":
            bpo_terms.add(node)
        elif ontology.namespace[node] == "cellular_component":
            cco_terms.add(node)
        else:
            raise(ValueError,"%s has no namespace" % node)
    return (mfo_terms, bpo_terms, cco_terms)

def benchmark_species_split(benchmark_path,taxonid):
    '''
    This splits the benchmark file by species(taxon id)
    benchmark_path should be ontology-specific
    by Ashley 05/25/2016

    '''
    out = open('%s_%s.txt' % (os.path.splitext(benchmark_path)[0],str(taxonid)),'w')
    with open(benchmark_path,'r') as bp:
        for line in bp:
            if line[1:5]==str(taxonid):
                out.write('%s' % line)
    out.close()
            
