# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:24:16 2016

@author: nzhou
"""

#!/usr/bin/env python3

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

import re
import sys
import os
import gzip
from collections import defaultdict
from Ontology.IO import OboIO
import matplotlib.pyplot as plt
import numpy as np


pr_field = re.compile("^PR=[0,1]\.[0-9][0-9];$")
rc_field = re.compile("^RC=[0,1]\.[0-9][0-9]$")
# Fix to add HP ontology 2014-1-9
go_field = re.compile("^(GO|HP):[0-9]{5,7}$")
# Fix to add EFI targets 2014-1-9
target_field = re.compile("^(T|EFI)[0-9]{5,20}$")
confidence_field = re.compile("^[0,1]\.[0-9][0-9]$")

# Legal states: the CAFA prediction records fields, and their order. KEYWORDS and ACCURACY are
# optional
legal_states1 = ["author","model","keywords","accuracy","go_prediction","end"]
legal_states2 = ["author","model","keywords","go_prediction","end"]
legal_states3 = ["author","model","go_prediction","end"]
legal_states4 = ["go_prediction"]

legal_keywords = [
"sequence alignment", "sequence-profile alignment", "profile-profile alignment", "phylogeny",
"sequence properties",
"physicochemical properties", "predicted properties", "protein interactions", "gene expression",
"mass spectrometry",
"genetic interactions", "protein structure", "literature", "genomic context", "synteny", 
"structure alignment",
"comparative model", "predicted protein structure", "de novo prediction", "machine learning", 
"genome environment", 
"operon", "ortholog", "paralog", "homolog", "hidden markov model", "clinical data", "genetic data", 
"natural language processing", "other functional information"
]
class GOPred:
    """
    A class for reading and storing CAFA GO predictions
    self.data is a dictionary
       key: protein ID
       value: [{'term':go_term_1, 'threshold': threshold_1},...,{'term':go_term_n, 'threshold': threshold_n}]
    in self.read(pred_path) , pred_path should be a handle
    Edited by Ashley: 05/26/2016
    """
    def __init__(self):
        self.author = None
        self.model = None
        self.keywords = []
        self.data = defaultdict(list)

    def _author_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 2:
            correct = False
            errmsg = "AUTHOR: invalid number of fields. Should be 2"
        elif fields[0] != "AUTHOR":
            correct = False
            errmsg = "AUTHOR: First field should be AUTHOR"
        else:
            self.author = fields[1]
        return correct, errmsg

    def _model_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 2:
            correct = False
            errmsg = "MODEL: invalid number of fields. Should be 2"
        elif fields[0] != "MODEL":
            correct = False
            errmsg = "MODEL: First field should be MODEL"
        elif len(fields[1]) != 1 or not fields[1].isdigit():
            correct = False
            errmsg = "MODEL: second field should be single digit."
        else:
            self.model = int(fields[1])
        return correct, errmsg

    def _keywords_check(self,inrec):
        correct = True
        errmsg = None
        if inrec[:8] != "KEYWORDS":
            correct = False
            errmsg = "KEYWORDS: first field should be KEYWORDS"
        else:
            keywords = [i.strip().lower() for i in inrec[8:].split(",")]
            for keyword in keywords:
                #stupid full stop 
                if keyword[-1] == ".":
                    keyword = keyword[:-1]
                if keyword not in legal_keywords:
                    correct = False
                    errmsg = "KEYWORDS: illegal keyword %s" % keyword
                    break
                else:
                    self.keywords.append(keyword)
        return correct, errmsg

    def _accuracy_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 4:
            correct = False
            errmsg = "ACCURACY: error in number of fields. Should be 4"
        elif fields[0] != "ACCURACY":
            correct = False
            errmsg = "ACCURACY: first field should be 'ACCURACY'"
        elif not fields[1].isdigit() or len(fields[1]) != 1:
            correct = False
            errmsg = "ACCURACY: second field should be a single digit"
        elif not pr_field.match(fields[2]):
            correct = False
            errmsg = "ACCURACY: error in PR field"
        elif not rc_field.match(fields[3]):
            correct = False
            errmsg = "ACCURACY: error in RC field"
        return correct, errmsg


    def _go_prediction_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 3:
            correct = False
            errmsg = "GO prediction: wrong number of fields. Should be 3"
        elif not target_field.match(fields[0]):
            correct = False
            errmsg = "GO prediction: error in first (Target ID) field"
        elif not go_field.match(fields[1]):
            correct = False
            errmsg = "GO prediction: error in second (GO ID) field"
        elif not confidence_field.match(fields[2]):
            correct = False
            errmsg = "GO prediction: error in third (confidence) field"
        elif float(fields[2]) > 1.0:
            correct = False
            errmsg = "GO prediction: error in third (confidence) field. Cannot be > 1.0"
        else:
            self.data[fields[0]].append({'term': fields[1], 'confidence': float(fields[2])})
        return correct, errmsg

    def _end_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 1:
            correct = False
            errmsg = "END: wrong number of fields. Should be 1"
        elif fields[0] != "END":
            correct = False
            errmsg = "END: record should include the word END only"
        return correct, errmsg

    def _handle_error(self,correct, errmsg, inrec):
        if not correct:
            print (inrec)
            raise ValueError(errmsg)


    def read(self, pred_path):
        visited_states = []
        n_accuracy = 0
        first_prediction = True
        first_accuracy = True
        first_keywords = True
        n_models = 0
        #inline = open(pred_path)
        for inline in pred_path: 
            # gzipped files are in bytes. Need to convert to utf-8
            if type(inline) is bytes:
                inline = inline.decode("utf-8")
            inrec = [i.strip() for i in inline.split()]
            field1 = inrec[0]
            # Check which field type (state) we are in
            if field1 == "AUTHOR":
                state = "author"
            elif field1 == "MODEL":
                state = "model"
            elif field1 == "KEYWORDS":
                state = "keywords"
            elif field1 == "ACCURACY":
                state = "accuracy"
            elif field1 == "END":
                state = "end"
            else: #default to prediction state
                state = "go_prediction"
            # Check for errors according to state
            if state == "author":
                correct,errmsg = self._author_check(inline)
                self._handle_error(correct, errmsg,inline)
                visited_states.append(state)
            elif state == "model":
                n_models += 1
                n_accuracy = 0
                if n_models > 3:
                    raise ValueError("Too many models. Only up to 3 allowed")
                correct,errmsg = self._model_check(inline)
                self._handle_error(correct, errmsg,inline)
                if n_models == 1:
                    visited_states.append(state)
            elif state == "keywords":
                if first_keywords:
                    visited_states.append(state)
                    first_keywords = False
                correct, errmsg = self._keywords_check(inline)
                self._handle_error(correct, errmsg,inline)
            elif state == "accuracy":
                if first_accuracy:
                    visited_states.append(state)
                    first_accuracy = False
                n_accuracy += 1
                if n_accuracy > 3:
                    self._handle_error(False, "ACCURACY: too many ACCURACY records")
                else:
                    correct, errmsg = self._accuracy_check(inline)
            elif state == "go_prediction":
                correct, errmsg = self._go_prediction_check(inline)
                self._handle_error(correct, errmsg,inline)
                if first_prediction:
                    visited_states.append(state)
                    first_prediction = False
            elif state == "end":
                correct, errmsg = self._end_check(inline)
                self._handle_error(correct, errmsg,inline)
                visited_states.append(state)
            # End file forloop
        if (visited_states != legal_states1 and
            visited_states != legal_states2 and
            visited_states != legal_states3 and
            visited_states != legal_states4):
            print (visited_states)
            print ("file not formatted according to CAFA specs")
            print ("Check whether all these record types are in your file")
            print ("Check whether all these record types are in your file in the correct order")
            print ("AUTHOR, MODEL, KEYWORDS, ACCURACY (optional), predictions, END")
            raise ValueError

class PrecRec:
    """
    Precision-recall calculations for CAFA
    """
    def __init__(self):
        self.ancestors = {}

        # Key: protein; Value: set of true terms from CAFA benchmark
        self.true_base_terms = defaultdict(set)

        # Key: protein; value: set of terms from CAFA benchmark, and their ancestors
        self.true_terms = defaultdict(set)

        # Key: protein; Value: set of predicted terms 
        self.predicted_base_terms = defaultdict(set)

        # Key: protein; Value: set of predicted terms and their ancestors 
        self.predicted_terms = defaultdict(set)

    def read_ancestors(self,ancestors_path):
        # Call this method first, populates self.ancestors
        # Read GO ancestors file generated with obo2ancestors
        # File format: 
        # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
        with open(ancestors_path) as ancestors_input:
            for inline in ancestors_input:
                inrec = inline.strip().split('\t')
                term = inrec[0]
                if len(inrec) == 1:
                    self.ancestors[term] = set({})
                else:
                    term_ancestors = inrec[1]
                    self.ancestors[term] = set(term_ancestors.split(','))

    def read_benchmark(self, benchmark_path):
        # Call this second. 
        # Read the benchmark file that contains the true annotations
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self._true_base_terms[protein].add(term)

    def propagate_true_terms(self):
        # Call this method third. After this method is run, 
        # we have all the true terms for this protein
        for protein in self.true_base_terms:
            for term in self.true_base_terms[protein]:
                try:
                    ancestors = self.get_ancestors(term)
                except KeyError:
                    sys.stderr.write("not found %s\n" % term) 
                self.true_terms[protein].add(term) # Add the base term
                self.true_terms[protein] |= ancestors # Add all ancestors

    def get_predicted_terms(self, gopred):
        for i in gopred.data:
            self.predicted_base_terms.add(i[0])

    def propagate_predicted_terms(self):
        for protein in self.predicted_base_terms:
            for term in self.predicted_base_terms[protein]:
                try:
                    ancestors = self.get_ancestors(term)
                except KeyError:
                    sys.stderr.write("not found %s\n" % term) 
                self.predicted_terms[protein].add(term) # Add the base term
                self.predicted_terms[protein] |= ancestors # Add all ancestors

    def get_ancestors(self,term):
        return self.ancestors[term]
    
    def term_precision_recall(self, protein, term):
        # Precision-recall for a single term for a single protein
        # term: single GO term to be propagated
        # true_terms: set of true terms

        # pred_terms: a set of all ancestors of this term, including 
        # the term itself.
        if term in self.ancestors:
            pred_terms = self.get_ancestors(term)
        else:
            pred_terms = set([term])
            sys.stderr.write("No ancestors for %s\n" % term)
        pred_terms.add(term)
        prot_true_terms = self.true_terms[protein]
        if not prot_true_terms:
            precision, recall = None, None
        else:
            # TP / (TP+FP)
            precision = len(pred_terms & prot_true_terms) / len(pred_terms)
            recall = len(pred_terms & prot_true_terms) / len(prot_true_terms)
            if precision > 1 or recall > 1:
                sys.stderr.write("Precision=%.2f\n" % precision)
                sys.stderr.write("Recall=%.2f\n" % recall)
                sys.stderr.write("Term=%s\n" % term)
                sys.stderr.write("Protein=%s\n" % protein)
                sys.stderr.write("Predicted terms=%s\n" % str(pred_terms))
                sys.stderr.write("True terms=%s\n" % str(prot_true_terms))
                raise ValueError

            # TP / (TP+FN)
        if prot_true_terms and precision > 1 and recall >1:
            pass
            sys.stderr.write("*********\n")
            sys.stderr.write("%s\t%s\t%s\n" % (protein, precision, recall))
            sys.stderr.flush()
        return (precision, recall)

def f_1(precision, recall):
    return 2 * (precision * recall) / (precision +  recall)

def f_max(prec_rec_vector):
    f1_max = 0.
    for recall, precision in prec_rec_vector:
        cur_f1  = f_1(precision, recall)
        if cur_f1 > f1_max:
            f1_max = cur_f1
    return f1_max


def get_predicted_terms(prediction_path, benchmark_path, ancestors_path, go_path, namespace):
    
    prediction = GOPred()
    if os.path.splitext(prediction_path)[1] == ".gz":
        inpred = gzip.open(prediction_path)
    else:
        inpred = open(prediction_path)
    prediction.read(inpred)
    inpred.close()

    prec_rec = PrecRec()
    prec_rec.read_ancestors(ancestors_path)
    prec_rec.read_benchmark(benchmark_path)
    prec_rec.propagate_true_terms()

def precision_recall(prediction_path, benchmark_path, ancestors_path, go_path, namespace):
    # accepts a GOPred instantiation
    # does precision / recall calculation
    # TODO: take care of threshold!!
    done_proteins = set({})
    prec = defaultdict(float)
    rec = defaultdict(float)
    mprot = defaultdict(set)
    zprot = defaultdict(int)
    done_predictions = {}
    prec_rec_vector = []

    prec_rec = PrecRec()
    prec_rec.read_ancestors(ancestors_path)
    prec_rec.read_benchmark(benchmark_path)
    prec_rec.propagate_true_terms()


    go_graph = OboIO.OboReader(open(go_path)).read()

    prediction = GOPred()
    if os.path.splitext(prediction_path)[1] == ".gz":
        inpred = gzip.open(prediction_path)
    else:
        inpred = open(prediction_path)
    prediction.read(inpred)
    inpred.close()
    nprot = len(prec_rec.true_terms)
    for threshold in [i*0.01 for i in range(1,101)]:
        threshold_s = "%.2f" % threshold
        # loop over all proteins for which we have a true prediction
        # First get all the predicted terms and their ancestors.
        for protein in prediction.data:
            for u in prediction.data[protein]:
                term = u['term']
                if go_graph.get_namespace(term) != namespace:
                    continue
                confidence = u['confidence']
                if confidence > threshold:
                    precision, recall = prec_rec.term_precision_recall(protein, term)
                    # done_proteins.add(protein)
                    if precision is not None:
                        prec[(threshold_s,protein)] += precision # was +=
                        rec[(threshold_s,protein)] += recall # was +=
                        mprot[threshold_s].add(protein)
                        zprot[threshold_s] += 1
        if threshold_s in mprot:
            prec_sum = 0.
            rec_sum = 0.
            prec_denom = 0.
            for p in mprot[threshold_s]:
                prec_sum += prec[(threshold_s,p)] 
                prec_denom += zprot[threshold_s]
            for p in prec_rec.true_terms:
                rec_sum += rec[(threshold_s,p)]
                
            #precision = prec_sum / zprot[threshold_s]
            precision = prec_sum / prec_denom
            recall = rec_sum / nprot
            # recall on X axis, precision on Y axis
            prec_rec_vector.append((recall, precision)) 
            
            sys.stderr.write("%.2f\t%.2f\t%s\t%d\t%d\n" % (recall, precision, threshold_s, len(mprot[threshold_s]),nprot))
            # sys.stderr.write("%s\t%d\n" % (threshold_s, len(mprot[threshold_s])))
    return prec_rec_vector

def prediction_ontology_split_write(pred_path, obo_path):
    """
    Separate the prediction file into the different ontologies
    """
    
    all_pred = GOPred()
    all_pred.read(pred_path)
    go_graph = OboIO.OboReader(open(obo_path)).read()
    mfo_out = open("%s_MFO.txt" % os.path.splitext(pred_path),"w")
    bpo_out = open("%s_BPO.txt" % os.path.splitext(pred_path),"w")
    cco_out = open("%s_CCO.txt" % os.path.splitext(pred_path),"w")

    for u in all_pred.data:
        if go_graph.get_namespace(go_graph.get_term(u['term'])) == 'molecular_function': 
            mfo_out.write("%s\t%.2f\n" % (u['term'], u['confidence']))
        elif go_graph.get_namespace(go_graph.get_term(u['term'])) == 'biological_process': 
            bpo_out.write("%s\t%.2f\n" % (u['term'], u['confidence']))
        elif go_graph.get_namespace(go_graph.get_term(u['term'])) == 'cellular_component': 
            cco_out.write("%s\t%.2f\n" % (u['term'], u['confidence']))
        else:
            raise ValueError ("Term %s not found in any ontology" % u['term'])
            
    
    mfo_out.close()
    bpo_out.close()
    cco_out.close()

def go_ontology_ancestors_split_write(obo_path):
    """
    Input: an OBO file
    Output: 3 files with ancestors
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

def stub(prediction_path, benchmark_path, ancestors_path):
    # Just a stub to test stuff
    mypred = GOPred()
    mypred.read(prediction_path)
    prec_rec_vector = []
    prec_rec = PrecRec()
    # Read ancestors for a given ontology (mfo, bpo, or cco)
    prec_rec.read_ancestors(ancestors_path)
    # Read benchmark for a given ontology (mfo, bpo, or cco)
    prec_rec.read_benchmark(benchmark_path)
    prec_rec.propagate_true_terms()
    return mypred, prec_rec


if __name__ == '__main__':
    if len(sys.argv) != 6:
        raise ValueError("Usage: ./GOPred.py prediction_path benchmark_path ancestors_path go_path namespace")
    prv = precision_recall(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    print(prv[:5])
    print("Fmax", f_max(prv))
    plt.plot([i[0] for i in prv],[i[1] for i in prv],'ro',ls='-')
    plt.axis([0, 1, 0, 1])
    plt.show()
