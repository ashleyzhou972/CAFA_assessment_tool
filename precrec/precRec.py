# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
New CAFA precision recall assessment

Created on Thu Apr 14 17:14:30 2016
dependency: - 
@author: Ashley Zhou
last updated: March 21, 2017
"""

import sys
from collections import defaultdict
import numpy
import os
from Ontology.IO import OboIO


legal_species = [
"all",
"eukarya",
"prokarya",
'HELPY',
 'ECOLI',
 'RAT',
 'DANRE',
 'SULSO',
 'DROME',
 'PSEPK',
 'STRPN',
 'PSEAE',
 'BACSU',
 'MYCGE',
 'HUMAN',
 'METJA',
 'DICDI',
 'YEAST',
 'SCHPO',
 'ARATH',
 'XENLA',
 'MOUSE',
 'PSESM',
 'SALTY',
 'CANAX',
 'SALCH']
legal_types = ["type1","type2","all"]
legal_subtypes = ["easy","hard"]
#easy and hard are only in NK benchmarks!!!!
legal_modes = ["full", "partial"]
root_terms = ['GO:0008150','GO:0005575','GO:0003674']


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
    
    
def go_ontology_ancestors_split_write(obo_path):
    """
    Input: an OBO file
    Output: 3 files with ancestors
    by Dr. Friedberg
    """
    
    obo_bpo_out = open("%s_ancestors_bpo.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_cco_out = open("%s_ancestors_cco.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_mfo_out = open("%s_ancestors_mfo.txt" % (os.path.splitext(obo_path)[0]),"w")
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
    return([len(bpo_terms),len(cco_terms),len(mfo_terms)])
    
    
    
class benchmark:
    def __init__(self,ancestor_path,benchmark_path):
        '''
        Initialize the benchmark.
        
        Input: 
        benchmark_path is ontology specific
        ancestor_path is ontology specific
        '''
        
        #Key: protein
        #Value: set of benchmark leaf terms
        self.ancestors = defaultdict(set)
        # Read GO ancestors file generated with go_ontology_ancestors_split_write()
        # File format: 
        # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
        with open(ancestor_path) as ancestors_input:
            for inline in ancestors_input:
                inrec = inline.strip().split('\t')
                term = inrec[0]
                if len(inrec) == 1:
                    self.ancestors[term] = set({})
                else:
                    term_ancestors = inrec[1]
                    self.ancestors[term] = set(term_ancestors.split(','))
                
        self.true_base_terms = defaultdict(set)
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self.true_base_terms[protein].add(term)
                
    def propagate(self):
        '''
        Progate Benchmark terms.
        '''
        
        #Key: protein
        #Value: set of benchmark propagated terms
        self.true_terms = defaultdict(set)
        for protein in self.true_base_terms:
            for term in self.true_base_terms[protein]:
                try:
                    
                    ancestors = self.ancestors[term].difference(root_terms)    
                #delete root term in self.true_terms
                #modified on 20170203
                except KeyError:
                    sys.stderr.write("%s not found \n" % term) 
                self.true_terms[protein].add(term)
                self.true_terms[protein] |= ancestors
        

def read_benchmark(namespace, species, types, fullbenchmarkfolder, obopath):
    '''
    Read Benchmark.
    
    Input:
    namespace 
    species
    types
    fullbenchmarkfolder
    obopath
    
    Output:
    bench
    obocountDict
    '''
    #Ancestor files here are precomputed
    #To get the ancestor files, use preprocess.py (go_ontology_ancestors_split_write) DOES NOT EXIST
    legal_types = ["type1","type2","typex"]
    legal_subtypes = ["easy","hard"]
    legal_namespace = ["bpo","mfo","cco","hpo"] 
    #fullbenchmarkfolder = './precrec/benchmark/'
    if namespace not in legal_namespace:
        sys.stderr.write("Namespace not accepted, choose from 'bpo', 'cco', 'mfo' and 'hpo'\n")
    elif (species not in legal_species) and (species not in legal_subtypes):
        sys.stderr.write('Species not accepted')
    elif types not in legal_types:
        sys.stderr.write('Type not accepted, choose from "type1","type2" and "typex"\n')
    else:
        matchname = namespace+'_'+species+'_'+types+'.txt'
    #generate ancestor files
    obocounts = go_ontology_ancestors_split_write(obopath)
    obocountDict={'bpo':obocounts[0],'cco':obocounts[1],'mfo':obocounts[2]}
    #ontology-specific calculations
    if namespace == 'bpo':
        full_benchmark_path = fullbenchmarkfolder+'/groundtruth/'+'leafonly_BPO.txt'
        ancestor_path = os.path.splitext(obopath)[0]+"_ancestors_bpo.txt"
    elif namespace =='cco':
        full_benchmark_path = fullbenchmarkfolder+'/groundtruth/'+'leafonly_CCO.txt'
        ancestor_path = os.path.splitext(obopath)[0]+"_ancestors_cco.txt"
    elif namespace == 'mfo':
        full_benchmark_path = fullbenchmarkfolder+'/groundtruth/'+'leafonly_MFO.txt'
        ancestor_path = os.path.splitext(obopath)[0]+"_ancestors_mfo.txt"
    benchmarkListPath = fullbenchmarkfolder+'/lists/'+matchname
    if os.path.isfile(benchmarkListPath) and os.path.getsize(benchmarkListPath)>0:
        handle = open(fullbenchmarkfolder+'/lists/'+matchname, 'r')
        prots  = set()
        for line in handle:
            prots.add(line.strip())
        handle.close()
        tempfilename = 'temp_%s_%s_%s.txt' % (namespace, species,types)
        tempfile = open(fullbenchmarkfolder+'/'+tempfilename ,'w')
        for line in open(full_benchmark_path,'r'):
            prot = line.split('\t')[0]
            if prot in prots:
                tempfile.write(line)
        tempfile.close()
        bench = benchmark(ancestor_path,tempfile.name)
        bench.propagate()
        os.remove(tempfile.name)
    else:
        print('Benchmark set is empty.\n')
        bench=None
    return(bench, obocountDict)

class result:    
#WHY DOES THIS EXIST -> THIS IS GOOD ONCE WE ADD OTHER TYPES BUT OVERALL NEEDS REWORK TO APPROACH THIS STYLE

    def __init__(self):
        ''' Initalize the result '''
        
        #type include: precision/recall, weighted pr, ru/mi, weighted ru/mi
        self.type = ''
        self.precision = []
        self.recall = []
        self.opt = int
        #opt is the optimized value including fmax, smin, wfmax, etc
        self.thres = float
        #thres is the threshold value that gives the optimized value
        self.author = ''
        self.model = ''
        self.keywords = ''
        self.taxon = ''
        self.ontology = ''
        self.mode = ''
        self.TYPE = ''
        self.coverage = float
        #NK or LK
        
        
    def read_from_GOPred(self,GOPred):
        ''' Read in GOPred object '''
        
        self.author = GOPred.author
        self.model = GOPred.model
        self.taxon = GOPred.taxon
        self.keywords = GOPred.keywords
        

class PrecREC: #Can we rename to PRRC ? Just a thought
    '''
    New code by Ashley
    updated: 12/01/2016
    A class for doing precision recall calculations
    '''
    def __init__(self, benchmark, os_pred_path,obocounts):
        '''
        Initalize.
        
        Input:
        benchmark -- instance of the benchmark class
        os_pred_path -- !ontology-specific! prediction file separated in GOPred (without headers)
        obocounts -- total number of terms in the ontology
        '''
        
        #Initialize variables
        self.exist                            = True
        self.ancestors                        = benchmark.ancestors
        self.true_terms                       = benchmark.true_terms
        self.obocount                         = obocounts
        #flag corresponding if the program has been run, allowing print functions
        self.ran                              = False
        #Set of all obsolete terms found
        self.obsolete                         = set()
        #count_above_threshold is the number of proteins with at least one term above threshold
        self.count_above_threshold            = defaultdict()
        #count_predictions_in_benchmark is the number of predicted proteins in this file that are in the benchmark file (for coverage)
        self.count_predictions_in_benchmark   = 0
        self.count_true_terms                 = len(self.true_terms)
        #self.data is the dictionary for all predicted terms
        self.data                             = defaultdict(list)
        #self.predicted_bench is the dictionary for all predicted terms that are benchmarks
        self.predicted_bench                  = defaultdict(defaultdict)
        #key:protein
        #value: list of dictionaries
        #key: GO term
        #value: tuple(confidence, Boolean) Boolean is whether in self.true_terms     
        
        #Read in prediction file        
        if os.path.getsize(os_pred_path) > 0:
            for inrec in open(os_pred_path,'r'):
                fields = [i.strip() for i in inrec.split()]
                self.data[fields[0]].append({'term': fields[1], 'confidence': float(fields[2])})
            
            #Propagated prediction
            for protein in self.data:
                if self.true_terms[protein]:
                    '''
                    benchmark.true_terms[protein] not an empty set
                    The protein is in the benchmark file
                    i.e. gained experimental annota
                    '''
                    self.count_predictions_in_benchmark += 1
                    for tc in self.data[protein]:
                        try:
                            ancterms = self.ancestors[tc['term']].difference(root_terms)
                        except KeyError:
                            # Add unknown term to the obsolete set
                            self.obsolete.add(tc['term'])
                            continue
                        #For each term
                        if tc['term'] in self.predicted_bench[protein]:
                            # Term already exists, update confidence
                            self.update_confidence(protein,tc)
                        else:
                            # Add term to self.predicted_bench
                            # Add confidence and compare with self.true_terms
                            #Regardless of comparision, propagate
                            if tc['term'] not in root_terms:
                                self.predicted_bench[protein][tc['term']]=self.compare(protein, tc)
                                for ancterm in ancterms:
                                    #Make a new TC with ancestors
                                    newtc = {'term':ancterm,'confidence':tc['confidence']}
                                    if ancterm in self.predicted_bench[protein]:
                                        # Term already exists, update confidence
                                        self.update_confidence(protein, newtc)
                                    else:
                                        # Add term to self.predicted_bench
                                        self.predicted_bench[protein][ancterm]=self.compare(protein, newtc)
                                        
            if self.count_predictions_in_benchmark==0:
                self.exist = False
                print("No protein in this predicted set became a benchmark\n")
        else:
            self.exist=False
            print('No prediction made in this ontology.\n')
            
            
    def coverage(self):
        '''
        Determine the coverage.
    
        Finds the value of the coverage-> BETTER EXPLANATION?
        
        Output:
        The coverage as a float
        '''
        
        return float(self.count_predictions_in_benchmark)/self.count_true_terms


    def update_confidence(self, protein, tc):
        '''
        Update Confidence for given protein and propagate.
        
        This function compares the confidence value in tc to the confidence in self.predicted
        If tc is larger, than it overwrites the confidence in self.predicted
        
        Input:
        protein -- chosen protein 
        tc -- a dictionary format as: {'confidence':0.57,'term':'GO:006644'}
        
        Output:
        None
        '''
        
        #Defined for readablity
        confidence = tc['confidence'] 
        
        if confidence > self.predicted_bench[protein][tc['term']][0]:
            #Update the confidence
            self.predicted_bench[protein][tc['term']][0] = confidence
            #Propagate changes if necessary
            for ancterm in self.ancestors[tc['term']].difference(root_terms):
                if confidence > self.predicted_bench[protein][ancterm][0]:
                    #Update the confidence
                    self.predicted_bench[protein][ancterm][0] = confidence
                    
                    
    def compare(self, protein, tc):
        '''
        Check if tc['term'] is a True term.
        
        This function compares if tc['term'] is in self.true_terms
        
        Input:
        protein -- chosen protein to for this comparision
        tc -- dictionary with {'confidence':0.57,'term':'GO:006644'}
        
        Output:
        A list containing:
        'confidence'
        Boolean Value
        '''
        
        if tc['term'] in self.true_terms[protein]:
            return [tc['confidence'],True]
        else:
            return [tc['confidence'],False]
            
        
    def getObsolete(self):
        ''' Get obsolete terms used by the prediction team. '''
        
        return(self.obsolete)


    def term_precision_recall(self, threshold, protein):
        '''
           Calculate the precision recall of a single protein.
           
           Input:
           threshold -- the value for PRRC threshold
           protein -- the chosen protein
        '''        
        
        #Initalize Variables
        TP = 0      # True positive
        count = 0   # Count how many terms are above the threshold
        
        if threshold == 0:
            # At threshold 0, every term (propagated) in the ontology is considered predicted
            # so TP is all the true terms in the benchmark
            # recall = 1
            # precision is TP over all terms in the ontology
            TP = float(len(self.true_terms[protein]))
            count = self.obocount        
        else:
            for term in self.predicted_bench[protein]:
                # If it is above the threshold, increment the count
                if self.predicted_bench[protein][term][0] >= threshold:
                    count += 1
                    # If it is actually True, increment TP
                    if self.predicted_bench[protein][term][1] :
                        TP += 1
        # Find PR: TP / (TP + FP)
        try:
            precision = TP / count 
        except ZeroDivisionError:
            precision = None
        # Find RC: TP (TP + FN)
        recall = TP/len(self.true_terms[protein])
        # Safe becuase protein is in the benchmark file
        return (precision, recall, count , TP)

                        
    def precision_recall(self, threshold, mode):
        '''
        Calculate the overall (average) PRRC of the team.
        
        Input:
        threshold -- The threshold for PRRC
        mode -- Which mode to run in, options below:
        'full'   :  recall is averaged over all benchmark proteins, unpredicted have 0 recall
        'partial':  recall is averaged over all benchmark proteins that are predicted
        
        Output:
        precision --
        recall -- 
        '''
        
        # Initialize Variables
        PR = 0.0
        RC = 0.0
        self.count_above_threshold[threshold] = 0

        for protein in self.predicted_bench:
            pr,rc = self.term_precision_recall(threshold,protein)[0:2]
            if pr is not None:
                PR += pr
                self.count_above_threshold[threshold] += 1
            if rc is not None:
                RC += rc   
                
        if mode=='partial':
            try:
                recall = RC/self.count_predictions_in_benchmark
            except ZeroDivisionError:
                recall = 0
                print("No protein in this predicted set became benchmarks\n")
                
        elif mode=='full':
            try:
                recall = RC/self.count_true_terms
            except ZeroDivisionError:
                recall = 0
                print("No protein in this benchmark set\n")
                
        try:
            precision = PR/self.count_above_threshold[threshold]   
        except ZeroDivisionError:
            precision = None
            print("No prediction is made above the %.2f threshold\n" % threshold)
           
        #PRRC has run
        self.ran = True
        return (precision, recall)
    

    def Fmax_output(self, mode):
        '''
        Compute PRRC for every threshold and Fmax value.
        
        Input:
        mode -- The evaluation mode used ('full' or 'partial')
        
        Output:
        A list containing:
        PR              -- A list containing the precison values
        RC              -- A list containing the recall values
        fmax            -- the maximun f value found
        fmax_threshold  -- The threshold that corresponds to fmax    
        '''
        
        # Intialize Variables
        fmax = 0.0
        fmax_threshold = 0.0
        PR = []
        RC = []
        
        # Run over all threshold values from 0 to 1, two signifigant digits
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            # Run PRRC on given threshold
            pr,rc = self.precision_recall(threshold, mode)
            if pr is None:
                # No prediction above this threshold 
                break
            else:
                PR.append(pr)
                RC.append(rc)
                # Find the F-value for this particular threshold
                try:
                    f = (2*pr*rc)/(pr+rc)
                except ZeroDivisionError:
                    f = None
                    
            if f is not None and f >= fmax: ###########QUESTION##############
                fmax = f
                fmax_threshold = threshold
        #Have found the Fmax at this point       
        return ([PR, RC, fmax, fmax_threshold])
        
        
    def printNumProteins(self, threshold):
        '''
        Print to console the counts if PRRC has been run.
        
        Input:
        threshold --The threshold value for count_above_threshold
        '''    
        
        if self.ran is True:
            print('number of benchmark proteins: %s\n'% self.count_true_terms)
            #Those with not-None recall: 
            #(predicted protein that are in benchmark, only one species per prediction file!! )
            print ('number of predicted proteins that are in the benchmark file: %s\n' % self.count_predictions_in_benchmark)
            #Those with not-None precision:
            print('number of proteins with at least one term above threshold: %s\n' % self.count_above_threshold[threshold] )
            
        else:
            print ("Run precision_recall(%s) first\n" % str(threshold))
            
        
    def printConfidence(self, output_path):
        '''
        Print confidence and True/False to a file.
        
        Input:
        output_path -- Where to save the file
        '''
        
        protindex = 0
        out = open(output_path, 'w')
        for prot in self.predicted_bench:
            protindex += 1
            for term in self.predicted_bench[prot]:
                out.write("%s\t%s\t%s\n" % (str(protindex), self.predicted_bench[prot][term][0], self.predicted_bench[prot][term][1]))
        out.close()
        


