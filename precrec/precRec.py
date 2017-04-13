# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
new CAFA precision recall assessment

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




legal_species = ["all",
            "eukarya",
            "prokarya",
            "ARATH",
            "DANRE",
            "DICDI",
            "DROME",
            "ECOLI",
            "HUMAN",
            "MOUSE",
            "PSEAE",
            "RAT",
            "SCHPO",
            "YEAST"]
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
    
    
    
class benchmark:
    def __init__(self,ancestor_path,benchmark_path):
        '''
        Here benchmark_path should be ontology specific
        ancestor_path should also be ontology specific
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
    #Ancestor files here are precomputed
    #To get the ancestor files, use preprocess.py (go_ontology_ancestors_split_write)
    #last updated 12/12/2016
    #TODO: hpo ancestor file?
    #TODO: changed relative path to within CAFAAssess, need to change the whole package!
    legal_species = ["all",
            "eukarya",
            "prokarya",
            "ARATH",
            "DANRE",
            "DICDI",
            "DROME",
            "ECOLI",
            "HUMAN",
            "MOUSE",
            "PSEAE",
            "RAT",
            "SCHPO",
            "YEAST"]
    legal_types = ["type1","type2","typex"]
    legal_subtypes = ["easy","hard"]
    legal_namespace = ["bpo","mfo","cco","hpo"] 
    #fullbenchmarkfolder = './precrec/benchmark/'
    if namespace not in legal_namespace:
        sys.stderr.write("Namespace not accepted, choose from 'bpo', 'cco', 'mfo' and 'hpo'\n")
    elif (species not in legal_species) and (species not in legal_subtypes):
        sys.stderr.write('Species not accepted, choose from "easy", "hard", "eukarya","prokarya","ARATH","DANRE","DICDI","DROME","ECOLI","HUMAN","MOUSE","PSEAE","RAT","SCHPO" and"YEAST"\n')
    elif types not in legal_types:
        sys.stderr.write('Type not accepted, choose from "type1","type2" and "typex"\n')
    else:
        matchname = namespace+'_'+species+'_'+types+'.txt'
    #generate ancestor files
    go_ontology_ancestors_split_write(obopath)
    #ontology-specific calculations
    if namespace == 'bpo':
        full_benchmark_path = fullbenchmarkfolder+'leafonly_BPO.txt'
        ancestor_path = os.path.splitext(obopath)[0]+"_ancestors_bpo.txt"
    elif namespace =='cco':
        full_benchmark_path = fullbenchmarkfolder+'leafonly_CCO.txt'
        ancestor_path = os.path.splitext(obopath)[0]+"_ancestors_cco.txt"
    elif namespace == 'mfo':
        full_benchmark_path = fullbenchmarkfolder+'leafonly_MFO.txt'
        ancestor_path = os.path.splitext(obopath)[0]+"_ancestors_mfo.txt"
    handle = open(fullbenchmarkfolder+'/lists/'+matchname, 'r')
    prots  = set()
    for line in handle:
        prots.add(line.strip())
    handle.close()
    tempfilename = 'temp_%s_%s_%s.txt' % (namespace, species,types)
    tempfile = open(fullbenchmarkfolder+tempfilename ,'w')
    for line in open(full_benchmark_path,'r'):
        prot = line.split('\t')[0]
        if prot in prots:
            tempfile.write(line)
    tempfile.close()
    bench = benchmark(ancestor_path,tempfile.name)
    bench.propagate()
    os.remove(tempfile.name)
    return bench

class result:    
    def __init__(self):
        self.type = ''
        #type include: precision/recall, weighted pr, ru/mi, weighted ru/mi
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
        self.author = GOPred.author
        self.model = GOPred.model
        self.taxon = GOPred.taxon
        self.keywords = GOPred.keywords
        

class PrecREC:
    '''
    New code by Ashley
    updated: 12/01/2016
    A class for doing precision recall calculations
    '''
    def __init__(self, benchmark, os_pred_path):
        '''
        benchmark is an instance of the benchmark class
        os_pred_path should be an !ontology-specific! prediction file separated in GOPred (without headers)
        countb is the number of predicted proteins in this file that are in the benchmark file (for coverage)
        counta is the number of proteins with at least one term above threshold
        '''
        self.exist = True
        self.ancestors = benchmark.ancestors
        self.true_terms = benchmark.true_terms
        self.obsolete = set()
        self.counta = defaultdict()
        self.countb = 0
        self.countc = len(self.true_terms)
        self.data = defaultdict(list)
        self.predicted_bench = defaultdict(defaultdict)
        
        #predicted_base_terms is the same as self.data
        #No need creating another dictionary
        
        #self.data is the dictionary for all predicted terms
        #self.predicted_bench is the dictionary for all predicted terms that are benchmarks
        
        #Now propogate the predicted terms
        #key:protein
        #value: list of dictionaries
        #key: GO term
        #value: tuple(confidence, True/False) whether in true terms or not
        #take the largest confidence
        #Take care of obsolete terms as well
        
        #Read in prediction file        
        if os.path.getsize(os_pred_path) > 0:
            for inrec in open(os_pred_path,'r'):
                fields = [i.strip() for i in inrec.split()]
                self.data[fields[0]].append({'term': fields[1], 'confidence': float(fields[2])})
            
            #Propagated prediction
            for prot in self.data:
                if benchmark.true_terms[prot]:
                    '''
                    benchmark.true_terms[prot] not an empty set
                    The protein is in the benchmark file
                    i.e. gained experimental annotation
                    '''
                    self.countb += 1
                    for tc in self.data[prot]:
                        try:
                            ancterms = self.ancestors[tc['term']].difference(root_terms)
                            #delete root terms
                            #This only delete terms in ancestors, not if the prediction 
                            #itself contains root terms!
                            #modified on 20170203
                        except KeyError:
                            #sys.stderr.write("%s not found\n" % tc['term'])
                            self.obsolete.add(tc['term'])
                            continue
                        if tc['term'] in self.predicted_bench[prot]:
                            #This term has already been added
                            #maybe as an ancestor of other terms
                            #update confidence with the maximum one
                            #propagate and update all ancestor confidence
                            self.__update_confidence__(prot,tc)
                        else:
                            #add this term to self.predicted_bench
                            #add confidence, and compare with self.true_terms
                            #if term in self.true_terms, True
                            #else False
                            #No matter true or false, propagate
                            if tc['term'] not in root_terms:
                                #Delete root terms!
                                #Modified on 20170217
                                self.predicted_bench[prot][tc['term']]=self.__compare__(prot,tc)
                                for ancterm in ancterms:
                                    newtc = {'term':ancterm,'confidence':tc['confidence']}
                                    if ancterm in self.predicted_bench[prot]:
                                        self.__update_confidence__(prot,newtc)
                                    else:
                                        self.predicted_bench[prot][ancterm]=self.__compare__(prot,newtc)
        else:
            self.exist=False
            print('No prediction made in this ontology.\n')

    def __update_confidence__(self,prot,tc):
        '''
        tc is a dictionary with {'confidence':0.57,'term':'GO:006644'}
        prot is a protein
        This function compares the confidence value in tc, and if it's larger than
        the confidence that's been added for term in self.predicted
        we update that confidence
        It also updates all propagated terms of tc
        '''
        if tc['confidence']>self.predicted_bench[prot][tc['term']][0]:
            self.predicted_bench[prot][tc['term']][0]=tc['confidence']
            for ancterm in self.ancestors[tc['term']].difference(root_terms):
                #added difference 20170219
                if tc['confidence']>self.predicted_bench[prot][ancterm][0]:
                    self.predicted_bench[prot][ancterm][0]=tc['confidence']
                    
                    
    def __compare__(self,prot,tc):
        '''
        tc is a dictionary with {'confidence':0.57,'term':'GO:006644'}
        prot is a protein
        This function compares if tc['term'] is in self.true_terms
        returns a list ["confidence","True/False"]
        '''
        if tc['term'] in self.true_terms[prot]:
            return [tc['confidence'],True]
        else:
            return [tc['confidence'],False]
            
        
    def getObsolete(self):
        '''
        return all obsolete terms used by the prediction team
        '''
        return(self.obsolete)

    def term_precision_recall(self,threshold,protein):
        #This function calculates the precision recall of a single protein
        TP = 0.0
        count = 0
        #count is to count how many terms are above the threshold
        for term in self.predicted_bench[protein]:
            if self.predicted_bench[protein][term][0]>=threshold:
                #greater or equal to
                count+=1
                #print(term)
                if self.predicted_bench[protein][term][1] :
                    TP+=1
        try:
            precision = TP/count
            #print(count)
            #print(TP)
        except ZeroDivisionError:
            precision=None
        if threshold==0:
            #At threshold 0, everything is considered predicted, therefore guarantee perfect recall
            #updated 20170405
            recall=1
        else:
            recall = TP/len(self.true_terms[protein])
        #recall should not have zerodivision problem
        #since if self.predicted_bench[protein] is not None
        #This protein is in the benchmark file
        #i.e. gained experimental annotation
        #len(self.true_terms[protein]) should not be 0
        return (precision, recall, count , TP)

                        
    def precision_recall(self,threshold,mode):
        '''
        this calculates the overall (average) precision recall of the team, given a threshold,
        For one prediction file
        mode = 'full': recall is averaged over all benchmark proteins, those not predicted have 0 recall
        mode = 'partial': recall is averaged over all benchmark proteins that are predicted
        '''
        prec = float(0)
        self.counta[threshold] = 0
        rec = float(0)
        for prot in self.predicted_bench:
            a,b = self.term_precision_recall(threshold,prot)[0:2]
            if a is not None:
                prec +=a
                self.counta[threshold]+=1
            if b is not None:
                rec +=b      
        if mode=='partial':
            try:
                recall = rec/self.countb
            except ZeroDivisionError:
                recall = 0
                print("No protein in this predicted set became benchmarks\n")
        elif mode=='full':
            try:
                recall = rec/self.countc
            except:
                recall = 0
                print("No protein in this benchmark set\n")
        try:
            precision = prec/self.counta[threshold]   
        except ZeroDivisionError:
            precision = 0
            print("No prediction is made above the %.2f threshold\n" % threshold)
        return (precision, recall)
    
    def getNumProteins(self,threshold,mode):
         
        '''
        run precision_recall first
        '''         
        print('number of benchmark proteins: %s\n'% len(self.true_terms))
        #Those with not-None recall: 
        #(predicted protein that are in benchmark, only one species per prediction file!! )
        print ('number of predicted proteins that are in the benchmark file: %s\n' % self.countb)
        #Those with not-None precision:
        try:
            print('number of proteins with at least one term above threshold: %s\n' % self.counta[threshold] )
        except KeyError:
            sys.stderr.write("Run precision_recall(%s) first\n" % str(threshold))
            
    def Fmax_output(self,mode):
        '''
        returns the fmax value AND outputs the precision-recall values for each threshold
        This computes precision and recall for every threshold
        mode can be 'full' or 'partial'
        full mode
        '''
        fmax = 0
        f_thres = 0.00
        pre = []
        rec = []
        for thres in numpy.arange(0.00,1.01,0.01,float):
            thres = numpy.around(thres,decimals=2)
            a,b = self.precision_recall(thres,mode)
            if a==0:
                #No prediction above this threshold 
                break
            else:
                pre.append(a)
                rec.append(b)
                f = 2*a*b/(a+b)
            if f>=fmax:
                fmax = f
                f_thres = thres
        print("Number of points in the P-R plot is %s\n" % len(pre))
        #coverage = float(self.countb)/sum([i!=set() for i in self.true_terms.values()])
        coverage = float(self.countb)/self.countc
        return ([pre,rec,fmax,f_thres,coverage])
    
    def printConfidence(self,output_path):
        '''
        print confidence and True/False to a file
        '''
        protindex = 0
        out = open(output_path,'w')
        for prot in self.predicted_bench:
            protindex += 1
            for term in self.predicted_bench[prot]:
                out.write("%s\t%s\t%s\n" % (str(protindex), self.predicted_bench[prot][term][0],self.predicted_bench[prot][term][1]))
        out.close()
        


