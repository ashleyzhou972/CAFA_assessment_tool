# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
 CAFA module
 
"""

import sys
import math
import collections
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

######################################BENCHMARK START#################################
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
    
    obo_bpo_out = open("%s_ancestors_bpo.txt" % (os.path.splitext(obo_path)[0]), "w")
    obo_cco_out = open("%s_ancestors_cco.txt" % (os.path.splitext(obo_path)[0]), "w")
    obo_mfo_out = open("%s_ancestors_mfo.txt" % (os.path.splitext(obo_path)[0]), "w")
    obo_parser = OboIO.OboReader(open(obo_path))
    go = obo_parser.read()
    mfo_terms, bpo_terms, cco_terms = go_ontology_split(go)
    for term in mfo_terms:
        ancestors = go.get_ancestors(term)
        obo_mfo_out.write("%s\t%s\n" % (term, ",".join(ancestors)))
    for term in bpo_terms:
        ancestors = go.get_ancestors(term)
        obo_bpo_out.write("%s\t%s\n" % (term, ",".join(ancestors)))
    for term in cco_terms:
        ancestors = go.get_ancestors(term)
        obo_cco_out.write("%s\t%s\n" % (term, ",".join(ancestors)))

    obo_mfo_out.close()
    obo_bpo_out.close()
    obo_cco_out.close()
    return([len(bpo_terms), len(cco_terms), len(mfo_terms)])
    
    
def read_benchmark(namespace, species, types, fullbenchmarkfolder, obopath):
    '''
    Read Benchmark.
    
    Input:
    namespace --
    species --
    types --
    fullbenchmarkfolder--
    obopath --
    
    Output:
    bench --
    obocountDict --
    '''
    # Ancestor files here are precomputed
    legal_types = ["type1","type2","typex"]
    legal_subtypes = ["easy","hard"]
    legal_namespace = ["bpo","mfo","cco","hpo"] 
    # fullbenchmarkfolder = './precrec/benchmark/'
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
    obocountDict = {'bpo':obocounts[0],'cco':obocounts[1],'mfo':obocounts[2]}
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
        proteins  = set()
        for line in handle:
            proteins.add(line.strip())
        handle.close()
        tempfilename = 'temp_%s_%s_%s.txt' % (namespace, species, types)
        tempfile = open(fullbenchmarkfolder+'/'+tempfilename , 'w')
        for line in open(full_benchmark_path,'r'):
            protein = line.split('\t')[0]
            if protein in proteins:
                tempfile.write(line)
        tempfile.close()
        bench = benchmark(ancestor_path, tempfile.name)
        bench.propagate()
        os.remove(tempfile.name)
    else:
        print('Benchmark set is empty.\n')
        bench = None
    return(bench, obocountDict)   
    
########################END OF FUNCTIONS OUTSIDE OF A CLASS ##############################
    
class benchmark:
    '''
    
    '''
    
    
    def __init__(self,ancestor_path,benchmark_path):
        '''
        Initialize the benchmark.
        
        Input: 
        benchmark_path -- ontology specific file location
        ancestor_path -- ontology specific file location
        '''
        
        # Key: protein
        # Value: set of benchmark leaf terms
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
        
        self.true_terms = defaultdict(set)
        # Key: protein
        # Value: set of benchmark propagated terms
        
        for protein in self.true_base_terms:
            for term in self.true_base_terms[protein]:
                try:
                    ancestors = self.ancestors[term].difference(root_terms)    
                
                except KeyError:
                    sys.stderr.write("%s not found \n" % term) 
                self.true_terms[protein].add(term)
                self.true_terms[protein] |= ancestors
        

#####################################BENCHMARK END ############################

###############################################################################
#################################### START OF INFO ############################
###############################################################################

class Info:
    ''' Holds all the stuff we need, allows for just importing once '''
    
   
    def __init__(self, benchmark, os_prediction_path, obocounts):
        '''
        Initalize.
        
        Input:
        benchmark -- instance of the benchmark class
        os_prediction_path -- !ontology-specific! prediction file 
                              separated in GOPred (without headers)
        obocounts -- total number of terms in the ontology
        '''
        
        # Initialize variables
        self.exist                            = True
        self.ancestors                        = benchmark.ancestors
        self.true_terms                       = benchmark.true_terms
        self.obocount                         = obocounts
        # Set of all obsolete terms found
        self.obsolete                         = set()
        # count_above_threshold is the number of proteins 
        # with at least one term above threshold
        self.count_above_threshold            = defaultdict()
        # count_predictions_in_benchmark is the number of predicted proteins 
        # in this file that are in the benchmark file (for coverage)
        self.count_predictions_in_benchmark   = 0
        self.count_true_terms                 = len(self.true_terms)
        # self.data is the dictionary for all predicted terms
        self.data                             = defaultdict(list)
        # self.predicted_bench is the dictionary 
        # for all predicted terms that are benchmarks
        self.predicted_bench                  = defaultdict(defaultdict)
        # key:protein
        # value: list of dictionaries
        # key: GO term
        # value: tuple(confidence, Boolean) Boolean is whether in self.true_terms     
        
        # Read in prediction file        
        if os.path.getsize(os_prediction_path) > 0:
            for inrec in open(os_prediction_path,'r'):
                fields = [i.strip() for i in inrec.split()]
                self.data[fields[0]].append({'term': fields[1], 'confidence': float(fields[2])})
            
            # Propagated prediction
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
                                self.predicted_bench[protein][tc['term']] = self.compare(protein, tc)
                                for ancterm in ancterms:
                                    #Make a new TC with ancestors
                                    newtc = {'term':ancterm,'confidence':tc['confidence']}
                                    if ancterm in self.predicted_bench[protein]:
                                        # Term already exists, update confidence
                                        self.update_confidence(protein, newtc)
                                    else:
                                        # Add term to self.predicted_bench
                                        self.predicted_bench[protein][ancterm] = self.compare(protein, newtc)
                                        
            if self.count_predictions_in_benchmark==0:
                self.exist = False
                print("No protein in this predicted set became a benchmark\n")
        else:
            self.exist=False
            print('No prediction made in this ontology.\n')
            
            
    def coverage(self):
        ''' Determine the coverage. '''
        
        return float(self.count_predictions_in_benchmark)/self.count_true_terms  
       
     
    def update_confidence(self, protein, tc):
        '''
        Update Confidence for given protein and propagate.
        
        This function compares the confidence value in tc to the confidence in self.predicted
        If tc is larger, than it overwrites the confidence in self.predicted
        
        Input:
        protein -- chosen protein 
        tc -- a dictionary format as: {'confidence':0.57,'term':'GO:006644'}
        '''
        
        # Defined for readablity
        confidence = tc['confidence'] 
        term = tc['term']
        
        if confidence > self.predicted_bench[protein][term][0]:
            # Update the confidence
            self.predicted_bench[protein][term][0] = confidence
            # Propagate changes if necessary
            for ancterm in self.ancestors[term].difference(root_terms):
                if confidence > self.predicted_bench[protein][ancterm][0]:
                    # Update the confidence
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
        
        
    def check(self, tool, mode):
        '''
        Effective main function, call this to get results
        
        Inputs:
        tool -- {Fmax, WFmax, Smin, NSmin, ALL}
        mode -- {full, partial, both}
        
        Output:
        result --  A result object with all needed info as fields
        '''
        r = Result()
        #### FOR NOW ####
        if(tool == "Fmax"):
            f = Fmax()
            return f.output(self, mode)
        elif(tool == "WFmax"):
            WFmax.output(self, mode)
        elif(tool == "Smin"):
            Smin.output(self, mode)
        elif(tool == "NSmin"):
            NSmin.output(self, mode)
        else:
            #THROW ERROR
            1
        
        return r #Results
        
############################################# END OF INFO ###########################################################

class Result:
    ''' Stores results in a common format '''   
    def __init():
        ''' State all variables needed '''

        #Fmax
        Fmax      = 0.0
        PR        = []
        RC        = []
        Threshold = 0.0
        
        #
        WFmax = 0.0
        Smin  = 0.0
        NSmin = 0.0
        
        #IC Value
        IC = 0 #Likely an array/list
        
        
        #####Add over info that you would want to output #####

class IC:
    '''
    Information Content
    
    Some code in this section taken / based on https://github.com/ppypp/debias/blob/master/lib/debias.py
    '''
        
    def __init__(info):
        '''
        Make local copy of needed data 
        
        Copy info and create data dictionary
        '''        
        data = dict()
        counter = 1
        #for every annotation in our set ###########FIGURE OUT WHERE THAT WOULD BE##############
        for annotation in something: #GAF in the other code 
            id = "annotation" + str(counter)
            data[id] = annotation
            counter += 1
        #Data is setup
            
            
        
    def GO_term_frequency(self):
        ''' Count fequency for each Go term '''
        
        frequency = dict()
        data = self.data
        for annotation in data:
            if data[annotation]['GO_ID'] in frequency:
                frequency[data[annotation]['GO_ID']] += 1
            else:
                frequency[data[annotation]['GO_ID']] = 1
        return frequency    
        
        
    def PL_IC(self, threshold):
        '''
        Calculate Phillip Lord Iimformation Content
        
        Input:
        threshold --
        '''
        
        data = self.data
        
        go_terms = []
        results = dict()
        ic = dict()
        
        for annotation in data:
            go_terms.append(data[annotation]["GO_ID"])
            
        # Makes an object that contains the count of each Go_term
        GO_term_to_PL_info=collections.Counter(go_terms)
                
        for term in GO_term_to_PL_info:
            ic_term = -math.log(GO_term_to_PL_info[term] / len(go_terms), 2)
            ic['term'] = ic_term
              
        for attnid in data:
            annotation=data[attnid]
            if GO_term_to_PL_info[annotation["GO_ID"]]>=threshold:
                results[attnid]=data[attnid]
        
        return results
    
    
##################################################################################################
##################################################################################################
##################################################################################################    
    
class Fmax:
    '''
    F maximum
    '''

    def f(self, precision, recall):
        ''' Calculate F function '''
        
        try:
            f = (2*precision*recall)/(precision+recall)
        except ZeroDivisionError:
            f = None
        return f


    def output(self, info, mode):
        ''' 
        Calculate the Fmax 
        
        Input:
        info --
        mode -- 
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
            pr, rc = self.PRRC_average(info, threshold, mode)
            if pr is None:
                # No prediction above this threshold 
                break
            else:
                PR.append(pr)
                RC.append(rc)
                # Find the F-value for this particular threshold
                try:
                    f = self.f(pr, rc)
                except ZeroDivisionError:
                    f = None
                    
            if f is not None and f >= fmax: ###########QUESTION##############
                fmax = f
                fmax_threshold = threshold
        #Have found the Fmax at this point       
        return ([PR, RC, fmax, fmax_threshold])
        
     
    def PRRC(self, info, threshold, protein):
        '''
        Calculate the PRRC of a single protein
        
        Input:
        info --
        threshold --
        protein --
        '''
        
        # Initalize Variables
        TP = 0.0      # True positive
        count = 0   # Count how many terms are above the threshold
        TT_length = len(info.true_terms[protein]) # Number of True terms
        
        if(threshold == 0):
            TP = TT_length
            count = info.obocount
        else:
            # For every term related to the protein
            for term in info.predicted_bench[protein]:
             # If it is above the threshold, increment the count
                if info.predicted_bench[protein][term][0] >= threshold:
                    count += 1
                    # If it is actually True, increment TP
                    if info.predicted_bench[protein][term][1] :
                        TP += 1
        # Find PR: TP / (TP + FP)
        try:
            precision = TP / count 
        except ZeroDivisionError:
            precision = None
        # Find RC: TP (TP + FN)
        recall = TP/TT_length
        return (precision,recall)
            
            
    def PRRC_average(self, info, threshold, mode):
        '''
        Calculate the overall PRRC of file
        
        Input:
        info --
        threshold --
        mode -- 
        '''
        
        # Initialize Variables
        PR = 0.0
        RC = 0.0
        info.count_above_threshold[threshold] = 0

        for protein in info.predicted_bench:
            pr,rc = self.PRRC(info, threshold, protein)
            if pr is not None:
                PR += pr
                info.count_above_threshold[threshold] += 1
            if rc is not None:
                RC += rc   
                
        if mode=='partial':
            try:
                recall = RC/info.count_predictions_in_benchmark
            except ZeroDivisionError:
                recall = 0
                print("No protein in this predicted set became benchmarks\n")
                
        elif mode=='full':
            try:
                recall = RC/info.count_true_terms
            except ZeroDivisionError:
                recall = 0
                print("No protein in this benchmark set\n")
                
        try:
            precision = PR/info.count_above_threshold[threshold]   
        except ZeroDivisionError:
            precision = None
            print("No prediction is made above the %.2f threshold\n" % threshold)
           
        return (precision, recall)        
                    
    
#####################################################################################################    
    
class WFmax(Fmax): # Will use parts of Fmax just using IC values
    '''
    Weighted F maximun
    '''
    def ouput(self, info, mode):
        '''
        Calculate WFmax
        '''
        # Intialize Variables
        wfmax = 0.0
        wfmax_threshold = 0.0
        WPR = []
        WRC = []
        
        # Run over all threshold values from 0 to 1, two signifigant digits
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            # Run PRRC on given threshold
            wpr, wrc = self.PRRC_average(info, threshold, mode)
            if wpr is None:
                # No prediction above this threshold 
                break
            else:
                WPR.append(wpr)
                WRC.append(wrc)
                # Find the F-value for this particular threshold
                try:
                    wf = self.f(wpr, wrc)
                except ZeroDivisionError:
                    wf = None
                    
            if wf is not None and wf >= wfmax: ###########QUESTION##############
                wfmax = wf
                wfmax_threshold = threshold
        #Have found the Fmax at this point       
        return ([WPR, WRC, wfmax, wfmax_threshold])
        
        
    def WPRRC_average(self, info, threshold, mode):
        '''
        Calculate the overall PRRC of file
        
        Input:
        info --
        threshold --
        mode -- 
        '''
        
        # Initialize Variables
        WPR = 0.0
        WRC = 0.0
        info.count_above_threshold[threshold] = 0

        for protein in info.predicted_bench:
            wpr, wrc = self.WPRRC(info, threshold, protein)
            if wpr is not None:
                WPR += wpr
                info.count_above_threshold[threshold] += 1
            if wrc is not None:
                WRC += wrc   
                
        if mode == 'partial':
            try:
                recall = WRC/info.count_predictions_in_benchmark
            except ZeroDivisionError:
                recall = 0
                print("No protein in this predicted set became benchmarks\n")
                
        elif mode == 'full':
            try:
                recall = WRC/info.count_true_terms
            except ZeroDivisionError:
                recall = 0
                print("No protein in this benchmark set\n")
                
        try:
            precision = WPR/info.count_above_threshold[threshold]   
        except ZeroDivisionError:
            precision = None
            print("No prediction is made above the %.2f threshold\n" % threshold)
           
        return (precision, recall) 
        
        
    def WPRRC(self, info, threshold, protein):
        '''
        Calculate the PRRC of a single protein
        
        Input:
        info --
        threshold --
        protein --
        '''
        
        # Initalize Variables
        total = 0.0        
        TP_total = 0.0      # True positive IC sum
        
        if(threshold == 0):
            # This is all predicted terms
            5 #temp code ##############################################################################################
        else:
            # For every term related to the protein
            for term in info.predicted_bench[protein]:
                if info.predicted_bench[protein][term][0] >= threshold:
                    #Add IC value to total
                    total += 0 #WHEREEVER WE HAVE THAT VALUE
                    # If it is actually True, add its IC to TP_total
                    if info.predicted_bench[protein][term][1] :
                        TP_total += 0 #WHEREEVER WE HAVE THAT VALUE
        
                        
        # Find PR: TP / (TP + FP)
        try:
            precision = TP_total / total 
        except ZeroDivisionError:
            precision = None
        # Find RC: TP / (TP + FN)
        recall = TP_total/5 ##We need the FN terms as well, how?###############
        
        return (precision,recall)
    
    
#########################################################################################################################################    
    
class Smin:
    '''
    S minimum
    '''
    def ru(T, P):
        '''
        Calculate Remaining Uncertainity
        
        Input:
        T -- Truth
        P -- Prediction
        '''
        
        total = 0.0
        # Sum the IC values of every element in T not in P
        for term in T:
            if term not in P:
                total += term.ic_value ################################################### MAKE SURE CORRECT FORMAT ONCE DTERMINED
        return total 
        
        
    def mi(T, P):
        '''
        Calculate Misinformation
        
        Input:
        T -- Truth
        P -- Prediction
        '''
        
        total = 0.0
        # Sum the IC values of every element in P not in T
        for term in P:
            if term not in T:
                total += term.ic_value ################################################## MAKE SURE CORRECT FORMAT ONCE DTERMINED
        return total                
        
        
    
    def s(k, ru, mi):
        ''' Semantic Distance '''
        
        s = (ru^k + mi^k)^(1/k)
        return s
        

    def s_average(info, k, threshold, mode):
        '''
        Semantic Distance 
        
        At a particular threshold, averaged over all proteins
        
        Input:
        info --
        k --
        threshold --
        mode --
        '''
        
        
        

    def output(self, info, k, mode):
        '''
        Calculate Smin
        
        k = 2 by convention
        '''
    
        # Intialize Variables
        smin = 0.0
        smin_threshold = 0.0
        RU = []
        MI = []
        
        # Run over all threshold values from 0 to 1, two signifigant digits
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            # Run S on given threshold
            ru, mi = self.s_average(info, k, threshold, mode)
            RU.append(ru)
            MI.append(mi)
            # Find the S-value for this particular threshold
            s = self.s(k, ru, mi)
            if s is not None and s <= smin: ###########QUESTION##############
                smin = s
                smin_threshold = threshold
        # Have found the Smin at this point       
        return ([RU, MI, smin, smin_threshold])
    
    
    ################################## END OF SMIN ###############################
    
class NSmin(Smin):
    '''
    Normalized S minimum
    '''
    
