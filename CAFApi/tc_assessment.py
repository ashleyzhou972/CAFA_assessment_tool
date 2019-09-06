# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 16:02:45 2017
term-centric evaluations
@author: nzhou
"""

#There is a difference between if a file was created but empty (this species is predicted, but the term was not predicted)
#And when there is no file, the species was not predicted

import sys
sys.path.append('../')
import os
os.chdir('/home/nzhou/git/newCAFAAssess/termcentric')
#This is temporary, don't forget to remove
from Ontology.IO import OboIO
import os
from collections import defaultdict
root_terms = ['GO:0008150','GO:0005575','GO:0003674']
import numpy 
from sklearn import metrics


def go_ontology_ancestors_write(obo_path):
    """
    Input: an OBO file
    Output: One Ancestor file
    """
    terms = set({})
    out = open("%s_ancestors_all.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_parser = OboIO.OboReader(open(obo_path))
    go = obo_parser.read()
    for node in go.get_ids():
        terms.add(node)        
    for term in terms:
        ancestors = go.get_ancestors(term)
        out.write("%s\t%s\n" % (term,",".join(ancestors)))
    out.close()




def read_prediction_and_propagate(prediction_path, goterm, ancestor_path, outfolder):
    #This reads the protein-centric prediction file and extracts records that predicts the goterm of interest
    #as well as its child terms
    #We propagate all predictions first and then filter out all other terms
    taxon = os.path.splitext(os.path.basename(prediction_path))[0].split('_')[2]
    ancestors = defaultdict(set)
    # Read GO ancestors file generated with go_ontology_ancestors_write()
    # File format: 
    # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
    with open(ancestor_path) as ancestors_input:
        for inline in ancestors_input:
            inrec = inline.strip().split('\t')
            term = inrec[0]
            if len(inrec) == 1:
                ancestors[term] = set({})
            else:
                term_ancestors = inrec[1]
                ancestors[term] = set(term_ancestors.split(','))
    #below reads prediction file
    data = defaultdict(list)        
    obsolete = set()        
    predicted = dict()    
    #should be dictionary of with protein as key and confidence as value
    if os.path.isfile(prediction_path):
        if os.path.getsize(prediction_path) > 0:
            exist = True
            for inrec in open(prediction_path,'r'):
                if inrec.startswith('T'):
                    #predictions
                    fields = [i.strip() for i in inrec.split()]
                    data[fields[0]].append({'term': fields[1], 'confidence': float(fields[2])})
                
            #Propagated prediction
            for prot in data:
                for tc in data[prot]:
                    try:
                        ancterms = ancestors[tc['term']].difference(root_terms)
                        #delete root terms
                        #This only delete terms in ancestors, not if the prediction 
                        #itself contains root terms!
                        #modified on 20170203
                    except KeyError:
                        #sys.stderr.write("%s not found\n" % tc['term'])
                        obsolete.add(tc['term'])
                        continue
                    ancterms.add(tc['term'])
                    #include itself in the set of terms
                    if goterm in ancterms:
                        if prot in predicted.keys():
                            #this protein already exists
                            if tc['confidence']>predicted[prot]:
                                predicted[prot]=tc['confidence']
                        else:
                            predicted[prot]=tc['confidence']
            del data
            outfilename = os.path.join(outfolder, "TC_"+os.path.splitext(os.path.basename(prediction_path))[0]+"_"+goterm+".txt")
            with open(outfilename,'w') as outhandle:        
                for prot in predicted:
                    outhandle.write("%s\t%s\n" % (prot, predicted[prot]))
        else:
            exist=False
            print('No prediction made on this term %s\n' % goterm)
    else:
        exist = False
        print('No prediction made for %s\n' % taxon)
    return(exist, outfilename)            



def checkPropagation(file1,file2):
    #TODO
    set1 =set()
    set2 = set()
    with open(file1,'r') as f1:
        for line in f1:
            prot = line.strip().split()[0]
            set1.add(prot)
    
    with open(file2,'r') as f2:
        for line in f2:
            prot = line.strip().split()[0]
            set2.add(prot)
    
def checkUniqueProtein(file1):
    set1 =set()
    linecount = 0
    with open(file1,'r') as f1:
        for line in f1:
            linecount+=1
            prot = line.strip().split()[0]
            set1.add(prot)
    if len(set1)==linecount:
        return(True)
    else:
        return(False)

def read_benchmark(groundtruth):
    #groundtruth file should have two column: 1.cgd id of the protein 2.whether it has function :goterm (T or F)
    y_true = []
    with open(groundtruth, 'r') as gt:
        for line in gt:
            y = line.strip().split()[1]
            if y=='T':
                y_true.append(1)
            else:
                y_true.append(0)
    y_true = numpy.array(y_true)
    return(y_true)
    
def readMappingFile(mappingfile):
    cafaiddict = dict()
    #the cafa3 mapping file between CAFA ID and cgd id
    with open(mappingfile,'r') as mf:
        for line in mf:
            cafaid = line.strip().split('\t')[0]
            cgd_id = line.strip().split('\t')[1]
            cafaiddict[cgd_id] = cafaid
    return(cafaiddict)
    
    
def readTCfile(TC_file, cafaiddict, groundtruth):
    #TC_FriedbergLab_1_237561_GO:0008150.txt should have two columns, first column is protein (CAFAid), 
    #second column is confidence
    #return two dictionaries, one for TPs, one for FPs
    #Check for duplicate proteins before using this!
    #one protein can only appear once!
    tc_dict = dict()
    y_score = []
    if os.path.isfile(TC_file):
        if os.path.getsize(TC_file):
            with open(TC_file,'r') as tc:
                for line in tc:
                    if line.startswith('AUTHOR') or line.startswith('MODEL') or line.startswith('KEYWORDS') or line.startswith('END'):
                        continue  #This is temporary
                        #TODO check AUTHOR and MODEL agrees with file name
                        #TODO save KEYWORDS
                    fields = line.strip().split()
                    tc_dict[fields[0]] = fields[1]
            with open(groundtruth,'r') as gt:
                for line in gt:
                    cgd = line.strip().split()[0]
                    cafaid = cafaiddict[cgd]
                    if cafaid not in tc_dict.keys():
                        #this groundtruth protein is not predicted
                        score =  0.0
                    else:
                        score = round(float(tc_dict[cafaid]),2)
                    y_score.append(score) 
            y_score = numpy.array(y_score)
            return(y_score)
        else:
            print('TC file exists but is empty.')
            return(None)
    else:
        print('TC file does not exist.')
        return(None)
                
def getAUC(y_true, y_score):
    fpr, tpr, thresholds = metrics.roc_curve(y_true,y_score)
    #plt.plot(fpr,tpr)
    auc = metrics.roc_auc_score(y_true,y_score)
    return(fpr,tpr,thresholds, auc)

"""
updated 20180607
use precision-recall as additional metric
"""

def getPR(y_true, y_score):
    pr, rc, thres = metrics.precision_recall_curve(y_true, y_score)
    #print(pr)
    #print(rc)
    #print(thres)
    numpy.seterr(divide='ignore', invalid='ignore')
    f1 = 2*(pr*rc)/(pr+rc) #this F1 disregards threshold
    #updated 20181130
    ap = metrics.average_precision_score(y_true, y_score)
    #f1_from_package = metrics.f1_score(y_true, y_score)
    f1_from_package=None
    f1_max = max(f1)
    max_thres = numpy.nanargmax(f1)
    f1_0 = f1[0]  #The first in the array corresponds to when threshold is 0
    #where all truth are considered predicted
    return(ap, f1_0, f1_max, f1_from_package, max_thres, pr, rc, thres, f1)



resultfile = '/home/nzhou/git/CAFA3_termcentric_assessment/candidaResults/TC_CaoLab7_3_237561_GO:0042710_results.txt'


def read_resultfile(resultfile):
    with open(resultfile,'r') as rf:
        auc = rf.readline().strip()
        line2 = rf.readline()
        line3 = rf.readline()
        fpr = line2.strip().split()
        tpr = line3.strip().split()
    return(fpr,tpr)

        
        
def plotMultiple(title,listofResults,smooth):
    '''
    supply lists of precision+recall+name lists
    '''
    fontP = FontProperties()
    fontP.set_size('small')
    num = len(listofResults)
    if num>1:
        pal=sns.color_palette("Paired", num)
        colors=pal.as_hex()
        for j,i in enumerate(listofResults):
            if j==num-2 or j==num-1:
                linetype = '--'
            else:
                linetype = '-'
            if smooth=='Y':
                ax = plt.subplot()
                precision = curveSmooth(i)[0][1:]
                recall = curveSmooth(i)[1][1:]
                ax.plot(recall,precision,linetype,color=colors[j],label=i.author+':\nF=%s C=%s'%(i.opt,i.coverage)) 
                ax.plot(i.recall[int(i.thres*100)],i.precision[int(i.thres*100)],'o',color=colors[j])
            elif smooth=='N':
                ax = plt.subplot()
                ax.plot(i.recall,i.precision,linetype,color=colors[j],label=i.author+':\nF=%s C=%s'%(i.opt,i.coverage))
                ax.plot(i.recall[int(i.thres*100)],i.precision[int(i.thres*100)],'o',color=colors[j])
        plt.axis([0,1,0,1])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        plt.yticks(numpy.arange(0,1,0.1))
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title(title)
        if title == None:
            figurename = './plots/Combined_plot.png'
        else:
            figurename = './plots/'+title+'.png'
        
        plt.savefig(figurename,dpi=200)
        plt.close()
    
