# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 15:02:03 2016
This creates the command-line interface of the precision-recall calculations in
CAFA assessment tool
@author: Ashley Zhou
last updated 12/29/2016
"""

import argparse
import tests
from precrec.precRec import PrecREC,read_benchmark
from precrec.GOPred import GOPred
import precrec.preprocess as pp
import matplotlib.pyplot as plt
import numpy
import os

'''
arguments supplied to main:
updated on 07/12/2016
    ontology: 'BPO', 'MFO', 'CCO', 'HPO'
    team number: e.g. 117
    species: use taxon ID e.g. 9606 for Homo Sapien
    model: which model of the submission, e.g. 1,2 or 3 
    saveplot path: the path where the PR curve plot is saved (should be an absolute path??)
    #use default go.obo: True or supply obopath (should be an absolute path??)
    #use default benchmark: True or supply benchmark path (Ontology-specific!!) (should be an absolute path)
    #submission folder: folder containing all submissions with each team in separate  folders titled by team number
'''

def get_namespace_index(namespace):
    '''
    copied from confidence.py 07/15/2016
    '''
    num = None
    if namespace=='BPO' or namespace=='bpo':
        num = 0
    elif namespace=='MFO' or namespace=='mfo':
        num = 1
    elif namespace=='CCO' or namespace=='cco':
        num =2
    else:
        raise ValueError("name space not found, check prediction files")
        print(namespace)
    return num

'''
def compare_curves(numberOfTeams,teamIDs):
    #teamIDs should be a string separated by comas(,)
    #e.g. compare_curves(3,'117,118,119')
    #TODO
    

'''
'''
def team_name_converter(team_number):
    
    #convert team number to team name for plotting purposes
    #07/18/2016: waiting for input from actual CAFA3
    
    return "teamName"
'''
def test_output_pr(precREC, filename):
    #This function outputs the specific precision and recall and number of GO terms for each protein
    #The above output is saved in a file 
    #precRec should be a precREC object
    outfile = open(filename,'w')
    for thres in numpy.linspace(0.00,1.00,101):
        thres=numpy.around(thres,decimals=2)
        for prot in precREC.predicted_bench:
            p,r,count,TP = precREC.term_precision_recall(thres,prot)
            outfile.write('%s:%s:%s\t%s\t%s\t%s\n' % (prot,thres,p,r,TP,count))
    outfile.close()
    
    
def test_output_terms(precREC, prot, filename):
    #This function outputs the GO terms associated with a certain protein, 
    #predicted, propagated and with confidence values
    #precREC is a precREC object
    #prot is a protein CAFA ID
    outfile = open(filename,'w')
    for i in precREC.predicted_bench[prot]:
        outfile.write('%s\t%s\n' % (i,precREC.predicted_bench[prot][i]))
    outfile.close()
    #print(precREC.term_precision_recall(0.94,prot))
    

def taxon_name_converter(taxonID):
    #convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
    taxonTable = {'10116':'RAT','9606':'HUMAN','3702':'ARATH','7955':'DANRE','44689':'DICDI','7227':'DROME','83333':'ECOLI','10090':'MOUSE','208963':'PSEAE','4896':'SCHPO','4932':'YEAST'}
    return taxonTable[taxonID]    


    #after readBenchmark, the benchmark sets are ontology-specific, but do not distinguish between NK and LK
    #Get NK/LK benchmark sets, select the groundtruths from lists, from Supplementary_data/data/benchmark
    #list_dir is the directory path where the CAFA2 lists of benchmark protein IDs are with different criteria
    #The criteria are: species/all, type1/type2/typex (type1 is NK, type2 is LK )

    
    
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
    #12202016 deleted input needs for team name, ONTOLOGY and taxon, as can be found in file name
    #parse everything from filename without input    
    #parser.add_argument('ontology',help='Input ontology',choices=['BPO','MFO','CCO'])
    #parser.add_argument('team 1',help = 'Input team number',type=int)
    #parser.add_argument('taxon', help= 'Input taxon ID, this will only be used to name the plot', type=int)
    #If it's all species combined, enter 0 as taxon ID(07/18/2016)
    #parser.add_argument('model',help = 'Input model number', choices=['1','2','3'])
    #parser.add_argument('filetype', type=int,help='If prediction file is ontology-specific, enter 1; if prediction file is raw (from team submission), enter 2')
    parser.add_argument('file',type=open,
                        help='Input the path of the prediction file. Filename should follow CAFA formats')
    #CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    #example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    #If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)                  
    
    #parser.add_argument('plotfile', help='Input path+filename to save the PR plot')

    parser.add_argument('type',help = 'Input evaluation type: No Knowledge or Limited Knowledge', choices=['type1','type2','all'])    
    parser.add_argument('teamnum',help='Input team number',type=int)
    parser.add_argument('mode', help = 'Input the evaluation mode: full or partial', choices = ['full','partial'])
    args = parser.parse_args()
    
    print('Evaluating %s.\n' % args.file.name)
    

    '''
    #TODO multiple curves on one plot
    plt.plot(fm[1],fm[0])
    plt.axis([0,1,0,1])
    plt.yticks(numpy.arange(0,1,0.1))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(str(args.team)+" "+str(args.taxon)+" "+args.ontology)
    plt.savefig(args.plotfile,dpi=200)
    plt.close()
    
    print('fmax value for this prediction is: %s.\n' % fm[2])
    print('PR plot is saved to %s.\n' % args.plotfile)
    '''
    
    
    #official code:
    #first split the prediction file into three ontologies
    all_pred = GOPred()
    pred_path = args.file
    obo_path = './precrec/go_20130615-termdb.obo'
    all_pred.read_and_split_and_write(obo_path,pred_path)
    info = [all_pred.author,all_pred.model,all_pred.keywords,all_pred.taxon]
    print('author: %s\n' % info[0])
    print('model: %s\n' % info[1])
    print('keywords: %s\n' % info[2])
    print('species:%s\n' % info[3])
    
    #parse file name 
    namefields = os.path.basename(pred_path.name).split('.')[0].split('_')
        
    #read benchmark for three ontologies
    for onto in ['bpo','cco','mfo']:
        print('ontology: %s\n' % onto)
        b = read_benchmark(onto, taxon_name_converter(namefields[2]),args.type,'./precrec/benchmark/',obo_path)
        path = os.path.splitext(pred_path.name)[0]+'_'+onto.upper()+'.txt'
        c = PrecREC(b,path)
        #test_output_pr(c,'./TianLab_10090_2_pr_'+onto+'.txt')
        #test_output_terms(c,'T100900003687','./TianLab_10090_2_terms_'+onto+'.txt')
        if c.exist:
            fm = c.Fmax_output(args.mode)
            print('fmax: %s\n' % fm[2])
            print('threshold giving fmax: %s\n' % fm[3])
            print('coverage: %s\n' % fm[4])
            yx = tests.read_Cafa2_sheet(onto, info[3],args.teamnum,info[1],args.type,args.mode)
            if yx[1]==round(fm[2],3):
                print('Fmax Match!\n')
                if yx[0]==round(fm[4],2):
                    print('Coverage Match!\n')
                    if yx[2]==fm[3]:
                        print('Threshold Match!\n')
                '''
            else:
                #test_output_pr(c,'./TianLab/'+str(info[0])+'_'+str(info[3])+'_'+str(info[1])+'_pr_'+onto+'.txt')
                test_output_terms(c,'T'+str(info[3])+'00013457',str(info[0])+'_'+str(info[3])+'_'+str(info[1])+'_terms_'+onto+'.txt')
                '''

    #Below are codes for testing
    #os.chdir('./CAFAAssess')
    '''
    obo_path = './precrec/go_20130615-termdb.obo'
    #all_pred = GOPred()
    #all_pred.read_and_split(obo_path,pred_path)
    #info = [all_pred.author,all_pred.model,all_pred.keywords,all_pred.taxon]
    #b=read_benchmark('bpo','RAT',args.type,'./precrec/benchmark/',obo_path)
    b=read_benchmark('bpo','RAT','type1','./precrec/benchmark/',obo_path)
    #pred_path = open('/home/nzhou/git/CAFAAssess/TianLab/TianLab_2_10116.txt','r')
    #pred_path = args.file
    #bpo_path = os.path.splitext(pred_path.name)[0]+'_BPO.txt'
    #bpo_path = './TU_1_10116_CCO.txt'
    bpo_path = './Doegroup_1_10116_BPO.txt'
    c = PrecREC(b, bpo_path)
    
    for thresnumpy.arange(0.00,1.01,0.01,float)
    fm = c.Fmax_output('full')
    print('fmax: %s\n' % fm[2])
    print('threshold giving fmax: %s\n' % fm[3])
    print('coverage: %s\n' % fm[4])
   

    
    outfile = open('./doegroup_pr_seq.txt','w')
    for thres in numpy.linspace(0.00,1.00,101):
        for prot in c.predicted_bench:
            p,r,count,TP = c.term_precision_recall(thres,prot)
            outfile.write('%s:%s:%s\t%s\t%s\t%s\n' % (prot,thres,p,r,TP,count))
    #fm = c.Fmax_output()
    print(fm[2],fm[3])
   '''
   
   
   
'''   
    outfile1 = open('./doe_5839_predicted_propagated_conf.txt','w')
    for i in c.predicted_bench['T101160005839']:
        outfile1.write('%s\t%s\n' % (i,c.predicted_bench['T101160005839'][i]))
    outfile1.close()
'''
'''    
2*fm[1][11]*fm[0][11]/(fm[0][11]+fm[1][11])
benchmarkprot = 'T101160001453'
nonbenchmarkprot = 'T101160007359'

for thres in numpy.linspace(0.01,0.99,99):
    pr,rc = c.term_precision_recall(thres,benchmarkprot)
    print(pr,rc,thres)
    
#On benchmarkprot, the CAFA participant predicted 1217 GO terms, 13 of them are true    


prec = float(0)
for prot in c.predicted:
    a=c.term_precision_recall(0.68,prot)[0]
    if a is not None:
        print(prot)
        prec +=a
        
        
c.term_precision_recall(0.68,'T101160006596')
c.predicted['T101160006596']['GO:0008150']
'''
'''
a=open('doe_5839_predicted_propagated.txt','w')
for i in c.predicted_bench['T101160005839'].keys():
    a.write('%s\n' % i )
a.close()
f = open('doe_5839_true_propagated.txt','w')
for i in b.true_terms['T101160005839']:
    f.write('%s\n' % i)
f.close()
'''    
