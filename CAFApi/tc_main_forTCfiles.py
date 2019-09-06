# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 15:47:06 2017
term_centric_main
This reads TC prediction directly
@author: nzhou
"""

import tc_assessment as tca
import argparse
import os
import numpy as np

if __name__=='__main__':
    '''
    groundtruth = "/home/nzhou/git/CAFA3_annotation/termcentric/candida_0042710.txt"
    obo_path = "/home/nzhou/git/CAFA3_annotation/go_20170429.obo"
    taxon = "237561"
    goterm = "GO:0042710"
                    
    prediction_path = '/home/nzhou/git/CAFA3_submissions/prot/ZhuLab2/all_submissions_ZhuLab2/ZhuLab2_1_237561.txt'     
    ancestor_path ='/home/nzhou/git/CAFA3_annotation/go_20170429_ancestors_all.txt' 
    outfolder = '/home/nzhou/git/CAFA3_submissions/tc/'
    #go_ontology_ancestors_write(obo_path)
    
    file1 = '/home/nzhou/git/CAFA3_submissions/tc/TC_ZhuLab2_1_237561_GO:0042710.txt'
    file2 = '/home/nzhou/git/CAFA3_submissions/tc/TC_ZhuLab2_1_237561_0042710.txt'
    a= read_prediction_and_propagate(prediction_path, goterm, ancestor_path, taxon, outfolder)  
    lenPos, lenNeg, benchset = read_benchmark(groundtruth)    
    mappingfile = '/home/nzhou/git/CAFA3_annotation/Mapping files/mapping.237561.map'
    tpd, fpd = readTCfile(file1,benchset)
    auc = getAUC(tpd, fpd, lenPos, lenNeg)
    '''
    def extant_file(x):
        if not os.path.isfile(x):
            raise argparse.ArgumentTypeError("{0} does not exist".format(x))
        else:
            return(x)
    
    parser = argparse.ArgumentParser(description='CAFA3 term-centric evaluations', )
    parser.add_argument('file',type=extant_file,
                       help='Input term-centric prediction file. Filename should follow CAFA formats.')
    parser.add_argument('goterm', type = str, help = 'GO term of interest')                   
    parser.add_argument('-t','--t' ,dest='taxon', help = 'Species to evaluate')
    parser.add_argument('-b','--b',dest='benchfile', help = 'Input the path to the benchmark file')
    #parser.add_argument('-out','--out', dest='ofolder', help = 'Output folder for the generated TC file')
    parser.add_argument('-r','--r', dest='resultfolder', help = 'Output folder for result')
    parser.add_argument('-m','--m', dest='metric', type=str, 
                        help = 'Evaluation metric: default is AUROC',
                        default = 'ROC', choices = ['ROC', 'PR'])
    args = parser.parse_args()
    #First generate TC file from protein-centric submission file
    tc_file = args.file
    exist=True
    if exist:
        y_true = tca.read_benchmark(args.benchfile)
        #print(len(y_true))
        #print(list(y_true))
        mappingfile = '/home/nzhou/git/CAFA3_annotation/Mapping files/mapping.%s.map' % args.taxon
        cafaiddict = tca.readMappingFile(mappingfile)
        y_score = tca.readTCfile(tc_file,cafaiddict, args.benchfile)
        #print(len(y_score))
        #print(list(y_score))
        if y_score is not None:
            if args.metric=='ROC':
                fpr, tpr, thres, auc = tca.getAUC(y_true, y_score)
                outfilename = os.path.splitext(os.path.basename(tc_file))[0]+'_ROC_results.txt'
                outhandle = open(os.path.join(args.resultfolder,outfilename), 'w')
                outhandle.write('AUC:%s\n' % auc)
                outhandle.write('FPR\tTPR\tThreshold\n')
                for index, item in enumerate(fpr):
                    outhandle.write('%s\t%s\t%s\n' % (item, tpr[index], thres[index]))
                outhandle.close()
            else:
                #updated 20181226
                #only two metrics allowed, ROC or PR
                #updated 20181130
                ap, f1_0, f1_max, f1_from_package, max_thres, pr, rc, thres, f1 = tca.getPR(y_true, y_score)
                #add dummy value to end of thres
                thres=np.append(thres,1)
                outfilename = os.path.splitext(os.path.basename(tc_file))[0]+'_PR_results.txt'
                outhandle = open(os.path.join(args.resultfolder,outfilename), 'w')
                outhandle.write('AP:%s\n' % ap)
                outhandle.write('F1(0):%s\n' % f1_0)
                #outhandle.write('F1 from package:%s\n' % f1_from_package)
                outhandle.write('Fmax:%s\t%s\n' % (f1_max, max_thres))
                outhandle.write('Precision\tRecall\tThreshold\tf1\n')
                for index, item in enumerate(pr):
                    outhandle.write('%s\t%s\t%s\t%s\n' % (item, rc[index], thres[index],f1[index]))
                outhandle.close()
            
            
    

