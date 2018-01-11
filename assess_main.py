#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from precrec.precRec import PrecREC,read_benchmark
from precrec.GOPred import GOPred
import os
import sys
import errno    
import gc
import yaml





def get_namespace_index(namespace):
    '''
    convert namespace into indices
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

def taxon_name_converter(taxonID):
    #convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
    taxonTable = {'10116':'RAT','9606':'HUMAN','3702':'ARATH','7955':'DANRE','44689':'DICDI',
    '7227':'DROME','83333':'ECOLI','10090':'MOUSE','208963':'PSEAE',
    '237561':'CANAX','559292':'YEAST','284812':'SCHPO','8355':'XENLA','224308':'BACSU',
    '99287':'SALTY','243232':'METJA','321314':'SALCH','160488':'PSEPK','223283':'PSESM',
    '85962':'HELPY','243273':'MYCGE','170187':'STRPN','273057':'SULSO','all':'all','prokarya':'prokarya','eukarya':'eukarya'}
    return taxonTable[taxonID]    
 



def typeConverter(oldType):
    if oldType=='type1':
        newType = 'NK'
    elif oldType == 'type2':
        newType = 'LK'
    elif oldType == 'all':
        newType = 'All'
    return(newType)

def extant_file(x):
    if not os.path.isfile(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    else:
        return(open(x,'r'))
        
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
            
def read_config():
    parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
    
    parser.add_argument('config_stream',type=extant_file, help='Configuration file')
    # CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    # example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    # If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)                  
    args = parser.parse_args()
    # Load config file to dictionary
    try:
        config_dict = yaml.load(args.config_stream)['assess']
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()
        
    obo_path = config_dict['obo']
    bfolder = config_dict['benchmark']
    results_folder = config_dict['results']
    f = config_dict['file']
    return(obo_path, bfolder, results_folder, f)


if __name__=='__main__':
    obo_path,benchmarkFolder, resultsFolder,f = read_config()
    mkdir_p(resultsFolder)
    mkdir_p(resultsFolder+'/pr_rc/')
    print('\nEvaluating %s.\n' % f)
    
    all_pred = GOPred()
    pred_path = open(f,'r')
    all_pred.read_and_split_and_write(obo_path,pred_path)
    info = [all_pred.author,all_pred.model,all_pred.keywords,all_pred.taxon]
    #clear memory
    del all_pred
    gc.collect()
    print('AUTHOR: %s\n' % info[0])
    print('MODEL: %s\n' % info[1])
    print('KEYWORDS: %s\n' % info[2][0])
    print('Species:%s\n' % info[3])
    author = info[0]
    taxon = info[3]
    model = info[1]
    keywords = info[2][0]        
    resulthandle= open(resultsFolder+"/%s_results.txt" % (os.path.basename(f).split('.')[0]),'w')
    prhandle = open(resultsFolder+"/pr_rc/%s_prrc.txt" % (os.path.basename(f).split('.')[0]),'w')
    resulthandle.write('AUTHOR:%s\n' % info[0])
    resulthandle.write('MODEL: %s\n' % info[1]) 
    resulthandle.write('KEYWORDS: %s\n' % info[2][0])  
    resulthandle.write('Species:%s\n' % info[3])
    resulthandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ('Ontology','Type','Mode','Fmax','Threshold','Coverage'))
    for onto in ['bpo','cco','mfo']:
        path = os.path.splitext(pred_path.name)[0]+'_'+onto.upper()+'.txt'
        print('ontology: %s\n' % onto)
        for Type in ['type1','type2']:
            print('benchmark type:%s\n' % typeConverter(Type))
            b,obocountDict = read_benchmark(onto, taxon_name_converter(taxon),Type,benchmarkFolder,obo_path)
            if b==None:
                sys.stderr.write('No benchmark is available for the input species and type')
            c = PrecREC(b,path,obocountDict[onto])
            if c.exist:
                for mode in ['partial', 'full']:
                    print('mode:%s\n' % mode)
                    fm = c.Fmax_output(mode)
                    precision = fm[0]
                    recall = fm[1]
                    opt = fm[2]
                    thres = fm[3]
                    coverage = c.coverage()
                    #fm.append(os.path.splitext(os.path.basename(pred_path.name))[0])
                    #print(fm)
                    print('fmax: %s\n' % opt)
                    print('threshold giving fmax: %s\n' % thres)
                    print('coverage: %s\n' % coverage)
                    resulthandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % (onto,typeConverter(Type),mode,opt,thres,coverage))
                    prhandle.write('>%s\t%s\t%s\n |' %  (onto,typeConverter(Type),mode))
                    prhandle.write(" ".join([str(i) for i in precision]) +'\n')
                    prhandle.write(" ".join([str(i) for i in recall]) +'\n') 
            del c
            gc.collect()
    resulthandle.close()
    prhandle.close()
