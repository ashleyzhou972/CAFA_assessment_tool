# -*- coding: utf-8 -*-
"""
Command-line interface of the precision-recall calculations in
CAFA assessment tool
This is a simple version:
- Supplied prediction file must be ontology specific: FriedbergLab_1_9606_BPO.txt
- Can only supply one file at a time
- No plots 
@author: Ashley Zhou
last updated 04/06/2017
"""

import argparse
from precrec.precRec import PrecREC,read_benchmark,result
from precrec.GOPred import GOPred
import os
import sys
import errno    

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


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



if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
    
    def extant_file(x):
        if not os.path.isfile(x):
            raise argparse.ArgumentTypeError("{0} does not exist".format(x))
        else:
            return(open(x,'r'))
        
    parser.add_argument('file',type=extant_file,
                        help='Input prediction file. Filename should follow CAFA formats. Accepts more than one predictions.',
                        nargs = '+',)
    #CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    #example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    #If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)                  
    

    
    parser.add_argument('-t','--t',dest='type',help = 'Input evaluation type: No Knowledge or Limited Knowledge', choices=['type1','type2','all'],required=True)    
    parser.add_argument('-o','--o', dest= 'obo_path',help = 'Input the obo file path',default = './precrec/go_20130615-termdb.obo')
    parser.add_argument('-m','--m',dest='mode', help = 'Input the evaluation mode: full or partial', choices = ['full','partial'],required = True)
    parser.add_argument('-b','--b',dest='bfolder', help = 'Input the path to the benchmark folder, default CAFA2 benchmarks provided', default = './precrec/benchmark/')
    parser.add_argument('-title', dest='title',help = 'Input title of combined plot, if multiple prediction files are supplied',default = ' ')
    parser.add_argument('-s', dest = 'smooth', help='Option to have the P-R curves smoothed. Enter "Y" or "N". Default is "N". Recommended if plotting multiple curves', default = 'N', choices=['Y','N'], action='store')
    args = parser.parse_args()
    mkdir_p('./results/')
    
    num = len(args.file)
    #print('Number of predictions supplied: %s\n' % num)
    f=args.file
    print('Evaluating %s.\n' % f.name)
    resulthandle = open("./results/%s_%s_%s_results.txt" % (os.path.basename(f.name).split('.')[0],args.mode,args.type),'w')
    #first split the prediction file into three ontologies
    #all_pred = GOPred()
    pred_path = f
    obo_path = args.obo_path
    benchmarkFolder = args.bfolder
    #all_pred.read_and_split_and_write(obo_path,pred_path)
    resulthandle.write('%s:\t%s\t%s\t%s\n' % ('Ontology','Fmax','Threshold','Coverage'))
    #parse file name 
    #namefields = os.path.basename(pred_path.name).split('.')[0].split('_')
        
    onto = os.path.splitext(pred_path.name)[0].split('_')[3]
    taxon = os.path.splitext(pred_path.name)[0].split('_')[2]
    #res.read_from_GOPred(all_pred)
    print('ontology: %s\n' % onto)
    b,obocountDict = read_benchmark(onto, taxon_name_converter(taxon),args.type,benchmarkFolder,obo_path)
    if b==None:
        sys.stderr.write('No benchmark is available for the input species and type')
    path = os.path.splitext(pred_path.name)[0]+'_'+onto.upper()+'.txt'
    c = PrecREC(b,path,obocountDict[onto])
    if c.exist:
        fm = c.Fmax_output(args.mode)
        precision=fm[0]
        print(precision)
        recall=fm[1]
        print(recall)
        opt = fm[2]
        thres = fm[3]
        coverage = fm[4]
        #fm.append(os.path.splitext(os.path.basename(pred_path.name))[0])
        #print(fm)
        print('fmax: %s\n' % opt)
        print('threshold giving fmax: %s\n' % thres)
        print('coverage: %s\n' % coverage)
        resulthandle.write('%s:\t%s\t%s\t%s\n' % (onto,opt,thres,coverage))
        resulthandle.write('%s:\t%s\n') % (onto, precision)
        resulthandle.write('%s:\t%s\n') % (onto, recall)
        resulthandle.close()
