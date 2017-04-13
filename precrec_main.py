# -*- coding: utf-8 -*-
"""
Command-line interface of the precision-recall calculations in
CAFA assessment tool
@author: Ashley Zhou
last updated 04/06/2017
"""

import argparse
from precrec.precRec import PrecREC,read_benchmark,result
from precrec.GOPred import GOPred
import numpy
import os
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import errno    
import seaborn as sns

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
    taxonTable = {'10116':'RAT','9606':'HUMAN','3702':'ARATH','7955':'DANRE','44689':'DICDI','7227':'DROME','83333':'ECOLI','10090':'MOUSE','208963':'PSEAE','4896':'SCHPO','4932':'YEAST'}
    return taxonTable[taxonID]    


def curveSmooth(result):
    #This function removes a p-r pair if there exists another p-r pair that's greater in both precision and recall
    #precision and recall should both be lists of the same length
    precision = []
    recall = []
    for i in range(len(result.precision)):
        remove = False
        for j in range(i):
            if result.precision[i]<result.precision[j] and result.recall[i]<result.recall[j]:
                remove = True
                break
        if not remove:
            precision.append(result.precision[i])
            recall.append(result.recall[i])
    return([precision,recall])

def plotSingle(result,smooth):
    '''
    an result object
    pay attention to order of precision-recall
    recall is on x-axis, but supplied second
    '''
    if smooth=='Y':
        precision = curveSmooth(result)[0]
        recall = curveSmooth(result)[1]
        ax = plt.subplot()
        ax.plot(recall,precision,'-g',label=result.author+': fmax='+ '%.3f'%result.opt)
        ax.plot(result.recall[int(result.thres*100)],result.precision[int(result.thres*100)],'gD')
    elif smooth=='N':
        ax = plt.subplot()
        ax.plot(result.recall,result.precision,'-g',label=result.author+': fmax='+ '%.3f'%result.opt)
        ax.plot(result.recall[int(result.thres*100)],result.precision[int(result.thres*100)],'gD')
    plt.axis([0,1,0,1])
    plt.yticks(numpy.arange(0,1,0.1))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='best')
    figuretitle = result.author+" "+str(result.model)+" "+str(result.taxon)+ " " + result.ontology+' '+ 'mode:'+result.mode+' '+'type:'+result.TYPE
    plt.title(figuretitle)
    figurename = './plots/'+result.author+"_"+str(result.model)+"_"+str(result.taxon)+ "_" + result.ontology+'_'+ result.mode+'_'+result.TYPE+'.png'
    plt.savefig(figurename,dpi=200)
    plt.close()

def plotMultiple(title,listofResults,smooth):
    '''
    supply lists of precision+recall+name lists
    '''
    num = len(listofResults)
    pal=sns.color_palette("Set2", num)
    colors=pal.as_hex()
    for j,i in enumerate(listofResults):
        if smooth=='Y':
            ax = plt.subplot()
            precision = curveSmooth(i)[0]
            recall = curveSmooth(i)[1]
            ax.plot(recall,precision,'-',color=colors[j],label=i.author+': fmax='+ '%.3f'%i.opt) 
            ax.plot(i.recall[int(i.thres*100)],i.precision[int(i.thres*100)],'o',color=colors[j])
        elif smooth=='N':
            ax = plt.subplot()
            ax.plot(i.recall,i.precision,'-',color=colors[j],label=i.author+': fmax='+ '%.3f'%i.opt )
            ax.plot(i.recall[int(i.thres*100)],i.precision[int(i.thres*100)],'o',color=colors[j])
    plt.axis([0,1,0,1])
    plt.yticks(numpy.arange(0,1,0.1))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='best')
    plt.title(title)
    if title == None:
        figurename = './plots/Combined_plot.png'
    else:
        figurename = './plots/'+title+'.png'
    
    plt.savefig(figurename,dpi=200)
    plt.close()
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
    parser.add_argument('file',type=open,
                        help='Input prediction file. Filename should follow CAFA formats. Accepts more than one predictions.',
                        nargs = '+')
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
    mkdir_p('./plots/')
    mkdir_p('./results/')
    
    num = len(args.file)
    resultBPO = []
    resultCCO = []
    resultMFO = []
    for f in args.file:
        print('Evaluating %s.\n' % f.name)
        resulthandle = open("./results/%s_results.txt" % os.path.basename(f.name),'w')
        #first split the prediction file into three ontologies
        all_pred = GOPred()
        pred_path = f
        obo_path = args.obo_path
        benchmarkFolder = args.bfolder
        all_pred.read_and_split_and_write(obo_path,pred_path)
        info = [all_pred.author,all_pred.model,all_pred.keywords,all_pred.taxon]
        print('AUTHOR: %s\n' % info[0])
        resulthandle.write('AUTHOR:%s\n' % info[0])
        print('MODEL: %s\n' % info[1])
        resulthandle.write('MODEL: %s\n' % info[1])
        print('KEYWORDS: %s\n' % info[2][0])
        resulthandle.write('KEYWORDS: %s\n' % info[2][0])
        print('Species:%s\n' % info[3])
        resulthandle.write('Species:%s\n' % info[3])
        print('benchmark type:%s\n' % typeConverter(args.type)) 
        resulthandle.write('benchmark type:%s\n' % typeConverter(args.type))
        print('mode:%s\n' % args.mode)
        resulthandle.write('mode:%s\n' % args.mode)
        resulthandle.write('%s:\t%s\t%s\t%s\n' % ('Ontology','Fmax','Threshold','Coverage'))
        #parse file name 
        #namefields = os.path.basename(pred_path.name).split('.')[0].split('_')
            
        #read benchmark for three ontologies
        for onto in ['bpo','cco','mfo']:
            res = result()
            res.read_from_GOPred(all_pred)
            res.mode = args.mode
            res.TYPE = typeConverter(args.type)
            print('ontology: %s\n' % onto)
            res.ontology = onto
            b = read_benchmark(onto, taxon_name_converter(res.taxon),args.type,benchmarkFolder,obo_path)
            path = os.path.splitext(pred_path.name)[0]+'_'+onto.upper()+'.txt'
            c = PrecREC(b,path)
            if c.exist:
                fm = c.Fmax_output(args.mode)
                res.precision = fm[0]
                res.recall = fm[1]
                res.opt = fm[2]
                res.thres = fm[3]
                res.coverage = fm[4]
                #fm.append(os.path.splitext(os.path.basename(pred_path.name))[0])
                #print(fm)
                print('fmax: %s\n' % res.opt)
                print('threshold giving fmax: %s\n' % res.thres)
                print('coverage: %s\n' % res.coverage)
                plotSingle(res,args.smooth)
                resulthandle.write('%s:\t%s\t%s\t%s\n' % (onto,res.opt,res.thres,res.coverage))
            if onto=='bpo':
                resultBPO.append(res)
            elif onto=='cco':
                resultCCO.append(res)
            elif onto=='mfo':
                resultMFO.append(res)

        resulthandle.close()
    print(num)
    if num>1:
        plotMultiple(args.title+'_BPO', resultBPO,args.smooth)
        plotMultiple(args.title+'_CCO', resultCCO,args.smooth)
        plotMultiple(args.title+'_MFO', resultMFO,args.smooth)
