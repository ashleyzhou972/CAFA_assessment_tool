#!/home/nzhou/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:32:54 2017
This script provides top ten labs from the rankings and output their respective prediction files
@author: nzhou
"""

import os
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy
import sys
from matplotlib.font_manager import FontProperties
import yaml
import argparse
import errno



def typeConverter(oldType):
    if oldType=='type1':
        newType = 'NK'
    elif oldType == 'type2':
        newType = 'LK'
    elif oldType == 'all':
        newType = 'All'
    return(newType)


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
    def read_info(self,onto,Type,mode,method):
        fields = method.split('_')
        self.author=fields[0]
        self.model=int(fields[1])
        self.taxon=fields[2]
        self.ontology=onto
        self.mode=mode
        self.TYPE=Type
        self.method = method
    def read_prrc(self,pr,rc):
        self.precision=pr
        self.recall=rc
    def calculate_fmax(self):
        fmax = 0
        f_thres = 0.00
        for i in range(101):
            try:
                a = float(self.precision[i])
                b = float(self.recall[i])
            except IndexError:
                print('cutoff')
                break
            try:
                f = 2*a*b/(a+b)
            except ZeroDivisionError:
                f = None
            if f!=None and f>=fmax:
                fmax = f
                f_thres = numpy.around(i*0.01,decimals=2)
        self.opt = numpy.around(fmax, decimals=5)
        self.thres = f_thres
        print('fmax is %s\n' % fmax)
        print('thres is %s\n' % f_thres)
        
    def check_fmax(self,onto,Type,mode,method, results_folder):
        #This functino checks the fmax calculated from the saved pr_rc values
        #vs the onces calculated in results
        result_file = os.path.join(results_folder,'%s_results.txt' % method)
        with open(result_file) as f:
            for line in f:
                if line.startswith(onto):
                    if line.split('|')[0].split()[1]==Type and line.split('|')[0].split()[2]==mode:
                        fmax = line.split('|')[1].split()[0]
                        print(fmax)
                        nfmax = numpy.around(float(fmax), decimals=5)
                        thres = line.split('|')[1].split()[1]
                        break
        if '%.5f' % nfmax == '%.5f' % self.opt  and str(thres)==str(self.thres):
            return(True)
        else:
            print('calculated fmax: %s, result fmax: %s\n' % (self.opt, nfmax))
            print('calculated thres: %s, result thres: %s\n' % (self.thres, thres))
            return(False)           
            
    def getCoverage(self,onto,Type,mode,method,results_folder):
        result_file = os.path.join(results_folder,'%s_results.txt' % method)
        with open(result_file) as f:
            for line in f:
                if line.startswith(onto):
                    if line.split('|')[0].split()[1]==Type and line.split('|')[0].split()[2]==mode:
                        coverage = line.split('|')[1].split()[2]
                        self.coverage=str(numpy.around(float(coverage),decimals=2))
                        break            
        
        
        
def getprrc(onto,Type,mode,prrcfolder,method, results_folder):
    #Type need to be 'NK' and 'LK'
    filename = os.path.join(prrcfolder,'%s_prrc.txt' % (method))
    r = result()
    r.read_info(onto,Type,mode,method)
    with open(filename,'r') as f:
        read = False
        for line in f:
            if line.startswith('>'):
                fields = line.strip()[1:].split('\t')
                if fields[0]==onto and fields[1]==Type and fields[2]==mode:
                    read=True
                    #print(line)
                    pr = f.next().split('|')[1].strip().split()
                   # print(pr)
                    rc = f.next().strip().split()
                    #print(rc)
                    r.read_prrc(pr,rc)
                    r.calculate_fmax()
                    if r.check_fmax(onto,Type,mode,method, results_folder):
                        r.getCoverage(onto,Type,mode,method, results_folder)                
                    else:
                        print("check results\n")
                    break
                #both pr and rc are lists
        if not read:
            print('result not found for %s %s %s %s' % (method, onto, Type, mode))
    return(r)
    


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
    
    
def plotMultiple(title,listofResults,smooth):
    '''
    supply lists of precision+recall+name lists
    '''
    fontP = FontProperties()
    fontP.set_size('small')
    num = len(listofResults)
    pal=sns.color_palette("Paired", num)
    colors=pal.as_hex()
    for j,i in enumerate(listofResults):
        linetype = '-'
        if smooth=='Y':
            ax = plt.subplot()
            precision = curveSmooth(i)[0][1:]
            recall = curveSmooth(i)[1][1:]
            ax.plot(recall,precision,linetype,color=colors[j],label=i.method+':\nF=%s C=%s'%(i.opt,i.coverage)) 
            ax.plot(i.recall[int(i.thres*100)],i.precision[int(i.thres*100)],'o',color=colors[j])
        elif smooth=='N':
            ax = plt.subplot()
            ax.plot(i.recall,i.precision,linetype,color=colors[j],label=i.method+':\nF=%s C=%s'%(i.opt,i.coverage))
            ax.plot(i.recall[int(i.thres*100)],i.precision[int(i.thres*100)],'o',color=colors[j])
    plt.axis([0,1,0,1])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.yticks(numpy.arange(0,1,0.1))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(title)
    figurename = os.path.join('./plots/',title)       
    plt.savefig(figurename,dpi=200)
    plt.close()
        
        
        
def extant_file(x):
    if not os.path.isfile(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    else:
        return(open(x,'r')) 
        
        
        
def read_config():
    parser = argparse.ArgumentParser(description='Precision-Recall Curves plot', )
    

    parser.add_argument('config_stream',type=extant_file, help='Configuration file')
    #CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    #example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    #If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)                  
    args = parser.parse_args()
    try:
        config_dict = yaml.load(args.config_stream)['plot']
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()
    Num_files = len(config_dict)-3
    results_folder = config_dict['results']
    title = config_dict['title']
    smooth = config_dict['smooth']
    methods = set()
    for i in xrange(Num_files):
        keyname = 'file'+str(i+1)
        methods.add(config_dict[keyname])
    return(results_folder, title,smooth, methods)   


def check_existence(results_folder, methods):
    re = True
    if os.path.isdir(results_folder):
        for method in methods:
            file_result = os.path.join(results_folder,'%s_results.txt' % method)
            if not os.path.isfile(file_result):
                sys.stderr.write('file %s not found\n' % file_result)
                re = False
                break
        prrc_folder = results_folder + '/pr_rc/'
        if os.path.isdir(prrc_folder):
            #check if files exist
            for method in methods:
                file_prrc = os.path.join(prrc_folder,'%s_prrc.txt' % method)
                if not os.path.isfile(file_prrc):
                    re = False
                    sys.stderr.write('file %s not found\n' % file_prrc)
                    break
        else:
            sys.stderr.write('directory %s not found\n' % prrc_folder)
            re = False
    else:
        sys.stderr.write('directory %s not found\n' % results_folder)
        re = False
    return(re)
    
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
     
if __name__=='__main__':
    results_folder, title, smooth, methods = read_config()
    if check_existence(results_folder, methods):
        prrcfolder = results_folder+'/pr_rc/'
        for onto in ['bpo','cco','mfo']:
            for Type in ['LK','NK']:
                for mode in ['partial','full']:
                    specific_title = '%s_%s_%s_%s_fmax.png' % (title,onto,Type,mode)
                    print('\nPlotting %s\n' % title)
                    mkdir_p('./plots')
                    result_list = []
                    for method in methods:
                        res = getprrc(onto, Type,mode, prrcfolder,method, results_folder)
                        result_list.append(res)
                    plotMultiple(specific_title,result_list,smooth)
    else:
        sys.exit()
 
