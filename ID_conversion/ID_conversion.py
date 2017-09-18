#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Not available for 237561, 208963 and 7227

import os
import sys

#below are all 23 cafa3 species
taxons = ['85962',
 '83333',
 '10116',
 '7955',
 '273057',
 '7227',
 '160488',
 '170187',
 '208963',
 '224308',
 '243273',
 '9606',
 '243232',
 '44689',
 '559292',
 '284812',
 '3702',
 '8355',
 '10090',
 '223283',
 '99287',
 '237561',
 '321314']

#ID mapping is not available for 237561, 208963 and 7227, because the inital CAFA3 target sequences were not based on UniProt
reduced_taxons = ['85962',
 '83333',
 '10116',
 '7955',
 '273057',
 '160488',
 '170187',
 '224308',
 '243273',
 '9606',
 '243232',
 '44689',
 '559292',
 '284812',
 '3702',
 '8355',
 '10090',
 '223283',
 '99287',
 '321314']
 
def __read_target_mapping__(taxon, targetfolder):
    #target dict maps between cafaid and uniprot gene name
    targetdict = dict()
    if taxon in reduced_taxons:
        if taxon in ['208963','7227','237561']:
            filename = 'mapping.'+str(taxon)+'.map'
        else:
            filename = 'sp_species.'+str(taxon)+'.map'
        handle = open(targetfolder+filename,'r')
        for line in handle:
            fields = line.strip().split('\t')
            name = fields[1]
            cafaid = fields[0]
            targetdict[cafaid]=name
        handle.close()
    else:
        print('%s is not a CAFA3 species' % taxon)
    return(targetdict)

def __uniprot_mapping__(taxon, uniprot_ac_to_id_folder):
    #convert between uniprot gene name and uniprot accession
    #ac to id files for all CAFA3 species are available
    #a dictionary is created
    folder = uniprot_ac_to_id_folder
    uniprotdict = dict()
    mapping = False
    if taxon in reduced_taxons:
        if str(taxon) in ['10090','10116','284812','3702','44689','559292','7227','7955','83333','9606']:
            filename = 'uniprot_ac_to_id_'+taxon+'.map'
            mapping = True
        else:
            filename = 'uniprot_ac_to_id_'+taxon+'.tab'
        
        with open(os.path.join(folder,filename),'r') as f:
            f.readline()
            for line in f:
                if mapping:
                    name = line.strip().split()[2]
                else:
                    name = line.strip().split()[1]
                accession= line.strip().split()[0]
                if name not in uniprotdict.keys():            
                    uniprotdict[name] = accession
                else:
                    print("Repeated uniprot gene name %s\t" % line)
    else:
        print('%s is not a CAFA3 species' % taxon)
    return(uniprotdict)
    
 

def cafaid_to_uniprot(taxon, cafaids):
    #cafaids should be a set/list of CAFA IDs
    #return a dictionary with the cafaids as keys and uniprot ac as values
    targetfolder = './CAFA_mapping/'
    targetdict = __read_target_mapping__(taxon, targetfolder)
    uniprotfolder = './uniprot_mapping/'
    uniprotdict = __uniprot_mapping__(taxon, uniprotfolder)   
    uniprotids_dict = dict()
    for cafaid in cafaids:
        uniprotac = uniprotdict[targetdict[cafaid]] 
        uniprotids_dict[cafaid] = uniprotac
    return(uniprotids_dict)


def uniprotac_to_cafaid(taxon, uniprotacs):
    #uniprotacs should be a set/list of UniProt Accession IDs
    #returns a dictionary
    targetfolder = './CAFA_mapping/'
    targetdict = __read_target_mapping__(taxon, targetfolder)
    uniprotfolder = './uniprot_mapping/'
    uniprotdict = __uniprot_mapping__(taxon, uniprotfolder)   
    targetdict_reverse = {v: k for k, v in targetdict.items()}
    uniprotdict_reverse = {v: k for k, v in uniprotdict.items()}
    cafaids_dict = dict()
    for uniprotac in uniprotacs:
        try:
            cafaid = targetdict_reverse[uniprotdict_reverse[uniprotac]]
            cafaids_dict[uniprotac] = cafaid
        except KeyError:
            sys.stderr.write('%s  not in CAFA3 target.\n' % (uniprotac))
    return(cafaids_dict)


if __name__=='__main__':
    accessions = set()
    file_with_uniprot_ac_list = sys.argv[1]
    taxon = sys.argv[2]
    outputfile = sys.argv[3]
    with open(file_with_uniprot_ac_list,'r') as f:
        for line in f:
            ac = line.strip()
            accessions.add(ac)
    cafaiddict = uniprotac_to_cafaid(taxon,accessions)
    with open(outputfile,'w') as w:
        for key in cafaiddict:
            w.write('%s\t%s\n' % (key,cafaiddict[key]))
