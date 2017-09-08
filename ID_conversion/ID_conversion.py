# -*- coding: utf-8 -*-

#Not available for 208963 and 7227

import os

def read_target_mapping(taxon, targetfolder):
    #target dict maps between cafaid and uniprot gene name
    targetdict = dict()
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
    return(targetdict)

def uniprot_mapping(taxon, uniprot_ac_to_id_folder):
    #convert between uniprot gene name and uniprot accession
    #ac to id files for all CAFA3 species are available
    #a dictionary is created
    folder = uniprot_ac_to_id_folder
    if str(taxon) in ['10090','10116','284812','3702','44689','559292','7227','7955','83333','9606']:
        filename = 'uniprot_ac_to_id_'+taxon+'.map'
    else:
        filename = 'uniprot_ac_to_id_'+taxon+'.tab'
    uniprotdict = dict()
    with open(os.path.join(folder,filename),'r') as f:
        f.readline()
        for line in f:
            name = line.strip().split()[2]
            accession= line.strip().split()[0]
            if name not in uniprotdict.keys():            
                uniprotdict[name] = accession
            else:
                print("Repeated uniprot gene name %s\t" % line)
    return(uniprotdict)
    
 

def cafaid_to_uniprot(taxon, cafaids):
    #cafaids should be a set/list of CAFA IDs
    #return a dictionary with the cafaids as keys and uniprot ac as values
    targetfolder = './CAFA_mapping/'
    targetdict = read_target_mapping(taxon, targetfolder)
    uniprotfolder = './uniprot_mapping/'
    uniprotdict = uniprot_mapping(taxon, uniprotfolder)   
    uniprotids_dict = dict()
    for cafaid in cafaids:
        uniprotac = uniprotdict[targetdict[cafaid]] 
        uniprotids_dict[cafaid] = uniprotac
    return(uniprotids_dict)


def uniprotac_to_cafaid(taxon, uniprotacs):
    #uniprotacs should be a set/list of UniProt Accession IDs
    #returns a dictionary
    targetfolder = './CAFA_mapping/'
    targetdict = read_target_mapping(taxon, targetfolder)
    uniprotfolder = './uniprot_mapping/'
    uniprotdict = uniprot_mapping(taxon, uniprotfolder)   
    targetdict_reverse = {v: k for k, v in targetdict.items()}
    uniprotdict_reverse = {v: k for k, v in uniprotdict.items()}
    cafaids_dict = dict()
    for uniprotac in uniprotacs:
        cafaid = targetdict_reverse[uniprotdict_reverse[uniprotac]]
        cafaids_dict[uniprotac] = cafaid
    return(cafaids_dict)

