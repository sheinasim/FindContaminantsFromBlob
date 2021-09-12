#!/usr/bin/python

import json
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("blobdir", help="The blob directory containing bestsumorder and identifiers files in .json format")
args = parser.parse_args()

with open(args.blobdir + '/identifiers.json', 'r') as f:
        identifiers = json.loads(f.read()) 

with open(args.blobdir + '/bestsumorder_phylum.json', 'r') as f:
        bestsumorder_phylum = json.loads(f.read())

with open(args.blobdir + '/bestsumorder_species.json', 'r') as f:
		bestsumorder_species = json.loads(f.read())

with open(args.blobdir + '/bestsumorder_class.json', 'r') as f:
		bestsumorder_class = json.loads(f.read())

with open(args.blobdir + '/bestsumorder_order.json', 'r') as f:
		bestsumorder_order = json.loads(f.read())

with open(args.blobdir + '/bestsumorder_family.json', 'r') as f:
		bestsumorder_family = json.loads(f.read())

with open(args.blobdir + '/bestsumorder_genus.json', 'r') as f:
		bestsumorder_genus = json.loads(f.read())		

taxa = {"phylum" : bestsumorder_phylum,  
        "class" : bestsumorder_class, 
        "order" : bestsumorder_order, 
        "family" : bestsumorder_family, 
        "genus": bestsumorder_genus,
        "species" : bestsumorder_species}

def findKeepersandContaminants(identifiersFile, bestsumorderFile):
	df_id = pd.json_normalize(identifiersFile, record_path=['values'])
	df_id = df_id.set_axis(["contig"], axis=1)
	df_bsc_values = pd.json_normalize(bestsumorderFile, record_path=['values'])
	df_bsc_values = df_bsc_values.set_axis(["value"], axis=1)
	df_bsc_keys = pd.json_normalize(bestsumorderFile, record_path=['keys'])
	df_bsc_keys = df_bsc_keys.set_axis(["taxa"], axis=1)
	df_bsc_keys['value'] = df_bsc_keys.index
	df_bsc = df_bsc_values.join(df_bsc_keys, lsuffix='_value', rsuffix='_key', on = 'value')['taxa']
	df_id = df_id.join(df_bsc)
	keepers = df_id[df_id['taxa'].isin(['Arthropoda', 'no-hit', 'undef'])]
	contaminants = df_id[~df_id['taxa'].isin(['Arthropoda', 'no-hit', 'undef'])]
	return contaminants, keepers

def assignContaminantstoSpecies(identifiersFile, taxa_dict):
    df_id = pd.json_normalize(identifiersFile, record_path=['values'])
    df_id = df_id.set_axis(["contig"], axis=1)
    for taxon, file in taxa_dict.items():
        for key in file:
            values = pd.json_normalize(file, record_path=['values'])
            values = values.set_axis(['value'], axis=1)
            keys = pd.json_normalize(file, record_path=['keys'])
            keys = keys.set_axis([taxon], axis=1)
            keys['value'] = keys.index
            df_vk = values.join(keys, lsuffix='_value', rsuffix='_key', on = 'value')[taxon]
        df_id = df_id.join(df_vk)
    return df_id

speciesTbl = assignContaminantstoSpecies(identifiers, taxa) 

prefix = args.blobdir.split('/')[-1]

contaminants_outfile = prefix +'_contaminants.txt'
keepers_outfile = prefix +'_keepers.txt'
species_outfile = prefix + '_taxa.tsv'

contaminantTigs, keeperTigs = findKeepersandContaminants(identifiers, bestsumorder_phylum)

contaminantTigs['contig'].to_csv(contaminants_outfile, sep='\t', index=False, header=False)
keeperTigs['contig'].to_csv(keepers_outfile, sep='\t', index=False, header=False)
speciesTbl.to_csv(species_outfile, sep='\t', index=False, header=True)
