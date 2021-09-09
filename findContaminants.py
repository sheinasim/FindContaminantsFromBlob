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

def findKeepersandContaminants(identifiersFile, bestsumorderFile):
	df_id = pd.json_normalize(identifiersFile, record_path=['values'])
	df_id = df_id.set_axis(["contig"], axis=1)
	df_bsc_values = pd.json_normalize(bestsumorderFile, record_path=['values'])
	df_bsc_values = df_bsc_values.set_axis(["value"], axis=1)
	df_bsc_keys = pd.json_normalize(bestsumorderFile, record_path=['keys'])
	df_bsc_keys = df_bsc_keys.set_axis(["taxa"], axis=1)
	df_bsc_keys['value'] = df_bsc_keys.index
	df_bsc = pd.merge(df_bsc_values, df_bsc_keys, on ='value')
	df_id = df_id.join(df_bsc)
	keepers = df_id[df_id['taxa'].isin(['Arthropoda', 'no-hit', 'undef'])]
	contaminants = df_id[~df_id['taxa'].isin(['Arthropoda', 'no-hit', 'undef'])]
	return contaminants, keepers

def assignContaminantstoSpecies(identifiersFile, bestsumorderFile):
	df_id = pd.json_normalize(identifiersFile, record_path=['values'])
	df_id = df_id.set_axis(["contig"], axis=1)
	df_bsc_values = pd.json_normalize(bestsumorderFile, record_path=['values'])
	df_bsc_values = df_bsc_values.set_axis(["value"], axis=1)
	df_bsc_keys = pd.json_normalize(bestsumorderFile, record_path=['keys'])
	df_bsc_keys = df_bsc_keys.set_axis(["taxa"], axis=1)
	df_bsc_keys['value'] = df_bsc_keys.index
	df_bsc = pd.merge(df_bsc_values, df_bsc_keys, on ='value')
	df_id = df_id.join(df_bsc)
	return df_id

speciesTbl = assignContaminantstoSpecies(identifiers, bestsumorder_species)

prefix = args.blobdir.split('/')[-1]

contaminants_outfile = prefix +'_contaminants.txt'
keepers_outfile = prefix +'_keepers.txt'
species_outfile = prefix + '_species.tsv'

contaminantTigs, keeperTigs = findKeepersandContaminants(identifiers, bestsumorder_phylum)

contaminantTigs['contig'].to_csv(contaminants_outfile, sep='\t', index=False, header=False)
keeperTigs['contig'].to_csv(keepers_outfile, sep='\t', index=False, header=False)
speciesTbl[['contig', 'taxa']].to_csv(species_outfile, sep='\t', index=False, header=True)
