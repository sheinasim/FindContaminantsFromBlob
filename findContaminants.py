#!/usr/bin/python

import json
import pandas as pd
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("blobdir", help="The blob directory containing bestsumorder and identifiers files in .json format")
parser.add_argument("fasta", help=".fasta BLASTed and Diamond BLASTed")
args = parser.parse_args()

df_contig_lengths = pd.DataFrame(columns=('contig', 'length'))
for record in SeqIO.parse(args.fasta, "fasta"):
	df_contig_lengths = df_contig_lengths.append({'contig' : record.name, 'length' : len(record.seq)}, ignore_index = True)

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

def assignContaminantstoSpecies(identifiersFile, taxa_dict, ctgLengths):
    df_id = ctgLengths
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

speciesTbl = assignContaminantstoSpecies(identifiers, taxa, df_contig_lengths) 
keeperTigs = speciesTbl[speciesTbl['phylum'].isin(['Arthropoda', 'no-hit', 'undef'])]
contaminantTigs = speciesTbl[~speciesTbl['phylum'].isin(['Arthropoda', 'no-hit', 'undef'])]


prefix = args.blobdir.split('/')[-1]

contaminants_outfile = prefix +'_contaminants.tsv'
keepers_outfile = prefix +'_keepers.tsv'
species_outfile = prefix + '_taxa.tsv'

contaminantTigs[['contig', 'length']].to_csv(contaminants_outfile, sep='\t', index=False, header=False)
keeperTigs[['contig', 'length']].to_csv(keepers_outfile, sep='\t', index=False, header=False)
speciesTbl.to_csv(species_outfile, sep='\t', index=False, header=True)

print("Find Contaminants from BLOBs is finished.")
