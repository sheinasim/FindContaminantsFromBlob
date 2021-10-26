#!/usr/bin/python

import json
import pandas as pd
import argparse
import numpy as np, scipy.stats as st
import os
import re
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("blobdir", help="The blob directory containing bestsumorder and identifiers files in .json format")
parser.add_argument("fasta", help=".fasta BLASTed and Diamond BLASTed")
parser.add_argument("buscoTable", help="The full table output of BUSCO (usually called: full_table.tsv for BUSCO5 or full_table_{outprefix}.tsv for BUSCO3)")
args = parser.parse_args()

df_buscos = pd.read_csv(args.buscoTable, sep='\t', comment='#')
df_buscos.columns = ["busco id", "status", "record", "start", "end", "strand", "score", "Length", "orthoDB url", "description"]
df_buscos = df_buscos[["busco id", "status", "record"]]

df_record_lengths = pd.DataFrame(columns=('record', 'length'))
for record in SeqIO.parse(args.fasta, "fasta"):
	df_record_lengths = df_record_lengths.append({'record' : record.name, 'length' : len(record.seq)}, ignore_index = True)

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

path = args.blobdir

for filename in os.listdir(path):
    if re.match(r".*reads_cov\.json", filename):
         with open(os.path.join(path, filename), 'r') as f:
                readCov = json.loads(f.read())
                
coverage = pd.json_normalize(readCov, record_path=['values'])

def assignContaminantstoSpecies(taxa_dict, ctgLengths, cov):
    df_id = ctgLengths
    df_id = pd.merge(left = df_id, right = cov, left_index=True, right_index=True)
    df_id = df_id.rename(columns={0: "coverage"})
    for taxon, file in taxa_dict.items():
        for key in file:
            values = pd.json_normalize(file, record_path=['values'])
            values = values.set_axis(['value'], axis=1)
            keys = pd.json_normalize(file, record_path=['keys'])
            keys = keys.set_axis([taxon], axis=1)
            keys['value'] = keys.index
            df_vk = values.join(keys, lsuffix='_value', rsuffix='_key', on = 'value')[taxon]
        df_id = df_id.join(df_vk)
    low,high = st.t.interval(0.95, len(df_id["coverage"])-1, loc=np.mean(df_id["coverage"]), scale=st.sem(df_id["coverage"]))
    lessthanhalflow = low/2
    morethantwicehigh = high*2
    conditions = [
        (df_id['coverage'] <= lessthanhalflow),
        (df_id['coverage'] > lessthanhalflow) & (df_id['coverage'] <= morethantwicehigh),
        (df_id['coverage'] > morethantwicehigh)
        ]
    values = ['low', 'medium', 'high']
    df_id['coverage class'] = np.select(conditions, values)
    return df_id

def addBUSCOs(df, buscos):
    df_complete = buscos[buscos["status"] == "Complete"]
    numbercomplete = df_complete["record"].value_counts().to_frame("complete BUSCOs")
    numbercomplete["record"] = numbercomplete.index
    df_duplicate = buscos[buscos["status"] == "Duplicated"]
    numberduplicated = df_duplicate["record"].value_counts().to_frame("duplicated BUSCOs")
    numberduplicated["record"] = numberduplicated.index
    df_busco_cols = pd.merge(left = df, right = numbercomplete, how = "left", left_on = "record", right_on = "record")
    df_busco_cols = pd.merge(left = df_busco_cols, right = numberduplicated, how = "left", left_on = "record", right_on = "record")
    df_busco_cols["complete BUSCOs"] = df_busco_cols["complete BUSCOs"].replace(np.nan, 0)
    df_busco_cols["complete BUSCOs"] = df_busco_cols["complete BUSCOs"].astype(int)
    df_busco_cols["duplicated BUSCOs"] = df_busco_cols["duplicated BUSCOs"].replace(np.nan, 0)
    df_busco_cols["duplicated BUSCOs"] = df_busco_cols["duplicated BUSCOs"].astype(int)
    return df_busco_cols

out1 = assignContaminantstoSpecies(taxa, df_record_lengths, coverage)
out2 = addBUSCOs(out1, df_buscos)

speciesTbl = out2[["record", "length", "coverage", "coverage class", "complete BUSCOs", "duplicated BUSCOs", "phylum", "class", "order", "family", "genus", "species"]]
keeperTigs = speciesTbl[speciesTbl['phylum'].isin(['Arthropoda', 'no-hit', 'undef'])]
contaminantTigs = speciesTbl[~speciesTbl['phylum'].isin(['Arthropoda', 'no-hit', 'undef'])]

prefix = args.blobdir.split('/')[-1]

contaminants_outfile = prefix +'_contaminants.tsv'
keepers_outfile = prefix +'_keepers.tsv'
species_outfile = prefix + '_all.tsv'

contaminantTigs.to_csv(contaminants_outfile, sep='\t', index=False, header=True)
keeperTigs.to_csv(keepers_outfile, sep='\t', index=False, header=True)
speciesTbl.to_csv(species_outfile, sep='\t', index=False, header=True)

print("Find Contaminants from BLOBs is finished.")
