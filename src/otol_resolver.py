# %% Load packages

import pandas as pd
import numpy as np
import zipfile
import glob
import os
import sys
import time
import shlex
import subprocess
from tqdm import tqdm
from tqdm import tqdm_notebook
from opentree import OT
import json
from pandas import json_normalize



# %% defininbg display options

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100)

# We deactivate the iloc warning see https://stackoverflow.com/a/20627316
pd.options.mode.chained_assignment = None  # default='warn'



# %% specifying input path

metadata_table_path = '../data/in/Inmuno_metadata.tsv'



org_column_header = 'Species'#'Sp_filtered'

# timer is started
start_time = time.time()


# %% Resolving the taxon information from the input file

# the metadata table is loaded 

samples_metadata = pd.read_csv(metadata_table_path, sep='\t')


# %% The sp. string are removed as they can be problematic for the otol resolving

samples_metadata[org_column_header] = samples_metadata[org_column_header].str.replace(' sp.' , '')


# Now we want to get the taxonomic information for each of the samples
# so we want to extract the species information from the metadata file
# %%
samples_metadata[org_column_header].dropna(inplace = True)
samples_metadata[org_column_header] = samples_metadata[org_column_header].str.lower()
species = samples_metadata[org_column_header].unique()
len_species = len(species)

print("%s unique species have been selected from the metadata table." % len_species )
# %%

species_tnrs_matched = OT.tnrs_match(species, context_name=None, do_approximate_matching=True, include_suppressed=False)



# %%

with open('species.json', 'w') as out:
    sf = json.dumps(species_tnrs_matched.response_dict, indent=2, sort_keys=True)
    out.write('{}\n'.format(sf))

# %%
with open("species.json") as tmpfile:
        jsondic = json.loads(tmpfile.read())

json_normalize(jsondic)
# %%

df_species_tnrs_matched = json_normalize(jsondic,
               record_path=['results', 'matches']
               )



df_species_tnrs_unmatched = json_normalize(jsondic,
               record_path=['unmatched_names']
               )
# %%

df_species_tnrs_matched.info()


# %%
len(df_species_tnrs_matched['taxon.ott_id'].unique())
# %%


# We then want to match with the accepted name instead of the synonym in case both are present. 
# We thus order by matched_name and then by is_synonym status prior to returning the first row.

df_species_tnrs_matched.sort_values(['search_string', 'is_synonym'], axis = 0, inplace = True)
df_species_tnrs_matched_unique = df_species_tnrs_matched.drop_duplicates('search_string', keep = 'first')

# both df are finally merged
merged_df = pd.merge(samples_metadata, df_species_tnrs_matched_unique, how='left', left_on=org_column_header, right_on='search_string', indicator=True)


# %%
# We'll stop here for now because this is what we need at the moment.
# (to check !!)

merged_df.replace([r"^ND", r"^NA", r"^nd", r'^\s*$'], np.nan, regex=True, inplace=True)


merged_df.to_csv('../data/out/Inmuno_metadata_otoled.tsv', index = None, sep = '\t')


# Actually we also want the NCBI IDs for Redu.
# Since these are retrieved at the previous OTOL fetching step we'll use these for now


merged_df["ncbi_id"] = merged_df["taxon.tax_sources"].str[0]


import re 
merged_df['ncbi_id'].replace(re.compile('gbif:.*'), np.nan, inplace=True)


# Here we want to extract the digits of the NCBI code only and format them for the Redu format e.g. 
# 1027072|Banisteriopsis muricata

merged_df['ncbi_id'] = merged_df['ncbi_id'].str.extract('(\d+)')

merged_df['ncbi_id_redu'] = merged_df['ncbi_id'] + '|' + merged_df['taxon.name']

# We now format the table so it is ReDu compatible see requirements here https://docs.google.com/spreadsheets/d/1v71bnUd8fiXX51zuZIUAvYETWmpwFQj-M3mu4CNsHBU/edit#gid=1747621311

# We set the Massive ID to the uploaded dataset here https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b753bf1e39cb4875bdf3b786e747bc15
#merged_df['MassiveID'] = 'MSV000087728'

# filename
# column is already existing

# SampleType

#conditions = [
#    (merged_df['sample_type'] == 'sample'),
#    (merged_df['sample_type'] == 'BK'),
#    (merged_df['sample_type'] == 'QC')]
#choices = ['plant', 'blank_extraction', 'blank_QC']

#merged_df['SampleType'] = np.select(conditions, choices)

# SampleTypeSub1

#merged_df['SampleTypeSub1'] = 'not specified'

# NCBITaxonomy
# We already have this column so we'll rename it in place.
# needd to understand how we deal with the QC and BK entries here ...

merged_df.rename(columns={'ncbi_id_redu': 'NCBITaxonomy'}, inplace=True)


# YearOfAnalysis
#merged_df['YearOfAnalysis'] = '2017'

#merged_df['SampleCollectionMethod'] = 'liquid, solid phase extraction (C18)'
#merged_df['SampleExtractionMethod'] = 'ethyl acetate (100%)'
#merged_df['InternalStandardsUsed'] = 'none'
#merged_df['MassSpectrometer'] = 'Q Exactive|MS:1001911'
#merged_df['IonizationSourceAndPolarity'] = 'electrospray ionization (positive)'
#merged_df['ChromatographyAndPhase'] = 'reverse phase (C18)'

# Now we subset the required columns

#cols_to_keep = ['MassiveID',
#'filename', 'SampleType', 'SampleTypeSub1', 'NCBITaxonomy', 'YearOfAnalysis',
#'SampleCollectionMethod','SampleExtractionMethod','InternalStandardsUsed','MassSpectrometer','IonizationSourceAndPolarity','ChromatographyAndPhase' ]

cols_to_keep = ['Family', 'Genus','Species', 'matched_name', 'taxon.name', 'taxon.ott_id', 'taxon.synonyms', 'NCBITaxonomy' ]

redu_df = merged_df[cols_to_keep]

redu_df.to_csv('../data/out/Inmuno_metadata_redued.tsv', index = None, sep = '\t')


# %%
#Now we want to retrieve the upper taxa lineage for all the samples

# Firsst we retrieve a list of unique ott_ids

# Here when checking the columns datatype we observe that the ott_ids are as float.
# We need to keep them as int
# displaying the datatypes 
display(merged_df.dtypes) 

# converting 'ott_ids' from float to int (check the astype('Int64') whic will work while the astype('int') won't see https://stackoverflow.com/a/54194908)
merged_df['taxon.ott_id'] = merged_df['taxon.ott_id'].astype('Int64')
  

# However, we then need to put them back to 
merged_df['taxon.ott_id']
ott_list = list(merged_df['taxon.ott_id'].dropna().astype('int'))

#ott_list = ott_list[0:10]

# %%

taxon_info = []

for i in ott_list:
    query = OT.taxon_info(i, include_lineage=True)
    taxon_info.append(query)

# %%


tl = []

for i in taxon_info:
    with open('taxon_info.json', 'w') as out:
        tl.append(i.response_dict)
        yo = json.dumps(tl)
        out.write('{}\n'.format(yo))

# %%

with open("taxon_info.json") as tmpfile:
        jsondic = json.loads(tmpfile.read())

df = json_normalize(jsondic)


# %%

df_tax_lineage = json_normalize(jsondic,
               record_path=['lineage'],
               meta = ['ott_id', 'unique_name'],
               record_prefix='sub_',
               errors='ignore'
               )
# %%
# This keeps the last occurence of each ott_id / sub_rank grouping https://stackoverflow.com/a/41886945

df_tax_lineage_filtered = df_tax_lineage.groupby(['ott_id', 'sub_rank'], as_index=False).last()
# %%
#Here we pivot long to wide to get the taxonomy

df_tax_lineage_filtered_flat = df_tax_lineage_filtered.pivot(index='ott_id', columns='sub_rank', values='sub_name')

# %%
# Here we actually also want the lowertaxon (species usually) name

df_tax_lineage_filtered_flat = pd.merge(df_tax_lineage_filtered_flat, df_tax_lineage_filtered[['ott_id', 'unique_name']], how='left', on='ott_id', )

#Despite the left join ott_id are duplicated 

df_tax_lineage_filtered_flat.drop_duplicates(subset = ['ott_id', 'unique_name'], inplace = True)

# %%
# we keep the fields of interest

df_tax_lineage_filtered_flat[['ott_id', 'kingdom', 'phylum',
                              'class', 'order', 'family', 'genus', 'unique_name']]



# We now rename our columns of interest

renaming_dict = {'kingdom': 'query_otol_kingdom',
                 'phylum': 'query_otol_phylum',
                 'class': 'query_otol_class',
                 'order': 'query_otol_order',
                 'family': 'query_otol_family',
                 'genus': 'query_otol_genus',
                 'unique_name': 'query_otol_species'}


df_tax_lineage_filtered_flat.rename(columns=renaming_dict, inplace=True)

# We select columns of interest 

cols_to_keep = ['ott_id',
                'query_otol_kingdom',
                'query_otol_phylum',
                'query_otol_class',
                'query_otol_order',
                'query_otol_family',
                'query_otol_genus',
                'query_otol_species']

df_tax_lineage_filtered_flat = df_tax_lineage_filtered_flat[cols_to_keep]


# We merge this back with the samplemetadata only if we have an ott.id in the merged df 

samples_metadata = pd.merge(merged_df[pd.notnull(merged_df['taxon.ott_id'])], df_tax_lineage_filtered_flat, how='left', left_on='taxon.ott_id', right_on='ott_id' )

cols_to_keep = ['Family', 'Genus','Species', 'Part', 'matched_name', 'taxon.synonyms', 'NCBITaxonomy', 'ott_id', 'query_otol_kingdom', 'query_otol_phylum', 'query_otol_class', 'query_otol_order','query_otol_family', 'query_otol_genus', 'query_otol_species']

samples_metadata = samples_metadata[cols_to_keep]

samples_metadata.to_csv('../data/out/Inmuno_metadata_final.tsv', index = None, sep = '\t')