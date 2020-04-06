#!/usr/bin/env python
# coding: utf-8

# ### Analysis of Plasma proteomes based on MS1 and MS2 information
#
# This jupyter notebook allows to analyze proteome discover output using exiting data from identified plasma proteins to re-score.

# In[167]:
import math

import pandas as pd
import seaborn as sns
import re
import numpy as np
import matplotlib.pyplot as plt

# Load the public plasma proteome:
#
# **Gen**: Gen accession number
# **Description**: Protein name
# **Log_Conc**: Concentration Log2-based
# **Zscore**: Zscore per proteins
# **Pvlue**: P value per proteins based on Zscore (probability in a normal distribution)

# In[168]:


public_plasma_proteome = pd.read_csv("data/public-plasma-proteome.csv")
public_plasma_proteome.head()

# In[169]:


sns.distplot(public_plasma_proteome.Log_Conc)

# Load internal database of plasma proteins:
#
# **Uniprot_Accession**: Uniprot accession
# **Pvalue**: P value based in the Zscore

# In[170]:


inhouse_plasma_proteome = pd.read_csv("data/custom-plasma-proteome.csv")
inhouse_plasma_proteome.head()

# In[171]:


dataset = pd.read_csv("data/sample-example.csv")
dataset.head()

# In[172]:


dataset = pd.merge(dataset, inhouse_plasma_proteome, on="UniprotAccession", how='outer')
dataset = pd.merge(dataset, public_plasma_proteome, on="Gene", how='outer')
dataset.head()


# In[173]:


def compute_probability(row, sample):
    if np.isnan(row['Pvalue']):
        pvalue = 0
    else:
        pvalue = row['Pvalue']

    if np.isnan(row['PPvalue']):
        ppvalue = 0
    else:
        ppvalue = row['PPvalue']

    abundance = pvalue + ppvalue - (pvalue * ppvalue)
    if pvalue == 0 and ppvalue == 0:
        abundance = 0.009

    if row['Found[' + sample + ']'] == 'High':
        abundance = 1 * row['Abundance[' + sample + ']']
    if row['Found[' + sample + ']'] == 'Peak Found':
        abundance = abundance * row['Abundance[' + sample + ']']
    if row['Globulins'] == 1:
        abundance = 1 * row['Abundance[' + sample + ']']

    return abundance


def is_decoy(row):
    decoy = 0
    if (row['UniprotAccession'] == "Q14573"):
        print("protein")
    if np.isnan(row['Pvalue']) and np.isnan(row['PPvalue']):
        decoy = 1
    if row['Found[' + sample + ']'] == 'High' or row['Globulins'] == 1:
        decoy = 0
    if row['Found[' + sample + ']'] == 'Not Found':
        decoy = math.nan
    return decoy

samples = []
for column in list(dataset.columns.values):
    m = re.search(r"\[([A-Za-z0-9_]+)\]", column)
    if m is not None and (len(m.group(1))):
        samples.append(m.group(1))

samples = np.unique(np.array(samples))
for sample in samples:
    dataset['AbundanceRecall[' + sample + ']'] = dataset.apply(lambda row: compute_probability(row, sample), axis=1)
    dataset['Decoy[' + sample + ']'] = dataset.apply(lambda row: is_decoy(row), axis=1)
dataset.head()

# In[174]:


for sample in samples:
    dataset = dataset.sort_values(by=['AbundanceRecall[' + sample + ']'], ascending=False)
    fdr_values = []
    decoy_count = 0
    target_count = 0
    for index, row in dataset.iterrows():
        if row['Decoy[' + sample + ']'] == 1:
            decoy_count += 1
        elif row['Decoy[' + sample + ']'] == 0:
            target_count += 1
        fdr = float(decoy_count) / float(decoy_count + target_count)
        if math.isnan(row['Decoy[' + sample + ']']):
            fdr = math.nan
        fdr_values.append(fdr)
    dataset['FDR[' + sample + ']'] = fdr_values
dataset.head(10)


# ## Distribution of targets vs decoy proteins
#
# The following plots show the distribution targets (proteins in plasma databases) vs decoys (proteins not present in plasma databases).

# In[178]:


def plot_targets_vs_decoy(fdr_thershold: float = 1.0, bins=20, hist=True):
    fig, axs = plt.subplots(int(len(samples) / 2), 2,
                            figsize=(15, 35))  # adjust the geometry based on your number of columns to plot
    for ax, sample in zip(axs.flatten(), samples):
        filtered_dataset = dataset.loc[(dataset['FDR[' + sample + ']'] < fdr_thershold)]
        filtered_dataset = filtered_dataset.dropna(subset=['FDR[' + sample + ']'])
        targets = filtered_dataset.loc[filtered_dataset['Decoy[' + sample + ']'] == 0]
        decoys = filtered_dataset.loc[filtered_dataset['Decoy[' + sample + ']'] == 1]

        ms2 = filtered_dataset.loc[filtered_dataset['Found[' + sample + ']'] == 'High']
        ms1 = filtered_dataset.loc[filtered_dataset['Found[' + sample + ']'] == 'Peak Found']

        targets = targets['AbundanceRecall[' + sample + ']'].to_list()
        decoys = decoys['AbundanceRecall[' + sample + ']'].to_list()
        ms2 = ms2['AbundanceRecall[' + sample + ']'].to_list()
        ms1 = ms1['AbundanceRecall[' + sample + ']'].to_list()

        # targets = sorted(i for i in targets if i >= 1)
        # decoys = sorted(i for i in decoys if i >= 1)
        #
        # ms2 = sorted(i for i in ms2 if i >= 1)
        # ms1 = sorted(i for i in ms1 if i >= 1)

        targets = np.log2(targets)
        targets = sorted(i for i in targets)

        decoys = np.log2(decoys)
        ms2 = np.log2(ms2)
        ms1 = np.log2(ms1)
        print(targets)
        sns.distplot(targets, bins=bins, kde=False, label='Target Proteins', hist=hist, ax=ax)
        sns.distplot(decoys, bins=bins, kde=False, label='Decoy Proteins', hist=hist, ax=ax)
        sns.distplot(ms2, bins=bins, kde=True, label='MS2 Proteins', hist=False, ax=ax)
        sns.distplot(ms1, bins=bins, kde=True, label='MS1 Proteins', hist=False, ax=ax)

        ax.legend(labels=['Target Proteins', 'Decoy Proteins', "MS2 Proteins", "MS1 Proteins"])

        ax.set_title('Sample Accession -- ' + sample + " FDR Thershold: " + str(fdr_thershold))

plot_targets_vs_decoy(fdr_thershold=0.05, bins=100)

# In[ ]:




