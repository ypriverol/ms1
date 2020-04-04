import pandas as pd
import seaborn as sns
import re
import numpy as np

public_plasma_proteome = pd.read_csv("data/public-plasma-proteome.csv")
public_plasma_proteome.head()
sns.distplot(public_plasma_proteome.Log_Conc)

inhouse_plasma_proteome = pd.read_csv("data/custom-plasma-proteome.csv")
inhouse_plasma_proteome.head()

dataset = pd.read_csv("data/sample-example.csv")
dataset.head()

dataset = pd.merge(dataset, inhouse_plasma_proteome, on="UniprotAccession", how='outer')
dataset = pd.merge(dataset, public_plasma_proteome, on="Gene", how='outer')
dataset.head()


def compute_probability(row, sample):
    abundance = row['Pvalue'] + row['PPvalue'] - (row['Pvalue'] * row['PPvalue'])
    if np.isnan(row['Pvalue']) or np.isnan(row['PPvalue']):
        abundance = 0.009

    if row['Found[' + sample + ']'] == 'High':
        abundance = 1 * row['Abundance[' + sample + ']']
    if row['Found[' + sample + ']'] == 'Peak Found':
        abundance = abundance * row['Abundance[' + sample + ']']

    return abundance


samples = []
for column in list(dataset.columns.values):
    m = re.search(r"\[([A-Za-z0-9_]+)\]", column)
    if m is not None and (len(m.group(1))):
        samples.append(m.group(1))

samples = np.unique(np.array(samples))
for sample in samples:
    dataset['AbundaceRecall[' + sample + ']'] = dataset.apply(lambda row: compute_probability(row, sample), axis=1)
dataset.head()
