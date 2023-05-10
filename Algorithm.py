#This is the machine learning algorthim that is going to be used to classify the input molecules

from rdkit import Chem
from rdkit.Chem import *
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import PandasTools
from rdkit.Chem import RDKFingerprint
import numpy as np
import pandas as pd
import math
from sklearn.ensemble import RandomForestClassifier

test_table = pd.read_csv('/Users/liam/GitHub/CLPS 0950/Untitled/Module6 Test Repository/CLPS0950_FinalProject/Test SMILES sequences - Classes.csv')
test_table = test_table.dropna(axis = 1)
PandasTools.AddMoleculeColumnToFrame(test_table, 'SMILES sequence', 'Molecule')

hash_list = []
for mol in test_table['Molecule']:
  bi = {}
  fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512, bitInfo=bi)
  fp_arr = np.zeros((1,))
  DataStructs.ConvertToNumpyArray(fingerprint, fp_arr)
  fp_list = fp_arr.tolist()
  hash = 0
  for entry in fp_list:
    hash += entry * math.pow(2, fp_list.index(entry))
  hash_list.append(hash)
test_table['fingerprint'] = hash_list
test_table = test_table.drop(index = 0)

model = RandomForestClassifier(n_estimators = 200, max_depth = 7)

test_table['weight'] = [Descriptors.MolWt(mol) for mol in test_table['Molecule']]
test_table['count'] = [Descriptors.NumValenceElectrons(mol) for mol in test_table['Molecule']]
test_table['charge'] = [Descriptors.MinAbsPartialCharge(mol) for mol in test_table['Molecule']]
test_table['energy'] = [Descriptors.MinAbsEStateIndex(mol) for mol in test_table['Molecule']]
test_table['bonds'] = [Descriptors.NumRotatableBonds(mol) for mol in test_table['Molecule']]
test_table['kappa'] = [Descriptors.Kappa1(mol) for mol in test_table['Molecule']]
test_table['chi'] = [Descriptors.Chi0(mol) for mol in test_table['Molecule']]

features = test_table[['fingerprint', 'weight', 'count', 'charge', 'energy', 'bonds', 'kappa', 'chi']]
target = test_table['Class']
model.fit(features, target)

test = Chem.MolFromSmiles('C/C=C(C)/CC(C(C)CC)CCC')
df = pd.DataFrame()
fingerprint = AllChem.GetMorganFingerprintAsBitVect(test, radius = 2, nBits = 512)
fp_arr = np.zeros((1,))
DataStructs.ConvertToNumpyArray(fingerprint, fp_arr)
fp_list = fp_arr.tolist()
hash = 0
for entry in fp_list:
  hash += entry * math.pow(2, fp_list.index(entry))
df['fingerprint'] = [hash]
df['weight'] = [Descriptors.MolWt(test)]
df['count'] = [Descriptors.NumValenceElectrons(test)]
df['charge'] = [Descriptors.MinAbsPartialCharge(test)]
df['energy'] = [Descriptors.MinAbsEStateIndex(test)]
df['bonds'] = [Descriptors.NumRotatableBonds(test)]
df['kappa'] = [Descriptors.Kappa1(test)]
df['chi'] = [Descriptors.Chi0(test)]
group = model.predict(df)
print(group)