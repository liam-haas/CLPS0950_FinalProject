#This is the machine learning algorthim that is going to be used to classify the input molecules

#Import hub
from rdkit import Chem
from rdkit.Chem import *
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import PandasTools
import numpy as np
import pandas as pd
import math
from sklearn.ensemble import RandomForestClassifier

#Reading in data we generated for training
test_table = pd.read_csv('/Users/liam/GitHub/CLPS 0950/Untitled/Module6 Test Repository/CLPS0950_FinalProject/Test SMILES sequences - Classes.csv')
test_table = test_table.dropna(axis = 1)
test_table = test_table[test_table.Class != 3]
PandasTools.AddMoleculeColumnToFrame(test_table, 'SMILES sequence', 'Molecule')

#Converting the SMILE strings to single number descriptors
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

#Adding features to the data
test_table['weight'] = [Descriptors.MolWt(mol) for mol in test_table['Molecule']]
test_table['count'] = [Descriptors.NumValenceElectrons(mol) for mol in test_table['Molecule']]
test_table['charge'] = [Descriptors.MinAbsPartialCharge(mol) for mol in test_table['Molecule']]
test_table['energy'] = [Descriptors.MinAbsEStateIndex(mol) for mol in test_table['Molecule']]
test_table['bonds'] = [Descriptors.NumRotatableBonds(mol) for mol in test_table['Molecule']]
test_table['kappa'] = [Descriptors.Kappa1(mol) for mol in test_table['Molecule']]
test_table['chi'] = [Descriptors.Chi0(mol) for mol in test_table['Molecule']]

#Instantiating the model
model = RandomForestClassifier(n_estimators = 200, max_depth = 7)

#Fitting the model
features = test_table[['fingerprint', 'weight', 'count', 'charge', 'energy', 'bonds', 'kappa', 'chi']]
target = test_table['Class']
model.fit(features, target)

#Running a test on a new molecule
test = Chem.MolFromSmiles('CCCC=CCCC')
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