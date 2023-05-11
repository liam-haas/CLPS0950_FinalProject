print('Hello')

print('Hi: {}'.format(8))

import rdkit
from rdkit import Chem
from rdkit.Chem import *
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np

mol = Chem.MolFromSmiles('CCCC')

Draw.ShowMol(mol)