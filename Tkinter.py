import tkinter as tk
from PIL import ImageTk, Image
from rdkit.Chem import Draw
from tkinter import messagebox
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
import pandastable

class ChemApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.geometry('500x500')
        self.root.title('Title Page')
        self.root.configure(bg = '#ADD8E6')
        
        self.label = tk.Label(self.root, text = 'Organic Chemistry Helper!',
                              font = ('Serif', 25), bg = '#ADD8E6')
        self.label.pack(pady = 10)

        self.img = ImageTk.PhotoImage(
            Image.open('/Users/liam/GitHub/CLPS 0950/Untitled/Module6 Test Repository/CLPS0950_FinalProject/Benzene.jpeg'))
        self.benzene = tk.Label(self.root, image = self.img)
        self.benzene.pack(pady = 10)

        self.startB = tk.Button(self.root, text = 'Click Here to Begin',
                                font = ('Serif', 15), bg = '#ADD8E6',
                                command = self.text_page)
        self.startB.pack(pady = 10)

        self.root.mainloop()

    def text_page(self):
        self.page1 = tk.Toplevel()
        self.page1.geometry('500x500')
        self.page1.title('Page One')
        self.page1.configure(bg = '#ADD8E6')

        self.description = tk.Label(self.page1, text = 'This interface has been designed to predict the product\n of a simple organic chemistry reaction. So far only the\n hydrobromination of alkenes is supported. To get started please \nenter a SMILE string in the box below.',
                                    font = ('Serif', 15), bg = '#ADD8E6')
        self.description.pack(pady = 10)

        self.SMILE = tk.Entry(self.page1, bg = 'white', fg = 'black')
        self.SMILE.bind('<KeyPress>', self.SMILEMOL)
        self.SMILE.pack()

    def SMILEMOL(self, event):
        if event.keysym == 'Return':
            self.Mol = Chem.MolFromSmiles(self.SMILE.get())
            Draw.MolToFile(self.Mol, 'input SMILE.png')
            self.ask = tk.Label(self.page1, text = 'Just to confirm, is this the molecule you input?',
                                 bg = '#ADD8E6', font = ('Serif', 15))
            self.ask.pack(pady = 10)
            self.img2 = ImageTk.PhotoImage(Image.open('/Users/liam/GitHub/CLPS 0950/Untitled/Module6 Test Repository/CLPS0950_FinalProject/input SMILE.png'))
            self.Disp = tk.Label(self.page1, image = self.img2, height = 250)
            self.Disp.pack()

            self.buttonframe = tk.Frame(self.page1)
            self.buttonframe.columnconfigure(0, weight = 1)
            self.buttonframe.columnconfigure(1, weight = 1)
            
            self.yesbutton = tk.Button(self.buttonframe, text = 'Yes', font = ('Serif', 15),
                                     bg = '#ADD8E6', command = self.descriptor_page)
            self.yesbutton.grid(row = 0, column = 0, sticky = 'ew')
            self.nobutton = tk.Button(self.buttonframe, text = 'No', font = ('Serif', 15),
                                      bg = '#ADD8E6', command = self.message_box)
            self.nobutton.grid(row = 0, column = 1, sticky = 'ew')

            self.buttonframe.pack(pady = 10)
            
    
    def descriptor_page(self):
        self.page2 = tk.Toplevel()
        self.page2.geometry('500x500')
        self.page2.title('Page Two')
        self.page2.configure(bg = '#ADD8E6')

        test_table = pd.read_csv('/Users/liam/GitHub/CLPS 0950/Untitled/Module6 Test Repository/CLPS0950_FinalProject/Test SMILES sequences - Classes.csv')
        test_table = test_table.dropna(axis = 1)
        PandasTools.AddMoleculeColumnToFrame(test_table, 'SMILES sequence', 'Molecule')
        hash_list = []
        for mol in test_table['Molecule']:
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)
            fp_arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(fingerprint, fp_arr)
            fp_list = fp_arr.tolist()
            hash = 0
            for entry in fp_list:
                hash += entry * math.pow(2, fp_list.index(entry))
            hash_list.append(hash)
        test_table['fingerprint'] = hash_list
        test_table = test_table.drop(index = 0)
        test_table = test_table[test_table.Class != 3]
        test_table['weight'] = [Descriptors.MolWt(mol) for mol in test_table['Molecule']]
        test_table['count'] = [Descriptors.NumValenceElectrons(mol) for mol in test_table['Molecule']]
        test_table['charge'] = [Descriptors.MinAbsPartialCharge(mol) for mol in test_table['Molecule']]
        test_table['energy'] = [Descriptors.MinAbsEStateIndex(mol) for mol in test_table['Molecule']]
        test_table['bonds'] = [Descriptors.NumRotatableBonds(mol) for mol in test_table['Molecule']]
        test_table['kappa'] = [Descriptors.Kappa1(mol) for mol in test_table['Molecule']]
        test_table['chi'] = [Descriptors.Chi0(mol) for mol in test_table['Molecule']]
        model = RandomForestClassifier(n_estimators = 200, max_depth = 7)
        features = test_table[['fingerprint', 'weight', 'count', 'charge', 'energy', 'bonds', 'kappa', 'chi']]
        target = test_table['Class']
        model.fit(features, target)

        input = Chem.MolFromSmiles(self.SMILE.get())
        df = pd.DataFrame()
        Infingerprint = AllChem.GetMorganFingerprintAsBitVect(input, radius = 2, nBits = 512)
        Infp_arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(Infingerprint, Infp_arr)
        Infp_list = Infp_arr.tolist()
        Inhash = 0
        for entry in Infp_list:
            Inhash += entry * math.pow(2, Infp_list.index(entry))
        df['fingerprint'] = [Inhash]
        df['weight'] = [Descriptors.MolWt(input)]
        df['count'] = [Descriptors.NumValenceElectrons(input)]
        df['charge'] = [Descriptors.MinAbsPartialCharge(input)]
        df['energy'] = [Descriptors.MinAbsEStateIndex(input)]
        df['bonds'] = [Descriptors.NumRotatableBonds(input)]
        df['kappa'] = [Descriptors.Kappa1(input)]
        df['chi'] = [Descriptors.Chi0(input)]
        self.group = model.predict(df)

        self.datatitle = tk.Label(self.page2, text = 'Here are some of the properties used in classifying your molecule.',
                                  bg = '#ADD8E6', font = ('Serif', 15))
        self.datatitle.pack(pady = 10)
        
        self.dataframe = tk.Frame(self.page2)
        self.dataframe.columnconfigure(0, weight = 1)
        self.dataframe.columnconfigure(1, weight = 1)
        self.dataframe.columnconfigure(2, weight = 1)
        self.dataframe.columnconfigure(3, weight = 1)
        self.dataframe.configure(bg = '#ADD8E6')

        self.fingerprint = tk.Label(self.dataframe, text = 'Molecular Fingerprint',
                                    bg = 'white', font = ('Serif', 10))
        self.fingerprint.grid(row = 0, column = 0, sticky = 's')
        self.weight = tk.Label(self.dataframe, text = 'Molecular Weight',
                                    bg = 'white', font = ('Serif', 10))
        self.weight.grid(row = 0, column = 1, sticky = 's')
        self.valence = tk.Label(self.dataframe, text = 'Valence Electron Count',
                                    bg = 'white', font = ('Serif', 10))
        self.valence.grid(row = 0, column = 2, sticky = 's')
        self.rotate = tk.Label(self.dataframe, text = 'Rotatable Bond Number',
                                    bg = 'white', font = ('Serif', 10))
        self.rotate.grid(row = 0, column = 3, sticky = 's')
        self.fingerdata = tk.Label(self.dataframe, text = '{}'.format(df.loc[0, 'fingerprint']),
                                   bg = 'white', font = ('Serif', 10))
        self.fingerdata.grid(row = 1, column = 0, sticky = 'n')
        self.weightdata = tk.Label(self.dataframe, text = '{}'.format(df.loc[0, 'weight']),
                                   bg = 'white', font = ('Serif', 10))
        self.weightdata.grid(row = 1, column = 1, sticky = 'n')
        self.valencedata = tk.Label(self.dataframe, text = '{}'.format(df.loc[0, 'count']),
                                   bg = 'white', font = ('Serif', 10))
        self.valencedata.grid(row = 1, column = 2, sticky = 'n')
        self.rotatedata = tk.Label(self.dataframe, text = '{}'.format(df.loc[0, 'bonds']),
                                   bg = 'white', font = ('Serif', 10))
        self.rotatedata.grid(row = 1, column = 3, sticky = 'n')

        self.dataframe.pack(pady = 10)

        self.disp = tk.Label(self.page2, text = 'This molecule is part of group {}'.format(self.group),
                             bg = '#ADD8E6', font = ('Serif', 20))
        self.disp.pack(pady = 10)

        self.reactbutton = tk.Button(self.page2, text = 'REACT WITH HBr',
                                    font = ('Serif', 20), command = self.reaction_page,
                                    fg = 'red', height = 12, width = 20)
        self.reactbutton.pack(pady = 5, fill = 'x')

    def message_box(self):
        messagebox.showinfo(title = 'You clicked No', message = 'Please rerun the code and start over')

    def reaction_page(self):
        pass


ChemApp()