from __future__ import print_function
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
from rdkit.Chem.Draw import IPythonConsole

class ChemApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.geometry('500x500')
        self.root.title('Title Page')
        self.root.configure(bg = '#ADD8E6')
        
        self.label = tk.Label(self.root, text = 'Organic Chemistry Helper!',
                              font = ('Serif', 25), bg = '#ADD8E6')
        self.label.pack(pady = 10)

        self.img = ImageTk.PhotoImage(Image.open('Benzene.jpeg'))
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
            self.img2 = ImageTk.PhotoImage(Image.open('input SMILE.png'))
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

        test_table = pd.read_csv('Test SMILES sequences - Classes.csv')
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
        self.weightdata = tk.Label(self.dataframe, text = '{:.2f}'.format(df.loc[0, 'weight']),
                                   bg = 'white', font = ('Serif', 10))
        self.weightdata.grid(row = 1, column = 1, sticky = 'n')
        self.valencedata = tk.Label(self.dataframe, text = '{}'.format(df.loc[0, 'count']),
                                   bg = 'white', font = ('Serif', 10))
        self.valencedata.grid(row = 1, column = 2, sticky = 'n')
        self.rotatedata = tk.Label(self.dataframe, text = '{}'.format(df.loc[0, 'bonds']),
                                   bg = 'white', font = ('Serif', 10))
        self.rotatedata.grid(row = 1, column = 3, sticky = 'n')

        self.dataframe.pack(pady = 10)

        if self.group == 0:
            self.disp = tk.Label(self.page2, text = 'This molecule is just a standard alkene',
                             bg = '#ADD8E6', font = ('Serif', 20))
            self.disp.pack(pady = 10)
        elif self.group == 1:
            self.disp = tk.Label(self.page2, text = 'This molecule has branching carbon chains',
                             bg = '#ADD8E6', font = ('Serif', 20))
            self.disp.pack(pady = 10)
        elif self.group == 2:
            self.disp = tk.Label(self.page2, text = 'This molecule has a benzene ring',
                             bg = '#ADD8E6', font = ('Serif', 20))
            self.disp.pack(pady = 10)
        else:
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
        self.page3 = tk.Toplevel()
        self.page3.geometry('500x500')
        self.page3.title('Page Three')
        self.page3.configure(bg = '#ADD8E6')
        
        IPythonConsole.ipython_useSVG=False
        string = self.SMILE.get()

        def phenyl_detector(string):
            str_as_list = []
            str_as_list[:0] = string
            length = len(str_as_list)
            b = str_as_list.index('1')
            a = b-1  
            if str_as_list[a] == 'C' and str_as_list[a+1] == '1' and str_as_list[a+2] == '=' and str_as_list[a+3] == 'C' and str_as_list[a+4] == 'C' and str_as_list[a+5] == '=' and str_as_list[a+6] == 'C' and str_as_list[a+7] == 'C' and str_as_list[a+8] == '=' and str_as_list[a+9] == 'C' and str_as_list[a+10] == '1':
                str_as_list.pop(a+10)
                str_as_list.pop(a+9)
                str_as_list.pop(a+8)
                str_as_list.pop(a+7)
                str_as_list.pop(a+6)
                str_as_list.pop(a+5)
                str_as_list.pop(a+4)
                str_as_list.pop(a+3)
                str_as_list.pop(a+2)
                str_as_list.pop(a+1)
                str_as_list[a] = 'Ph'
                new_string = ''.join(str_as_list)
                return new_string
            else:
                return string
            
        def sub_detector(string):
            i = string.index('=')
            length = len(string)
            if length >= 11:
                if '1' in string:
                    string = phenyl_detector(string)
                    i = string.index('=')
                    length = len(string)
                else:
                    string
            else:
                string
  
            if i >= length/2:
                str_as_list = []
                str_as_list[:0] = string
                str_as_list.reverse()
                for a in range(len(str_as_list)):
                    if str_as_list[a] == '(':
                        str_as_list[a] = ')'
                    elif str_as_list[a] == ')':
                        str_as_list[a] = '('
                    elif str_as_list[a] == 'P':
                        str_as_list[a] = 'h'
                    elif str_as_list[a] == 'h':
                        str_as_list[a] = 'P'
                new_string = ''.join(str_as_list)
                string = new_string

            j = string.index('=')
            if j == 1:
                sub = 'right'
                return sub

            if string[j-2] == 'C':
                if string[j+2] == 'C':
                    sub = 'equal'
                    return sub
                elif string[j+2] == '(':
                    sub = 'right'
                    return sub
                elif string[j+2] == 'P':
                    sub = 'right'
                    return sub
            if string[j+2] == 'C':
                if string[j-2] == ')':
                    sub = 'left'
                    return sub
                elif string[j-2] == 'h':
                    sub = 'left'
                    return sub
            if string[j-2] == 'h':
                sub = 'left'
                return sub
            if string[j+2] == 'P':
                sub = 'right'
                return sub
            if string[j-2] == ')' and string[j-2] == '(':
                if string[j-3] == 'C' and string[j+3] != 'C':
                    sub = 'left'
                    return sub
                elif string[j-3] != 'C' and string[j-3] == 'C':
                    sub = 'right'
                    return sub
                elif string[j-3] == 'C' and string[j+3] == 'C':
                    if string[j-4] == 'C' and string[j+4] != 'C':
                        sub = 'left'
                        return sub
                    elif string[j-4] != 'C' and string[j+4] == 'C':
                        sub = 'right'
                        return sub
                    elif string[j-4] == 'C' and string[j+4] == 'C':
                        sub = 'equal'
                        return sub

        if string == 'C=C':
           self.response = tk.Label(self.page3, text = 'No rxn: 1º carbocations are too unstable.',
                                    bg = '#ADD8E6', font = ('Serial', 20))
           self.response.pack(pady = 20)
    
        if '=' not in string:
            self.response = tk.Label(self.page3, text = 'No alkene found!',
                                     bg = '#ADD8E6', font = ('Serial', 20))
            self.response.pack(pady = 20)
    
        if '/' or '\\' in string:
            counter = 0
            for a in range(len(string)):  
                if string[a] == '/':
                    counter += 1
                elif string[a] == '\\':
                    counter += 1
          
            for a in range(len(string)-counter):
                str_as_list = []
                str_as_list[:0] = string
                if str_as_list[a] == '/':
                    str_as_list.pop(a)
                elif str_as_list[a] == '\\':
                    str_as_list.pop(a)
                new_string = ''.join(str_as_list)
                string = new_string
    
        i = string.index('=')
        length = len(string)
    
        if i >= length/2:
            str_as_list = []
            str_as_list[:0] = string
            str_as_list.reverse()
            for a in range(len(str_as_list)):
                if str_as_list[a] == '(':
                    str_as_list[a] = ')'
                elif str_as_list[a] == ')':
                    str_as_list[a] = '('
                elif str_as_list[a] == 'P':
                    str_as_list[a] = 'h'
                elif str_as_list[a] == 'h':
                    str_as_list[a] = 'P'
            new_string = ''.join(str_as_list)
            string = new_string
    
        a = Chem.MolFromSmiles(string)
        if sub_detector(string) == 'equal':
            rxn1 = AllChem.ReactionFromSmarts('[#6:1]=[#6:2]>>[#6:1]([Br])[#6:2]')
            rxn2 = AllChem.ReactionFromSmarts('[#6:1]=[#6:2]>>[#6:1][#6:2]([Br])')
            pdt1 = rxn1.RunReactants((a, ))[0][0]
            pdt2 = rxn2.RunReactants((a, ))[0][0]        
            Draw.MolToFile(pdt1, "product1.png")
            Draw.MolToFile(pdt2, "product2.png")
            self.response = tk.Label(self.page3, text = 'Reacting with HBr will generate two products shown below',
                                     bg = '#ADD8E6', font = ('Serif', 15))
            self.response.pack(pady = 10)
            self.pdt1 = ImageTk.PhotoImage(Image.open('product1.png'))
            self.pdt1Disp = tk.Label(self.page3, image = self.pdt1, bg = '#ADD8E6', height = 200)
            self.pdt1Disp.pack(fill = 'x', pady = 5)
            self.pdt2 = ImageTk.PhotoImage(Image.open('product2.png'))
            self.pdt2Disp = tk.Label(self.page3, image = self.pdt2, bg = '#ADD8E6', height = 250)
            self.pdt2Disp.pack(fill = 'x', pady = 5)
        elif sub_detector(string) == 'right':
            rxn = AllChem.ReactionFromSmarts('[#6:1]=[C:2]>>[#6:1][C:2][Br]')
            pdt = rxn.RunReactants((a, ))[0][0]
            Draw.MolToFile(pdt, "product.png")
            self.response = tk.Label(self.page3, text = 'Here is the product of the reaction',
                                     bg = '#ADD8E6', font = ('Serif', 20))
            self.response.pack(pady = 10)
            self.pdt = ImageTk.PhotoImage(Image.open('product.png'))
            self.pdtDisp = tk.Label(self.page3, image = self.pdt, bg = '#ADD8E6')
            self.pdtDisp.pack(fill = 'x')
        elif sub_detector(string) == 'left':
            rxn = AllChem.ReactionFromSmarts('[CH:1]=[C:2]>>[CH]([Br])[C:2]')
            pdt = rxn.RunReactants((a, ))[0][0]
            Draw.MolToFile(pdt, "product.png")
            self.response = tk.Label(self.page3, text = 'Here is the product of the reaction',
                                     bg = '#ADD8E6', font = ('Serif', 20))
            self.response.pack(pady = 10)
            self.pdt = ImageTk.PhotoImage(Image.open('product.png'))
            self.pdtDisp = tk.Label(self.page3, image = self.pdt, bg = '#ADD8E6')
            self.pdtDisp.pack(fill = 'x')
        else:
            self.response = tk.Label(self.page3, text = 'Please enter a valid SMILES sequence',
                                     bg = '#ADD8E6', font = ('Serif', 20))
            self.response.pack(pady = 10)


ChemApp()