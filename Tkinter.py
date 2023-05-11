import tkinter as tk
from PIL import ImageTk, Image
from rdkit import Chem
from rdkit.Chem import Draw

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
        self.page1 = tk.Tk()
        self.page1.geometry('500x500')
        self.page1.title('Page One')
        self.page1.configure(bg = '#ADD8E6')

        self.description = tk.Label(self.page1, text = 'This interface has been designed to predict the product\n of a simple organic chemistry reaction. So far only the hydrobromination\n of alkenes is supported. To get started please \nenter a SMILE string in the box below.',
                                    font = ('Serif', 15), bg = '#ADD8E6')
        self.description.pack(pady = 10)

        self.SMILE = tk.Entry(self.page1, bg = 'white', fg = 'black')
        self.SMILE.bind('<KeyPress>', self.SMILEMOL)
        self.SMILE.pack()

    def SMILEMOL(self, event):
        if event.keysym == 'Return':
            self.Mol = Chem.MolFromSmiles(self.SMILE.get())
            self.Molimg = Draw.MolToFile(self.Mol, 'input SMILE')
            self.ask = tk.Label(self.page1, text = 'Just to confirm, is this the molecule you input?',
                                 bg = '#ADD8E6', font = ('Serif', 15))
            self.ask.pack(pady = 10)
            self.Disp = tk.Label(self.page1, image = ImageTk.PhotoImage(Image.open('/Users/liam/GitHub/CLPS 0950/Untitled/Module6 Test Repository/CLPS0950_FinalProject/input SMILE.png')))
            self.Disp.pack()
            

ChemApp()