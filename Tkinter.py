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
        pass

    def message_box(self):
        pass

ChemApp()