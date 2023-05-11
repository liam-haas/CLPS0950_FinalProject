import tkinter as tk
from PIL import ImageTk, Image

class ChemApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.geometry('500x500')
        self.root.title('Alkene Hydrobromination Predictor')
        self.root.configure(bg = '#ADD8E6')
        self.label = tk.Label(self.root, text = 'Organic Chemistry Helper!',
                              font = ('Serif', 25), bg = '#ADD8E6')
        self.label.pack(pady = 10)
        self.root.mainloop()

ChemApp()