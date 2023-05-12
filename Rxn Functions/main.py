from tkinter import *
from tkinter import ttk

main = Tk()
main.geometry("750x250")

def display_text():
   global entry
   string = entry.get()
   label.configure(text=string)

#Initialize a Label to display the User Input
label=Label(main, text="", font=("Courier 22 bold"))
label.pack()

#Create an Entry widget to accept User Input
entry= Entry(main, width= 40)
entry.focus_set()
entry.pack()

#Create a Button to validate Entry Widget
ttk.Button(main, text= "Submit",width= 20, command= display_text).pack(pady=20)

canvas = Canvas(main, width = 300, height = 300)      
canvas.pack()      
img = PhotoImage(file="product.png")      
canvas.create_image(20,20, anchor=NW, image=img)

main = Tk()
main.mainloop()