"""
from Tkinter import Tk, RIGHT, BOTH, RAISED, TOP
from ttk import Frame, Button, Style, Entry


class GUI(Frame):
  
    def __init__(self, parent, series):
        Frame.__init__(self, parent)   
        self.series = series
        self.parent = parent
        self.initUI()

    def OnClickOk(self):
        print self.textbox.get()
        self.series.updatePeriod(float(self.textbox.get()))

    def OnClickClose(self):
        self.parent.destroy()
        
    def initUI(self):
        self.parent.title("Particle Filter")
        self.style = Style()
        self.style.theme_use("default")
        
        frame = Frame(self, relief=RAISED, borderwidth=1)
        frame.pack(fill=BOTH, expand=True)
        
        self.pack(fill=BOTH, expand=True)
        self.textbox = Entry(self)
        self.textbox.pack(side=TOP)
        closeButton = Button(self, text="Close", command=self.OnClickClose)
        closeButton.pack(side=RIGHT, padx=5, pady=5)
        okButton = Button(self, text="OK", command=self.OnClickOk)
        okButton.pack(side=RIGHT)
"""