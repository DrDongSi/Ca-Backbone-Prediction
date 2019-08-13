from tkinter import *
from tkinter import filedialog
from email_validator import validate_email, EmailNotValidError
from tkinter import messagebox
from prediction.prediction import run_predictions
import os
# NOTE: MUST USE PYTHON3+

    


class Window(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master
        
        welcomelabel = Label(self.master, font='bold', background='white', text='Welcome to Automated Cryo-EM Modeler!')
        welcomelabel.pack(side=TOP, pady=10)
        
        # *** Frame Input ***
        self.inputframe = Frame(self.master)
        self.inputframe.pack(fill=X)
        self.inputframe.configure(background='white')
        # Label
        inputlabel = Label(self.inputframe, background='white', text='Choose a input directory', width=20, justify=LEFT, anchor='w')
        inputlabel.pack(side=LEFT, padx=5, pady=5)
        # Button       
        inputbutton = Button(self.inputframe, text='Select Input', width=13, command=self.__selectinput)
        inputbutton.pack(side=LEFT, padx=5)
        
        # *** Frame Output ***
        self.outputframe = Frame(self.master)
        self.outputframe.pack(fill=X)
        self.outputframe.configure(background='white')
        # Label
        outputlabel = Label(self.outputframe, background='white', text='Choose a output directory', width=20, justify=LEFT, anchor='w')
        outputlabel.pack(side=LEFT, padx=5, pady=5)
        # Button
        outputbutton = Button(self.outputframe, text='Select Output', width=13, command=self.__selectoutput)
        outputbutton.pack(side=LEFT, padx=5)
    
        # *** Frame Email ***
        self.emailframe = Frame(self.master)
        self.emailframe.pack(fill=X)
        self.emailframe.configure(background='white')
        # Label
        emaillabel = Label(self.emailframe, background='white', text='Enter email: ', width=20, justify=LEFT, anchor='w')
        emaillabel.pack(side=LEFT, padx=5, pady=5)
        # Button
        self.emailentry = Entry(self.emailframe, bd=5, width=30)
        self.emailentry.pack(side=LEFT, padx=5)


        # *** Frame Optional ***
        self.otherframe = Frame(self.master)
        self.otherframe.pack(fill=X)
        self.otherframe.configure(background='white')
        # Label
        otherlabel = Label(self.otherframe, font='bold', background='white', text='Optional Options', width=20, justify=LEFT, anchor='w')
        otherlabel.pack(side=LEFT, padx=5, pady=5)


        # *** Frame Threshold ***
        self.thresholdframe = Frame(self.master)
        self.thresholdframe.pack(fill=X)
        self.thresholdframe.configure(background='white')
        # Label
        thresholdlabel = Label(self.thresholdframe, background='white', text='Choose threshold file', width=20, justify=LEFT, anchor='w')
        thresholdlabel.pack(side=LEFT, padx=5)
        # Button
        thresholdbutton = Button(self.thresholdframe, text='Select File', width=13, command=self.__selectthreshold)
        thresholdbutton.pack(side=LEFT, padx=5)


        # *** Frame Hidedust ***
        self.hidedustframe = Frame(self.master)
        self.hidedustframe.pack(fill=X)
        self.hidedustframe.configure(background='white')
        # Label
        hidedustlabel = Label(self.hidedustframe, background='white', text='Choose hide dust file', width=20, justify=LEFT, anchor='w')
        hidedustlabel.pack(side=LEFT, padx=5)
        # Button
        hidedustbutton = Button(self.hidedustframe, text='Select File', width=13, command=self.__selectthreshold)
        hidedustbutton.pack(side=LEFT, padx=5)


        # *** Frame Check Existing ***
        self.cexistingframe = Frame(self.master)
        self.cexistingframe.pack(fill=X)
        self.cexistingframe.configure(background='white')
        # Label
        cexistinglabel = Label(self.cexistingframe, background='white', text='Check Existing', width=20, justify=LEFT, anchor='w')
        cexistinglabel.pack(side=LEFT, padx=5)
        # Check box
        self.cexistingcheck = BooleanVar()
        cexistingcheckbutton = Checkbutton(self.cexistingframe, background='white', variable=self.cexistingcheck, onvalue=True, offvalue=False)
        cexistingcheckbutton.pack(side=LEFT, padx=5)


        # *** Frame Debug ***
        self.debugframe = Frame(self.master)
        self.debugframe.pack(fill=X)
        self.debugframe.configure(background='white')
        # Label
        debuglabel = Label(self.debugframe, background='white', text='Debug Mode', width=20, justify=LEFT, anchor='w')
        debuglabel.pack(side=LEFT, padx=5)
        # Check box
        self.debugcheck = BooleanVar()
        debugcheckbutton = Checkbutton(self.debugframe, background='white', variable=self.debugcheck, onvalue=True, offvalue=False)
        debugcheckbutton.pack(side=LEFT, padx=5)

        # *** Frame Skip ***
        self.skipframe = Frame(self.master)
        self.skipframe.pack(fill=X)
        self.skipframe.configure(background='white')
        # Label
        skiplabel = Label(self.skipframe, background='white', text='Skip pipeline steps', width=20, justify=LEFT, anchor='w')
        skiplabel.pack(side=LEFT, padx=5)
        # Option Menu
        choices = [0,1,2,3,4,5,6,7]
        self.tkvar = IntVar()
        self.tkvar.set(0)
        skipmenu = OptionMenu(self.skipframe, self.tkvar, *choices)
        skipmenu.pack(side=LEFT, padx=5)


        # *** Verify Frame ***
        self.verifyframe = Frame(self.master)
        self.verifyframe.pack(fill=X)
        self.verifyframe.configure(background='white')
        # Label
        verifylabel = Label(self.verifyframe, font='bold', background='white', text='Verify Above', width=20, justify=LEFT, anchor='w')
        verifylabel.pack(side=LEFT, padx=5)

        # *** Run Frame ***
        self.runframe = Frame(self.master)
        self.runframe.pack(fill=X)
        self.runframe.configure(background='white')
        # Button
        runbutton = Button(self.runframe, text='RUN!', width=13, command=self.__run)
        runbutton.pack(padx=5)


         # *** Graphic Frame ***
        self.graphicframe = Frame(self.master)
        self.graphicframe.pack(fill=X)
        self.graphicframe.configure(background='white')
        # Photo
        self.numframes = 100
        self.numframestep = 10
        self.numframeindices = self.numframes // self.numframestep
        self.graphicimage = [PhotoImage(file=os.getcwd().replace(os.sep, '/') + '/gui/protein.gif', format='gif -index %i' %(i)) for i in range(0, self.numframes, self.numframestep)]
        self.graphiclabel = Label(self.graphicframe, image=self.graphicimage[0], background='white')
        self.graphiclabel.pack()

        self.inputdirectory = ''
        self.outputdirectory = ''
        self.thresholdfilename = None
        self.hidedustfilename = None
        self.inputdirlabel = Label(self.inputframe)
        self.outputdirlabel = Label(self.outputframe)
        self.threshlabel = Label(self.thresholdframe)

    def update(self, ind):
        self.frame = self.graphicimage[ind]
        ind += 1
        ind %= self.numframeindices
        self.graphiclabel.configure(image=self.frame)
        milliseconds = 300
        self.master.after(milliseconds, self.update, ind)
   
    def __selectthreshold(self):
        self.thresholdfilename =  filedialog.askopenfilename(initialdir = "/",title = "Select file")
        if self.thresholdfilename is None:
            self.thresholdfilename = None
        self.threshlabel.destroy()
        if self.thresholdfilename is not None:
            self.threshlabel = Label(self.thresholdframe, background='white', text='...' + self.thresholdfilename[-20:], width=20, justify=LEFT, anchor='w')
        else:
            self.threshlabel = Label(self.thresholdframe, background='white', text='...', width=20, justify=LEFT, anchor='w')
        self.threshlabel.pack(side=LEFT, padx=5)
 
    def __selecthidedust(self):
        self.hidedustfilename =  filedialog.askopenfilename(initialdir = "/",title = "Select file")
        if self.hidedustfilename is '':
            self.hidedustfilename = None

    def __selectinput(self):
        self.inputdirectory = filedialog.askdirectory()
        self.inputdirlabel.destroy()
        self.inputdirlabel = Label(self.inputframe, background='white', text='...' + self.inputdirectory[-20:], width=20, justify=LEFT, anchor='w')
        self.inputdirlabel.pack(side=LEFT, padx=5, pady=5)
        
    def __selectoutput(self):
        self.outputdirectory = filedialog.askdirectory()
        self.outputdirlabel.destroy()
        self.outputdirlabel = Label(self.outputframe, background='white', text='...' + self.outputdirectory[-20:], width=20, justify=LEFT, anchor='w')
        self.outputdirlabel.pack(side=LEFT, padx=5, pady=5)

    def __validate(self):
        try:
            validate_email(self.emailentry.get())
        except EmailNotValidError as e:
            messagebox.showinfo('Missing!', 'Please enter a valid email address')
            return False

        if self.inputdirectory is '':
            messagebox.showinfo('Missing!', 'Please select a valid input folder')
            return False

        if self.outputdirectory is '':
            messagebox.showinfo('Missing!', 'Please select a valid output folder')
            return False

        return True

    def __run(self):
        if self.__validate():
            self.inputdirectory += '/' if self.inputdirectory[-1] != '/' else ''
            self.outputdirectory += '/' if self.outputdirectory[-1] != '/' else ''
            run_predictions(self.inputdirectory, self.outputdirectory, self.thresholdfilename, self.tkvar.get(), self.cexistingcheck.get(), None, self.debugcheck.get())
        

def run_gui():
    root = Tk()
    app = Window(root)
    root.title('Automated Cryo-EM Modeler (ACEM)')
    root.geometry('500x600')
    root.resizable(0,0)
    root.configure(background='white')
    root.after(0, app.update, 0)
    root.mainloop()