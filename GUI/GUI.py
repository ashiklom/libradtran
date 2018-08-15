#! /usr/bin/env python
#
# 
# This file is part of libRadtran.
# Copyright (c) 2010 by Jonas Kylling.
# 
# libRadtran is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# libRadtran is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with libRadtran.  If not, see <http://www.gnu.org/licenses/>.

"""
GUI for uvspec

To run the GUI you need
    python 2.5 or 2.6. Older versions might work
    wxPython 2.8

Stdout and stderr are redirected to .UVSPEC_ERR. This can be turned on and off by setting redirect to True or False.

Plotting
The GUI should be able to plot most output file made with uvspec. Special care
is needed when plotting output files with the umu, zout and phi options. The GUI
uses the number of values in umu, phi and zout to parse the output file. 
If you want to plot an output file from a previous run, you have to open the inputfile
your created the output file with(or type in the same number of values in umu, phi
and zout).

Known Bugs
Resizing needed for scrollbars and expandable_on_off to behave properly

"""

import os
import sys
import re
import time
import threading
import subprocess
import itertools
import pickle

# Fix path
p = os.path.realpath(__file__)
os.chdir(os.path.dirname(p))

import variables

try:
    import wx
    import wx.lib.scrolledpanel
    import Plot
except:
    print "\nError occured while importing wx!"
    print "Quitting..."
    raise

try:
    if int("".join(wx.__version__.split("."))) < 2890:
        # Hopefully wxPython will keep the 2.8.9.0
        # format of their versionnumber.
        print "Your version of wx is %s. Version 2.8.9.0 or higher is recommended." % wx.__version__
        print "The GUI might not work as expected."
except:
    print "Your version of wx is %s. Version 2.8.9.0 or higher is recommended." % wx.__version__

try:
    import widgets
    import create
    # if variables.enable_plotting:
    #     import plotting
    # else:
    #     plotting = widgets.DummyClass()
    import MakeExampleDialog
    from MakeExampleDialog import MakeExampleDialogResult
    
except:
    print "\nError occured while importing widgets or create!"
    raise

import wx.lib.platebtn

def pr(path): # Returns path to resources
    global libRadtran_GUI_path
    return os.path.join(libRadtran_GUI_path, path)

def pe(path): #Returns path to examples
    global libRadtran_examples_path
    return os.path.join(libRadtran_examples_path, path)

widgets.pr = pr
widgets.pe = pe
create.pr = pr
create.pe = pe


vars_widgets = vars(widgets).values()

#if sys.platform == "darwin":
#    ClearButton = wx.lib.platebtn.PlateButton
#else:
ClearButton = wx.Button

def GetDocumentationFunction(name):
    def tmp(event):
        GetDocumentation(name)
    return tmp

def Debug(name, num):
    def tmp(event):
        print "%s: %i" % (name, num)
        event.Skip()
    return tmp

def OnClear(obj):
    def tmp(event):
        obj.Clear()
        event.Skip()
    return tmp

def SplashScreen():
    """
    Display a splashscreen
    """
    try:
        if not os.path.exists(pr("resources/splash.png")):
            print "Could not find splash.png"
            return
        image = wx.Image(pr("resources/splash.png"), wx.BITMAP_TYPE_ANY)
        image = image.ConvertToBitmap()
        splashstyle = wx.SPLASH_CENTRE_ON_SCREEN|wx.SPLASH_NO_TIMEOUT
        style = wx.SIMPLE_BORDER#|wx.STAY_ON_TOP
        SplashWindow = wx.SplashScreen(image, splashstyle, -1, parent=None)
        return SplashWindow
    except Exception, error:
        print "Could not create splash screen.'"
        print error

def getInput():
    data=[]
    for option in Options:
	if Options[option].IsChanged() and  Options[option].IsSet():
		value = Options[option].GetWriteValue()
		data.extend(value)
    return data

FileDialog = widgets.FileDialog
ErrorMessage = widgets.ErrorMessage


class Start(wx.Notebook):
    """The startup panel
    
    Events for the buttons in this class are bound in the main class.
    """
    def __init__(self, parent):
        wx.Notebook.__init__(self, parent, )

        Home = wx.Panel(self, )
        Examples = wx.Panel(self)
        
        self.WelcomeText = wx.StaticText(Home, -1,
                                         variables.start_welcome_text)

        self.New  = wx.Button(Home, -1, "New")
        self.Open = wx.Button(Home, -1, "Open")
        self.Try = wx.Button(Examples, -1, "Try")
        self.NewExample =  wx.Button(Examples, -1, "New example",)
        self.EditExample = wx.Button(Examples, -1, "Edit example")

        self.Examples = widgets.ExampleList(Examples,
                                            create.example_list)

        Logo = widgets.Logo(Home, -1, wx.Bitmap(pr("resources/esa.gif")))

        main_sizer = wx.BoxSizer(wx.HORIZONTAL)
        common_sizer = wx.BoxSizer(wx.VERTICAL)
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        example_sizer = wx.BoxSizer(wx.VERTICAL)
        example_button_sizer = wx.BoxSizer(wx.HORIZONTAL)

        button_sizer.Add(self.New, 0, flag=wx.TOP|wx.RIGHT|wx.CENTER,
                         border=20)
        button_sizer.Add(self.Open, 0, wx.TOP|wx.CENTER, border=20)
        
        common_sizer.Add(self.WelcomeText, 0, wx.CENTER)
        common_sizer.Add(button_sizer, 0, wx.CENTER)
        common_sizer.Add((-1, -1), 1)
        common_sizer.Add(Logo, 0, wx.ALIGN_RIGHT)

        example_button_sizer.Add((-1, -1), 1)
        example_button_sizer.Add(self.NewExample, 0,
                                 wx.CENTER|wx.RIGHT, border=15)
        example_button_sizer.Add(self.Try, 0,
                                 wx.CENTER|wx.RIGHT, border=15)
        example_button_sizer.Add(self.EditExample, 0, wx.CENTER)
        example_button_sizer.Add((-1, -1), 1)

        example_sizer.Add(self.Examples, 1, wx.CENTER|wx.EXPAND)
        example_sizer.Add(example_button_sizer, 0,
                          wx.CENTER|wx.TOP|wx.BOTTOM|wx.EXPAND, border=5)
        
        Home.SetSizer(common_sizer)
        Examples.SetSizer(example_sizer)
        
        self.AddPage(Home, "Home")
        self.AddPage(Examples, "Examples")


class Panel(wx.Panel):
    def __init__(self, parent, sizer_flag=wx.VERTICAL):
        """ Optional sizer_flag, default is wx.VERTICAL """
        wx.Panel.__init__(self, parent, style=variables.option_style)
        self.Sizer = wx.BoxSizer(sizer_flag)

        self.is_page_panel = True


class NoteBook(wx.Notebook):
    """
    Various options are added to sizers in wx.Panel (done in CreatePanel).
    The panels are added to a Page, and the pages are added to this notebook
    """
    def __init__(self, parent):
        wx.Notebook.__init__(self, parent, -1, style=wx.NB_TOP)

        self.PageNum = {}
        self.Pages = {}

        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGING, self.OnPageChange)

    def AddTabs(self, option_list, out):
        """ A page should alway be added by calling this function """
        for page in option_list:
            self.AddTab(page, out)

    def AddTab(self, create_dict, out):
        max_items = variables.max_column_items # Maximum number of options per column
        
        main_panel = widgets.ScrolledWindowFix(self, style=wx.VSCROLL)
        main_panel.Sizer = wx.BoxSizer(wx.HORIZONTAL)

        label = create_dict["title"]
        
        options = create_dict["options"]
        options_len = len(options)

        if options_len > max_items: # Divide the options into two panels
            times, left = divmod(options_len, 2)
            panel, option_dict = self.CreatePanel(main_panel,
                                                  options[:times + left])
            out.update(option_dict)
            panel2, option_dict = self.CreatePanel(main_panel,
                                                   options[times + left:])

            main_panel.Sizer.Add(panel, 1, wx.RIGHT, border=10)
            main_panel.Sizer.Add(panel2, 1,)
        else:
            panel, option_dict = self.CreatePanel(main_panel, options)
            
            main_panel.Sizer.Add(panel, 1, wx.EXPAND)
            
	self.InputButton= wx.Button(panel,-1,'Show input file')
	self.InputButton.Bind(wx.EVT_BUTTON,self.showInput)
	panel.Sizer.Add(self.InputButton, 0, flag=wx.ALL|wx.CENTER, border=5)
        out.update(option_dict)
        
        main_panel.SetSizer(main_panel.Sizer)
        main_panel.SetAutoLayout(True)

        self.AddPage(main_panel, label)

    def AddPage(self, panel, label):
        self.PageNum[label] = self.GetPageCount()
        self.Pages[label] = panel
        wx.Notebook.AddPage(self, panel, label)
    
    def CreatePanel(self, parent, options):
            out = {}
            panel = Panel(parent)
            
            #panel.SetBackgroundColour("RED")

            panel.Show(False)

            largest_string = 0
            for option in options:
                s = 0
                if not "hide_name" in option:
                    s += panel.GetTextExtent(option["name"])[0]
                if "plot" in option and not option["plot"] is None:
                    s += variables.plot_button_size
                largest_string = max(largest_string, s)

            for option in options:
                #print option["name"], "plot" in option and \
                #    not option["plot"] is None, largest_string
                # Various layout and initialization
                container = wx.Panel(panel, style=variables.option_style)
		container.parent=panel
                sizer = wx.BoxSizer(wx.HORIZONTAL)
           
                # Debug
                #print option, option["name"]
     

                if "class" in option:
                    obj = option["class"](container, **option["args"])
                else:
                    obj = widgets.Option(container,
                                         largest_string=largest_string,
                                         option=option,
					 non_unique=option.non_unique)

                sizer.Add(obj, 1, flag=wx.EXPAND|wx.CENTER)
                container.SetSizer(sizer)

                # Documentation setup
                
                if not "class" in option:
                    GetDocumentation = GetDocumentationFunction(option["name"])

                    for i in container, obj, obj.name_obj:
                        if not i: continue
                        i.Bind(wx.EVT_LEFT_DOWN, GetDocumentation)
                        i.Bind(wx.EVT_LEFT_DCLICK, GetDocumentation)
                
                if not "class" in option and option["help"]:
                    obj.SetToolTipString(option["help"])

                # Class specific sizing
                if type(obj) in (widgets.OpenButton,
                                 #widgets.SaveRunButtons,
                                 widgets.Text, 
                                 ):
                    panel.Sizer.Add(container, 0,
                                    flag=wx.ALL|wx.EXPAND|wx.CENTER, border=5)
                elif type(obj) in (widgets.Html,
                                    MakeExampleDialog.DescriptionImage):
                    panel.Sizer.Add(container, 1,
                                    wx.EXPAND|wx.ALL|wx.CENTER, border=5)
                else:
                    panel.Sizer.Add(container, 0,
                                    flag=wx.ALL|wx.CENTER|wx.EXPAND, border=5)
                    
                out[option["name"]] = obj

                if "external_event_bindings" in option:
                    for binding in option["external_event_bindings"]:
                        obj.Bind(binding["type"], binding["func"])

            panel.SetSizer(panel.Sizer)
            return panel, out

    def OnPageChange(self, event):
        n = event.GetSelection()
        o = event.GetOldSelection()
        n = self.GetPage(n)
        for w in n.GetChildren():
            w.Show(True)
        if o >= 0:
            o = self.GetPage(o)
            #o.Show(False)

        n.SetScrollbars(1, 4, 1, 100)
        event.Skip()        

    def showInput(self, event):
	InputDialog(self)

class InputDialog(wx.Frame):
    def __init__(self, parent):
	wx.Frame.__init__(self, parent, title="Input file",
                          style=wx.FRAME_TOOL_WINDOW|wx.FRAME_FLOAT_ON_PARENT |
                          wx.DEFAULT_FRAME_STYLE,
                          pos=variables.documentation_window_position)
	inputDialog = ''.join( inp+'\n' for inp in getInput() )
	sw=wx.ScrolledWindow(self, -1)
	text=wx.StaticText(sw,-1,inputDialog)
	sw.SetScrollbars(1,1,text.GetSize()[0]+1,text.GetSize()[1]+1)
	self.Show(True)
	

class GeneralPage(widgets.ScrolledWindowFix):
    def __init__(self, parent, mainframe):
        widgets.ScrolledWindowFix.__init__(self, parent, style=wx.VSCROLL)
        
        self.mainframe = mainframe
        self.panel = Panel(self)

        # Great for debug!
        #self.panel.SetBackgroundColour("RED")

        Sizer = wx.BoxSizer(wx.VERTICAL)
        
        # Add buttons
        self.save_buttons = widgets.SaveButtons(self.panel)
        self.run_buttons = widgets.SaveButtons(self.panel, lbl1="Save and Run",
                                               lbl2="Run")

        for obj in self.save_buttons, self.run_buttons:
            Sizer.Add(obj, 0, wx.EXPAND|wx.ALL^wx.BOTTOM, 5)

        # Add HTML
        self.html_info = widgets.Html(self.panel,
                                      filename=pr("resources/how_to.html"))
        Sizer.Add(self.html_info, 1, wx.EXPAND|wx.ALL^wx.BOTTOM, 5)

        # Add back button
        self.back_button = wx.Button(self.panel, -1, "Back")
        Sizer.Add(self.back_button, 0, wx.CENTER|wx.TOP|wx.BOTTOM, border=5)

        self.panel.SetSizer(Sizer)

        MainSizer = wx.BoxSizer(wx.HORIZONTAL)
        MainSizer.Add(self.panel, 1, wx.EXPAND)

        self.SetSizer(MainSizer)
        self.SetAutoLayout(True)

        self.save_buttons.Bind(wx.EVT_BUTTON, self.OnSave)
        self.run_buttons.Bind(wx.EVT_BUTTON, self.OnRun)
        self.Bind(wx.EVT_BUTTON, self.OnBack, self.back_button)

    def OnSave(self, event):
        self.mainframe.OnSave(event)

    def OnRun(self, event):
        self.mainframe.OnRun(event)

    def OnBack(self, event):
        self.mainframe.OnNew()

    def SetInputFile(self, *args):
        self.save_buttons.SetValue(*args)
        
    def SetOutputFile(self, *args):
        self.run_buttons.SetValue(*args)

    def EnableSaveButtons(self, value=True):
        self.save_buttons.Enable(value)

    def EnableRunButtons(self, value=True):
        self.run_buttons.Enable(value)

    def GetInputFile(self):
        return self.save_buttons.GetValue()

    def GetRunFile(self):
        return self.run_buttons.GetValue()

    def Clear(self, event):
        for obj in self.save_buttons, self.run_buttons:
            obj.Clear()

class PlotPage(widgets.ScrolledWindowFix):
    def __init__(self, parent, mainframe):
        widgets.ScrolledWindowFix.__init__(self, parent, style=wx.VSCROLL)
        
        self.mainframe = mainframe
        self.panel = Panel(self)

        # Great for debug!
        #self.panel.SetBackgroundColour("RED")

        Sizer = wx.BoxSizer(wx.VERTICAL)
        
        # Add buttons
        self.buttons = []
        self.texts = []
        button_data = (("Plot 2D", Plot.Plot2D),
                       ("Plot map", Plot.PlotMap),
                       ("Plot block", Plot.PlotBlock))
        for title, plot_type in button_data:
            tmp_sizer = wx.BoxSizer(wx.HORIZONTAL)
            button = widgets.SaveButtons(self.panel, lbl1="Open and Plot",
                                          lbl2="Plot", flags=wx.FD_OPEN)
            text = wx.StaticText(self.panel, -1, title)
            tmp_sizer.Add(text, 0, wx.RIGHT, 5)
            tmp_sizer.Add(button, 1)
            Sizer.Add(tmp_sizer, 0, wx.EXPAND|wx.ALL^wx.BOTTOM, 5)
            self.buttons.append(button)
            self.texts.append(text)

        # Fix Layout
        self.panel.SetSizer(Sizer)

        MainSizer = wx.BoxSizer(wx.HORIZONTAL)
        MainSizer.Add(self.panel, 1, wx.EXPAND)

        self.SetSizer(MainSizer)
        self.SetAutoLayout(True)

        # Bind events
        self.plot2D_buttons, self.plotMap_buttons, \
            self.plotBlock_buttons = self.buttons
        self.plot2D_buttons.Bind(wx.EVT_BUTTON, self.OnPlot2D)
        self.plotMap_buttons.Bind(wx.EVT_BUTTON, self.OnPlotMap)
        self.plotBlock_buttons.Bind(wx.EVT_BUTTON, self.OnPlotBlock)

    def OnPlot2D(self, event):
        fname = self.plot2D_buttons.GetValue()
	column_names = Options['output_user'].GetWriteValue()[0].split()[1:]
	if not column_names: 
		column_names = ["wavelength",
                               "edir",
                               "edn",
                               "eup",
                               "uavgdir", 
                               "uavgdn",
                               "uavgup"]
        plotter = Plot.Plot2D(column_names=(column_names) )
        
        plotter.ShowPlotWindow(self.mainframe, fname)

    def OnPlotMap(self, event):
        fname = self.plotMap_buttons.GetValue()
        plotter = Plot.PlotMap()
        plotter.ShowPlotWindow(self.mainframe, fname)
    
    def OnPlotBlock(self, event):
        fname = self.plotBlock_buttons.GetValue()
        plotter = Plot.PlotBlock(phi=Options["phi"].GetValue(),
                                 umu=Options["umu"].GetValue())
        plotter.ShowPlotWindow(self.mainframe, fname)

    def SetValue(self, value):
        self.plot2D_buttons.SetValue(value)
        #self.plotMap_buttons.SetValue(value)
        self.plotBlock_buttons.SetValue(value)


class Main(wx.Frame):
    """Main Frame
    """
    global Options
    Options = {"uvspec": widgets.UvspecOptionDummy()}    # Options contains all options
    #plotting.Options = Options
    global Other
    Other = {}
    edit_example = False
    Saved = True # Start as saved, changes if anything is changed
    def __init__(self, quick_start=False):
        self.ExampleThings = [{
                "title": "Example",
                "options": [
                    {"name":"Example name",
                     "class": widgets.Text,
                     "args": {
                            },
                     "external_event_bindings": ({
                                "type":wx.EVT_TEXT,
                                "func":self.OnExampleName,
                                },),
                     "parents":[],
                     "non_parents":[],
                     "non_parent_exceptions":[],
                     "childs":[],
                     },
                    {"name":"ExampleDescriptionImage",
                     "hide_name": True,
                     "class": MakeExampleDialog.DescriptionImage,
                     "args": {},
                     "parents":[],
                     "non_parents":[],
                     "non_parent_exceptions":[],
                     "childs":[],
                     }
                    ],
                },]
                 
        #-----------------------------------------------------------------------
        create_list = create.create_list # Edit create.py when adding new options
        #-----------------------------------------------------------------------
        wx.Frame.__init__(self, None,
                          title=variables.main_window_title,
                          size=variables.main_window_size)

        # Fix
        self.error_txt = [] #most errors should be collected here before they are displayed,
                                # it is better with only one long error message than ten small ones

        global tmp_panel
        tmp_panel = wx.Panel(self)
        tmp_panel.Show(False)
        
        # Edit contains all the pages
        self.Edit = NoteBook(self)
        if not quick_start:
            # Create the General tab
            self.General = GeneralPage(self.Edit, self)
            self.Edit.AddPage(self.General, "General")
            # Create all the other tabs
            self.Edit.AddTabs(create_list, Options)
            self.PlotPage = PlotPage(self.Edit, self)
            self.Edit.AddPage(self.PlotPage, "Plot")
            # Page for examples
            self.Edit.AddTabs(self.ExampleThings, Other)
            self.Edit.RemovePage(self.Edit.PageNum["Example"])

        # Only show Edit after either OnNew or OnOpenExistingINP have been run
        self.Edit.Show(False)

        self.Start = Start(self)
        
        self.HorSizer = wx.BoxSizer(wx.VERTICAL)
        self.HorSizer.Add(self.Start, 1,
                          wx.EXPAND|wx.BOTTOM|wx.RIGHT|wx.LEFT, border=5)
        
        self.SetSizer(self.HorSizer)

        self.Center(True)
        self.SetSizeHints(*variables.main_window_size_hints)

        if variables.documentation_window:
            documentation_window = widgets.DocWindow = \
                widgets.DocumentationWindow(self, doc_dir=pr("resources/html_doc"))
            global GetDocumentation
            GetDocumentation = documentation_window.GetDocumentation

        self.Layout()
        self.Show(True)

        widgets.Options = Options

        # Bind some special events
        self.Bind(wx.EVT_BUTTON, self.OnNewInp, self.Start.New)
        self.Bind(wx.EVT_BUTTON, self.OnOpenExistingINP, self.Start.Open)
        self.Bind(wx.EVT_BUTTON, self.OnExample, self.Start.Try)
        self.Bind(wx.EVT_BUTTON, self.OnNewExample, self.Start.NewExample)
        self.Bind(wx.EVT_BUTTON, self.OnEditExample, self.Start.EditExample)
        self.Bind(wx.EVT_CLOSE, self.OnExit, self)
        
    def OpenFile(self, fname):
        """
        Reads a file and set options inside
        """
        self.SetOptions(self.ReadFile(fname))
        self.General.SetInputFile(fname)
        widgets.Saved = True

    def OnNew(self):
        if self.GuardianAngel():
            self.error_txt = []
            self.Clear()
            self.HorSizer.Detach(self.Edit) # Remove the start panel
            self.Edit.Show(False)

            if self.edit_example: # Remove the example tab
                self.General.EnableSaveButtons()
                self.new_example = False
                self.Edit.RemovePage(self.Edit.PageNum["Example"])
                self.edit_example = False
                
            self.HorSizer.Insert(0, self.Start, 1, wx.EXPAND)
            self.Start.Show(True)
        
            self.HorSizer.Layout()
            self.Refresh()
            
    def Clear(self):
        for option in Options:
	    try:
        	Options[option].ClearAll()
	    except:
		Options[option].Clear()
        for special in Other:
            Other[special].Clear()
        
    def OnNewInp(self, event):
        """When the New Button is pressed, self.Main.Start.New"""
        tmp = FileDialog(ext_format="*.INP",
                         flags=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if tmp:
            self.InitEdit()
            self.General.SetInputFile(tmp)
            
    def OnOpen(self, event):
        """ Check if anything have to be saved """
        if not self.GuardianAngel():
            return
        tmp = FileDialog(ext_format="*.INP")
        if tmp:
            self.Clear()
            self.OpenFile(tmp)

    def OnOpenExistingINP(self, event):
        """
        When open is pressed on the start page.
        """
	wildcard = "Input files (*.inp; *.INP)|*.inp;*.INP|" \
         "All files (*)|*"
        tmp = FileDialog(ext_format=wildcard)
        if tmp:
            self.InitEdit()
            self.OpenFile(tmp)
            
    def OnExample(self, event):
        """
        When a user tries an example.
        """
        example = self.Start.Examples.GetCurrentPage()
        self.InitEdit()
        self.OpenFile(example.example_path)

    def OnNewExample(self, event):
        """
        When a user wants to make an example.
        """
        self.General.EnableSaveButtons(False)
        self.EditExample()

    def OnEditExample(self, event):
        """
        When the user wants to edit an existing example.
        """
        example = self.Start.Examples.GetCurrentPage()
        self.OpenFile(example.example_path)
        Other["Example name"].SetValue(example.title)
        Other["ExampleDescriptionImage"].SetDescription(example.txt)
        Other["ExampleDescriptionImage"].SetImage(example.img_path)
        self.EditExample()

    def OnExampleName(self, event):
        event.Skip()
        if self.edit_example:
            title = Other["Example name"].GetValue()
            example_inp = os.path.join(pe("examples/GUI"),
                                       title.lower().replace(" ", "_"),
                                       "inp.inp")
            self.General.SetInputFile(example_inp)

    def EditExample(self):
        self.edit_example = True
        self.Edit.AddPage(self.Edit.Pages["Example"], "Example")
        self.InitEdit()
        
    def InitEdit(self):
        """
        Removes the start panel and inserts the Edit panel
        Add menus.
        """
        # Fix new layout
        self.HorSizer.Detach(self.Start)
        self.Start.Show(False)
        
        self.HorSizer.Insert(0, self.Edit, 1,
                             wx.ALL|wx.EXPAND, border=5)
        self.Edit.Show(True)
        
        self.HorSizer.Layout()
        self.Refresh() # To prevent some ugly redraw errors
        
        self.AddMenus()

        for page in create.create_list :
            for opt in page["options"]:
                if "uvspec" in opt["parents"]:
                    Options["uvspec"].dependencies.append(opt["name"])

        for opt in Options["uvspec"].dependencies:
            Options[opt].CanChange()
        # for opt in Options.values():
        #     for child in opt.childs:
        #         try:
        #             Options[child].Enable(False)
        #         except KeyError:
        #             pass

        widgets.Saved = True

    def AddMenus(self):
        menu_bar = wx.MenuBar()

        for title, items in (("File",
                              (("Open", self.OnOpen),
                               ("Save", self.OnSave),
                               ("Save as...", self.OnSaveAs),
                               (), # Separator
                               ("Make example...", self.OnMakeExample),
                               (),
                               ("Quit", self.OnExit),
                               ),
                              ),
                             ):
            tmp = wx.Menu()
            for i in items:
                if not i:
                    tmp.AppendSeparator()
                else:
                    t, handler = i
                    a = tmp.Append(-1, t)
                    self.Bind(wx.EVT_MENU, handler, a)
            menu_bar.Append(tmp, title)

        self.SetMenuBar(menu_bar)

    def OnMakeExample(self, event):
        self.MakeExample()

    def OnSaveAs(self, event):
        tmp = FileDialog(".INP", flags=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if tmp:
            self.Save(tmp)

    def OnExit(self, event):
        if self.GuardianAngel():
            self.Destroy()

    def GuardianAngel(self):
        """ Protect the user from losing unsaved work."""
        if widgets.Saved:
            return True
        else:
            Dialog = wx.MessageDialog(None, "Save Input file?", "Save",
                                      wx.YES_NO|wx.CANCEL|wx.YES_DEFAULT)
            answer = Dialog.ShowModal()
            if answer == wx.ID_YES:
                save_file = self.General.GetInputFile()
                if not save_file:
                    save_file = FileDialog(format=".INP",
                                           flags=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
                    self.General.SetInputFile(save_file)
                self.Save(save_file)
                widgets.Saved = True
                Dialog.Destroy()
            elif answer == wx.ID_CANCEL:
                Dialog.Destroy()
                return False
        return True
    
    def OnSave(self, event):
        if self.edit_example:
            self.SaveExample()
            return
        tmp = self.General.GetInputFile()
        if not tmp:
            tmp = FileDialog(".INP",
                             flags=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
            if not tmp: return
            self.General.SetInputFile(tmp)
        self.Save(tmp)
        widgets.Saved = True
        
    def OnRun(self, event):
        out_file = self.General.GetRunFile() #Other["Output File"].GetValue()
        tmp_file = os.path.abspath(".tmp_UVSPEC.INP")
        tmp_save = self.Saved
        self.Save(tmp_file) # Save a temporary INP file with the latest modifications
        self.Saved = tmp_save # Set it to the previous value, since self.Saved sets it to True

        win = widgets.RunUvspec(self, tmp_file, out_file)        
        #Other["Plot_Out"].SetValue(out_file)
        self.PlotPage.SetValue(out_file)

    def env_check(self):
        # Not in use
        """Check if the environment variable and uvspec exists
        
        If the environment varaible "LIBRADTRAN_DATA_FILES" don't exist and uvspec 
        can't be found in the path the run button will be disabled.
        """
        if not os.getenv("LIBRADTRAN_DATA_FILES"):
            warning = wx.MessageDialog(None, 
            """The Environmentvariabel LIBRADTRAN_DATA_FILES was not found!\nYou can still edit .INP files and plot .OUT files, but you can not run UVSPEC! It is strongly recommended to set the environment variabel and restart the GUI!""",
            "Warning!", # Title
            wx.OK)      # Style
            warning.Center(True)
            warning.ShowModal()
            # Disable run Button
            Other["Output File"].Enable(False)
        try:
            subprocess.Popen(("uvspec", "-v"),
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE)
        except OSError:
            warning = wx.MessageDialog(self, 
            """Could not find UVSPEC in your PATH!\nYou can still edit .INP files and plot .OUT files, but you can not run UVSPEC! It is strongly recommended to change your PATH and restart the GUI!""",
            "Warning!", # Title
            wx.OK)      # Style
            warning.Center(True)
            warning.ShowModal()
            # Disable run Button
            Other["Output File"].Enable(False)

    def MakeExample(self):
        res = MakeExampleDialogResult(self)
        if not res: return
        title, img, description = res
        print res
        # Ensure that the directory does not collide
        example_dir = os.path.join(pe("examples/GUI"),
                                    title.lower().replace(" ", "_"))
        if os.path.exists(example_dir):
            msg = "Example directory already exists. Pleace choose a new name of your example!"
            ErrorMessage(msg)

        # Create the example (directory and necessary files)
        try:
            os.mkdir(example_dir)
            # Save the image
            img.SaveFile(os.path.join(example_dir, "big.jpg"),
                         wx.BITMAP_TYPE_JPEG)
            img.SetSize((60, 60))
            img.SaveFile(os.path.join(example_dir, "thumb.jpg"),
                         wx.BITMAP_TYPE_JPEG)
            # Save the text
            f = open(os.path.join(example_dir, "info.txt"), "w")
            f.write("Title: %s\n" % title)
            f.write(description)
            f.close()
            # Save the input file
            self.Save(os.path.join(example_dir, "inp.inp"))
        except (OSError, IOError), e:
            msg = ("Error encountered. Tried to make example in the directory %s" % \
                repr(example_dir)) + str(e)
            ErrorMessage(msg)

    def Save(self, fname):
        data = ["" + "#Generated by libRadtran GUI(%s)\n" % variables.version]
	data.extend(getInput())
        try:
            print "Saving input file to '%s'." % (fname)
            print repr(data)
            f = open(fname, "w")
            f.write("\n".join(data).strip())
	    f.write("\n")
            f.close()
        except IOError:
            msg = "Could not save input file to '%s'." % (fname)
            #print msg
            ErrorMessage(msg)

    def SaveExample(self):
        title = Other["Example name"].GetValue()
        if not title:
            msg = "Please specify a name for your example!"
            ErrorMessage(msg)
            return
        description, img = Other["ExampleDescriptionImage"].GetValue()
        example_dir = os.path.join(pe("examples/GUI"),
                                    title.lower().replace(" ", "_"))
        if os.path.exists(example_dir):
            if not os.path.isdir(example_dir):
                msg = "Exampledirectory exists, but is not a directory (%s). Please choose a new name of your example!" % repr(example_dir)
                ErrorMessage(msg)
        else:
            os.mkdir(example_dir)
            
        try:
            # Save the image
            img.SaveFile(os.path.join(example_dir, "big.jpg"),
                         wx.BITMAP_TYPE_JPEG)
            img.SetSize((60, 60))
            img.SaveFile(os.path.join(example_dir, "thumb.jpg"),
                         wx.BITMAP_TYPE_JPEG)
            # Save the text
            f = open(os.path.join(example_dir, "info.txt"), "w")
            f.write("Title: %s\n" % title)
            f.write(description)
            f.close()
            # Save the input file
            self.Save(os.path.join(example_dir, "inp.inp"))
        except (OSError, IOError), e:
            msg = ("Error encountered. Tried to make example in the directory %s" % \
                repr(example_dir)) + str(e)
            ErrorMessage(msg)

    def ReadFile(self, path, includes=[]):
        try:
            f = open(path, "r")
            data = f.readlines()
            f.close()
        except IOError:
            msg = "Could not open '%s'." % (path)
            self.error_txt.append(msg)
            #print " "*len(includes) + msg
            return

        data = [line.split() for line in data]
        line_nos = [(path, i) for i in xrange(1, len(data)+1)]

        # Remove lines with comments and merge continuous lines
        line = 0
        while line < len(data):
            # Skip empty lines
            if not data[line]:
                data.pop(line)
                line_nos.pop(line)
                continue
            # Skip comments
            if data[line][0].startswith("#"):
                data.pop(line)
                line_nos.pop(line)
                continue
            # Remove comments from the line
            elif [True for word in data[line] if (word.find("#") != -1)]:
                tmp = []
                for word in data[line]:
                    pos = word.find("#")
                    if pos != -1:
                        if word != "#":
                            tmp.append(word[:pos])
                        break
                    else:
                        tmp.append(word)
                data[line] = tmp
            # continous line
            elif data[line][-1].endswith("\\"):
                data[line][-1] = data[pos][-1][:-1] # remove the \
                    # if the \ was preceded and continued by whitespace
                if data[line][-1] == "":
                    data[line].pop()
                    line_nos.pop()
                data[line].extend(data[line + 1])
            else:
                line += 1

        # Get the includes and include them
        buff = 0
        this_includes = []
        for line in xrange(len(data)):
            if (data[line + buff][0] == "include" and \
                    len(data[line + buff]) == 2):
                this_includes.append(data.pop(line + buff)[1])
                line_nos.pop(line + buff)
                buff -= 1

        for include_path in this_includes:
            if not os.path.exists(path):
                msg = "Include '%s' in '%s' does not exist." % (include_path, path)
                self.error_txt.append(msg)
                #print " " * len(includes) + msg
            # If the file has been included before.
            elif (this_includes + includes).count(include_path) != 1:
                msg = "Include %s included more than one time." + \
                    " One time in %s. Please fix this " % (include_path, path)
                self.error_txt.append(msg)
                #print " " * len(includes) + msg
            else:
                include_data = self.ReadFile(include_path,
                                             includes + this_includes)
                # The include file might contains errors and return None.
                if include_data:
                    data.extend(include_data[0])
                    line_nos.extend(include_data[1])
        
        # Only Display error message in the most outer scope.
        #if self.error_txt and not includes:
        #    ErrorMessage(self.error_txt)
        #    self.error_txt = []

        return data, line_nos

    def SetOptions(self, from_file):
        if not from_file and self.error_txt:
            ErrorMessage(self.error_txt)
            return

        data, line_nos = from_file
        for line, pos in zip(data, line_nos):
            name = line[0]
            try:
                obj = Options[name]
                obj.AddValue(line[1:] if line[1:] else [])
            except KeyError:
                error = "Unknown option %s in file '%s' line %i." % \
                    (repr(name), pos[0], pos[1])
                self.error_txt.append(error)
            except widgets.SetValueError, e:
                self.error_txt.append("Illegal value %s for option '%s'"
                                      " in file '%s' line %i." % \
                                          (repr(e.value), e.name, pos[0], pos[1]))
        if self.error_txt:
            ErrorMessage(self.error_txt)

def main(GUI_path="", examples_path=".."):
    global libRadtran_GUI_path
    libRadtran_GUI_path = GUI_path
    global libRadtran_examples_path
    libRadtran_examples_path = examples_path

    create.create_example_list()
    
    app = wx.App(redirect=variables.redirect, filename=variables.redirect_file)
    print time.strftime("%Y-%m-%d-%H-%M-%S")
    splash = SplashScreen()
    frame = Main(quick_start=False)
    app.SetTopWindow(frame)
    if splash: splash.Close() # close directly after all GUI elements are initialized
    #frame.env_check() # have to be called after the Splashscreen is closed
                               # because the MessageDialog and Splashscreen 
                               # will block each other
    app.MainLoop()

if __name__ == "__main__":
    main()
