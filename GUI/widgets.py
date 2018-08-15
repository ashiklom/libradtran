#!/usr/bin/env python
#
#  widgets.py
#  GUI
#
# This file is part of libRadtran.
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

#
"""
This file contains "widgets" used by the GUI.
"""

import sys
sys.path.append("../src_py")

import wx
import os
import sys
import time
import itertools
import threading
import variables
import subprocess_fix as subprocess
import webbrowser

# """
# Fallback to old plotting could be added if wxmple is not installed.
# """
# if variables.enable_plotting:
#     import wxmpl


#from wx.lib import plot
from wx import grid
import wx.html
import wx.lib.customtreectrl
import wx.lib.scrolledpanel
import wx.lib.buttons
import wx.lib.platebtn

import option_definition as OptionDefinition
from GUI_definition import *

import Plot

import FloatSpin # Made by Andrea Gavana

# Quote from "http://www.wxpython.org/what.php":
# """
# wxPython is a cross-platform toolkit. This means that the same
# program will run on multiple platforms without modification.
# """
# Why do I then need the following code?
#try:
#    ClearButton = wx.lib.platebtn.PlateButton
#    raise NotImplementedError, "How do I test if PlateButton is avaiable?"
#except NotImplementedError:
ClearButton = wx.Button

Saved = True

def pr(path): return path
def pe(path): return path

def AttachFunction(attachFunc, orgFunc):
    """ Attach an 'on event' function to another function. """
    def newFunc(*args, **kwargs):
        r = orgFunc(*args, **kwargs)
        attachFunc(None)
        return r
    return newFunc

def FileDialog(ext_format="*", flags=wx.FD_OPEN, defaultDir="",
               defaultFile=""):
    """A simple file dialog
    
    FileDialog(format="") -> path
    
    The function takes an optional format as argument.
    Creates a wx.FileDialog and returns the file that have been chosen.
    """
    if flags&wx.FD_OPEN: title = "Open file"
    else: title = "Save file"
    Dialog = wx.FileDialog(None, title, defaultDir=defaultDir,
                           defaultFile=defaultFile, wildcard=ext_format,
                           style=flags)
    if Dialog.ShowModal() == wx.ID_OK:
        filename=Dialog.GetFilename()
        dirname=Dialog.GetDirectory()
        Dialog.Destroy()
        return os.path.join(dirname, filename)
    Dialog.Destroy()

def OnOrOff(value):
    """ Return on if value is true, else off
    Used by ExpandableOnOff
    """
    if value: return "on"
    else: return "off"

def style_test(style):
    while True:
        yield style

def ErrorMessage(message="", title="Error", style=wx.OK|wx.CENTRE):
    if not isinstance(message, basestring):
        message = "\n".join(message)
    msg = wx.MessageDialog(None, message=message, caption=title, style=style)
    msg.ShowModal()

class DummyClass():
    def __init__(self, *args, **kwargs):
        return None
    def __set__(self, *args, **kwargs):
        pass
    def __getattr__(self, *args, **kwargs):
        return self

global Options
Options = {}


class SetValueError(Exception):
    def __init__(self, name, value):
        self.name = name
        self.value = value


class Widget(wx.Panel):
    def __init__(self, parent, style=variables.option_style,
                 optional=False):
        wx.Panel.__init__(self, parent, style=style)

        self._parent = parent
        self._changed = False
        self._edit = False
        self._has_default = False
        self._default = None
        self._optional = optional

    def IsOptional(self):
        return self._optional

    def IsSet(self):
        """
        1. Value has been changed
         1.1 User has changed the value
         1.2 Value changed by file
        And
        2. I'm enabled
        """
        if self.IsChanged() and self.IsEnabled():
            return True
        #elif self._has_default:
        #    return True
        else:
            return False

    def OnEdit(self,event):
	self._edit = True
	self._parent.OnEdit()
	self.Touch()

    def IsEdited(self):
        """
        User has changed the value since last call to Clear
        """
        return self._edit
    def IsChanged(self):
        """
        Value has been changed either by user or by file.
        Returns True if the option should be written to file.
        """
        if self._changed or self._edit:
            return True
        else:
            return False
    def GetValue(self):
        raise NotImplementedError
    def SetValue(self, *args):
        raise NotImplementedError
    def Touch(self):
        self._changed = True
    def TouchValue(self, *args):
        self.touch()
        self.SetValue(*args)
    def Clear(self):
        """
        Clear an option.
        That is reset to default value, if default exists, else a
        "normal" value (e.g. 0 or "").
        IsEdited and IsChanged will return False after Clear.
        IsSet may return True if the option has a default value.
        """
        self._changed = False
        self._edit = False
        if self._has_default:
            self.SetValue(self._default)
        else:
            self.SetValue(None)
	try:	self._Clear()
	except: pass

def addInput(inp, parent):
    sizer_args = [0, wx.CENTER|wx.RIGHT]
    dt = inp.__class__
    kwargs = {"optional":inp.optional}
    if dt == FileInput:
        obj = OpenButton(parent, **kwargs)
        parent.file_obj = obj
        sizer_args[0] = 2
    elif dt == FloatInput:
        obj = Spinner(parent, spintype=float,
                      default=inp.default,
                      valid_range=inp.valid_range, **kwargs)
    elif dt == IntegerInput:
        obj = Spinner(parent, spintype=int,
                      default=inp.default,
                      valid_range=inp.valid_range, **kwargs)
    elif dt == ListInput or dt == IntegerListInput:
        obj = Choice(parent, default=inp.default,
                     choices=inp.valid_range,optional=inp.optional,logical_file=inp.logical_file)
    elif dt == TextInput:
        obj = Text(parent, default=inp.default, **kwargs)
        sizer_args[0] = 2
    elif dt == BooleanInput:
        obj = OnOff(parent, **kwargs)
    elif dt == VariableNumberOfLinesInput:
        obj = VariableLines(parent, inp.valid_range, **kwargs)
    else:
        raise ValueError, "Never reached!"

    return obj, sizer_args

class UvspecOptionDummy():
    """
    Class used as the parent uvspec.
    (Currently only useful for consistency)
    """
    def __init__(self, *args):
        self.name = "uvspec"
        self.dependencies = []
    def IsSet(self, *args):
        return True
    def IsChanged(self, *args):
        """
        Never write to file. (After all we do no implement GetValue)
        """
        return False
    def IsEdited(self, *args):
        return True
    def Clear(self, *args):
        pass
    def CanChange(self):
        pass


class Option(wx.Panel):
    def __init__(self, parent, largest_string, option, border=0,
                 colour=variables.option_background_colour,
                 style=variables.option_style,non_unique=False,is_child=False):
        wx.Panel.__init__(self, parent, style=style)
        
        plot = option.plot

        can_plot = not plot is None

        self._set = False
	self._mandatory = False
        self.continious_update = hasattr(option, "continious_update") and \
            option.continious_update

        self.name = option.name
        self.dependencies = option.dependencies
        self.canEnable = option.canEnable
        self.isMandatory = option.isMandatory
        self.inputs = []

        self.sub_names = []

        self.file_obj = None

        self.sizer = wx.BoxSizer(wx.HORIZONTAL)

	#NEW: unique option, main_sizer and varibles for allocation of new option class
       	self.main_sizer=wx.BoxSizer(wx.VERTICAL) 
	self.option=option
	self.largest_string=largest_string
	self.parent=parent
	self.non_unique=non_unique
	self._is_child=is_child

        # Add name
        name_str = self.name
        border = largest_string - self.GetTextExtent(name_str)[0] + 5
        border -= (variables.plot_button_size if can_plot else 0)
        self.name_obj = wx.StaticText(self, -1, name_str)
        self.sizer.Add(self.name_obj, 0, wx.RIGHT|wx.CENTER, border=border)

	self.font = wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
	self.fontIsSet = wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.BOLD)

        # Add plotting
        if variables.enable_plotting and can_plot:
            plot_type = {"map":Plot.PlotMap, "2D":Plot.Plot2D}[plot["plot_type"]]
            self.plot_button = Plot.PlotButton(self, filegetter=self,
                                               plot_type=Plot.Plot2D,
                                               optional_args=plot["optional_args"])
            self.sizer.Add(self.plot_button, 0, wx.RIGHT|wx.CENTER, 5)

        gui_inputs = option.gui_inputs

        # Add inputs
        # tokens = option.tokens
        # # TODO:
        # if isinstance(option, OptionDefinition.SignedFloats):
        #     tokens = [OptionDefinition.addToken("", str)]
        # if tokens is None or tokens == []:
        #     obj, sizer_args = addInput(None, self, self.name_obj)
        #     sizer.Add(obj, *sizer_args)            
        #     self.inputs.append(obj)
        if len(gui_inputs) == 1:
            obj, sizer_args = addInput(gui_inputs[0], self)
            self.inputs.append(obj)
            self.sizer.Add(obj, *sizer_args)
        else:
            small_font = wx.Font(10, wx.FONTFAMILY_DEFAULT,
                                 wx.FONTSTYLE_NORMAL,
                                 wx.FONTWEIGHT_NORMAL)

            for inp in gui_inputs:
                tmp_sizer = wx.BoxSizer(wx.VERTICAL)
                
                tmp_name = wx.StaticText(self, -1,
                                         inp.name.split(".")[-1])
                self.sub_names.append(tmp_name)
                tmp_name.SetFont(small_font)
                tmp_sizer.Add(tmp_name, 0, wx.CENTER)

                obj, sizer_args = addInput(inp, self)

                self.inputs.append(obj)
                tmp_sizer.Add(obj, 1, wx.EXPAND|wx.CENTER)

                self.sizer.Add(tmp_sizer, *sizer_args, border=3)

	if gui_inputs[0].__class__ == BooleanInput:  #Quick fix for boolean input; ugly
		self.Clearbutton = self.inputs[0]
		self.booleanInput = True
	else:
	        # Add reset/Clear button
		self.booleanInput = False
        	self.Clearbutton = ClearButton(self, -1, "Reset")
	        self.sizer.Add(self.Clearbutton, 0,
        	          flag=wx.CENTER|wx.LEFT|wx.ALIGN_RIGHT, border=5)
		if not self._is_child: self.Clearbutton.Show(False)
        self.Clearbutton.Bind(wx.EVT_BUTTON, self.onClear)

        # Button for adding new option line, if option is not unique
	self.Childbutton = wx.Button(self, -1, "Add new line")
	self.sizer.Add(self.Childbutton, 0, flag=wx.CENTER|wx.LEFT|wx.ALIGN_RIGHT, border=5)
	self.Childbutton.Bind(wx.EVT_BUTTON, self.AddChild)
	if not self.non_unique:
		self.Childbutton.Show(False)
	self.childs = []

        # Add tooltips
        help = option.help
        if help:
            for obj in self, self.name_obj, self.Clearbutton:
                self.SetToolTipString(help)
            for obj in self.inputs:
                obj.SetToolTipString(help)
            for obj in self.sub_names:
                obj.SetToolTipString(help)

        self.Enable(False)

	#NEW: unique option, main_sizer
	self.main_sizer.Add(self.sizer, 1, wx.EXPAND|wx.CENTER)
	self.SetSizer(self.main_sizer)


    def OnEdit(self):
        self.Clearbutton.Show(True)
	self.OnChange()

    def onClear(self, event):
        self.Clear()
	self.OnChange()

    def is_set_function(self, name):
        try:
            return Options[name].IsSet()
        except KeyError:
            return False

    def get_value_function(self, name):
        try:
            return Options[name].GetValue()
        except KeyError:
            return None

    def OnChange(self):
	self.IsSet()
	self.UpdateFont()
	self.AwakeDeps()
	self.Layout()

    def AddChild(self, event=None):
	child=Option(self, self.largest_string, self.option, is_child=True)
	child.Enable(True)
	self.main_sizer.Add(child, 1, wx.EXPAND|wx.CENTER)
	self.childs.append(child)
	self.main_sizer.Layout()
	self.Layout()

    def Layout(self):
	wx.Panel.Layout(self)
	self.parent.Layout()
	self.parent.parent.Layout()

    def CanChange(self):
        # Decide if the option should be mandatory
        self._mandatory = self.isMandatory(self.is_set_function, self.get_value_function)
	self.UpdateFont()
        # Decide if the option should be enabled
        self.Enable(self.canEnable(self.is_set_function,
                                   self.get_value_function))
        

    def AwakeDeps(self):
        for dep in self.dependencies:
            try:
                Options[dep].CanChange()
            except KeyError, e:
                # Useful for debug if you edit the src_py/* files
                #print "KeyError, %s has missing dependency %s." % (self.name, dep)
                pass

    def GetFileValue(self):
        if self.file_obj:
            return self.file_obj.GetValue()
        else:
            raise TypeError

    def ClearAll(self):
	self.Clear()
	for child in self.childs:
	    self.main_sizer.Remove(child)
	self.childs=[]

    def Clear(self):
        for obj in self.inputs:
            obj.Clear()
        if self.Clearbutton.__class__ == ClearButton:
            self.Clearbutton.Show(False)
	if self._is_child:
	    self.parent.main_sizer.Remove(self)

    def IsChanged(self):
        return any([obj.IsChanged() for obj in self.inputs])

    def UpdateFont(self):
	if self._set:		self.name_obj.SetFont(self.fontIsSet)
	else:			self.name_obj.SetFont(self.font)
	if self._mandatory:	self.name_obj.SetForegroundColour(wx.Colour(255, 0, 0))
        else:			self.name_obj.SetForegroundColour(wx.Colour(0, 0, 0))

    def IsSet(self):
        """ Returns True if an option is set """
	self._set=False
	if self.IsEnabled():
		if all([obj.IsSet() for obj in self.inputs if not obj._optional]):
			self._set=True
	return self._set

    def IsEdited(self):
        return any([obj.IsEdited() for obj in self.inputs])

    def Touch(self):
	for obj in self.inputs:
	    if not obj._optional:
		obj.Touch()
        self.OnEdit()

    def TouchValue(self, *args):
        self.SetValue(*args)
        self.Touch()

    def AddValue(self, *args):
        """
        Used when reading from file.
        Works like TouchValue, except for funny options like ...
        """
        self.TouchValue(*args)

    def SetValue(self, value):
        if isinstance(value, basestring):
            self.SetValue(value, value.split())
        else:
            l = len(self.inputs)
            if len(value) > l:
                self.SetValue(value[:l-1] + [" ".join(value[l-1:])])
            elif self.IsSet() and self.non_unique:
		self.AddChild()
		self.childs[-1].AddValue(value)
	    else:
                for i in xrange(len(value)):
		    if self.booleanInput: 	i_inp = i+1
		    else: 			i_inp = i
                    self.inputs[i_inp].SetValue(value[i])
		    self.inputs[i_inp].Touch()

    def GetValue(self):
        out = []
        for obj in self.inputs:
            if obj.IsOptional() and not obj.IsChanged(): continue
            out.append(obj.GetValue())
        return out

    def GetWriteValue(self):
	out= [self.name + " " + " ".join([str(i) for i in self.GetValue()])]
	for child in self.childs:
	    if child.IsChanged() and  child.IsSet():
		out[0] += '\n'+child.GetWriteValue()[0]
	return out
    
    def Enable(self, enable=True):
	for obj in self.inputs:
		obj.Enable(enable)
	for name in self.sub_names:
		name.Enable(enable)
	for child in self.childs:
		child.Enable(enable)
	self.Clearbutton.Enable(enable)
	self.Childbutton.Enable(enable)
	if enable:	self.UpdateFont()
	else:		self.name_obj.SetForegroundColour(wx.Colour(150,150,150))

class Spinner(Widget):
    """
    spintype is the type of the spinner, either float or int
    labels is a list of strings containing a label which will appear above a spinctrl.
    Leave labels empty if you want only one spinner.

    wx.SpinCtrl is not used for int because it misbehaves.
    """
    def __init__(self, parent, default, valid_range, spintype=float,
                 increment=1.0, **kwargs):
        Widget.__init__(self, parent, **kwargs)

        self.parent = parent

        self.spintype = spintype
        self._default = default

        value = valid_range[0] # This is the lower limit
        if not (self._default is None):
            self._has_default = True
            value = self._default


        self.min_val, self.max_val = valid_range

        increment = (valid_range[1] - valid_range[0])/100.
        increment = 1.0 if increment > 1 else increment

        Sizer = wx.BoxSizer(wx.HORIZONTAL)        
        
        if self.spintype == float:
            self.Spin = FloatSpin.FloatSpin(self, value=value,
                                            max_val=valid_range[1],
                                            min_val=valid_range[0],
                                            increment=increment)
        elif self.spintype == int or self.spintype == long:
            self.Spin = FloatSpin.FloatSpin(self, value=value,
                                            max_val=valid_range[1],
                                            min_val=valid_range[0],
                                            increment=1)
            self.Spin.SetFormat("{0:100.18g}")
            
        Sizer.Add(self.Spin, flag=wx.CENTER|wx.RIGHT)

        self.SetSizer(Sizer)

        if not self._has_default:
            self.SetBlank()
        
        self.Bind(wx.EVT_SPINCTRL, self.OnEdit)
        
    def SetValue(self, value):
        try:
            if value is None:
                self.SetBlank(True)
            else:
                self.SetBlank(False)
                self.Spin.SetValue(value)
        except ValueError:
            raise SetValueError(self.parent.name, value)
                
    def GetValue(self):
        return str(self.spintype(self.Spin.GetValue()))

    def Clear(self):
        self._changed = False
        self._edit = False
        if self._has_default:
            self.SetValue(self._default)
        else:
            self.SetBlank()
            #self.SetValue(0 if self.min_val < 0 else self.min_val)

    def SetBlank(self, isblank=True):
        #self.SetValue(None)        
        self.Spin.SetBlank(isblank)

    def SetToolTipString(self, *args, **kwargs):
        wx.Panel.SetToolTipString(self, *args, **kwargs)
        self.Spin.SetToolTipString(*args, **kwargs)


class Choice(Widget):
    """
    A list which the user can choose from
    """
    def __init__(self, parent, default=None, choices=[], optional=False, logical_file=False,
                 **kwargs):
        Widget.__init__(self, parent, optional=optional, **kwargs)

        if not default is None:
            self._has_default = True

        self.parent = parent

        self._default = default
        self.choices = choices
        if optional and not self.choices.count(''):
            self.choices.insert(0, "")

        Sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.Choice = wx.Choice(self, -1, choices=self.choices)
        
        #Sizer.Add(wx.Panel(self, -1,), 1, wx.CENTER)
        # I do not know why this one looks ugly
        Sizer.Add(self.Choice, 0, wx.CENTER)

	self.OpenButton=OpenButton(self)
	Sizer.Add(self.OpenButton,1,wx.EXPAND)
	self.OpenButton.Show(False)
	if logical_file:
		self.Choice.Append('user_defined')

        self.SetSizer(Sizer)

        self.SetValue(self._default)

        self.Bind(wx.EVT_CHOICE, self.OnEdit)

    def OnEdit(self,event=None):
	self._edit = True
	self._parent.OnEdit()
	self.Touch()
	if self.GetChoiceValue()=='user_defined':
		self.OpenButton.Show(True)
		self.Layout()
	else:
		self.OpenButton.Show(False)
		self.Layout()

    def SetValue(self, value):
        if value:
            if type(value) == list: value = value[0]
            try:
                self.Choice.SetSelection(self.choices.index(value.lower()))
	    except:
		try:
			self.Choice.SetStringSelection('user_defined')
			self.OpenButton.SetValue(value)
			self.OpenButton.Show(True)
        	except:
                	raise SetValueError(self.parent.name, value)
        else:
            self.Choice.SetSelection(0)

    def GetChoiceValue(self):
        #selected = self.choices[self.Choice.GetCurrentSelection()]
        selected = self.Choice.GetStringSelection()
        if selected != "None":
            return selected

    def GetValue(self):
	selected = self.GetChoiceValue()
	if selected == 'user_defined':
		return self.OpenButton.GetValue()
	elif selected != "None":
		return selected

    def SetToolTipString(self, *args, **kwargs):
        wx.Panel.SetToolTipString(self, *args, **kwargs)
        self.Choice.SetToolTipString(*args, **kwargs)

    def _Clear(self):
	self.OpenButton.Show(False)
	self.Layout()

    def Layout(self):
	Widget.Layout(self)
	self.parent.Layout()

class OnOff(Widget):
    """OnOff button"""
    def __init__(self, parent, default=None, out_value="", **kwargs):
        """
        write_value is the opposite of the default value, i.e. the logical value of
        the button when the state of the button will be written to the input file.
        """
        Widget.__init__(self, parent, **kwargs)
        
        if not default is None:
            self._has_default = True

        self.parent = parent

        self.out_value = out_value
        self._default = default
        self.write_value = not bool(default)

        Sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.Button = wx.ToggleButton(self, -1, "Off")

        Sizer.Add(self.Button, flag=wx.CENTER)
        self.SetValue(False)

        self.SetSizer(Sizer)
        
        self.Button.Bind(wx.EVT_TOGGLEBUTTON, self._OnEdit)
        
    def _OnEdit(self, event=None):
        self.OnEdit(event)
        self.ChangeLabel()
	if not self.Button.GetValue():
		self.parent.onClear(event)
        
    def ChangeLabel(self):
        """Fancy changing of the buttonlabel
        Makes it easier to see if the button is on or off.
        """
        if self.Button.GetValue():
            self.Button.SetLabel("Deactivate")
        else:
            self.Button.SetLabel("Activate")

    def SetValue(self, value):
        """
        Set the value of the button
        Any True value will turn the button to the self.write_value
        Any False value will turn the button to the opposite of the self.write_value
        """
        if value == "on":
            self.Button.SetValue(True)
        elif value == "off":
            self.Button.SetValue(False)
        elif value == "":
            self.Button.SetValue(True)
        else:
            if not value: # In case you want to turn the button off
                self.Button.SetValue(not self.write_value)
            else:
                self.Button.SetValue(self.write_value)
        self.ChangeLabel()

    def GetValue(self):
        value = self.Button.GetValue()
        if value == self.write_value:
            return self.out_value
        else:
            return ""

    def SetToolTipString(self, *args, **kwargs):
        wx.Panel.SetToolTipString(self, *args, **kwargs)
        self.Button.SetToolTipString(*args, **kwargs)


class OpenButton(Widget):
    """Let's the user open files and view and edit the path to the currently selected file

    ext_format restricts the matched files to files matching the pattern in ext_format.
    The default is *.* which should match any file on any platform (probably not).
    """
    def __init__(self, parent, default=None, buttonlabel="Open",
                 flags=wx.FD_OPEN, ext_format="*", **kwargs):
        Widget.__init__(self, parent, **kwargs)

        if not default is None:
            self._has_default = True

        self.parent = parent
        self._default = default

        self.ext_format = ext_format
        self.flags = flags

        Sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.FileInfo = wx.TextCtrl(self, -1)
        self.FileInfo.SetValue(default if not default is None else "")
        
        self.Button = wx.Button(self, -1, buttonlabel)
        
        Sizer.Add(self.FileInfo, 20, wx.CENTER|wx.RIGHT, border=10)
        Sizer.Add(self.Button, 0, flag=wx.CENTER)

        self.SetSizer(Sizer)

        self.FileInfo.Bind(wx.EVT_TEXT, self.OnEdit)

        self.Button.Bind(wx.EVT_BUTTON, self.OnEdit)
        self.Button.Bind(wx.EVT_BUTTON, self.onOpen)
        
    def SetValue(self, value):
        """Change the value of self.FileInfo"""
        self._set = True
        try:
            if value:
                if type(value) == str or type(value) == unicode:
                    self.FileInfo.ChangeValue(value)
                else:
                    self.FileInfo.ChangeValue(value[0])
            #self.Change() # Triggered by SetValue anyway.
            else:
                # Do not use SetValue, because it triggers self.Change
                self.FileInfo.ChangeValue("")
        except ValueError:
            raise SetValueError(self.parent.name, value)

    def GetValue(self):
        """Returns the path to the currently selected file"""
        val = self.FileInfo.GetValue()
        return val

    def onOpen(self, event):
        """FileDialog when the open button is pressed"""
        path = os.path.split(self.FileInfo.GetValue())
        val = FileDialog(ext_format=self.ext_format, flags=self.flags,
                         defaultDir=path[0], defaultFile=path[1])
        if val:
            self.SetValue(val)
            self.OnEdit(event)
            event.Skip()

    def disable(self):
        self.Button.Disable()
        self.FileInfo.Disable()

    def enable(self):
        self.Button.Enable()
        self.FileInfo.Disable()

    def SetToolTipString(self, *args, **kwargs):
        wx.Panel.SetToolTipString(self, *args, **kwargs)
        self.Button.SetToolTipString(*args, **kwargs)
        self.FileInfo.SetToolTipString(*args, **kwargs)


class Text(Widget):
    """ Allows the user to type in anything """
    def __init__(self, parent, default=None, value="", **kwargs):
        Widget.__init__(self, parent, **kwargs)

        if not default is None:
            self._has_default = True

        self.parent = parent
        self._default = default

        Sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.Text = wx.TextCtrl(self, -1, value)
        
        # Looks ugly without a large border
        Sizer.Add(self.Text, 1, flag=wx.CENTER|wx.RIGHT, border=20)

        self.SetSizer(Sizer)

        self.Bind(wx.EVT_TEXT, self.OnEdit)
        
    def SetValue(self, value):
        try:
            if value:
                if not isinstance(value, basestring):
                    value = " ".join(value)
                self.Text.ChangeValue(value)
            else:
                self.Text.ChangeValue("")
        except ValueError:
            raise SetValueError(self.parent.name, value)

    def GetValue(self):
        #if self.Has_Changed: # Hope this does not break anything 
        return self.Text.GetValue()

    def SetToolTipString(self, *args, **kwargs):
        wx.Panel.SetToolTipString(self, *args, **kwargs)
        self.Text.SetToolTipString(*args, **kwargs)


class VariableLines(Widget):
    """
    Widget for <option_name_I_do_not_remember> (The option which can occur on more than one line in the input file and list components).
    """
    def __init__(self, parent, choices=["1", "2", "3"], **kwargs):
        Widget.__init__(self, parent, **kwargs)
        
        self.parent = parent
        self.top_sizer = None
        self.choices = [""] + choices

        self.choice_ctrls = []
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.SetSizer(self.sizer)

        self.AddChoiceCtrl()

        self.Bind(wx.EVT_CHOICE, self._OnEdit)

    def AddChoiceCtrl(self, ind=0):
        obj = wx.Choice(self, -1, choices=self.choices)
        obj.SetSelection(ind)
        self.choice_ctrls.append(obj)
        self.sizer.Add(obj, 0, wx.CENTER|wx.BOTTOM, 3)
        
        self.UpdateLayout()
    
    def _OnEdit(self, event):
        self._edit = True
        self._changed = True
        self.parent.OnEdit()

        ind = [obj.GetSelection() != 0 for obj in self.choice_ctrls]

        if all(ind):
            self.AddChoiceCtrl()
        else:
            count = len([i for i in ind if i == 0])
            if count == 1:
                return
                
            new_ctrls = []
            for i, v in enumerate(ind):
                obj = self.choice_ctrls[i]
                if not v and count > 1:
                    self.sizer.Remove(obj)
                    del obj
                    count -= 1
                else:
                    new_ctrls.append(obj)
            self.choice_ctrls = new_ctrls
        self.UpdateLayout()

    def UpdateLayout(self):
        if not self.top_sizer:
            p = self.parent
            c = 0
            while not(hasattr(p, "is_page_panel")) and c < 20:
                p = p.GetParent()
                c += 1
            self.top_sizer = p.GetContainingSizer()
        if self.top_sizer:
            self.top_sizer.Layout()
    
    def Clear(self):
        self._changed = False
        for obj in self.choice_ctrls[:-1]:
            self.sizer.Remove(obj)
            del obj
        self.choice_ctrls = [self.choice_ctrls[-1]]

        self.UpdateLayout()

    def SetValue(self, a):
        self.Clear()
        if isinstance(a, basestring):
            a = [a]
        
        for i in a:
            self.AddChoiceCtrl(ind=self.choices.index(i))

    def GetValue(self):
        data = []
        for obj in self.choice_ctrls:
            i = obj.GetSelection()
            if i:
                data.append(self.choices[i])
        return data

    def AddValue(self, vals):
        # Error handling!
        if isinstance(vals, basestring):
            vals = [vals]
        try:
            inds = [self.choices.index(val) for val in vals]
        except ValueError:
            raise SetValueError(self.parent.name, val)
        
        if len(inds) == 0: return
        
        if self.choice_ctrls[-1].GetSelection() == 0:
            self.choice_ctrls[-1].SetSelection(inds[0])
            inds = inds[1:]

        for ind in inds:
            self.AddChoiceCtrl(ind)

        self.AddChoiceCtrl()


class TextInfo(wx.Panel):
    def __init__(self, parent, text="", border=0, filename=""):
        wx.Panel.__init__(self, parent)
        self.text = wx.StaticText(self, -1, text, style=wx.ALIGN_CENTRE)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.text, 1, wx.EXPAND)
        
        self.SetSizer(sizer)

        if filename:
            f = open(filename, "r")
            self.text.SetLabel(f.read())
            f.close()

class HtmlWindow(wx.html.HtmlWindow):
    def OnLinkClicked(self, url):
        webbrowser.open(url.GetHref(), new=2, autoraise=1)

class Html(wx.Panel):
    def __init__(self, parent, text="", border=0, filename="", **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)

        self.html = HtmlWindow(self)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.html, 1, wx.EXPAND)
        
        self.SetSizer(sizer)
        

        if filename:
            self.html.LoadFile(filename)
        if text:
            self.html.SetPage(text)
            
    def SetValue(self, value):
        self.html.SetPage(value)
    def GetValue(self):
        return None
    def Clear(self):
        pass
            
class SaveButtons(wx.Panel):
    """SaveButtons in general page
    
    The design of the savebuttons.
    """
    def __init__(self, parent, ext_format="*.*", lbl2="Save",
                 lbl1="Save As", flags=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT):
        wx.Panel.__init__(self, parent)
        
        self._ext_format = ext_format
        self._flags = flags

        Sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text = wx.TextCtrl(self)
        self.save_as_button = wx.Button(self, -1, lbl1)
        self.save_button = wx.Button(self, -1, lbl2)

        Sizer.Add(self.text, 1, wx.RIGHT, 5)
        Sizer.Add(self.save_as_button, 0, wx.RIGHT, 5)
        Sizer.Add(self.save_button, 0, wx.RIGHT, 5)

        self.SetSizer(Sizer)

        self.save_button.Bind(wx.EVT_BUTTON, self.OnSave)
        self.save_as_button.Bind(wx.EVT_BUTTON, self.OnSaveAs)

    def OnSave(self, event):
        if self.GetValue():
            event.Skip()
        else:
            self.OnSaveAs(event)

    def OnSaveAs(self, event):
        path = os.path.split(self.GetValue())
        val = FileDialog(ext_format=self._ext_format,
                         flags=self._flags,
                         defaultDir=path[0], defaultFile=path[1])
        if val:
            self.SetValue(val)
            event.Skip()
        
    def GetValue(self):
        return self.text.GetValue()

    def SetValue(self, value):
        self.text.SetValue(value)

    def Clear(self):
        self.SetValue("")

class RunUvspecProcess(threading.Thread):
    def __init__(self, window, input_file, output_file):
        threading.Thread.__init__(self)
        self.window = window
        self.input_file = input_file
        self.output_file = output_file
    def run(self):
        wx.CallAfter(self.window.AddToLog,
                     "Running uvspec.\nInput file: %s \nOutput file: %s\n" % \
                         (self.input_file, self.output_file))

        exit_value = 0

        process_in = open(self.input_file, "r")
        process_out = open(self.output_file, "w")
        process_err = open(self.output_file + ".err", "w")
        
        try:
            process = subprocess.Popen("uvspec", stdin=process_in,
                                    stdout=process_out, stderr=process_err)
        except OSError:
            wx.CallAfter(self.window.AddToLog, "\n!!!Could not find UVSPEC in your PATH!!!\n\nIt is strongly recommended to change your PATH and restart the GUI!\n")
            time.sleep(5)
            wx.CallAfter(self.window.UvspecFinished, "cancel")
            wx.CallAfter(self.window.Close)

        data = ""

        wx.CallAfter(self.window.AddToLog, "Uvspec has started.\n")
        while process.poll() == None: # None inidicates that uvspec is not yet finished
            time.sleep(0.1)
            if self.window.should_thread_die: # not sure if this is thread safe
                process.terminate()
                process_in.close()
                process_out.close()
                wx.CallAfter(self.window.UvspecFinished, "cancel")
                wx.CallAfter(self.window.Close)
                return
        
        process_err.close()
        process_err = open(process_err.name, "r")
    
        data = process_err.read()
        if data.upper().find("ERROR") != -1: # catch errors with return code 0
                exit_value = "Error"

        process_in.close()
        process_out.close()
        process_err.close()

        wx.CallAfter(self.window.AddToLog, data)
        
        if not exit_value or process.poll() != 0: # correcting uvspec returncodes
            exit_value = process.poll()
        wx.CallAfter(self.window.UvspecFinished, exit_value)
        
class RunUvspec(wx.Frame):
    def __init__(self, parent, input_file, output_file):
        size = 400, 200
        wx.Frame.__init__(self, parent,
                          style=wx.DEFAULT_FRAME_STYLE^wx.SYSTEM_MENU,
                          size=size)

        self.parent = parent
        self.input_file = input_file

        self.finished = False
        self.should_thread_die = False

        wx.SetCursor(wx.StockCursor(wx.CURSOR_WAIT))

        self.Log = wx.TextCtrl(self, -1, "",
                               style=wx.TE_MULTILINE|wx.TE_READONLY)
        self.Cancel = wx.Button(self, -1, "Stop run")

        self.PretendIAmWorking = wx.Gauge(self,)
        self.PretendIAmWorking.Pulse()

        MainSizer = wx.BoxSizer(wx.VERTICAL)
        HorSizer = wx.BoxSizer(wx.HORIZONTAL)

        MainSizer.Add(self.Log, 5, wx.CENTER|wx.EXPAND)
        HorSizer.Add(self.PretendIAmWorking, 2, wx.EXPAND|wx.ALL|wx.CENTER,
                     border=10)
        HorSizer.Add(self.Cancel, 1, wx.RIGHT|wx.CENTER, border=10)

        MainSizer.Add(HorSizer, 0, wx.EXPAND)

        self.SetSizer(MainSizer)
	self.Layout()

        self.Centre()
        self.SetSizeHints(*size)
        self.Show()

        
        self.parent.Enable(False)

        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Bind(wx.EVT_BUTTON, self.OnCancel)

        self.thread = RunUvspecProcess(self, input_file, output_file)
        self.thread.start()

        
    def OnClose(self, event):
        # veto if uvspec is still running
        if self.finished:
            self.parent.Enable(True)
            event.Skip()
    def OnCloseButton(self, event):
        self.Close()
    def OnCancel(self, event):
        self.should_thread_die = True
    def AddToLog(self, text):
        self.Log.AppendText(text)
        
    def UvspecFinished(self, exit_value):
        self.finished = True
        self.AddToLog("Uvspec finished with exit value %s\n" % exit_value)
        self.PretendIAmWorking.SetValue(100)
        wx.SetCursor(wx.StockCursor(wx.CURSOR_ARROW)) # standard cursor
        self.parent.Enable(True)
        del self.thread
        os.remove(self.input_file)

#        if exit_value == 0: # success
#            self.Close()
            
        self.Cancel.SetLabel("Close")
        self.Bind(wx.EVT_BUTTON, self.OnCloseButton)
        
#class ScrolledWindowFix(wx.lib.scrolledpanel.ScrolledPanel):
class ScrolledWindowFix(wx.ScrolledWindow):
    """
    Prevent the scroll windows rom jumping up and down, when a child past
    the window border is focused. 
    """
    
    def __init__(self, *args, **kwargs):
        #wx.lib.scrolledpanel.ScrolledPanel.__init__(self, *args, **kwargs)
        wx.ScrolledWindow.__init__(self, *args, **kwargs)
        self.EnableScrolling(False, True)
    
    def OnChildFocus(self, *args, **kwargs):
        """ Overwritten to prevent the jumping when a child is focused. """
        pass


class Example(wx.Panel):
    def __init__(self, parent, title, info, img_path, img,
                 example_path, example_dir):
        wx.Panel.__init__(self, parent, )#style=wx.SUNKEN_BORDER)

        self.example_dir = example_dir
        self.txt = info
        self.img_path = img_path
        self.title = title
        self.example_path = example_path

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        txt = wx.TextCtrl(self, -1, info, style=wx.TE_MULTILINE|wx.TE_READONLY)
        image = wx.StaticBitmap(self, -1, img)

        main_sizer.Add(image, 0, wx.CENTER|wx.BOTTOM, border=10)
        main_sizer.Add(txt, 1, wx.EXPAND)

        self.SetSizer(main_sizer)
        

class ExampleList(wx.Treebook):
    def __init__(self, parent, example_list):
        wx.Treebook.__init__(self, parent, -1, )

        img_list = wx.ImageList(*variables.example_thumb_size)
        img_id = 0

        for dirp, info_path, inp_path, img_path_thumb, img_path_big \
                in example_list:
            thumb = wx.Bitmap(img_path_thumb)
            thumb.SetSize((60, 60))
            img = wx.Bitmap(img_path_big)

            f = open(info_path)
            title = f.readline()[7:].strip() # Starts with 'Title: ', ends with a newline
            txt = f.read()
            f.close()

            obj = Example(self, title, txt, img_path_big, img, inp_path, dirp)

            self.AddPage(obj, title, imageId=img_id)

            img_list.Add(thumb)
            img_id += 1
            
        self.AssignImageList(img_list)


class Logo(wx.StaticBitmap):
    def __init__(self, *args, **kwargs):
        wx.StaticBitmap.__init__(self, *args, **kwargs)
        self.SetToolTipString("Visit esa.")
        self.Bind(wx.EVT_LEFT_DOWN, self.OnClick)

    def OnClick(self, event):
        webbrowser.open("http://www.esa.int")


class DocumentationWindow(wx.Frame):
    def __init__(self, parent, doc_dir="."):
        wx.Frame.__init__(self, parent, title="Documentation",
                          style=wx.FRAME_TOOL_WINDOW|wx.FRAME_FLOAT_ON_PARENT |
                          wx.DEFAULT_FRAME_STYLE,
                          pos=variables.documentation_window_position)

        self.display = wx.html.HtmlWindow(self)
        self.DocDir = doc_dir

        self.Show(False)

        self.Bind(wx.EVT_CLOSE, self.OnClose)
        
    def GetDocumentation(self, option):
        self.Show(True)
        self.SetTitle("UVspec Graphical Interface: " + option)
        p = os.path.join(self.DocDir, option + ".html")
        if os.path.exists(p):
           self.display.LoadFile(p)
        else:
           self.display.SetPage("<html><head></head><body>No documentation for %s.</body></html>" % repr(option))
    def OnClose(self, event):
        self.Show(False)


def test():
    print "Oooops! The tests probably does not work anymore..."
    app = wx.PySimpleApp()
    size = 600, 500
    frame = wx.Frame(None, -1, title="Widgets test", size=size)
    panel = wx.Panel(frame)

    debug_button = wx.Button(panel, -1, "Debug")
    
    GridSizer = wx.GridSizer(rows=11)
    
    for obj in [  # Add every new item to this list
                SpinTogether(panel, spinlabels=["loooong", "3"], increment=0.01),
                Spinner(panel, spinlabels=["loooong", "3"], increment=0.01),
                OpenButton(panel, filename="/", buttonlabel="Delete!"),
                OnOff(panel),
                Choice(panel, choices=["One", "Two", "Three"]),
                #InfoButton,
                ListAndFile(panel, choices=[unichr(0x03c0), unichr(0x2107)]),
                ExpandableOnOff(panel, choices=["excellent", "beautiful", "tremendous", "superb"]),
                Text(panel, value="Hello world!"),
                TextInfo(panel, text="Hei verden!"),
                Html(panel, text="<html><head></head><body><h1>Hallo!</h1>"+
                     "Hallo Welt!</body></html>"),
                #SaveRunButtons
        ]:
        t = wx.StaticText(panel, -1, obj.__class__.__name__)
        GridSizer.Add(t, 0, flag=wx.ALL|wx.EXPAND|wx.CENTER, border=5)
        GridSizer.Add(obj, 1, flag=wx.ALL|wx.EXPAND|wx.CENTER, border=5)
    GridSizer.Add(debug_button, 1)

    debug_button.Bind(wx.EVT_BUTTON, lambda x: GridSizer.Layout())
    
    panel.SetSizer(GridSizer)
    
    frame.Centre()
    frame.SetSizeHints(*size)
    frame.Show(True)
    app.MainLoop()

if __name__ == "__main__":
    test()


