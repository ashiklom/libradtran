import wx
import wx.wizard
import os
from widgets import FileDialog, Option
from widgets import pr

class OpenImage(wx.Panel):
    def __init__(self, parent,
                 default_image=pr("resources/default_example_image.png")):
        wx.Panel.__init__(self, parent, -1)
        
        self.pth = ""
        self.default_image = default_image

        self.ImageInfo = wx.StaticText(self, -1, "")
        self.Image = wx.StaticBitmap(self, -1, wx.Bitmap(default_image),
                                     size=(1, 1)) # Fixes sizing problems, somehow
        self.OpenButton = wx.Button(self, -1, "Open", size=(100, -1))

        MainSizer = wx.BoxSizer(wx.VERTICAL)
        MainSizer.Add(self.ImageInfo, 0)
        MainSizer.Add((-1, 5))
        MainSizer.Add(self.Image, 1, flag=wx.EXPAND|wx.CENTER)
        MainSizer.Add((-1, 10))
        MainSizer.Add(self.OpenButton, 0, flag=wx.CENTER|wx.BOTTOM, border=5)
        self.SetSizer(MainSizer)

        self.Bind(wx.EVT_BUTTON, self.OnOpen, self.OpenButton)

    def OnOpen(self, event):
        tmp = FileDialog(ext_format="*.png")
        if tmp:
            self.SetValue(tmp)
    def SetValue(self, pth):
        tmp_img = wx.Image(pth)
        w, h = tmp_img.GetWidth(), tmp_img.GetHeight()
        max_size = 500
        s = float(max_size)/max(w, h) if max(w, h) > max_size else 1
        w = w*s
        h = h*s
        tmp_img.Rescale(w, h, wx.IMAGE_QUALITY_HIGH)
        self.Image.SetBitmap(tmp_img.ConvertToBitmap())
        self.Image.CacheBestSize(self.Image.GetSize())
        self.ImageInfo.Label = pth
        self.pth = pth

    def GetValue(self):
        return self.pth
    def GetBitmap(self):
        return self.Image.GetBitmap()
    def Clear(self):
        self.pth = ""
        self.ImageInfo.Label = ""
        self.Image.SetBitmap(wx.Bitmap(self.default_image)) # Memory?

class DescriptionImage(wx.Panel):
    def __init__(self, parent, txt="", img=""):
        wx.Panel.__init__(self, parent, -1)

        self.Description = wx.TextCtrl(self, -1, txt, style=wx.TE_MULTILINE)
        self.Image = OpenImage(self)

        MainSizer = wx.BoxSizer(wx.VERTICAL)

        middle_sizer= wx.BoxSizer(wx.HORIZONTAL)
        for item, title in ((self.Description, "Description"),
                            (self.Image, "Image")):

            t_sizer = wx.BoxSizer(wx.VERTICAL)
            t_sizer.Add(wx.StaticText(self, -1, title), 0,
                        wx.CENTER|wx.BOTTOM|wx.CENTER, border=5)
            t_sizer.Add(item, 2, wx.EXPAND|wx.CENTER)

            middle_sizer.Add(t_sizer, 1, wx.EXPAND|wx.RIGHT, border=5)

        MainSizer.Add(middle_sizer, 2, wx.RIGHT|wx.LEFT|wx.EXPAND|wx.CENTER, 5)

        self.SetSizer(MainSizer)
        if img:
            self.Image.SetValue(img)
    def SetDescription(self, a):
        self.Description.SetValue(a)
    def SetImage(self, a):
        self.Image.SetValue(a)
    def GetValue(self):
        return self.Description.GetValue(), self.Image.GetBitmap()
    def Clear(self):
        self.Description.SetValue("")
        self.Image.Clear()

class MakeExampleDialog(wx.Dialog):
    """
    Wizard for making examples (is more of a dialouge).
    Not in use.
    """
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, -1, "Make an example file",
                           size=(800, 500))
        self.Panel = wx.Panel(self)
        
        self.done = False
        self.result = None

        MainSizer = wx.BoxSizer(wx.VERTICAL)

        self.TitleDescription = wx.TextCtrl(self.Panel, -1, "")
        self.Description = wx.TextCtrl(self.Panel, -1, "", style=wx.TE_MULTILINE)
        self.Image = OpenImage(self.Panel)
        self.DoneButton = wx.Button(self.Panel, -1, "Done")

        MainSizer.Add((-1, 10))
        top_sizer = wx.BoxSizer(wx.HORIZONTAL)
        top_sizer.Add(wx.StaticText(self.Panel, -1, "Name"), 0,
                      flag=wx.RIGHT|wx.RIGHT|wx.CENTER, border=5)
        top_sizer.Add(self.TitleDescription, 2,
                      flag=wx.CENTER|wx.EXPAND)
        MainSizer.Add(top_sizer, 0,
                      wx.EXPAND|wx.CENTER|wx.LEFT|wx.RIGHT, border=5)
        MainSizer.Add((-1, 10))

        middle_sizer = wx.BoxSizer(wx.HORIZONTAL)
        for item, title in ((self.Description, "Description"),
                            (self.Image, "Image")):

            t_sizer = wx.BoxSizer(wx.VERTICAL)
            t_sizer.Add(wx.StaticText(self.Panel, -1, title), 0,
                        wx.CENTER|wx.BOTTOM|wx.CENTER, border=5)
            t_sizer.Add(item, 2, wx.EXPAND|wx.CENTER)

            middle_sizer.Add(t_sizer, 1, wx.EXPAND|wx.RIGHT, border=5)

        MainSizer.Add(middle_sizer, 2, wx.RIGHT|wx.LEFT|wx.EXPAND|wx.CENTER, 5)
        MainSizer.Add((-1, 10))
        MainSizer.Add(self.DoneButton, 0,
                      flag=wx.CENTER|wx.EXPAND|wx.RIGHT|wx.LEFT, border=25)
        MainSizer.Add((-1, 10))

        self.Panel.SetSizer(MainSizer)

        tmp_sizer = wx.BoxSizer(wx.VERTICAL)
        tmp_sizer.Add(self.Panel, 1, wx.EXPAND)

        self.SetSizer(tmp_sizer)

        self.Center()

        self.Bind(wx.EVT_BUTTON, self.OnDone, self.DoneButton)
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def OnClose(self, event):
        if any((self.TitleDescription.GetValue(),
               self.Image.GetValue(),
               self.Description.GetValue())) and not self.done:
            d = wx.MessageDialog(self, "You have unsaved work.",
                                 style=wx.YES|wx.NO|wx.CENTRE|wx.ICON_EXCLAMATION,
                                 caption="Are you absolutely completly sure that you want to exit?")
            r = d.ShowModal()
            if r == wx.ID_NO:
                event.Veto()
                return
        event.Skip()
            
    def OnDone(self, event):
        name = self.TitleDescription.GetValue()
        image = self.Image.Image.GetBitmap()
        description = self.Description.GetValue()
        for v, e in ((name, "Please specify a name for your example."),
                     (image.IsOk(), "Please specify an image for your example."),
                     (description, "Pleace write a description for your example.")):
            if not v:
                d = wx.MessageDialog(self, e,
                                     style=wx.OK|wx.CENTRE|wx.ICON_ERROR,
                                     caption="Error!")
                d.ShowModal()
                return

        self.done = True
        self.result = name, image, description
        self.Close()
    def ShowModal(self):
        wx.Dialog.ShowModal(self)
        return self.result

def MakeExampleDialogResult(parent):
    m = MakeExampleDialog(parent)
    return m.ShowModal()

class MakeExampleFrame(wx.Frame):
    def __init__(self, ):
        wx.Frame.__init__(self, None, -1, "Example Frame")

        bs = wx.BoxSizer(wx.VERTICAL)
        bs.Add(wx.StaticText(self, -1, "Hello"), 0)
        bs.Add(DescriptionImage(self), 1, wx.EXPAND)
        self.SetSizer(bs)
        
        
        self.Center()
        self.Show()

def _test():
    app = wx.PySimpleApp()
    f = MakeExampleFrame()
    #print MakeExampleDialogResult(None)
    
    app.MainLoop()

if __name__ == "__main__": _test()
