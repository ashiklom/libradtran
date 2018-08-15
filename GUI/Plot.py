import sys

import wx
import widgets
import variables
import math

if variables.enable_plotting:
    import wxmpl
    import numpy as np
    import matplotlib.cm as cm

#widgets.ErrorMessage = widgets.widgets.ErrorMessage
#lambda x: (sys.stdout.write(x + "\n"), sys.stdout.flush())

class open_column_file(file):
    def readline(self, *args):
        l = file.readline(self, *args)
        l = self._parse_line(l)
        while not l:
            l = file.readline(self, *args)
            l = self._parse_line(l)
        return l
    def _parse_line(self, l):
        if not l:
            raise EOFError
        l = l.strip()
        return l[:l.find("#")].strip()

class Transformation:
    def __call__(self, x, ys):
        raise NotImplementedError

class Log(Transformation):
    def __call__(self, x, ys):
        return x, [[math.log(i) for i in y] for y in ys]

class LogLog(Transformation):
    def __call__(self, x, ys):
        return [math.log(i) for i in x], [[math.log(i) for i in y] for y in ys]

class Normalize(Transformation):
    """
    Normalize the ys so that they fit in the same plot frame
    """
    def __call__(self, x, ys):
        nys = []
        for y in ys:
            y = np.asarray(y)
            miny, maxy = min(y), max(y)
            d = max(y) - min(y)
            if d != 0:
                y = (y - miny)/d
            nys.append(y)
        return x, nys

Transformations = {"Log": Log(),
                   "Log-log": LogLog(),
                   "Normalize": Normalize(),
                   }

class PlotButton(wx.BitmapButton):
    def __init__(self, parent, filegetter,
                 plot_type, optional_args={}, **kwargs):
        self.parent = parent
        wx.BitmapButton.__init__(self, parent, -1, wx.ArtProvider.GetBitmap(
                wx.ART_MISSING_IMAGE, size=(10, 10)), **kwargs)
        self._plotter = plot_type(**optional_args)
        self._filegetter = filegetter
        self.Bind(wx.EVT_BUTTON, self.OnPlot)

    def OnPlot(self, event):
        #print type(self._filegetter)
        path = self._filegetter.GetFileValue()
        #print path
        if not path:
            msg = "Please choose a file to plot!"
            widgets.ErrorMessage(msg)
            return
        self._plotter.ShowPlotWindow(self.parent, path)


class PlotWindow(wx.Frame):
    def __init__(self, parent, path):
        wx.Frame.__init__(self, parent, -1, title=path)

        self.Panel = wx.Panel(self, -1)

        self.MainSizer = wx.BoxSizer(wx.VERTICAL)

        # Layout
        self.MainSizer.Add((-1, -1), 1)
        
        self.PlotButton = wx.Button(self.Panel, -1, "Plot")

        self.MainSizer.Add((-1, -1), 1)
        self.MainSizer.Add(self.PlotButton, 0,
                       wx.CENTER|wx.ALL^wx.TOP, border=5)
        self.Panel.SetSizer(self.MainSizer)


class Plot():
    def __init__(self):
        self.PlotOption = {}

    def ShowPlotWindow(self, parent, path):
        self._path = path
        self.frame = PlotWindow(parent, path)
        self.LayoutFrame()
        
        self.frame.MainSizer.Layout()
        self.frame.Center()

        self.frame.Show()
        self.frame.Bind(wx.EVT_BUTTON, self.OnPlot,
                        self.frame.PlotButton)
    def LayoutFrame(self):
        raise NotImplementedError
    def Plot(self, *args, **kwargs):
        raise NotImplementedError


class Plot2D(Plot):
    def __init__(self, column_names=["x", "y"]):
        Plot.__init__(self)
        self._l = len(column_names)
        self.column_names = column_names
        assert(self._l == len(column_names))

    def LayoutFrame(self):
        GridSizer = wx.GridSizer(len(self.column_names) + 1, 3)
            
        GridSizer.Add((-1, -1))

        self._transformations = "Log", "Log-log", "Normalize"

        for axis in ("x", "y"):
            GridSizer.Add(wx.StaticText(self.frame.Panel, -1, axis),
                          flag=wx.ALIGN_RIGHT)
            self.PlotOption[axis] = []
        
        for name in self.column_names:
            text = wx.StaticText(self.frame.Panel, -1, name)
            GridSizer.Add(text,
                          flag=wx.ALIGN_LEFT)
            for axis, cls in ("x", wx.RadioButton), ("y", wx.CheckBox):
                obj = cls(self.frame.Panel, -1)
                GridSizer.Add(obj, 1, flag=wx.ALIGN_RIGHT)
                self.PlotOption[axis].append(obj)

        self.PlotOption["x"][0].SetValue(True)
        self.frame.MainSizer.Insert(0, (-1, -1), 0, wx.TOP|wx.BOTTOM, 10)
        self.frame.MainSizer.Insert(1, GridSizer, 0, wx.EXPAND|wx.RIGHT|
                                    wx.LEFT, 5)

        # Transformations
        for name in self._transformations:
            obj = wx.CheckBox(self.frame.Panel, -1, name)
            self.PlotOption[name] = obj
            self.frame.MainSizer.Add(obj, 0, wx.EXPAND|wx.RIGHT|wx.LEFT,
                                 border=5)
            self.frame.MainSizer.Add((-1, -1), 0, wx.BOTTOM, border=5)

        # Check valid columns and disable if the file is too short
        try:
            f = open_column_file(self._path)
            l = f.readline()
            for i in xrange(len(l.split()), len(self.PlotOption["x"])):
                self.PlotOption["x"][i].Enable(False)
                self.PlotOption["y"][i].Enable(False)
            f.close()
                
        except (IOError, EOFError), e:
            msg = "IOError: %s" % e
            widgets.ErrorMessage(msg)
            self.frame.PlotButton.Enable(False)
            self.frame.SetTitle("Error: %s" % self._path)


        w, h = 0, 0
        for obj in self.frame.MainSizer.GetChildren():
            r = obj.CalcMin()
            w = max(r[0], w)
            if r[1] > 0: h += r[1]
            org = obj.GetWindow() if obj.IsWindow() else obj.GetSizer()

        # Adjust size for the title bar
        h += self.frame.GetSize()[1] - self.frame.GetClientSize()[1] \
            + 10

        size = (w, h)
        self.frame.SetSizeHints(*size)
        self.frame.SetSize(size)

    def GetData(self, path, xcol, ycols):
        """
        TODO: Rewrite with numpy
        Returns metadata and data
        """
        variables = [[] for i in xrange(self._l)]
        try:
            lcount = 0
            f = open(path)
            for l in f.readlines():
                lcount += 1
                if l.startswith("#"): continue
                l = l.split()
                for i in xrange(min(len(l), self._l)):
                    variables[i].append(float(l[i]))
            f.close()
            
        except (IOError, ValueError), e:
            msg = "Error in file %s at line %d: %s" % (path, lcount, e)
            widgets.ErrorMessage(msg)
        return (({"name":self.column_names[xcol]}, variables[xcol]),
                [({"name": self.column_names[i]}, variables[i]) \
                     for i in ycols])
        
    def OnPlot(self, event):
        errors = []
        # Which columns to plot
        for xcol, button in enumerate(self.PlotOption["x"]):
            if button.GetValue(): break
        ycols = []
        for i, button in enumerate(self.PlotOption["y"]):
            if button.GetValue(): ycols.append(i)

        if not ycols:
            msg = "Please specify a variabele for y!"
            widgets.ErrorMessage(msg)
            return 

        x, ys = self.GetData(self._path, xcol, ycols)

        # Transformations to be applied to the data
        transformations = []
        for tname in self._transformations:
            if self.PlotOption[tname].GetValue():
                transformations.append((tname, Transformations[tname]))

        # Apply transformations
        for tname, t in transformations:
            try:
                xn, ysn = t(x[1], [data for conf, data in ys])
                x = (x[0], xn)
                ys = [(prev[0], data) for prev, data in zip(ys, ysn)]
            except (ValueError, ArithmeticError), e:
                msg = "Error applying transformation %s: %s" % (tname, e)
                errors.append(msg)

        # Plot
        f = PlotFrame(self.frame, -1, title="%s (%d)" % (self._path,
                                                         len(x[1])))
        fig = f.get_figure()
        a = fig.add_subplot(1,1,1)

        xdata = x[1]
        for conf, ydata in ys:
            if len(ydata) != len(xdata):
                errors.append("Column %s and column %s have different sizes." % \
                                 (self.column_names[xcol], conf["name"])
                             )
            else:
                a.plot(xdata, ydata, label=conf["name"])
    
        a.set_xlabel(x[0]["name"])
        a.legend()
        a.axis("tight")

        try:
            f.draw()
        except OverflowError, e:
            msg = "Too much data: %s" % e
            errors.append(msg)
            f.Close()
        else:
            f.Show()
        if errors:
            widgets.ErrorMessage(errors)
            
    

class PlotMap(Plot):
    def LayoutFrame(self):
        self._colour_maps = [(name, cm.get_cmap(name)) \
                                 for name in ("gray", "spectral", "copper",
                                              "coolwarm", "ocean", "terrain",
                                              "bone", "prism", "pink",
                                              "rainbow", "binary", "brg",
                                              "afmhot", "autumn")
                             if cm.get_cmap(name) != None
                             ]
        self._choice = wx.Choice(self.frame.Panel, -1,
                                 choices=[name for name, map_t in \
                                              self._colour_maps])
        self.frame.MainSizer.Insert(1, self._choice, 0, wx.EXPAND|wx.ALL, 10)
    def getData(self, path):
        try:
            f = open(path, "r")
            X, Y, dx, dy = [int(i) for i in f.readline().split()]
            a = np.loadtxt(f)
            f.close()

            # Ensure indexing by zero
            a[:,0] -= 1
            a[:,1] -= 1
            b = np.zeros((Y, X))
            for x, y, z in a:
                x, y = int(x), int(y)
                b[x,y] = z

        except IOError, e:
            widgets.ErrorMessage("IOError for file %s: %s" % (path, e))
            return None
        except (ValueError, IndexError), e:
            widgets.ErrorMessage("Error parsing file %s: %s" % (path, e))
            return None
        else:
            return b
    def OnPlot(self, event):
        selected_map = self._colour_maps[self._choice.GetCurrentSelection()][1]

        Z = self.getData(self._path)
        if Z is None: return
        #cm = self._colour_maps()[1]

        # Plot
        f = PlotFrame(None, -1, title=self._path)
        fig = f.get_figure()
        a = fig.add_subplot(1,1,1)

        im = a.imshow(Z, interpolation='bilinear', cmap=selected_map,
                        origin='lower')

        try:
            f.draw()
        except OverflowError, e:
            msg = "Too much data: %s" % e
            f.Close()
            widgets.ErrorMessage(msg)
        else:
            f.Show()


class PlotBlock(Plot):
    def __init__(self, phi, umu):
        Plot.__init__(self)
        self._nphi = len(phi)
        self._numu = len(umu)

        self._top = ["lambda",
                     "edir",
                     "edn",
                     "eup",
                     "uavgdir",
                     "uavgdn",
                     "uavgup"]

        if self._numu > 0 and self._nphi == 0:
            self._top.append(self._numu)

    def ShowPlotWindow(self, parent, path):
        self._path = path
        self.frame = PlotWindow(parent, path)
        
        self.data = self.GetData(self._path)

        if not self.data:
            self.frame.Close()
            return

        GridSizer = wx.GridSizer(len(self._top) + 2, 3)
            
        GridSizer.Add((-1, -1))

        for axis in ("x", "y"):
            GridSizer.Add(wx.StaticText(self.frame.Panel, -1, axis),
                          flag=wx.ALIGN_RIGHT)
            self.PlotOption[axis] = []
        
        for name in self._top:
            text = wx.StaticText(self.frame.Panel, -1, name)
            GridSizer.Add(text,
                          flag=wx.ALIGN_LEFT)
            for axis, cls in ("x", wx.RadioButton), ("y", wx.CheckBox):
                obj = cls(self.frame.Panel, -1)
                GridSizer.Add(obj, 1, flag=wx.ALIGN_RIGHT)
                self.PlotOption[axis].append(obj)

        self.PlotOption["x"][0].SetValue(True)
        self.frame.MainSizer.Insert(0, (-1, -1), 0, wx.TOP|wx.BOTTOM, 10)
        self.frame.MainSizer.Insert(1, GridSizer, 0, wx.EXPAND|wx.RIGHT|
                                    wx.LEFT, 5)

        self.frame.MainSizer.Add((-1, -1), 1)
        # Add umu control
        tmp_sizer = wx.BoxSizer(wx.HORIZONTAL)
        text = wx.StaticText(self.frame.Panel, -1, "uu(lambda) ")
        self.lambda_choice = wx.Choice(self.frame.Panel, -1,
                                  choices=[str(l) for l in self.data["lambda"]])
        
        self.lambda_choice.SetToolTipString("Wavelength")

        tmp_sizer.Add(text, 0, wx.RIGHT, 5)
        tmp_sizer.Add(self.lambda_choice)

        self.frame.MainSizer.Add(tmp_sizer, 0, wx.CENTER|wx.ALL, 5)
        
        plot_uu_button = wx.Button(self.frame.Panel, -1, "Plot")
        self.frame.MainSizer.Add(plot_uu_button, 0, wx.CENTER|wx.BOTTOM, 10)
        
        # Layout
        self.SetPerfectSize()

        self.frame.MainSizer.Layout()
        self.frame.Center()

        self.frame.Show()
        self.frame.Bind(wx.EVT_BUTTON, self.OnPlot,
                        self.frame.PlotButton)
        plot_uu_button.Bind(wx.EVT_BUTTON, self.OnPlotUU, plot_uu_button)

    def SetPerfectSize(self):
        w, h = 0, 0
        for obj in self.frame.MainSizer.GetChildren():
            r = obj.CalcMin()
            w = max(r[0], w)
            if r[1] > 0: h += r[1]
            org = obj.GetWindow() if obj.IsWindow() else obj.GetSizer()

        # Adjust for the title bar
        h += self.frame.GetSize()[1] - self.frame.GetClientSize()[1] \
            + 10

        #print self.frame.GetSize()
        #print self.frame.GetClientSize()
        
        size = (w, h)
        self.frame.SetSizeHints(*size)
        self.frame.SetSize(size)

    def GetNUmu(self):
        return self.numu_ctrl.GetValue()
    def GetNPhi(self):
        return self.nphi_ctrl.GetValue()

    def GetData(self, path):
        errors = []
        try:
            f = open(path)
            lines = f.readlines()
            f.close()
        except IOError, e:
            msg = "Error reading file '%s'." % (path,)
            widgets.ErrorMessage(msg)
            return

        numu = -1
        nphi = -1

        umu = []
        phi = []
        top = [[] for i in xrange(len(self._top))]
        u0u = []
        uu = []
        u0u_block = []
        uu_block = []

        lcount = 1
        n = 0
        dn = 0
        block_size = 0
        first_block = True

        try:
            while n < len(lines):
                line = lines[n]
                if not line.strip() or line.startswith("#"):
                    n += 1
                    continue

                line = np.asarray([float(i) for i in lines[n].split()])

                if (block_size == 0 and dn == 0) or (block_size and \
                                                         (dn % block_size == 0)):
                    # In top
                    assert len(line) == len(top), "Wrong number of values!"
                    for a, b in zip(top, line):
                        a.append(b)
                elif (block_size == 0 and dn == 1) or block_size and \
                        (dn % block_size == 1):
                    # In phi or n * 1 blcok
                    if first_block:
                        nphi = len(line)
                        phi = line
                    elif first_block:
                        assert False, "Not implemented"
                else:
                    # In n*m blcok
                    # Guess end of block assuming wavelengths are above 2 nm
                    if first_block and (line[0] > 2.0 or len(line) != nphi + 2):
                        first_block = False
                        numu = len(umu)
                    
                        block_size = 2 + numu

                        u0u = [np.asarray(u0u_block)]
                        uu = [np.asarray(uu_block)]

                        uu_block = []
                        u0u_block = []

                        continue

                    assert len(line) == nphi + 2, "Wrong number of values."

                    if first_block:
                        umu.append(line[0])

                    u0u_block.append(line[1])
                    uu_block.append(line[2:])
                
                n += 1
                dn += 1
                if not first_block and dn % block_size == 0:
                    u0u.append(np.asarray(u0u_block))
                    uu.append(np.asarray(uu_block))

                    uu_block = []
                    u0u_block = []

        except ValueError, e:
            msg = "Value error: line %d: %s." % (n, repr(e))
            widgets.ErrorMessage(msg)
            self.frame.PlotButton.Enable(False)
            self.frame.SetTitle("Value error: %s" % self._path)
        except AssertionError, e:
            msg = "Parsing error: line %d: %s." % (n, e)
            widgets.ErrorMessage(msg)
            self.frame.PlotButton.Enable(False)
            self.frame.SetTitle("Parsing error: %s" % self._path)
        except IndexError, e:
            msg = "Parsing error: line %d: %s." % (n, e)
            widgets.ErrorMessage(msg)
            self.frame.PlotButton.Enable(False)
            self.frame.SetTitle("Parsing error: %s" % self._path)
        else:
            if first_block:
                u0u = [u0u_first]
                uu = [uu_first]

            top = np.asarray(top)
            u0u = np.asarray(u0u)
            uu = np.asarray(uu)
            
            d = dict(zip(self._top, top))
            d.update(dict((("umu", umu),
                           ("phi", phi),
                           ("u0u", u0u),
                           ("uu", uu)
                           )))
            return d

        return None

    def OnPlot(self, event):
        #r = self.GetData(self._path)
        #if r is None: return
        r = self.data

        xaxis = ""
        for i, button in enumerate(self.PlotOption["x"]):
            if button.GetValue(): xaxis = self._top[i]
    
        yaxis = []
        for i, button in enumerate(self.PlotOption["y"]):
            if button.GetValue(): yaxis.append(self._top[i])

        # Plot
        f = PlotFrame(None, -1, title=self._path)
        fig = f.get_figure()
        a = fig.add_subplot(1, 1, 1)
        for y in yaxis:
            #print r[y]
            a.plot(r[xaxis], r[y], label=y)

        a.legend()
    
        try:
            f.draw()
        except OverflowError, e:
            msg = "Too much data: %s." % e
            widgets.ErrorMessage(msg)
            f.Close()
        else:
            f.Show()
    
    def OnPlotUU(self, event):
        s = self.lambda_choice.GetCurrentSelection()

        f = PlotFrame(None, -1, title="%s (wavelength = %s)" % 
                      (self._path, str(self.data["lambda"][s])))
        fig = f.get_figure()
        a = fig.add_subplot(1, 1, 1)
        
        for i, phi in enumerate(self.data["phi"]):
            a.plot(self.data["umu"], self.data["uu"][s][:,i],
                   label="phi=%f" % phi)

        a.legend()
        try:
            f.draw()
        except OverflowError, e:
            msg = "Too much data: %s." % e
            widgets.ErrorMessage(msg)
            f.Close()
        else:
            f.Show()


class PlotFrame(wxmpl.PlotFrame):
    def __init__(self, parent, id_n, **kwargs):
        wxmpl.PlotFrame.__init__(self, parent, id_n, **kwargs)
        sizer = self.GetSizer()

        self.figures = []

        self.unzoomButton = wx.Button(self, -1, "unzoom")
        sizer.Add(self.unzoomButton, 0, wx.EXPAND|wx.TOP, 10)
        self.Fit()

        self.Bind(wx.EVT_BUTTON, self.onUnzoom, self.unzoomButton)
    def get_figure(self, **kwargs):
        fig = wxmpl.PlotFrame.get_figure(self, **kwargs)
        self.figures.append(fig)
        return fig
    def onUnzoom(self, event):
        for fig in self.figures:
            for a in fig.axes:
                a.axis("tight")
        self.draw()


def TerminalPlot(X, Y, title=""):
    i = 0
    for conf, ydata in Y:
        TerminalPlotHelper(X[1], ydata, title + " %i" % i)
        i += 1

def TerminalPlotHelper(X, Y, title):
    import curses, time
    excep = ""
    print title
    try:
        stdscr = curses.initscr()
        rows, cols = stdscr.getmaxyx()
        rows -= 2
        curses.endwin()

        scr = [[" " for a in xrange(cols)] for b in xrange(rows)]

        minx, maxx = min(X), max(X)
        miny, maxy = min(Y), max(Y)
        for x, y in zip(X, Y):
            x, y = float(x), float(y)
            x = min(int((x - minx)/(-minx + maxx)*cols), cols-1)
            y = min(int((y - miny)/(-miny + maxy)*rows), rows-1)
            scr[y][x] = "*"
            #stdscr.addch(y, x, "*")
            #print x, y
        #print scr
        scr.reverse()
        print "\n".join(["".join(i) for i in scr])
        sys.stdout.flush()
    finally:
        pass

class FilegetterTest():
    def __init__(self, data):
        self._it = iter(data)
    def GetFileValue(self):
        return self._it.next()

def _test():
    app = wx.PySimpleApp()
    f = wx.Frame(None)
    d = FilegetterTest(
        ["/Users/jonas/tmp.out.txt"] + \
        ["../examples/GUI/lbl_O2A/lbl_O2A.OUT"] +\
            ["../examples/UVSPEC_MC_ELEV.INP"] +\
            ["../examples/KTROP.UVSPEC"] + \
            ["../data/solar_flux/atlas3"]*3 + ["empty.txt"] \
            + ["", "asd", "../GUI/"]
        )
    
    b = PlotButton(f, d, Plot2D, {"column_names":["Long name",
                                                    "even loooooooooooooooonger",
                                                    "more"],
                                  },
                   pos=(10, 10)
                   )
    b = PlotButton(f, d, PlotMap,
                   pos=(10, 50)
                   )
    f.Center()
    f.Show()
    app.MainLoop()

def _test2():
    X = range(100)
    TerminalPlot(X, [x**2 for x in X])

def _test3():
    app = wx.PySimpleApp()
    
    p = PlotBlock(phi = [], umu = [])
    p.ShowPlotWindow(None, "/Users/jonas/libRadtran/trunk/examples/ANGRES_RADDIS_1_tmp.DAT")

    app.MainLoop()
    
if __name__ == "__main__": _test3()
