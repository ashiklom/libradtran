import sys
import os
sys.path.append("../src_py")

from writeLex import loadOptionsHelper
import variables

"""
The create_list is used by the GUI to create the layout.
The list is a list of dictionaries with these entries:
title: Title of the page in the NoteBook display.
options: A list of dictionaries, each dictionary in the list will make an option.

The entries in the dictionaries for options are:
name: Name of the option.
class: The classname of the option as a string.
args: Arguments for the class when it is initialized.

Please run create.py after editing the create_list to
check that there are no missing brackets of commas.
"""

def parseOption(opt):       
    parents = opt["parents"]
    if parents == "" or parents == []:
        opt["parents"] = ["uvspec"]

    return opt

def makePages(load_function):
    from_lex = load_function()
    pages = []

    mthree = variables.threedmystic_in_GUI

    for group in from_lex:
        try:
            name = group[0]["group"] # Please don't give me empty groups
        except AttributeError:
            name = group.group_name
	if name != "Special":
        	pages.append({"title":  name.title(),
                     "options": [parseOption(option) for option in group \
                                     if (mthree or not option.threedmystic) and len(option.gui_inputs) != 0 and option.showInGui]
                     })

    return pages

#from test_options import loadTestOptions
#create_list = makePages(loadTestOptions)
create_list = makePages(loadOptionsHelper)

def get_example_content(dirp):
    r = [dirp]
    for i in ("info.txt", "inp.inp", "thumb.jpg", "big.jpg"):
        p = os.path.join(dirp, i)
        if os.path.exists(p) and not os.path.isdir(p):
            r.append(p)
        else:
            return None
    return r

def create_example_list():
#    import shutil
    global example_list
    example_list = []
    examples_dir = pe("examples/GUI")
    dirs = os.listdir(examples_dir)
    dirs.sort()
    for dirp in dirs:
        dirp = os.path.join(examples_dir, dirp)
# CE: Commented out, inp.inp is now copied by the Makefile
#    for dirname in dirs:
#        dirp = os.path.join(examples_dir, dirname)
        if not os.path.isdir(dirp):
            continue
#	examplefile=os.path.join(dirp,dirname+'.INP')
#	destfile=os.path.join(dirp,'inp.inp')
#	try:
#	    shutil.copyfile(examplefile,destfile)
#	except Exception, e:
#	    print 'Error copying example file %s to %s' %(examplefile,destfile)
        r = get_example_content(dirp)
        if not r:
            print "Error fetching example in %s" % dirp
            continue
        example_list.append(r)

def _test():
    pages = {}
    for page in create_list:
        pages[page["title"]] = page
    pkeys = pages.keys()
    pkeys.sort()
    for page in pkeys:
        indent = 0
        options = pages[page]["options"]
        print pages[page]["title"] + " (%i)" % len(options)
        indent = 4
        for option in options:
            print " "*indent + option["name"] + \
                (" plot" if (not option["plot"] is None) else "")

if __name__ == "__main__": _test()
    
