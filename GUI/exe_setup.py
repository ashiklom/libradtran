from distutils.core import setup
import py2exe
import config_setup

OPTIONS = {'argv_emulation': True}

setup(
    windows=[
        {"script": config_setup.APP[0], 
         "icon_resources": [(0, "resources/windows.ico")]}],
    data_files= config_setup.DATA_FILES,
    **config_setup.META
)
