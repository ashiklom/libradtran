version = "Best version"
description = "GUI for uvspec."
name = "Uvspec GUI"
author = "Jonas Irgens Kylling"
author_email = "joirkyl@gmail.com"

license = "Opens Source of course!"

APP = ["GUI_main.py"]
DATA_FILES = [ ("resources",
		("resources/how_to.html", 
		"resources/splash.png", 
		"resources/GUI_DOC",
		"resources/Label.png",
		"resources/Selected.png",), ),
		]

META = {"name":name, "version":version, "description":description, "author":author, "author_email":author_email, "license":license}

PY2APP_OPTIONS = {"argv_emulation": True, "iconfile":"resources/mac.icns"}
