import wx

version = "1"

redirect = False
redirect_file = ".libRadtran_GUI.log"


documentation_window = True
tooltips = False

enable_plotting = True

# Main window
main_window_size = 1200, 750
main_window_size_hints = 800, 450
main_window_title = "GUI for libRadtran"

max_column_items = 20 # Maximum number of options in a single column

# Start panel
start_welcome_text = """Welcome to the uvspec Graphical User Interface.
        
To start a new input file, push the New button.\n\n To open an old input file, push the Open button.

Mandatory options are in red, the rest are optional.        
        
Find more information about uvspec and \n libRadtran at www.libradtran.org."""

# Widgets
option_style =  wx.DEFAULT_FRAME_STYLE #wx.BORDER_SIMPLE# is great for debugging!
option_background_colour = 255, -1, -1

documentation_window_position = 0, 30

plot_dpi = 96

# ColumnChoices
column_choices_y_size = 80

# ExampleList
example_thumb_size = 60, 60

maxint = 2**31 - 1 # sys.maxint, but not if python is 64 bits and wx 32 bits.

plot_button_size = 20

# Mysticthreed
threedmystic_in_GUI = True
