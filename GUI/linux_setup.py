from cx_Freeze import setup, Executable
                               
import config_setup


setup(
	executables = [Executable(config_setup.APP[0], base=None),],
	data_files = config_setup.DATA_FILES,
	**config_setup.META
)
  
 
