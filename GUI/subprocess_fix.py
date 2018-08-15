"""
Add a terminate method to python versions below 2.6.
Loosely based on the code in subprocess.py found in the 2.6 release of python.
http://www.python.org/download/releases/2.6.4/
"""
from subprocess import *

try:
   Popen.terminate
except AttributeError:
    import os
    import signal
    def terminate(self):
        """Terminate the process with SIGTERM
        """
        os.kill(self.pid, signal.SIGKILL)
    Popen.terminate = terminate