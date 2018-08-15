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

import sys
import os
from PIL import Image
from optparse import OptionParser

def resize(filename):
    img = Image.open(filename)
    thumb = img.copy()
    thumb.thumbnail((60, 60), Image.ANTIALIAS)
    img.thumbnail((200, 200), Image.ANTIALIAS)
    thumb.save(os.path.splitext(os.path.basename(filename))[0] + "_thumb.jpg")
    img.save(os.path.splitext(os.path.basename(filename))[0] + "_prev.jpg")

def main():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename", help="Image file to create images with the correct size from.")

    options, args = parser.parse_args()
    if options.filename and not args:
        filename = options.filename
    else:
        filename = args[0]
    resize(filename)

if __name__ == "__main__":
    main()
