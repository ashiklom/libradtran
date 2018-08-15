# sdoc2.awk
#                                                                      
# sdoc2 - Automatic source documentation                           
#                                                                      
# Author:                                                              
#   Bernhard Mayer                                                 
#   Lehrstuhl fuer Experimentelle Meteorologie
#   Ludwig-Maximilians-Universitaet Muenchen
#   Theresienstrasse 37
#   80333 Muenchen
#   Germany
#                                                                      
# --------------------------------------------------------------------
# Copyright (C) 2010 Bernhard Mayer                              
#                                                                      
# This program is free software; you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation; either version 1, or (at your option)  
# any later version.                                                   
#                                                                      
# This program is distributed in the hope that it will be useful,      
# but WITHOUT ANY WARRANTY; without even the implied warranty of       
# MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the         
# GNU General Public License for more details.                         
#                                                                      
# To obtain a copy of the GNU General Public License write to the      
# Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,   
# USA.                                                           
#---------------------------------------------------------------------


# sdoc2 is a tool for automatic documentation of source code. 
# It is a very very downgraded version of sdoc2 which only knows
# two markers to start and end a comment and extracts the text 
# between these two markers to stdout. 
# It may be used for documentation of C, Fortran and AWK code.
#
# Usage:
#   gawk -f sdoc.awk -vlanguage="c" <source file>
#
# Valid languages:
# Fortran, C, awk


# replace "_" or "\_" with "\_" (LaTeX doesn't like "_"
# while we don't like "\_" in the code) in \key{...}
function escape_underscore(key, string)
{
  n=split (string,a,"\\\\"key"[ ]*{");
  string=a[1];
  for (i=2;i<=n;i++) {
    m=split (a[i],b,"}");
    for (j=3;j<=m;j++)
      b[2]=b[2]"}"b[j];
    
    gsub("[\\\\]?_","\\_",b[1]);
    string=string"\\"key"{"b[1]"}"b[2];
  }
  return string;
}
# Remove all command strings from a comment  @fi@ @if@
function cleanup(string)
{

  # clear string, if only ****
  match(string, /\*+/);
  if (RLENGTH==length(string)) {
    string="";
    return string;
  }

  # 'copying' comment inclusion
  gsub("<lpdoc>","",string);
  gsub("</lpdoc>","",string);
#  gsub("@ci@","",string);
#  gsub("@ic@","",string);


  # 'overview' source inclusion
  gsub("@os@","",string);
  gsub("@so@","",string);

  # 'overview' comment inclusion
  gsub("@oc@","",string);
  gsub("@co@","",string);
  

  # 'functions' source inclusion
  gsub("@fs@","",string);
  gsub("@sf@","",string);

  # 'functions' comment inclusion
  gsub("@fc@","",string);
  gsub("@cf@","",string);

  # 'functions' inclusion
  gsub("@fi@","",string);
  gsub("@if@","",string);

  # 'definitions' source inclusion
  gsub("@ds@","",string);
  gsub("@sd@","",string);

  # 'definitions' comment inclusion
  gsub("@dc@","",string);
  gsub("@cd@","",string);


  # user defined sections:
  gsub(/@[0-9][0-9]c@/,"",string);
  gsub(/@c[0-9][0-9]@/,"",string);
  gsub(/@[0-9][0-9]_[0-9][0-9]c@/,"",string);
  gsub(/@c[0-9][0-9]_[0-9][0-9]@/,"",string);

  gsub(/@[0-9][0-9]s@/,"",string);
  gsub(/@s[0-9][0-9]@/,"",string);
  gsub(/@[0-9][0-9]_[0-9][0-9]s@/,"",string);
  gsub(/@s[0-9][0-9]_[0-9][0-9]@/,"",string);

  gsub(/@[0-9][0-9]i@/,"",string);
  gsub(/@i[0-9][0-9]@/,"",string);
  gsub(/@[0-9][0-9]_[0-9][0-9]i@/,"",string);
  gsub(/@i[0-9][0-9]_[0-9][0-9]@/,"",string);

  gsub(/@[0-9][0-9]n .+@/,"",string);
  gsub(/@[0-9][0-9]_[0-9][0-9]n .+@/,"",string);

  # newline character
  gsub("@nl@","\n",string);
  
  # remove blanks before @ifinfo etc.
  gsub(/^[ \t]+@ifinfo/,"@ifinfo",string);
  gsub(/^[ \t]+@end ifinfo/,"@end ifinfo",string);
  gsub(/^[ \t]+@iftex/,"@iftex",string);
  gsub(/^[ \t]+@end iftex/,"@end iftex",string);
                                                                       
  # replace "_" or "\_" with "\_" (LaTeX doesn't like "_"
  # while we don't like "\_" in the code) in 
  # \option and \code
  string = escape_underscore("option", string);
  string = escape_underscore("code", string);
  string = escape_underscore("fcode", string);
  string = escape_underscore("file", string);
  string = escape_underscore("parameter", string);

  return string;
}


BEGIN { 

  # set defaults 
  language="c";

  # C comments start with "/*" and end with "*/" 
  if(language=="c") {
    comment_start="\\/\\*";
    comment_end  ="\\*\\/";
  }

  # Fortran comments start with "*" and end with "\n" 
  if(language=="fortran") {
    comment_start="\\*"; 
    comment_end  ="\\n";
  }

  # awk comments start with "#" and end with "\n" 
  if(language=="awk") {
    comment_start="#"; 
    comment_end  ="\\n";
  }

  RS=comment_start; 
}


# if (record separator is comment_start), which means, we are in the source
# and need to do nothing
RS==comment_start {

  RS=comment_end;
  next;
}


# if (record separator is comment_end), which means, we are inside a comment
# and need to extract marked text 
RS==comment_end {

  # 'start of copying comment inclusion' marker
  if (match($0, "<lpdoc>")) {
    comment_inclusion=1;
  }

  # if 'comment inclusion'
  if (comment_inclusion) {

    # 'end of comment inclusion' marker
    if (match($0, "</lpdoc>"))
      comment_inclusion=0;

    comment=comment"\n"cleanup($0);
  }

  RS=comment_start;
  next;
}


END {
  print comment;
}
