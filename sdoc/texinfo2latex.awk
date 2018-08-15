# texinfo2latex.awk
#                                                                      
# texinfo2latex                                                   @ti@ 
#
# convert texinfo to latex
# ATTENTION: this is only a very poor converter which 
#            does only a small part of what it actually should do
#                                                                      
# Author:                                                              
#   Bernhard Mayer                                                @au@ 
#   Fraunhofer Institute for Atmospheric Environmental Research,  @ad@ 
#   82467 Garmisch-Partenkirchen,                                 @ad@ 
#   Germany                                                       @ad@ 
#                                                                      
# --------------------------------------------------------------------
# Copyright (C) 1997 Bernhard Mayer                               @ci@ 
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
# USA.                                                            @ic@ 
#---------------------------------------------------------------------


# replace texinfo sectioning command with corresponding latex command
function sectioning(texinfo, latex, string)
{
  texinfo = "@"texinfo;
  latex   = "\\"latex;


  if (gsub(texinfo, "", string) > 0)
    string=latex"{"string"}";

  return string;
}
  

# eat up line 
function eatline(command, string)
{
  if ($1==command)
    string="";

  return string;
}


# convert texinfo @item to latex \description[] \item 
function item(command, string, itemcommand)
{
  command = "@"command;

  if (gsub(command, "", string) > 0)
    string="\\item["itemcommand"{"string"}]";

  return string;
}


# convert texinfo source to latex source on a line by line basis
function texinfo2latex(str) 
{
  # latex does not like "_" and "&"
  gsub ("_","\\_", str);
  gsub ("\\&","\\\\&", str);

  # @asis does essentially nothing
  gsub ("@asis","", str);

  # ??? this has to be improved
  gsub ("@code","", str);

  # replace @strong with \underline
  gsub ("@strong","\\underline", str);

  # convert sectioning commands
  str = sectioning("unnumberedsubsubsec", "subsubsection*", str);
  str = sectioning("unnumberedsubsec", "subsection*", str);
  str = sectioning("unnumberedsec", "section*", str);
  str = sectioning("unnumbered", "chapter*", str);

  str = sectioning("subsubsection", "subsubsection", str);
  str = sectioning("subsection", "subsection", str);
  str = sectioning("section", "section", str);
  str = sectioning("chapter", "chapter", str);

  
  # nodes
  str = eatline("@node", str);

  # comments 
  gsub ("@c","%",str);

  # pagebreaks 
  gsub ("@page","\\newpage",str);


  # table start
  if (gsub(/@table/, "", str)) { 
    itemcommand = str;
    str="\\begin{description}";
  }

  # table end
  gsub(/@end[ \t]+table/,    "\\end{description}", str);

  # table item 
  str = item("item", str, itemcommand);


  # @example -> verbatim
  if (gsub("@example", "\\begin{verbatim}", str) > 0)
    inverbatim=1;
  
  if (gsub(/@end[ \t]+example/, "\\end{verbatim}", str) > 0)
    inverbatim=0;

  # in verbatim mode we prefer "_" instead of "\_"
  if (inverbatim) {
    gsub(/\\_/, "_", str);
    gsub(/\$\\backslash\$/, "\\", str);
  }

  # begin of function definition
  if (gsub("@deftypefun", "", str) > 0) {
    number1 = split(str,a,"(");
    str="\n\\noindent";

    str=str"{\\bf"a[1]"} (";
    
    for(i=2; i<=number1; i++)
      str=str""a[i];
  }
  
  # end of function definition
  gsub("@end deftypefun", "\n\n\\mbox{}\n", str);

  return str;
}


{
  print texinfo2latex($0);
}
  
