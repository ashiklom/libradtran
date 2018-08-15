# sdoc.awk
#                                                                      
# sdoc - Automatic source documentation                           @ti@ 
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


# @strong{sdoc} is a tool for automatic documentation of source code.      @oc@
# It may be used for documentation of C, Fortran and AWK code
# if the comments follow some simple rules. Extension to other
# languages should be straightforward. The generated documentation
# is written in texinfo. Comments which are to be included into
# the documentation may therefore contain any valid texinfo 
# command. An example for a sdoc generated file is the document which
# you are currently reading. 
#
# sdoc is a awk script which was tested with GNU awk 3.0.0
# (version 2.15 will @strong{NOT} work!).
#
# @strong{Usage:}
# @example
# gawk -f sdoc.awk -vcnf="filename" <source file> > <texinfo file> 
# @end example
#
# A configuration file, specified in the variable cnf on the command line
# may be used for some settings. Example:
# @example
# @lf@ ascii.cnf configuration @fl@ 
# @end example
#
#       "lang" is either "c", "fortran", or "awk", @nl@
#       "lev"  is either "document", "chapter", "section", or "subsection" @nl@
#       "topname" is the name of the chapter, section or subsection @nl@
#       "out"  is either "latex" or "texinfo" @nl@
#       "num"  is either "on" or "off"  @nl@
#       "page" is either "on" or "off"
#
# "lang" defaults to "c", "lev" to "document", "topname" to the filename
# without extension, "out" to "texinfo", "num" to "on", and "page" to "off" 
#
# @strong{level:} "document" means that a complete document will be created
# while "chapter", "section", or "subsection" will create only 
# part of a document.
#
# @strong{output:} if latex is chosen, the output of sdoc.awk has
# to be piped through 'gawk -f texinfo2latex.awk'
#
# @strong{numbering:} "on"   numbered sections,
#                     "off"  unnumbered sections
#
# Parts of other files may be included using @@lf@@ filename identifier @@fl@@ .
# The file filename will be parsed for text included between @@begin identifier@@
# and @@end identifier@@ . @@lf@@ ... @@fl@@ will be replaced by the found text. 
#  
# @end enumerate
#
#
# @strong{Examples} (example files are included in the package)
# @itemize @bullet
# @item C source:
#   @example
#   gawk -f sdoc.awk ascii.h ascii.c > ascii.texi 
#   @end example
# @item C source, more complicated:
#   @example
#   make uvbtexinfo
#   @end example
# will create a texinfo documentation for the uvb package which 
# consists of two libraries, ASCII and NUMERIC.
#   @example
#   make uvblatex
#   @end example 
# will create a LaTeX file with the same contents.
#
# @item Fortran:
#   @example
#   gawk -f sdoc.awk -vlanguage="fortran" absdep.f > absdep.texi 
#   @end example
# @item awk:
#   @example
#   gawk -f sdoc.awk -vlanguage="awk" sdoc.awk > sdoc.texi 
#   @end example
# or
#   @example
#   make sdoc 
#   @end example
# will create the documentation which you are currently reading.
# @end itemize
# @co@


# Remove all command strings from a comment  @fi@ @if@
function cleanup(string)
{
  # at first handle the 'lookupfile' feature
  string=lookupfile(string);

  # remove 'lookupfile' remainders
  gsub(/ @begin[ \t]+[^ ]+@ /,"",string);
  gsub(/ @end[ \t]+[^ ]+@ /,"",string);


  # clear string, if only ****
  match(string, /\*+/);
  if (RLENGTH==length(string)) {
    string="";
    return string;
  }

  # 'copying' comment inclusion
  gsub("@ci@","",string);
  gsub("@ic@","",string);


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
                                                                       

  return string;
}


# Print a header for the texinfo output @fi@ @if@
function print_texinfo_header(filename, infoname, year, author, address, title)
{
  if (level=="document") {
    print "\\input texinfo  @c -*-texinfo-*-";
    print "@c %**start of header";
    print "@setfilename "infoname;
    print "@settitle   "title;
    print "@c %**end of header";
    print "";
    print "";
    print "@c @setchapternewpage odd";
    print "";
    print "@ifinfo";
    print title;
    print "";
    print "Copyright @copyright{} "year" "author;
    print "@end ifinfo";
    print "";
    print "";
    print "@titlepage";
    print "@title "title;
    print "@author "author;
    print "";
    print "";
    print "@c  The following two commands";
    print "@c  start the copyright page.";
    print "@page";
    print "@vskip 0pt plus 1filll";
    print "Copyright @copyright{} "year" "author;
    print "";
    print "@end titlepage";
    print "";
    print "@ifinfo";
    print "@node Top, Copying, (dir), (dir)";
    print "@top sdoc";
    print "@end ifinfo";
    print "";
    print "@menu";
    print "* Copying::          Your rights and freedoms.";
    print "* Overview::         Overview of the library.";
    print "* Definitions::      Type and macro definitions.";
    print "* Functions::        Function definitions.";
    print "* Functions index::  Index of functions";
    print "@end menu";
    print "";
  } 
}


# Print a header for the latex output @fi@ @if@
function print_latex_header(filename, infoname, year, author, address, title)
{
  if (level=="document") {
    print "\\documentstyle[12pt]{report}";
    print "\\begin{document}";
    print "\\begin{titlepage}";
    print "\\title{"title"}";
    print "\\author{"author"}";
    print "\\end{titlepage}";
    print "\\maketitle";

    # ??? copyright note not considered up to now
    # print "Copyright @copyright{} "year" "author;
    # print "\\maketitle";
    # print "@c  The following two commands";
    # print "@c  start the copyright page.";
    # print "@page";
    # print "@vskip 0pt plus 1filll";
    # print "Copyright @copyright{} "year" "author;

    print "\\newpage";
    print "\\tableofcontents";
    print "\\newpage";
  }
}


# Print node and section  @fi@ @if@
function print_texinfo_node(nodename, title)
{

  if (level=="chapter") {
    print "@node "nodename;
    print chapter_header" "title;
    print "";
  }


  if (level=="section") {
    print "@node "nodename;
    print section_header" "title;
    print "";
  }

  if (level=="subsection") {
    print "@node "nodename;
    print subsection_header" "title;
    print "";
  }
}


# Print the texinfo 'Copying' chapter  @fi@ @if@
function print_texinfo_copying(string)
{
  if (!COPYING)
    return;

  if (level=="document")  {
    print "@node Copying, Overview, Top, Top";
    print chapter_header" Copying";

    print string;
    print "";
  
    if (pagebreak=="on")
      print "@page";
  }
}


# Print the texinfo 'Overview' chapter  @fi@ @if@
function print_texinfo_overview(string)
{
  if (!OVERVIEW)
    return;

  if (level=="document")  {
    print "@node Overview, Definitions, Copying, Top";
    print chapter_header" Overview";
  }


  if (level=="chapter" || level=="section" || level=="subsection")
# removed topname"." from node name, BM, April 26, 1998 
#    print "@node "topname".Overview, "topname".Definitions, "topname", "topname;
    print "@node Overview, Definitions, "topname", "topname;
    

  if (level=="chapter")
    print section_header" Overview";

  if (level=="section")
    print subsection_header" Overview";

  if (level=="subsection")
    print subsubsection_header" Overview";

  print string;
  print "";

  if (pagebreak=="on")
    print "@page";
}


# Print the texinfo 'Definitions' chapter  @fi@ @if@
function print_texinfo_definitions(string)
{
  if (!DEFINITIONS)
    return;

  if (level=="document")  {
    print "@node Definitions, Functions, Overview, Top";
    print chapter_header" Definitions";
  }

  if (level=="chapter" || level=="section" || level=="subsection")
# removed topname"." from node name, BM, April 26, 1998 
#  print "@node "topname".Definitions, "topname".Functions, "topname".Overview, "topname;
    print "@node Definitions, Functions, Overview, "topname;


  if (level=="chapter")
    print section_header" Definitions";

  if (level=="section")
    print subsection_header" Definitions";

  if (level=="subsection") 
    print subsubsection_header" Definitions";

  print string;
  print "";

  if (pagebreak=="on")
    print "@page";
}


# Print the texinfo 'Functions' chapter  @fi@ @if@
function print_texinfo_functions(string)
{
  if (!FUNCTIONS)
    return;

  if (level=="document")  {
    print "@node Functions, Functions index, Definitions, Top";
    print "@chapter Functions";
  }

  if (level=="chapter" || level=="section" || level=="subsection")
# removed topname"." from node name, BM, April 26, 1998 
#  print "@node "topname".Functions, , "topname".Definitions, "topname;
    print "@node Functions, , Definitions, "topname;

  if (level=="chapter")
    print section_header" Functions";

  if (level=="section")
    print subsection_header" Functions";

  if (level=="subsection")
    print subsubsection_header" Functions";

  print string;
  print "";

  if (pagebreak=="on")
    print "@page";
}


# Print a texinfo user defined chapter  @fi@ @if@
function print_texinfo_userchapter(chapter_name, chapter_text)
{
  if (level=="document")  {
    print "@node "chapter_name;
    print chapter_header" "chapter_name;
  }

  if (level=="chapter")  {
# removed topname"." from node name, BM, April 26, 1998 
#    print "@node "topname"."chapter_name;
    print "@node "chapter_name;
    print section_header" "chapter_name;
  }

  if (level=="section")  {
# removed topname"." from node name, BM, April 26, 1998 
#    print "@node "topname"."chapter_name;
    print "@node "chapter_name;
    print subsection_header" "chapter_name;
  }
  
  if (level=="subsection")  {
# removed topname"." from node name, BM, April 26, 1998 
#    print "@node "topname"."chapter_name;
    print "@node "chapter_name;
    print subsubsection_header" "chapter_name;
  }

  print chapter_text;
  print "";

  if (pagebreak=="on")
    print "@page";

}


# Print a texinfo user defined section  @fi@ @if@
function print_texinfo_usersection(section_name, section_text)
{
  if (level=="document")  {
    print "@node "section_name;
    print section_header" "section_name;
  }

  if (level=="chapter")  {
# removed topname"." from node name, BM, April 26, 1998 
#    print "@node "topname"."section_name;
    print "@node "section_name;
    print subsection_header" "section_name;
  }

  if (level=="section")  {
# removed topname"." from node name, BM, April 26, 1998 
#    print "@node "topname"."section_name;
    print "@node "section_name;
    print subsubsection_header" "section_name;
  }

  if (level=="subsection")  {
    print "@strong{"section_name"}";
  }

  print section_text;
  print "";

}

# Print a latex user defined chapter  @fi@ @if@
function print_latex_userchapter(chapter_name, chapter_text)
{
  # ??? to be done
  return;
}


# Print a latex user defined section  @fi@ @if@
function print_latex_usersection(section_name, section_text)
{
  # ??? to be done
  return;
}


# print a submenu (only if level=="section" or "subsection") @fi@ @if@
function print_texinfo_submenu()
{
  if(level=="chapter" || level=="section" || level=="subsection")  {
    print "@menu";
#    print "* "topname".Overview::         Overview of the module.";
#    print "* "topname".Definitions::      Definitions.";
#    print "* "topname".Functions::        Functions.";
    print "* Overview::         Overview of the module.";
    print "* Definitions::      Definitions.";
    print "* Functions::        Functions.";
    print "@end menu";
  }
}


# Print the end of the texinfo output (table of contents, index of functions, ...) @fi@ @if@
function print_texinfo_end()
{
  if (level=="document") {
    print "@node Functions index, , Functions, Top";
    print "@unnumbered Index of functions";
    print "@printindex fn";
    print "@contents";
    print "";
    print "@bye";
  }
}


# Print the end of the latex output @fi@ @if@
function print_latex_end()
{
  if (level=="document") {
    print "\\end{document}";
  }
}




# Replace @@lf@@ filename identifier @@fl@@ with text between            @fi@ 
# @@identifier@@ and @@end identifier@@  from a different file filename   @if@
function lookupfile(string,
		    rs_save, line_save, 
		    found, filename, 
		    start_identifier, end_identifier, 
		    replace,temp,tempnum,i)
{
  # save important things
  rs_save=RS;
  line_save=$0;


  while(match(string, /@lf@[ \t]+[^ ]+[ \t]+[^ ]+[ \t]+@fl@/))  {
    found=substr(string,RSTART,RLENGTH);
    split(found,a," ");
  
    filename=a[2];
    start_identifier = "@begin "a[3]"@";
    end_identifier   = "@end "a[3]"@";
    
    # read text from filename 
    RS = start_identifier;
    getline<filename;
    RS = end_identifier;
    getline replace < filename;
    close(filename);
    
    # important: replace '&' with '\&' because of special treatment 
    # of '&' in replacement strings
    gsub("&", "\\&", replace);
  
    # remove first comment_end and last comment_start 
    # because we don't like them here.
    # first
    sub(/[ ]+/,"",replace);
    sub(comment_end, "", replace);

    # last
    tempnum = split(replace,temp,comment_start);
    replace=temp[1];
    for (i=2; i<tempnum; i++)
      replace=replace""comment_start""temp[i];

    
    # replace command by the string found between start_identifier
    # and end_identifier in filename
    sub(found,replace,string);
  }  	
  
  # restore important things
  RS=rs_save;
  $0=line_save;

  return string;
}


# Read configuration file    @fi@ @if@
function read_config(cnf,
                     line, linecounter, comment)     # local variables
{
  # if no configuration file: return and use defaults
  if (cnf=="")
    return;


  while ((getline line < cnf) > 0)  {
    linecounter++;

    
    # get rid of comments 
    split(line,comment,"#");
    line=comment[1];

    split(line,token," ");

    if (tolower(token[1])=="language")
      language = tolower(token[2]);

    if (tolower(token[1])=="level")
      level = tolower(token[2]);

    if (tolower(token[1])=="output")
      output = tolower(token[2]);

    if (tolower(token[1])=="numbering")
      numbering = tolower(token[2]);

    if (tolower(token[1])=="pagebreak")
      pagebreak = tolower(token[2]);

    if (tolower(token[1])=="topname")  {
      split(line,string,"\\\"");
      topname = string[2];
    }

    if (tolower(token[1])=="title")  {
      split(line,string,"\\\"");
      title = string[2];
    }

    if (tolower(token[1])=="chaptername") {
      split(line,string,"\\\"");
      chapter_name[token[2]]=string[2];
    }

    if (tolower(token[1])=="sectionname") {
      split(line,string,"\\\"");
      section_name[token[2],token[3]]=string[2];
    }
  }

  # here should a range check be performed ???


  if (linecounter==0)
    print "ATTENTION: "cnf" empty or not existing";

  return;
}



BEGIN { 

  # set defaults 
  language="c";
  level="document";
  output="texinfo";
  numbering="on";
  pagebreak="off";
  topname="";
  title="";


  # read configuration file 
  read_config(cnf);

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


  # section headings 
  if (numbering=="on") {
    chapter_header       = "@chapter";
    section_header       = "@section";
    subsection_header    = "@subsection";
    subsubsection_header = "@subsubsection";
  }
  else {
    chapter_header       = "@unnumbered";
    section_header       = "@unnumberedsec";
    subsection_header    = "@unnumberedsubsec";
    subsubsection_header = "@unnumberedsubsubsec";
  }


  # how to indicate "verbatim" in Texinfo
  source_start_symbol = "@example";
  source_end_symbol   = "@end example";

  RS=comment_start; 
}


# if (record separator is comment_start), which means, we are in the source
RS==comment_start {

  if (first_overview_source_line==1) {
    overview = overview"\n"source_start_symbol;
    first_overview_source_line=0;
  }

  if (overview_source_inclusion)
    overview = overview"\n"cleanup($0);
  


  if (first_definition_source_line==1) {
    definitions = definitions"\n"source_start_symbol;
    first_definition_source_line=0;
  }

  if (definition_source_inclusion)
    definitions = definitions"\n"cleanup($0);
  


  if (first_function_source_line==1) {
    functions = functions"\n"source_start_symbol;
    first_function_source_line=0;
  }

  if (function_source_inclusion)
    functions = functions"\n"cleanup($0);
  

  # check if function is to be included in some chapter/section;
  # determine type of section header

  function_all_inclusion=function_inclusion;

  if (function_inclusion) {
    if (level=="document") 
      func_header=section_header;
    if (level=="chapter") 
      func_header=subsection_header;
    if (level=="section") 
      func_header=subsubsection_header;
    if (level=="subsection") 
      func_header=" ";
  }


  for (chapnum in chapter_function_inclusion) {
    function_all_inclusion += chapter_function_inclusion[chapnum];
    if (chapter_function_inclusion[chapnum]>0) {
      if (level=="document") 
        func_header=section_header;
      if (level=="chapter") 
        func_header=subsection_header;
      if (level=="section") 
        func_header=subsubsection_header;
      if (level=="subsection") 
        func_header=" ";
    }
  }


  for (secnum in section_function_inclusion) {
    function_all_inclusion += section_function_inclusion[secnum];
    if (section_function_inclusion[secnum]>0) {
      if (level=="document") 
        func_header=subsection_header;
      if (level=="chapter") 
        func_header=subsubsection_header;
      if (level=="section") 
        func_header=" ";
      if (level=="subsection") 
        func_header=" ";
    }
  }


  # if yes, interpret function documentation and create a string func_add,
  # which has to be added to the respective chapter/section

  if (function_all_inclusion) {

    if(!match($0, ")" ) )
      func_definition = func_definition"\n"cleanup($0);
    else {
      split($0,temp, ")");

      func_definition = temp[1]")\n";

      
#      if (func_definition="")    
#	func_definition = temp[1]")\n";
#      else 
#	func_definition = func_definition"\n"temp[1]")\n";

      # create single line with single space only
      gsub(/\n/,"",func_definition);
      gsub(/[ \t]+/," ",func_definition);

      # remove "$" from function definition if fortran
      if(language=="fortran")
	gsub(/\$/, "", func_definition);

   
      func_description=func_text;

      # split up function text
      split(func_text,temp,"Function:");      
      func_text=temp[2];

      split(func_text,temp,"Description:");      
      func_name = temp[1];
      func_text=temp[2];

      split(func_text,temp,"Parameters:");      
      if (temp[1]!="")
	func_description = temp[1];
      func_text=temp[2];

      split(func_text,temp,"Return value:");      
      func_parameters = temp[1];
      func_text=temp[2];

      split(func_text,temp,"Example:");      
      func_returnvalue = temp[1];
      func_text=temp[2];

      split(func_text,temp,"Files:");      
      func_example = temp[1];
      func_text=temp[2];

      split(func_text,temp,"Known bugs:");      
      func_files = temp[1];
      func_text  = temp[2];

      split(func_text,temp,"Author:");      
      func_bugs   = temp[1];
      func_author = temp[2];

      # reset parameters
      paranum=0;
      for(i in func_paraname) {
        delete func_paraname[i];
        delete func_paradesc[i];
      }

      # split parameters
      temp1num = split(func_parameters,temp1,"[\n]+");
      for(i=1; i<=temp1num; i++) {
        if(match(temp1[i],":")) {
	  split(temp1[i],temp2,":");
	  paranum++;
          func_paraname[paranum]=temp2[1];
          func_paradesc[paranum]=func_paradesc[paranum]" "temp2[2];
        }
        else {
          func_paradesc[paranum]=func_paradesc[paranum]" "temp1[i];
        }
      }

      # defaults
      match(func_author,/[ \t\n]*/);
      if (RLENGTH==length(func_author))
	func_author=author;

      match(func_files,/[ \t\n]*/);
      if (RLENGTH==length(func_files))
	func_files="none";

      match(func_bugs,/[ \t\n]*/);
      if (RLENGTH==length(func_bugs))
	func_bugs="none";
        
      func_add = func_add"\n  "func_header" Function "func_name"\n";  
      func_add = func_add"\n  @deftypefun "func_definition"\n\n""@end deftypefun\n";  

      func_add = func_add"  @table @asis\n";  
      func_add = func_add"  @item Description:\n"func_description"\n";
      func_add = func_add"  @item Parameters:\n  @table @code\n";

      for(i=1; i<=paranum; i++)
        func_add = func_add"  @item "func_paraname[i]"\n"func_paradesc[i]"\n";

      func_add = func_add"  @end table\n";  

      func_add = func_add"  @item Return value:\n"func_returnvalue"\n";
      func_add = func_add"  @item Example:\n"func_example"\n";  
      func_add = func_add"  @item Files:\n"func_files"\n";
      func_add = func_add"  @item Known bugs:\n"func_bugs"\n";
      func_add = func_add"  @item Author:\n"func_author"\n";

      func_add = func_add"  @end table\n";  


      # add generated text to chapter/section
      if (function_inclusion)
        functions=functions""func_add;

      for (chapnum in chapter_function_inclusion)
	if (chapter_function_inclusion[chapnum])
   	  chapter_text[chapnum] = chapter_text[chapnum]""func_add;

      for (secnum in section_function_inclusion)
	if (section_function_inclusion[secnum])
	  section_text[secnum] = section_text[secnum]""func_add;

      # reset func_add
      func_add = "";


      # reset function_inclusion
      function_inclusion=0;

      for (chapnum in chapter_function_inclusion)
	chapter_function_inclusion[chapnum]=0; 

      for (secnum in section_function_inclusion)
	section_function_inclusion[secnum]=0;


      func_name="";
      func_text="";

    }
  }


 

  for (chapnum in chapter) {
    if (first_chapter_source_line[chapnum]) {
      chapter_text[chapnum] = chapter_text[chapnum]"\n"source_start_symbol;
      first_chapter_source_line[chapnum]=0;
    } 

    if (chapter_source_inclusion[chapnum])
      chapter_text[chapnum] = chapter_text[chapnum]"\n"cleanup($0);
  }


  for (secnum in section) {
    if (first_section_source_line[secnum]) {
      section_text[secnum] = section_text[secnum]"\n"source_start_symbol;
      first_section_source_line[secnum]=0;
    } 

    if (section_source_inclusion[secnum])
      section_text[secnum] = section_text[secnum]"\n"cleanup($0);
  }


  RS=comment_end;
  next;
}


# if (record separator is comment_end), which means, we are inside a comment
RS==comment_end {

  # author
  if (gsub("@au@", ""))
    author = $0;

  # address
  if (gsub("@ad@", ""))
    address = address""$0;

  # year
  if (gsub("@yr@", ""))
    year = year""$0;

  # title
  # if (gsub("@ti@", ""))
  #   title = title""$0;
  # @ti@ is ignored, Bernhard, September 23, 1997


  # 'start of overview source inclusion' marker
  if (match($0, "@os@")) {
    OVERVIEW=1;
    overview_source_inclusion=1;
    first_overview_source_line=1;
  }

  # 'end of source inclusion' marker
  if (match($0, "@so@")) {
    overview_source_inclusion=0;
    overview = overview"\n"source_end_symbol;
  }


  # 'start of copying comment inclusion' marker
  if (match($0, "@ci@")) {
    COPYING=1;
    copying_comment_inclusion=1;
  }

  # if 'comment inclusion'
  if (copying_comment_inclusion) {

    # 'end of comment inclusion' marker
    if (match($0, "@ic@"))
      copying_comment_inclusion=0;

    copying=copying"\n"cleanup($0);
  }


  # 'start of overview comment inclusion' marker
  if (match($0, "@oc@")) {
    OVERVIEW=1;
    overview_comment_inclusion=1;
  }

  # if 'overview comment inclusion'
  if (overview_comment_inclusion) {

    # 'end of comment inclusion' marker
    if (match($0, "@co@"))
      overview_comment_inclusion=0;

    overview=overview"\n"cleanup($0);
  }




  # 'start of definitions source inclusion' marker
  if (match($0, "@ds@")) {
    DEFINITIONS=1;
    definition_source_inclusion=1;
    first_definition_source_line=1;
  }

  # 'end of source inclusion' marker
  if (match($0, "@sd@")) {
    definition_source_inclusion=0;
    definitions = definitions"\n"source_end_symbol;
  }


  # 'start of definition comment inclusion' marker
  if (match($0, "@dc@")) {
    DEFINITIONS=1;
    definition_comment_inclusion=1;
  }


  # if 'definition comment inclusion'
  if (definition_comment_inclusion) {

    # 'end of comment inclusion' marker
    if (match($0, "@cd@"))
      definition_comment_inclusion=0;

    definitions = definitions"\n"cleanup($0);
  }




  # 'start of function source inclusion' marker
  if (match($0, "@fs@")) {
    FUNCTIONS=1;
    function_source_inclusion=1;
    first_function_source_line=1;
  }

  # 'end of source inclusion' marker
  if (match($0, "@sf@")) {
    function_source_inclusion=0;
    functions = functions"\n"source_end_symbol;
  }


  # 'start of function comment inclusion' marker
  if (match($0, "@fc@")) {
    FUNCTIONS=1;
    function_comment_inclusion=1;
  }

  # if 'function comment inclusion'
  if (function_comment_inclusion) {

    # 'end of comment inclusion' marker
    if (match($0, "@cf@"))
      function_comment_inclusion=0;

    functions = functions"\n"cleanup($0);
  }


  # 'start of function inclusion' marker
  if (match($0, "@fi@")) {
    FUNCTIONS=1;
    function_inclusion=1;
    function_description=1;
  }
  
  # start function inclusion in user defined chapter
  if(match($0,/@[0-9][0-9]i@/)) {
    index1=substr($0,RSTART+1,2);
    chapter[index1]=1;

    chapter_function_inclusion[index1]=1;
    function_description=1;
  }
  
  # start function inclusion in user defined section
  if(match($0,/@[0-9][0-9]_[0-9][0-9]i@/)) {
    index1 = substr($0,RSTART+1,2);
    index2 = substr($0,RSTART+4,2);
    section[index1,index2]=1;

    section_function_inclusion[index1,index2]=1;
    function_description=1;
  }


  # if inside function description
  if (function_description) {

    # 'end of function inclusion' marker
    if (match($0, "@if@"))
      function_description=0;

    # end function inclusion in user defined chapter
    if(match($0,/@i[0-9][0-9]@/))
      function_description=0;

    # end function inclusion in user defined section
    if(match($0,/@i[0-9][0-9]_[0-9][0-9]@/))
      function_description=0;

    func_text = func_text"\n"cleanup($0);
  }


  # start comment inclusion in user defined chapter
  if(match($0,/@[0-9][0-9]c@/)) {
    index1=substr($0,RSTART+1,2);
    chapter[index1]=1;
    chapter_comment_inclusion[index1]=1;
  }

  # start source inclusion in user defined chapter
  if(match($0,/@[0-9][0-9]s@/)) {
    index1=substr($0,RSTART+1,2);
    chapter[index1]=1;
    chapter_source_inclusion[index1]=1;
    first_chapter_source_line[index1]=1;
  }

  # start comment inclusion in user defined section
  if(match($0,/@[0-9][0-9]_[0-9][0-9]c@/)) {
    index1 = substr($0,RSTART+1,2);
    index2 = substr($0,RSTART+4,2);
    section[index1,index2]=1;
    section_comment_inclusion[index1,index2]=1;
  }

  # start source inclusion in user defined section
  if(match($0,/@[0-9][0-9]_[0-9][0-9]s@/)) {
    index1 = substr($0,RSTART+1,2);
    index2 = substr($0,RSTART+4,2);
    section[index1,index2]=1;
    section_source_inclusion[index1,index2]=1;
    first_section_source_line[index1,index2]=1;
  }

  # include comment in user defined chapter
  for (chapnum in chapter_comment_inclusion)
    if (chapter_comment_inclusion[chapnum])  
      chapter_text[chapnum] = chapter_text[chapnum]"\n"cleanup($0);
      
  # include comment in user defined section
  for (secnum in section_comment_inclusion)
    if (section_comment_inclusion[secnum])  
      section_text[secnum] = section_text[secnum]"\n"cleanup($0);
      

  # end comment inclusion into user defined chapter
  if(match($0,/@c[0-9][0-9]@/)) {
    index1=substr($0,RSTART+2,2);
    chapter_comment_inclusion[index1]=0;
  }

  # end source inclusion into user defined chapter
  if(match($0,/@s[0-9][0-9]@/)) {
    index1=substr($0,RSTART+2,2);
    chapter_source_inclusion[index1]=0;
    chapter_text[index1] = chapter_text[index1]"\n"source_end_symbol;
  }

  # end comment inclusion into user defined section
  if(match($0,/@c[0-9][0-9]_[0-9][0-9]@/)) {
    index1 = substr($0,RSTART+2,2);
    index2 = substr($0,RSTART+5,2);
    section_comment_inclusion[index1,index2]=0;
  }

  # end source inclusion into user defined section
  if(match($0,/@s[0-9][0-9]_[0-9][0-9]@/)) {
    index1 = substr($0,RSTART+2,2);
    index2 = substr($0,RSTART+5,2);
    section_source_inclusion[index1,index2]=0;
    section_text[index1,index2] = section_text[index1,index2]"\n"source_end_symbol;
  }

  # chapter name
  # if(match($0,/@[0-9][0-9]n .+@/))
  #   chapter_name[substr($0,RSTART+1,2)]=substr($0,RSTART+5,RLENGTH-6);

  # section name
  # if(match($0,/@[0-9][0-9]_[0-9][0-9]n .+@/))
  #   section_name[substr($0,RSTART+1,2), substr($0,RSTART+4,2)]=substr($0,RSTART+7,RLENGTH-8);


  RS=comment_start;
  next;
}


END {

  # replace several "\n"s by just two "\n"s
  gsub(/\n\n+/,"\n\n", copying);
  gsub(/\n\n+/,"\n\n", overview);
  gsub(/\n\n+/,"\n\n", definitions);
  gsub(/\n\n+/,"\n\n", functions);

  # create the name for the .info file
  number = split (FILENAME,temp,".");
  if (number>1)
    basename=temp[1];
  else 
    basename=FILENAME;

  for (i=2; i<number; i++)
    basename = basename"."temp[i];
  
  infoname = basename".info";

  if (topname=="")
    topname=basename;  

  if (output=="texinfo") {
    print_texinfo_header(FILENAME, infoname, year, author, address, title);
    print_texinfo_node(topname, title);
    print_texinfo_submenu();
    print_texinfo_copying(copying);
    print_texinfo_overview(overview);
    print_texinfo_definitions(definitions);
    print_texinfo_functions(functions);


# removed the next 8 lines and replaced them with the following 36
# in order to fix the sort problem: 
#  user-defined chapters sections should now be sorted by their number

#    for(chapnum in chapter) {
#      print_texinfo_userchapter(chapter_name[chapnum], chapter_text[chapnum]);
#      for (secnum in section) {
#	split(secnum,sep_secnum,SUBSEP);
#	if(sep_secnum[1]==chapnum)
#          print_texinfo_usersection(section_name[sep_secnum[1],sep_secnum[2]], section_text[sep_secnum[1],sep_secnum[2]]);
#      }
#    }


    chapnum=0;
    while (chapnum+0<100) {
      if (chapnum<10) 
        chapnum="0"chapnum;

      if (chapter_text[chapnum]!="")
        print_texinfo_userchapter(chapter_name[chapnum], chapter_text[chapnum]);

      for (secnum in section) {
	split(secnum,sep_secnum,SUBSEP);

	if(sep_secnum[1]==chapnum) {
          temp_section_text[sep_secnum[2]] = section_text[sep_secnum[1],sep_secnum[2]];
          temp_section_name[sep_secnum[2]] = section_name[sep_secnum[1],sep_secnum[2]];
        } 
      }

      secnum=0;
      
      while (secnum+0<100) {
        if (secnum<10) 
          secnum="0"secnum;

	if (temp_section_text[secnum]!="")
          print_texinfo_usersection(temp_section_name[secnum], temp_section_text[secnum]);
	# reset text

	temp_section_text[secnum]="";
	secnum++;
      }
      chapnum++;
    }

# end replace, September 12, 1997 
   
    print_texinfo_end();
  }

  if (output=="latex") {
    print_latex_header(FILENAME, infoname, year, author, address, title);

    # latex does not like "\"
    gsub(/\\/, "$\\backslash$", copying);
    gsub(/\\/, "$\\backslash$", overview);
    gsub(/\\/, "$\\backslash$", definitions);
    gsub(/\\/, "$\\backslash$", functions);

    print_texinfo_node(topname, title);
    print_texinfo_copying(copying);
    print_texinfo_overview(overview);
    print_texinfo_definitions(definitions);
    print_texinfo_functions(functions);

    for(chapnum in chapter) {
      gsub(/\\/, "$\\backslash$", chapter_text[chapnum]);
      print_latex_userchapter(chapter_name[chapnum], chapter_text[chapnum]);
      for (secnum in section) {
	split(secnum,sep_secnum,SUBSEP);
	if(sep_secnum[1]==chapnum) {
          gsub(/\\/, "$\\backslash$", section_text[sep_secnum[1],sep_secnum[2]]);
          print_latex_usersection(section_name[sep_secnum[1],sep_secnum[2]], section_text[sep_secnum[1],sep_secnum[2]]);
          }
      }
    }


    print_latex_end();

  }
}
