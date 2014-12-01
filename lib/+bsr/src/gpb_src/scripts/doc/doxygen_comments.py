#!/usr/bin/python2

# Convert C++ commenting convention to one appropriate for automatically 
# generating documentation with doxygen.

import os;
import string;
import sys;

# Get list of lines that start valid top-level (class/function) comment blocks.
#
# Top-level comment blocks are valid description blocks provided that a line
# starting with/containing one of the bind symbols is encountered following
# them before a line starting with/containing one of the stop symbols.
def top_comment_locations(filename):
   # define bind and stop symbols
   bind_start_symbols   = ["#define"];
   bind_contain_symbols = ["{"];
   stop_start_symbols   = ["}", "#if", "namespace", "using"];
   stop_contain_symbols = [");", ") const;"];
   # open file
   file = open(filename,'r');
   # initialize line and position in file
   line = file.readline();
   curr_line = 0;
   curr_comment = -1;
   in_comment = 0;
   # initialize valid comment line list
   c_top_locs = [];
   # find top-level comment blocks
   while (len(line) > 0):
      # check for start of new top-level comment
      s = line.strip();
      if ((s.startswith("/*")) and (len(s) == 2)):
         curr_comment = curr_line;
         in_comment = 1;
      elif ((s.startswith("*/")) and (len(s) == 2)):
         in_comment = 0;
      elif ((curr_comment >= 0) and (in_comment == 0)):
         # check for stopping symbols
         pos_stop = len(s) + 1;
         if (len(s) == 0):
            pos_stop = 0;
         for sym in stop_start_symbols:
            if (s.startswith(sym)):
               pos_stop = 0;
         for sym in stop_contain_symbols:
            pos = s.find(sym);
            if ((pos >= 0) and (pos < pos_stop)):
               pos_stop = pos;
         # check for binding symbols
         pos_bind = len(s) + 1;
         for sym in bind_start_symbols:
            if (s.startswith(sym)):
               pos_bind = 0;
         for sym in bind_contain_symbols:
            pos = s.find(sym);
            if ((pos >= 0) and (pos < pos_bind)):
               pos_bind = pos;
         # check for binding with pure virtual method declaration
         if ((s.startswith("virtual")) and (not (s.startswith("virtual ~"))) and (s.endswith("= 0;"))):
            pos_bind = 0;
         # handle stop/binding case
         if (pos_bind < pos_stop):
            c_top_locs.append(curr_comment);
            curr_comment = -1;
         elif (pos_stop <= len(s)):
            curr_comment = -1;
      # update position in file
      line = file.readline();
      curr_line = curr_line + 1;
   # close file
   file.close();
   return c_top_locs;

# Get list of lines that contain valid parameter comments.
# 
# A comment is a valid parameter description if it occurs within class scope
# and is located after some non-comment text on the same line.
def param_comment_locations(filename):
   # define symbols that start class sections and start/end blocks
   section_symbols = ["public:", "private:", "proctected:"];
   class_end = "};";
   block_start = "{";
   block_end = "}";
   # open file
   file = open(filename,'r');
   # initialize line and position in file
   line = file.readline();
   curr_line = 0;
   block_history = []; # push 1 when entering class, 0 when encountering {
   # initialize valid comment line list
   c_param_locs = [];
   # find parameter comments
   while (len(line) > 0):
      # ignore top-level comment blocks
      s = line.strip();
      if ((not (s.startswith("/*"))) and (not (s.startswith("*")))):
         # get current block status
         in_class = 0;
         if (len(block_history) > 0):
            in_class = block_history[len(block_history)-1];
         # check if currently in a class
         if (in_class == 1):
            # handle case where current line ends the class
            if (s.startswith(class_end)):
               block_history.pop();
            else:
               # check if current line contains parameter comment
               if (s.find("/*") > 0):
                  c_param_locs.append(curr_line);
               # check if current line begins block
               if (s.find(block_start) >= 0):
                  block_history.append(0);
               # check if current line ends block
               if (s.find(block_end) >= 0):
                  block_history.pop();
         else:
            # check for start of class section
            found_start = 0;
            for sym in section_symbols:
               if (s.startswith(sym)):
                  found_start = 1;
            if (found_start):
               l = len(block_history);
               if (l > 0):
                  if (block_history[l-1] == 0):
                     block_history[l-1] = 1;
                  else:
                     block_history.append(1);
               else:
                  block_history.append(1);
      # update position in file
      line = file.readline();
      curr_line = curr_line + 1;
   # close file
   file.close();
   return c_param_locs;

# Convert comments
def convert_comments(filename):
   # find locations of comments to update
   c_top_locs   = top_comment_locations(filename);
   c_param_locs = param_comment_locations(filename);
   c_top_locs.reverse();
   c_param_locs.reverse();
   # initialize locations of next comment to update
   next_top_loc   = -1;
   next_param_loc = -1;
   if (len(c_top_locs) > 0):
      next_top_loc = c_top_locs.pop();
   if (len(c_param_locs) > 0):
      next_param_loc = c_param_locs.pop();
   # open file
   file = open(filename,'r');
   # intialize line and position in file
   line = file.readline();
   curr_comment = -1;
   curr_line = 0;
   # update lines with doxygen style comments
   while (len(line) > 0):
      # update top-level comments
      if (curr_line == next_top_loc):
         line = line.replace("/*", "/**", 1);
         if (len(c_top_locs) > 0):
            next_top_loc = c_top_locs.pop();
         else:
            next_top_loc = -1;
      # update parameter comments
      if (curr_line == next_param_loc):
         cut_ind = line.find("{");
         if (cut_ind < 0):
            cut_ind = len(line);
         line_left  = line[0:cut_ind];
         line_right = line[cut_ind:len(line)];
         line = (line_left.replace("/*", "/**<")) + line_right;
         if (len(c_param_locs) > 0):
            next_param_loc = c_param_locs.pop();
         else:
            next_param_loc = -1;
      # output line
      sys.stdout.write(line);
      # get next line
      line = file.readline();
      curr_line = curr_line + 1;
   # close file
   file.close();

# main program
if (len(sys.argv) != 2):
   print("Usage: doxygen_comments [filename]");
   raise Exception("invalid argument(s)");
else:
   # parse arguments
   filename = sys.argv[1];
   # convert comments
   convert_comments(filename);
