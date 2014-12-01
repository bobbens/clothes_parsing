#!/usr/bin/python2

# Automatically determine the list of object (.o) files required for linking 
# a specified application.

import os;
import string;
import sys;

# the path info class stores build path settings
class PathInfo:
   # constructor - initialize path information
   def __init__(self, bin_dir, depend_dir, include_dir, src_dir):
      self.bin_dir     = bin_dir;
      self.depend_dir  = depend_dir;
      self.include_dir = include_dir;
      self.src_dir     = src_dir;

   # print the path settings
   def display(self):
      print("bin dir     = " + self.bin_dir);
      print("depend dir  = " + self.depend_dir);
      print("include dir = " + self.include_dir);
      print("src dir     = " + self.src_dir);

# get an initial header list from a depend (.d) file
def hh_list_from_d(filename, path_to_root, path_info):
   # read the list from the file
   hh_list = [];
   path = path_to_root + path_info.include_dir + "/";
   file = open(filename,'r');
   line = file.readline().lstrip().rstrip(" \\\n");
   while (len(line) > 0):
      s_set = line.split(' ');
      for s in s_set:
         if ((s.startswith(path)) and (s.endswith(".hh"))):
            hh_list.append(s.replace(path, "", 1));
      line = file.readline().lstrip().rstrip(" \\\n");
   file.close();
   # only return unique items in the header list
   hh_list_unique = [];
   hh_seen = {};
   for hh in hh_list:
      if (not hh_seen.has_key(hh)):
         hh_list_unique.append(hh);
         hh_seen[hh] = 1;
   return hh_list_unique;

# get the source (.cc) containing the implementation for the given header (.hh) file
# return the empty string if there is no such source file
def cc_from_hh(hh):
   return hh[0:len(hh)-3] + ".cc";

# get list of source (.cc) files corresponding to given header (.hh) files
def cc_list_from_hh_list(hh_list):
   cc_list = [];
   for hh in hh_list:
      cc = cc_from_hh(hh);
      cc_list.append(cc);

# get the depend (.d) file for the given source (.cc) file
def d_from_cc(cc):
   return cc[0:len(cc)-3] + ".d";

# get the object (.o) file for the given source (.cc) file
def o_from_cc(cc):
   return cc[0:len(cc)-3] + ".o";

# check existance of a source (.cc) file
def cc_file_exists(cc, path_to_root, path_info):
   cc_filename = path_to_root + path_info.src_dir + "/" + cc;
   return os.path.exists(cc_filename); 

# get the relative filename and path to root for the given depend file
def d_get_filename(d, path_to_root, path_info):
   d_path_from_root = path_info.depend_dir + "/" + d;
   d_filename = path_to_root + d_path_from_root;
   n_slashes = d_path_from_root.count("/");
   d_path_to_root = "";
   for n in range(n_slashes):
      d_path_to_root = d_path_to_root + "../";
   return [d_filename, d_path_to_root];
   
# generate list of all object (.o) files that must be linked into an application 
# given a depend (.d) file describing its list of header file dependencies
def o_list_from_d(filename, path_to_root, path_info):
   # initialize dependency search 
   hh_pending = hh_list_from_d(filename, path_to_root, path_info);   # list of possible header dependencies
   hh_touched = {};                                                  # headers already processed or currently pending
   o_list     = [];                                                  # object file dependencies
   # add all pendinging items to touched dictionary
   for hh in hh_pending:
      hh_touched[hh] = 1;
   # perform dependency search
   while (len(hh_pending) > 0):
      hh = hh_pending.pop();
      cc = cc_from_hh(hh);
      if (cc_file_exists(cc, path_to_root, path_info)):
         # add object file dependency
         o = o_from_cc(cc);
         o_list.append(path_to_root + path_info.bin_dir + "/" + o);
         # process dependencies for that object file
         d = d_from_cc(cc);
         [d_filename, d_path_to_root] = d_get_filename(d, path_to_root, path_info);
         hh_list = hh_list_from_d(d_filename, d_path_to_root, path_info);
         for hh_new in hh_list:
            if (not hh_touched.has_key(hh_new)):
               hh_pending.append(hh_new);
               hh_touched[hh_new] = 1;
   # return list of dependencies
   o_list.sort();
   return o_list;

# main program
if (len(sys.argv) != 8):
   print("Usage: depend_o [target name] [filename.d] [path to root] [bin dir] [depend dir] [include dir] [src dir]");
   raise Exception("invalid argument(s)");
else:
   # parse arguments
   target_name    = sys.argv[1];
   filename       = sys.argv[2];
   path_to_root   = sys.argv[3];
   bin_dir        = sys.argv[4];
   depend_dir     = sys.argv[5];
   include_dir    = sys.argv[6];
   src_dir        = sys.argv[7];
   path_info = PathInfo(bin_dir, depend_dir, include_dir, src_dir);
   # compute dependencies 
   o_list = o_list_from_d(filename, path_to_root, path_info);
   # open output file, write target name
   filename_out = filename[0:len(filename)-2] + ".do";
   filename_o   = filename[0:len(filename)-2] + ".o";
   file = open(filename_out,'w');
   file.write(target_name + ": " + filename_o);
   # write dependencies
   n_o = len(o_list);
   if (n_o > 0):
      file.write(" \\\n");
   for n in range(n_o - 1):
      file.write("  " + o_list[n] + " \\\n");
   if (n_o > 0):
      file.write("  " + o_list[n_o-1]);
   # close output file
   file.close();
