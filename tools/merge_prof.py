import numpy as np
import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file_tag', type=str, help='common file tag')
parser.add_argument('x_col',  type=int, help='x column')
parser.add_argument('y_col',  type=int, help='y column')
parser.add_argument('-shift',  type=float, default=0.0, help='')
parser.add_argument('-scale',  type=float, default=1.0, help='')
parser.add_argument('-transpose',  type=int, default=0, help='')

args = parser.parse_args()
file_tag = args.file_tag
x_col = args.x_col
y_col = args.y_col
shift = args.shift
scale = args.scale
transpose = bool(args.transpose)

len_tag = len(file_tag)

files = {}
for file in os.listdir(path='.'):
   if file_tag in file:
      iter_ = file[len_tag:]
      try:
         files[int(iter_)] = file
      except:
         ''

# sort the dictionary
iters = list(files.keys())
iters.sort()
files = {i: files[i] for i in iters}

# parse headers and x-axis from a file
x_header = ""
y_header = ""
x_vals = []
with open(files[0], "r") as foo:
   lsplit = foo.readline().split()
   x_header = lsplit[x_col]
   y_header = lsplit[y_col]
   while True:
      line = foo.readline()
      if not line:
         break
      x_vals.append(float(line.split()[x_col]))

nx = len(x_vals)
ny = len(files)

y_vals = {}
for iter_ in iters:
   y_file = []
   with open(files[iter_], "r") as foo:
      foo.readline() # skip first
      for ix in range(nx):
         y_val = float(foo.readline().split()[y_col])
         y_val += shift
         y_val *= scale
         y_file.append(float(y_val))
   y_vals[iter_] = y_file
   #print(iter_, files[iter_])


# export_csv

if not transpose:
   file_out = y_header+".csv"
   print("Exporting to file: %s" % (file_out))
   with open(file_out, "w") as goo:
      goo.write("%15s" % (""))
      for iter_ in iters:
         goo.write(",%15d" % (iter_))
      goo.write("\n")
      for ix in range(nx):
         goo.write("%15.10f" % (x_vals[ix]))
         for iter_ in iters:
            goo.write(",%15.10f" % (y_vals[iter_][ix]))
         goo.write("\n")
      
if transpose:
   file_out = y_header+".csv"
   print("Exporting to file: %s" % (file_out))
   with open(file_out, "w") as goo:
      goo.write("%15s" % (""))
      for ix in range(nx):
         goo.write(",%15.10f" % (x_vals[ix]))
      goo.write(",\n")
      for iter_ in iters:
         goo.write("%15d" % (iter_))
         for ix in range(nx):
            goo.write(",%15.10f" % (y_vals[iter_][ix]))
         goo.write(",\n")
      
