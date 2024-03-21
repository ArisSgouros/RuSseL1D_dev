import numpy as np
import sys
import os

file_tag = sys.argv[1]
x_col = int(sys.argv[2])
y_col = int(sys.argv[3])

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
         y_val = foo.readline().split()[y_col]
         y_file.append(float(y_val))
   y_vals[iter_] = y_file
   #print(iter_, files[iter_])


# export_csv
file_out = y_header+".csv"
print("Exporting to file: %s" % (file_out))
with open(file_out, "w") as goo:
   goo.write("%15s," % (x_header))
   for iter_ in iters:
      goo.write("%15d," % (iter_))
   goo.write("\n")
   for ix in range(nx):
      goo.write("%15.7f," % (x_vals[ix]))
      for iter_ in iters:
         goo.write("%15.7f," % (y_vals[iter_][ix]))
      goo.write("\n")
      
      

