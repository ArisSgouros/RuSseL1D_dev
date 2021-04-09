import os
import sys
import math as m
import copy as cp
from scipy import optimize
#

#EXEC_NAME   = "/home/cjrevelas/fd_1d_tests/isolated/fd_1d_2020-01-14.exe"
EXEC_NAME      = "./fd_1d.exe"
LOG_NAME       = "scft.out.txt"
INPUT_FNAME    = "input.in.txt"
ENERGIES_FNAME = "energies.out.txt"

EXEC_LOG = True
NPROCS   = 1
#
# End of parameter section
#


LARGE_NUMBER = 10000000

g_sim_evals = 0

def set_lx(lx_value):
   try:
      line = os.popen("""grep '! lx'""" + " " + INPUT_FNAME).read()
      linesplit = line.split()
      linesplit[0] = "%10.9e             " % lx_value
      string = " ".join(linesplit)
      os.system("sed -i '/! lx/c\ "+string+"' " + INPUT_FNAME)
      print("   lx was set to ",lx_value)
   except:
      print("problem with setting lx: "+INPUT_FNAME)
      exit()
   return

def set_nx(nx_value):
   try:
      line = os.popen("""grep '! nx'""" + " " + INPUT_FNAME).read()
      linesplit = line.split()
      linesplit[0] = "%d             " % nx_value
      string = " ".join(linesplit)
      os.system("sed -i '/! nx/c\ "+string+"' " + INPUT_FNAME)
      print("   nx was set to ",nx_value)
   except:
      print("problem with setting nx: "+INPUT_FNAME)
      exit()
   return

def run_simulation():
   if EXEC_LOG:
      os.system(EXEC_NAME+" >> LOG")
   else:
      os.system(EXEC_NAME)
   return

def get_energies():
   try:
      energy_line = 2
      line = os.popen("sed -n " + str(energy_line)+"p " + ENERGIES_FNAME).read()
      line_split = line.split()
      energies = []
      energies.append(float(line_split[0]))
      energies.append(float(line_split[1]))
      energies.append(float(line_split[2]))
      energies.append(float(line_split[3]))
      energies.append(float(line_split[4]))
   except:
      print("problem with reading the " + ENERGIES_FNAME + " output file")
      exit()
   return energies

def get_grafting_density():
   try:
      line = os.popen("""grep "Grafting density     = " """ + LOG_NAME).read()
      graf_dens = float(line.split()[-1])
   except:
      print("problem with reading the grafting density from the output file")
      exit()
   return graf_dens

def get_grafted_initial_condition():
   try:
      line = os.popen("""grep "! initial value of" """ + INPUT_FNAME).read()
      q_init = float(line.split()[0])
   except:
      print("problem with reading the grafting density from the output file")
      exit()
   return q_init

def store_output(path):
   os.system("mkdir "+path)
   os.system("mv *.out.txt "+path)
   #os.system("mv *.out.bin "+path)
   os.system("cp *.in.txt "+path)
   return

def copy_last_field():
   field_filename = "field.out.bin"
   os.system("cp "+field_filename+" field.in.bin")
   return

def read_field(read_flag):
   if read_flag is False:
      os.system("sed -i '/! read field/c\ 0                             ! read field' "+INPUT_FNAME)
   else:
      os.system("sed -i '/! read field/c\ 1                             ! read field' "+INPUT_FNAME)
   return


d_x    = 0.1
l_min  = 10.0
l_max  = 200.0
l_step = 5.0 
n_step = int((l_max - l_min) / l_step) + 1

print(n_step)

foo = open("pylog.out","w")
foo.write(("%-20s "*9+"\n") % ("iter","lx","nx","gdens","term1","term1","term3","term4","total"))
for i_step in range(n_step):
   l_x_temp = l_min + i_step * l_step
   n_x = int(l_x_temp / d_x)

   if n_x % 2 != 0:
      n_x += 1

   l_x = n_x * d_x

   set_nx(n_x)
   set_lx(l_x)

   run_simulation()
   [term1, term2, term3, term4, etotal] = get_energies()
   gdens = get_grafting_density()
   foo.write(("%-20d %-20f %-20d %-20f %-20e %-20e %-20e %-20e %-20e\n") % (i_step, l_x, n_x, gdens, term1, term2, term3, term4, etotal))

foo.close()
