import os
import sys
import math as m
import copy as cp

yy_target = 0.001
tol_bisection = 0.00001
tol_secant = 0.00000001
x_init = 10000.0

EXEC_NAME   = "./fd_1d.exe"
LOG_NAME    = "scft.out.txt"
INPUT_FNAME = "input.in.txt"

NPROCS = 1

def set_grafted_initial_condition(q_init):
   try:
      str_q_init = str(q_init)

      os.system("sed -i '/! initial value of grafting point/c\ "+str(q_init)+"                             ! initial value of grafting point' "+INPUT_FNAME)
   except:
      print("problem with reading file: "+INPUT_FNAME)
      exit()
   return

def run_simulation():
   os.system(EXEC_NAME)
   return


def get_grafting_density():
   try:
      line = os.popen("""grep "Grafting density     = " """ + LOG_NAME).read()
      graf_dens = float(line.split()[-1])
   except:
      print("problem with reading the grafting density from the output file")
      exit()
   return graf_dens

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

#################################################################################
# This function performs a simulation
#################################################################################
def ff(x, target, copy_field):
   if copy_field:
      copy_last_field()
   set_grafted_initial_condition( x )
   run_simulation()
   return get_grafting_density() - target
#################################################################################
# Bisection minimization
#################################################################################
def bisection(ff, copy_field, xx_in, yy_in, yy_target, tol_bisection):

   xx = cp.deepcopy(xx_in)
   yy = cp.deepcopy(yy_in)

   ii = 0

   # Find the lower and upper bounds
   while True:
      ii += 1
      yy[1] = ff(xx[1], yy_target, copy_field)
      print("bisection proconditioning:", ii, xx[0], xx[1], yy[0], yy[1], abs(yy[1]) - tol_bisection)

      if yy[1] < 0.0:
         xx[1] *= 2.0
      else:
         break

   # Perform the bisection
   while True:
      ii += 1
   
      xx_next = 0.5 * (xx[0] + xx[1])
      yy_next = ff(xx_next, yy_target, copy_field)

      if yy_next < 0:
         xx[0] = xx_next
         yy[0] = yy_next
      else:
         xx[1] = xx_next
         yy[1] = yy_next

      #LOG.write(str(ii)+" "+str(xx[0])+" "+str(xx[1])+" "+str(yy[0])+" "+str(yy[1])+"\n")
      print("bisection details:", ii, xx[0], xx[1], yy[0], yy[1], abs(yy[1]) - tol_bisection)

      if (abs(yy_next) < tol_bisection):
         break

   return xx, yy

#################################################################################
# Secant minimization
#################################################################################
def secant(ff, copy_field, xx_in, yy_in, yy_target, tol_secant):

   xx = cp.deepcopy(xx_in)
   yy = cp.deepcopy(yy_in)

   ii = 0
   while True:
      ii += 1

      dx = -yy[1] * (xx[1] - xx[0]) / (yy[1] - yy[0])
      xx_next = xx[1] + dx

      xx[0] = xx[1]
      yy[0] = yy[1]

      xx[1] = xx_next


      yy[1] = ff(xx[1], yy_target, copy_field)

      print("secant details:", ii, xx[0], xx[1], yy[0], yy[1], abs(yy[1]) - tol_secant)

      if (abs(yy[1]) < tol_secant):
         break

   return xx, yy

xx = [0.0, x_init]
yy = [0, 0]

read_field(False)
copy_field = False

yy[0] = -yy_target
yy[1] = ff(xx[1], yy_target, copy_field)
read_field(True)
copy_field = True

xx, yy = bisection(ff, copy_field, xx, yy, yy_target, tol_bisection)
xx, yy = secant(ff, copy_field, xx, yy, yy_target, tol_secant)

quit()
