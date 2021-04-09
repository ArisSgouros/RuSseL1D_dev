import os
import sys
import math as m
import copy as cp
from scipy import optimize
#
# Set the parameters of the minimization
#
g_target_grafting_density = 0.002

# energy optimization
q_lo_try = m.pow(10,-60)
q_hi_try = m.pow(10,10)
q_lo_min = m.pow(10,-200)
q_hi_max = m.pow(10,50)
lx_lo = 90.0
lx_hi = 100.0
gss_tol = 1.e-3

# optimization for given lx
#q_lo_try = 0.001
#q_hi_try = 10000.0
#q_lo_min = m.pow(10,-170)
#q_hi_max = m.pow(10,170)

# convergence settings
g_bisection_y_tol   = 0.000000001
g_bisection_x_ratio = 1.0
g_secant_y_tol      = 0.0000001
g_secant_x_tol      = 0.0000001


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

def set_grafted_initial_condition(q_init):
   try:
      str_q_init = "%10.9e" % q_init

      os.system("sed -i '/! initial value of grafting point/c\ "+str_q_init+"              ! initial value of grafting point' "+INPUT_FNAME)
   except:
      print("problem with reading file: "+INPUT_FNAME)
      exit()
   return

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

def run_simulation():
   if EXEC_LOG:
      os.system(EXEC_NAME+" >> LOG")
   else:
      os.system(EXEC_NAME)
   return

def get_energy():
   try:
      energy_line = 2
      energy_column = 5
      line = os.popen("sed -n " + str(energy_line)+"p " + ENERGIES_FNAME).read()
      energy = float(line.split()[energy_column-1])
   except:
      print("problem with reading the " + ENERGIES_FNAME + " output file")
      exit()
   return energy

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

#################################################################################
# This function returns the grafting density for given grafting initial condition
#################################################################################
def get_gdens_for_given_qinit(x):

   global g_target_grafting_density

   if x <= 0:
      print("NEGATIVE Q_INIT!!!")
      quit()
      return 0.0 - g_target_grafting_density

   set_grafted_initial_condition( x )
   run_simulation()
   copy_last_field()

   res = get_grafting_density() - g_target_grafting_density

   global g_sim_evals
   g_sim_evals += 1

   print("         eval:",g_sim_evals, ", f(", x,") = ", res)
   return res

def get_bounds(function, args):

   q_lo = args["q_lo_try"]
   q_hi = args["q_hi_try"]
   q_hi_max = args["q_hi_max"]
   q_lo_min = args["q_lo_min"]

   (q_lo, q_hi) = (min(q_lo, q_hi), max(q_lo, q_hi))
   x_init = m.sqrt(q_hi*q_lo)
   scale = m.sqrt(q_hi/q_lo)

   #if abs(q_lo-q_hi) < 0.000001:
   if q_lo == q_hi:
      if q_hi >= q_hi_max:
         status = "LOWER_BOUND_EXCEEDED"
         print("         q_lo = q_hi >= q_hi_max")
         set_grafted_initial_condition( q_hi )
         return [],[], status
      elif q_lo <= q_lo_min:
         status = "HIGHER_BOUND_EXCEEDED"
         print("         q_lo = q_hi <= q_lo_min")
         set_grafted_initial_condition( q_lo )
         return [],[], status
      else:
         print("Problem in get_bounds: q_lo = q_hi != q_lo_min || != q_hi_max")
         exit()

   # find if the initial value is larger or smaller than the target value
   read_field(False)
   y_init = function(x_init)
   #read_field(False)
   read_field(True)

   x_range = [0, 0]
   y_range = [0, 0]

   is_lo_bound = False
   is_hi_bound = False

   # Generate the upper/lower bound
   if y_init < 0.0:
      x_range[0] = x_init
      y_range[0] = y_init
      is_lo_bound = True
   elif y_init > 0.0:
      x_range[1] = x_init
      y_range[1] = y_init
      is_hi_bound = True
      scale = 1.0 / scale

   # Trace the lower/upper bound
   xx = x_init
   while True:

      xx = xx * scale
      yy = function(xx)

      if is_lo_bound:
         if yy < 0.0:
            x_range[0] = xx
            y_range[0] = yy
         else:
            x_range[1] = xx
            y_range[1] = yy
            print("         x_range = ", x_range)
            break

      if is_hi_bound:
         if yy > 0.0:
            x_range[1] = xx
            y_range[1] = yy
         else:
            x_range[0] = xx
            y_range[0] = yy
            print("         x_range = ", x_range)
            break

      if xx >= q_hi_max:
         status = "LOWER_BOUND_EXCEEDED"
         return x_range, y_range, status
      if xx <= q_lo_min:
         status = "HIGHER_BOUND_EXCEEDED"
         return x_range, y_range, status

   status = "success"
   return [x_range, y_range, status]

def bisection_geom(function, x_range, y_range, x_ratio, y_tol):
   x_range = cp.deepcopy(x_range)
   y_range = cp.deepcopy(y_range)
   while True:

      ratio = x_range[1] / x_range[0] - 1
      if ratio < x_ratio:
         status = "x_ratio"
         print("         exiting bisection due to xhi/xlo = ", ratio, " < ", x_ratio )
         break

      xx_next = m.sqrt( (x_range[0] * x_range[1]) )
      yy_next = get_gdens_for_given_qinit(xx_next)
      if yy_next < 0:
         x_range[0] = xx_next
         y_range[0] = yy_next
      else:
         x_range[1] = xx_next
         y_range[1] = yy_next
      print("         x_range = ", x_range)
      if (abs(yy_next) < y_tol):
         status = "y_tol"
         break
   return [x_range, y_range, status]

def bisection_arith(function, x_range, y_range, y_tol):
   x_range = cp.deepcopy(x_range)
   y_range = cp.deepcopy(y_range)
   while True:
      xx_next = m.mean(x_range[0],x_range[1])
      yy_next = get_gdens_for_given_qinit(xx_next)
      if yy_next < 0:
         x_range[0] = xx_next
         y_range[0] = yy_next
      else:
         x_range[1] = xx_next
         y_range[1] = yy_next
      print("         x_range = ", x_range)
      if (abs(yy_next) < y_tol):
         status = "y_tol"
         break
   return [x_range, y_range]

def secant(get_gdens_for_given_qinit, xx_in, yy_in, x_tol, y_tol):
   xx = cp.deepcopy(xx_in)
   yy = cp.deepcopy(yy_in)
   while True:
      dx = -yy[1] * (xx[1] - xx[0]) / (yy[1] - yy[0])
      xx_next = xx[1] + dx
      xx[0] = xx[1]
      yy[0] = yy[1]
      xx[1] = xx_next
      yy[1] = get_gdens_for_given_qinit(xx[1])
      print("         x_range = ", xx)
      if (abs(yy[1]) < y_tol):
         print("         Secant converged: x_range = ", xx)
         status = "success"
         break
      if (abs(xx[1] - xx[0]) < x_tol):
         print("         Secant converged: x_range = ", xx)
         status = "success"
         break
   return [xx, yy, status]


def get_energy_and_qinit_for_given_lx(lx, args):

   set_lx(lx)

   global g_bisection_y_tol, g_bisection_x_ratio, g_secant_y_tol, g_secant_x_tol

   print("      Getting bounds")
   [x_range, y_range, status] = get_bounds(get_gdens_for_given_qinit, args)

   # In case the initial value is such that lx < lx_min give the energy a high value that increases linearly with decreasing lx
   if status == "LOWER_BOUND_EXCEEDED":
      print("     ", status)
      energy = LARGE_NUMBER-lx
   # In case the initial value is such that lx > lx_max give the energy a high value that increases linearly with increasing lx
   elif status == "HIGHER_BOUND_EXCEEDED":
      print("     ", status)
      energy = LARGE_NUMBER+lx
   else:
      print("      Geometric bisection")
      [x_range, y_range, status] = bisection_geom(get_gdens_for_given_qinit, x_range, y_range, g_bisection_x_ratio, g_bisection_y_tol)
      print("      Secant optimization")
      [x_range, y_range, status] = secant(get_gdens_for_given_qinit, x_range, y_range, g_secant_x_tol, g_secant_y_tol)
      energy = get_energy()

   q_final = get_grafted_initial_condition()
   return [energy, q_final]




def gss(f, args, a, b, tol=1e-5):

    print("\n*Golden section search*\n")

    g_lx_evals = 0

    foo = open("pylog.out","w")
    foo.write(("%-20s "*15+"\n") % ("lx_iter","sim_evals","lx","q_lo_try","q_hi_try","q_final","energy","la","lc","ld","lb","qa","qc","qd","qb"))
    format = "%-20d "*2+"%-20.9e "*13+"\n"


    invphi = (m.sqrt(5) - 1) / 2  # 1 / phi
    invphi2 = (3 - m.sqrt(5)) / 2  # 1 / phi^2

    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return (a, b)

    # Required steps to achieve tolerance
    n = int(m.ceil(m.log(tol / h) / m.log(invphi)))
    c = a + invphi2 * h
    d = a + invphi * h

    q_of_L = {}
    q_of_L['a'] = args["q_hi_try"]
    q_of_L['b'] = args["q_lo_try"]

    [yc, qc] = f(c,args)
    g_lx_evals += 1
    q_of_L['c'] = qc
    foo.write(format % (g_lx_evals,g_sim_evals,c,args["q_lo_try"],args["q_hi_try"],qc,yc, a,c,d,b,q_of_L['a'],q_of_L['c'],0.0        ,q_of_L['b']))
    foo.flush()
    print("   it=",g_lx_evals," sim",g_sim_evals," lx=",c," q_lo=",args["q_lo_try"]," q_hi=",args["q_hi_try"]," q_target=",qc," eng=",yc," qa=",q_of_L['a']," qc=",q_of_L['c']," qd=",0.00," qb=",q_of_L['b'])

    [yd, qd] = f(d,args)
    g_lx_evals += 1
    q_of_L['d'] = qd
    foo.write(format % (g_lx_evals,g_sim_evals,d,args["q_lo_try"],args["q_hi_try"],qd,yd, a,c,d,b,q_of_L['a'],q_of_L['c'],q_of_L['d'],q_of_L['b']))
    foo.flush()
    print("   it=",g_lx_evals," sim",g_sim_evals," lx=",d," q_lo=",args["q_lo_try"]," q_hi=",args["q_hi_try"]," q_target=",qd," eng=",yd," qa=",q_of_L['a']," qc=",q_of_L['c']," qd=",q_of_L['d']," qb=",q_of_L['b'])

    for k in range(n-1):
        if yc < yd:
            b = d
            d = c
            q_of_L['b'] = q_of_L['d']
            q_of_L['d'] = q_of_L['c']
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            args["q_lo_try"] = q_of_L['d']
            args["q_hi_try"] = q_of_L['a']
            [yc, qc] = f(c,args)
            g_lx_evals += 1
            q_of_L['c'] = qc
            foo.write(format % (g_lx_evals,g_sim_evals,c,args["q_lo_try"],args["q_hi_try"],qc,yc, a,c,d,b,q_of_L['a'],q_of_L['c'],q_of_L['d'],q_of_L['b']))
            foo.flush()
            print("   it=",g_lx_evals," sim",g_sim_evals," lx=",c," q_lo=",args["q_lo_try"]," q_hi=",args["q_hi_try"]," q_target=",qc," eng=",yc," qa=",q_of_L['a']," qc=",q_of_L['c']," qd=",q_of_L['d']," qb=",q_of_L['b'])
        else:
            a = c
            c = d
            q_of_L['a'] = q_of_L['c']
            q_of_L['c'] = q_of_L['d']
            yc = yd
            h = invphi * h
            d = a + invphi * h
            args["q_lo_try"] = q_of_L['b']
            args["q_hi_try"] = q_of_L['c']
            [yd, qd] = f(d,args)
            g_lx_evals += 1
            q_of_L['d'] = qd
            foo.write(format % (g_lx_evals,g_sim_evals,d,args["q_lo_try"],args["q_hi_try"],qd,yd, a,c,d,b,q_of_L['a'],q_of_L['c'],q_of_L['d'],q_of_L['b']))
            foo.flush()
            print("   it=",g_lx_evals," sim",g_sim_evals," lx=",d," q_lo=",args["q_lo_try"]," q_hi=",args["q_hi_try"]," q_target=",qd," eng=",yd," qa=",q_of_L['a']," qc=",q_of_L['c']," qd=",q_of_L['d']," qb=",q_of_L['b'])

    if yc < yd:
        return (a, d)
    else:
        return (c, b)

# Add the parameters of the optimized function into a dictionary
args = {}
args["q_lo_try"]=q_lo_try
args["q_hi_try"]=q_hi_try
args["q_lo_min"]=q_lo_min
args["q_hi_max"]=q_hi_max

#gss(get_energy_and_qinit_for_given_lx, args, lx_lo, lx_hi, gss_tol)

get_energy_and_qinit_for_given_lx(12,args)

#print(  get_energy_and_qinit_for_given_lx(8.0,args)  )
