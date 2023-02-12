import numpy as np
import matplotlib.pyplot as plt
from math import log

def convergence_plot(h,error,filename):
  
  plt.loglog(h, error,'-b',label = 'error',marker = 'o')
  plt.xlabel(r'$h$')
  plt.ylabel(r'$E(h)$')
  plt.savefig(filename, dpi=500)
  plt.clf()

def convergence_rate(error):
  for i in range(1,np.size(error)):
    covergence_rate = (log(error[i-1]) - log(error[i]))/log(2)
    print(covergence_rate)



#load h and error abs(p - pexact)
h, error = np.loadtxt("error.dat", unpack=True)
convergence_plot(h,error,"error.png")


#rate of convergence:
print("Rate of convergence")
convergence_rate(error)
