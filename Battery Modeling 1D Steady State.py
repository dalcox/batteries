
# coding: utf-8

# In[1]:

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
get_ipython().magic('matplotlib inline')


# In[2]:

# initializing constants
K = 1.
s = 1.
ao = 1.
io = 1.
L = 1.
n = 1
F = 96485
R = 8.314
T = 298
I = 1.


#graph analytically
X = np.linspace(0., L, 100)
Y=y = X/L
v = L*np.sqrt(ao*io*(n*F)/(R*T)*(K + s)/(K*s))
i2 = I * K/(K + s)*(1 + (s*(K**-1)*np.sinh(v*(1-y)) - np.sinh(v*y))/np.sinh(v))


# In[3]:

plt.plot(X, i2)


# In[14]:

#solve numerically
def simplebattfunc(i, x):
    i0, i1 = i
    di = i1
    d2i = ao*io*(n*F)/(R*T)*(-I/s + i0*(1/s + 1/K))
    return di, d2i

def battfunc(IV, x):
    i1, i2, V1, V2 = IV
    di2 = ao*io*(n*F)/(R*T)*(V1 - V2)
    di1 = -di2
    dV1 = -i1/s
    dV2 = -i2/K
    return di1, di2, dV1, dV2


# In[23]:

t = np.linspace(0., 1., 100)
batt = odeint(battfunc, [0, I, 0, 0.01], t)
plt.plot(t, batt)


# In[9]:

from scipy.optimize import fsolve

u1_0 = I
def objective(u2_0):
    """
    The thing we want to set equal to zero
    """
    U = odeint(simplebattfunc, [u1_0, u2_0], t)
    print(U[-1,0])
    return U[-1,0]

u2_0, = fsolve(objective, 0)
print(u2_0)

i = odeint(simplebattfunc, [u1_0, u2_0], t)
i1 = I - i[:,0]
plt.plot(t, i[:,0], label = 'i2 - ionic')
plt.plot(t, i1, label = 'i1 - electronic')
plt.legend(loc = 'best')


# In[ ]:



