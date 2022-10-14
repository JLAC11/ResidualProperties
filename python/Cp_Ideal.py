import scipy.integrate as integrate
from scipy.special import gamma
import numpy as np
# from sympy import *
import matplotlib.pyplot as plt

print(integrate.quad(lambda T: 29.862-0.009664*T+0.0001097*T**2-1.549E-7*T**3+1.561E-5*T**4, 698, 298.15)[0])
