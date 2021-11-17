import numpy as np
from myst_nb import glue
from pint import UnitRegistry
ureg = UnitRegistry()
ureg.default_format = '~L'

import sympy.physics.units as unit
from sympy.physics.units import speed_of_light, meter, gram, second, ampere
from sympy import *
init_printing(use_latex='mathjax')

import pandas as pd

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 5] #Plotgröße anpassen
import matplotlib
matplotlib.rcParams['text.usetex'] = True
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
matplotlib.rc('font', **font)

from IPython.display import display, Markdown, Latex, Math, HTML


def assign_meta_data(dic,symbol,entry):
    try:
        dic[symbol] = entry
    except:
        dic = dic | {symbol:entry}
    return dic

def assign_numpy_function(vals,equation,variable):
    f = lambdify(variable, equation.rhs.subs(vals))
    return assign_meta_data(np_functions,equation.lhs,f)

def show_numerical_value(vals,symbol):
    #display(float(vals[symbol]))
    #display(Latex(' ' + latex(symbol) + ' = '))
    text_template = ' %.2e'
    #display(Latex(' ' + latex(symbol) + '  =  \\textrm{' + (text_template % vals[symbol]) + '}' ) )
    try:
        val = vals[symbol].evalf()
    except:
        val = vals[symbol]

    #display(eq1)
    display(Latex('$ ' + latex(symbol) + ' =  ' + str(val) + units[symbol] + ' $'))
    
    
def show_legend_entry(sym_legend,symbol):
    display(Latex(' ' + latex(symbol) + '  ...  \\textrm{' + sym_legend[symbol] + '}' ) )

vals = {}
sym_legend = {}
np_functions = {}
