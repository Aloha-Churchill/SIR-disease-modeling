import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

df = pd.read_csv("CDPHE_COVID19_Daily_State_Statistics.csv")
I0 = df['Cases'][0]

def approximateLambda(n):
    l = []
    for i in range(1,n):
        I = df['Cases'][i]
        l.append(np.log(I/I0)/i)
    print(sum(l)/len(l))

def approximateExponential(days, l_approx):
    exp_fun = []
    for i in range(days):
        exp_fun.append(I0*math.exp(i*(l_approx)))

    #calculating r^2
    correlation_matrix = np.corrcoef(df['Cases'][0:n], exp_fun)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2
    print(r_squared)

    return exp_fun

n = 20
approximateLambda(n)
func = approximateExponential(n, 0.21)
plt.plot(np.arange(n), df['Cases'][0:n])
plt.plot(np.arange(n), func)

plt.show()