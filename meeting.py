import math
import numpy as np
import matplotlib.pyplot as plt

S0 = 1
I0 = 1.27*math.pow(10, -6)

beta1 = 0.5 #rate of contact between susceptible and infectives
beta2 = 1
gamma = 0.6 #constant rate at which infectives recover/die

#q = beta/gamma #contact ratio --> fraction of population that comes in contact with infected individual during infection
#rr = beta*S0/gamma #reproductive rate --> secondary infections caused by single primary infection approx 3-4 with covid at beginning of pandemic



#using euler's method to approximate solutions to S,I,R model
def getPoints(days, iterations):
    S_low = []
    S_high = []
    I = []
    R = []

    S_low_current = S0/2
    S_high_current = S0/2
    I_current = I0
    R_current = 0

    delta_t = days/iterations

    for i in range(iterations):
        S_new_low = S_low_current + -1*beta1*S_low_current*I_current*delta_t
        S_new_high = S_high_current + -1*beta2*S_high_current*I_current*delta_t
        I_new = I_current + (beta1*S_low_current*I_current + beta2*S_high_current*I_current - gamma*I_current)*delta_t
        R_new = R_current + gamma*I_current*delta_t

        S_low.append(S_new_low)
        S_high.append(S_new_high)
        I.append(I_new)
        R.append(R_new)

        S_low_current = S_new_low
        S_high_current = S_new_high
        I_current = I_new
        R_current = R_new

    return S_low, S_high ,I,R

"""
def maxInfectives():
    I_max = I0 + S0 - (1/q)*(1+np.log(q*S0))
    return(I_max)
"""


def main():
    days = 300
    iterations = 150

    S_l, S_h, I, R = getPoints(days, iterations)
    time = np.linspace(0, days, iterations)

    """
    plt.plot(time, S_l, label = "low susceptible")
    plt.plot(time, S_h, label = "high susceptible")
    #plt.plot(time, np.array(S_l)+np.array(S_h), label = "total susceptible")
    plt.plot(time, I, label = "infected")
    plt.plot(time, R, label = "recovered")
    """
    #plt.plot(time, np.array(S_l)*beta1 + np.array(S_h)*beta2, label = "Average rate of S-->I")
    plt.plot(time, np.array(S_l)/(np.array(S_l)+np.array(S_h))*beta1 + np.array(S_h)/(np.array(S_l)+np.array(S_h))*beta2)

    plt.xlabel("time (days)"), plt.ylabel("Percentage of Population (%)")
    plt.legend()
    plt.show()

main()