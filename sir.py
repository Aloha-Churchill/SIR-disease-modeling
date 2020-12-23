import math
import numpy as np
import matplotlib.pyplot as plt

S0 = 1
I0 = 1.27*math.pow(10, -6)

beta = 0.8 #rate of contact between susceptible and infectives
gamma = 0.5 #constant rate at which infectives recover/die

q = beta/gamma #contact ratio --> fraction of population that comes in contact with infected individual during infection
rr = beta*S0/gamma #reproductive rate --> secondary infections caused by single primary infection approx 3-4 with covid at beginning of pandemic



#using euler's method to approximate solutions to S,I,R model
def getPoints(days, iterations):
    S = []
    I = []
    R = []

    S_current = S0
    I_current = I0
    R_current = 0

    delta_t = days/iterations

    for i in range(iterations):
        S_new = S_current + -1*beta*S_current*I_current*delta_t
        I_new = I_current + (beta*S_current*I_current - gamma*I_current)*delta_t
        R_new = R_current + gamma*I_current*delta_t

        S.append(S_new)
        I.append(I_new)
        R.append(R_new)

        S_current = S_new
        I_current = I_new
        R_current = R_new

    return S,I,R

def maxInfectives():
    I_max = I0 + S0 - (1/q)*(1+np.log(q*S0))
    return(I_max)

def main():
    days = 150
    iterations = 150

    S, I, R = getPoints(days, iterations)
    time = np.linspace(0, days, iterations)

    plt.plot(time, S, label = "susceptible")
    plt.plot(time, I, label = "infected")
    plt.plot(time, R, label = "recovered")
    plt.xlabel("time (days)"), plt.ylabel("Percentage of Population (%)")
    plt.legend()
    plt.show()

main()