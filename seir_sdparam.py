import math
import numpy as np
import matplotlib.pyplot as plt

S0 = 999
E0 = 1
I0 = 0
R0 = 0
N = S0 + E0 + I0 + R0


S_norm = S0/N
E_norm = E0/N
I_norm = I0/N
R_norm = R0/N




beta = 0.8 #rate of contact between susceptible and infectives
gamma = 0.3 #constant rate at which infectives recover/die
sigma = 0.5 #rate at which incubated population becomes infected
rho = 1 #social distancing parameter

q = beta/gamma #contact ratio --> fraction of population that comes in contact with infected individual during infection
rr = beta*S0/gamma #reproductive rate --> secondary infections caused by single primary infection approx 3-4 with covid at beginning of pandemic


#using euler's method to approximate solutions to S,I,R model
def getPoints(days, iterations):
    S = []
    E = []
    I = []
    R = []

    S_current = S_norm
    I_current = I_norm
    E_current = E_norm
    R_current = R_norm

    delta_t = days/iterations

    for i in range(iterations):
        S_new = S_current + -1*rho*beta*S_current*I_current*delta_t
        E_new = E_current + (beta*rho*S_current*I_current -1*sigma*E_current)*delta_t
        I_new = I_current + (sigma*E_current - gamma*I_current)*delta_t
        R_new = R_current + gamma*I_current*delta_t

        S.append(S_new), E.append(E_new), I.append(I_new), R.append(R_new)

        S_current, E_current, I_current, R_current = S_new, E_new, I_new, R_new

    return S,E,I,R

def main():
    days = 300
    iterations = 150

    S,E, I, R = getPoints(days, iterations)
    time = np.linspace(0, days, iterations)

    plt.plot(time, S, label = "susceptible")
    plt.plot(time, E, label = "exposed")
    plt.plot(time, I, label = "infected")
    plt.plot(time, R, label = "recovered")
    plt.xlabel("time (days)"), plt.ylabel("# Individuals")
    plt.legend()
    plt.show()

main()