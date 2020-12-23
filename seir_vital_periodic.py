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


beta = 0.5 #rate of contact between susceptible and infectives
gamma = 0.4 #constant rate at which infectives recover/die
sigma = 0.5 #rate at which incubated population becomes infected
#rho = 0.7 #social distancing parameter

mu = 0.094 #birth rate
nu = 0.04 #death rate



q = beta/gamma #contact ratio --> fraction of population that comes in contact with infected individual during infection
rr = beta*S0/gamma #reproductive rate --> secondary infections caused by single primary infection approx 3-4 with covid at beginning of pandemic


#using euler's method to approximate solutions to S,I,R model
#for more precision, use runge-kutta methods
def getPoints(days, iterations, cycles):
    S = []
    E = []
    I = []
    R = []
    Population = []

    S_current = S0
    I_current = I0
    E_current = E0
    R_current = R0

    delta_t = days/iterations
    for j in range(cycles):
        I_current += 100
        for i in range(iterations):
            S_new = S_current + (mu*N + -1*beta*S_current*I_current*(1/N) - nu*S_current)*delta_t
            E_new = E_current + (beta*S_current*I_current*(1/N) - 1*sigma*E_current - nu*E_current)*delta_t
            I_new = I_current + (sigma*E_current - gamma*I_current - nu*I_current)*delta_t
            R_new = R_current + (gamma*I_current - nu*R_current)*delta_t
            Population_new = S_new + E_new + I_new + R_new

            S.append(S_new), E.append(E_new), I.append(I_new), R.append(R_new), Population.append(Population_new)

            S_current, E_current, I_current, R_current = S_new, E_new, I_new, R_new



    return S,E,I,R, Population

def main():
    days = 900
    iterations = 1000
    num_cycles = 2

    S,E, I, R, Population = getPoints(days, iterations, num_cycles)
    time = np.linspace(0, days, iterations*num_cycles)

    plt.plot(time, S, label = "susceptible")
    plt.plot(time, E, label = "exposed")
    plt.plot(time, I, label = "infected")
    plt.plot(time, R, label = "recovered")
    plt.plot(time, Population, label = "Population")
    plt.xlabel("time (days)"), plt.ylabel("# Individuals")
    plt.legend()
    plt.show()

main()