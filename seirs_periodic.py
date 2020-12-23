import math
import numpy as np
import matplotlib.pyplot as plt

S0 = 999
E0 = 1
I0 = 0
R0 = 0
N0 = S0 + E0 + I0 + R0


S_norm = S0/N0
E_norm = E0/N0
I_norm = I0/N0
R_norm = R0/N0


beta = 0.9 #rate of contact between susceptible and infectives
gamma = 0.3 #constant rate at which infectives recover/die
sigma = 0.7 #rate at which incubated population becomes infected
xi = 0.9 #rate at which recovered loose immunity

#mu = 0.001 #birth rate
#nu = 0 #death rate



q = beta/gamma #contact ratio --> fraction of population that comes in contact with infected individual during infection
rr = beta*S0/gamma #reproductive rate --> secondary infections caused by single primary infection approx 3-4 with covid at beginning of pandemic


#using euler's method to approximate solutions to S,I,R model
#for more precision, use runge-kutta methods
def getPoints(days, iterations):
    S = []
    E = []
    I = []
    R = []
    Population = []

    S_current = S0
    I_current = I0
    E_current = E0
    R_current = R0
    N_current = N0

    delta_t = days/iterations

    for i in range(iterations):
        S.append(S_current), E.append(E_current), I.append(I_current), R.append(R_current), Population.append(N_current)

        # WITH VITAL DYNAMICS
        """
        S_new = S_current + (-1*beta*S_current*I_current*(1/N_current) - nu*S_current)*delta_t
        E_new = E_current + (beta*S_current*I_current*(1/N_current) - 1*sigma*E_current - nu*E_current)*delta_t
        I_new = I_current + (sigma*E_current - gamma*I_current - nu*I_current)*delta_t
        R_new = R_current + (gamma*I_current - nu*R_current)*delta_t
        N_new = S_new + E_new + I_new + R_new
        """

        #W/O VITAL DYNAMICS
        S_new = S_current + (-1*beta*S_current*I_current*(1/N_current) + xi*R_current)*delta_t
        E_new = E_current + (beta*S_current*I_current*(1/N_current) - 1*sigma*E_current)*delta_t
        I_new = I_current + (sigma*E_current - gamma*I_current)*delta_t
        R_new = R_current + (gamma*I_current - xi*R_current)*delta_t
        N_new = S_new + E_new + I_new + R_new

        S_current, E_current, I_current, R_current, N_current = S_new, E_new, I_new, R_new, N_new

    return S,E,I,R, Population

def main():
    days = 300
    iterations = 2000

    S,E, I, R, Population = getPoints(days, iterations)
    time = np.linspace(0, days, iterations)

    plt.plot(time, S, label = "susceptible")
    plt.plot(time, E, label = "exposed")
    plt.plot(time, I, label = "infected")
    plt.plot(time, R, label = "recovered")
    plt.plot(time, Population, label = "Population")
    plt.xlabel("time (days)"), plt.ylabel("# Individuals")
    plt.legend()
    plt.show()

main()