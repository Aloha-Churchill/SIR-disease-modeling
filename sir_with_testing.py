import math
import numpy as np
import matplotlib.pyplot as plt

S0 = 999
I0 = 1
R0 = 0
N = S0 + I0 + R0

beta = 0.9 #rate of contact between susceptible and infectives
gamma = 0.1 #constant rate at which infectives recover/die

q = beta/gamma #contact ratio --> fraction of population that comes in contact with infected individual during infection
rr = beta*S0/gamma #reproductive rate --> secondary infections caused by single primary infection approx 3-4 with covid at beginning of pandemic


#using euler's method to approximate solutions to S,I,R model
def getPoints(days, iterations, S0, I0, R0, N, beta, gamma):
    S = []
    I = []
    R = []

    S_current = S0
    I_current = I0
    R_current = R0

    delta_t = days/iterations

    for i in range(iterations):
        S_new = S_current + -1*beta*S_current*I_current*delta_t*(1/N)
        I_new = I_current + (beta*S_current*I_current*(1/N) - gamma*I_current)*delta_t
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

def samplePoints(iterations, n, I):
    #this function returns what data (no. positive results) you get from randomly sampling n people on a given day
    #could also begin to factor in that tests are not perfect (get false +, false -)
    #distribution of positive results follows binomial distribution X~Binomial(n, p = I_curr/N)
    #using numpy.random.binomial, but don't know if this is correct approach

    positives = []
    for i in range(iterations):
        p = I[i]/N
        X = np.random.binomial(n, p)
        positives.append(X) #X*100 to show on graph

    return positives

def getApproximation(time, iterations, I):
    exp_fun = []
    for i in range(iterations):
        exp_fun.append(I0*math.exp(time[i]*(beta-gamma)))

    #calculating r^2
    correlation_matrix = np.corrcoef(I, exp_fun)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2
    print(r_squared)

    return exp_fun


def main():

    days = 15
    iterations = 30

    S, I, R = getPoints(days, iterations, S0, I0, R0, N, beta, gamma)
    time = np.linspace(0, days, iterations)

    approx_func = getApproximation(time, iterations, I)

    #plt.plot(time, S, label = "susceptible")
    plt.plot(time, I, label = "infected")
    plt.plot(time, R, label = "recovered")
    plt.plot(time, approx_func, label = "Estimated # of Infectives")


    test_results = samplePoints(iterations, 20, I)
    #as n->infinity, converges to differential curve of i
    #plt.plot(time, test_results, label = "positive tests")


    plt.xlabel("time (days)"), plt.ylabel("# Individuals ")
    plt.legend()
    plt.show()




main()