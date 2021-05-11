# SIR Disease Modeling

This repository contains code for modeling various versions of the SIR model for infectious diseases. 

**Introduction**
Important constants are q and R0. q is the contact ratio or fraction of population that comes in contact with infected individuals during periods when they are infectious. R0 is the reproductive ratio which is essentially the secondary infections caused by single primary infection.

**Variations on SIR Model - SEIR(S)**
E is the exposed/incubation period where an individual is infected but not contagious
In SEIR model, there is one additional parameter (sigma)  which is the rate at which exposed people become infectious
In SEIRS model, there is another parameter (xi) which is the rate that people lose immunity and go back into the susceptible category

**SEIR with vital dynamics**
Introducing vital dynamics (births and deaths) mu = birth rate and nu = death rate.

**Discrete SIR Model**
https://kingaa.github.io/short-course/stochsim/stochsim.html
sir_binomial.py

**Estimating Rates**
In practice, limited to data and difficult to estimate parameters (reproduction number, transmission rate)
Parameters can be biologically predicted (due to molecular biology of disease) but are there ways to derive it given data?

***Estimating Reproduction Number (R0)***
Most epidemics can be approximated by exponential function during initial phase of epidemic.
sir_with_testing.py

Approx 0.2 for lambda. Then using data from https://www.who.int/bulletin/online_first/20-255695.pdf, gamma is equal to approx 1/(period of infectiousness) where the period of infectiousness was approx 7.5 days. Using the equation for R0, I got the overall reproduction number to be around 3 for Boulder County. For COVID, it's estimated around 2-3.5.

