from sir_with_testing import getPoints
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

#def getPoints(days, iterations, S0, I0, R0, N, beta, gamma, rho)
days = 150
iterations = 150

S0 = 999
I0 = 1
R0 = 0
N = S0 + I0 + R0

beta = 0.5 #rate of contact between susceptible and infectives
gamma = 0.3 #constant rate at which infectives recover/die
rho = 0.9 #social distancing parameter btwn 0-1

#S, I, R = getPoints(days, iterations, S0, I0, R0, N, beta, gamma, rho)
time = np.linspace(0,days, iterations)

fig = go.Figure()
for step in np.arange(0,1,0.1):
    S, I, R = getPoints(days, iterations, S0, I0, R0, N, beta, gamma, step)
    fig.add_trace(
        go.Scatter(x= time, y = I)
    )


steps = []
for i in range(len(fig.data)):
    step = dict(method = "update",
                args=[{"visible": [False] * len(fig.data)},
                {"title": "Slider switched to step: " + str(i)}])
    step["args"][0]["visible"][i] = True
    steps.append(step)

sliders = [dict(
    active=10,
    currentvalue={"prefix": "Frequency: "},
    pad={"t": 50},
    steps=steps)]

fig.update_layout(
    sliders=sliders
)

fig.show()