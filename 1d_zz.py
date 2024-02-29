import numpy as np 
import matplotlib.pyplot as plt

def poisson_time(a, b, u):
    if b > 0:
        if a < 0:
            return np.sqrt(-np.log(u)*2.0/b) - a/b
        else:  # a > 0
            return np.sqrt((a/b)**2 - np.log(u)*2.0/b) - a/b
    elif b == 0:
        if a > 0:
            return -np.log(u)/a
        else:
            return float('inf')
    else:  # b < 0
        if a <= 0:
            return float('inf')
        elif -np.log(u) <= -a**2/b + a**2/(2*b):
            return -np.sqrt((a/b)**2 - np.log(u)*2.0/b) - a/b
        else:
            return float('inf')

def pdmp_1d(x, v, post_mean, post_var, T):
    t = 0.0
    trace = [(t, x)]
    tau = poisson_time((x-post_mean)*v/post_var, v**2/post_var, np.random.rand())
    while t < T:
        x += tau*v
        v *= -1  
        t += tau
        trace.append((t, x, v))
        tau = poisson_time((x-post_mean)*v/post_var, v**2/post_var, np.random.rand())
    return trace

# Initial conditions and parameters
x_initial = 0.0
v_initial = 1
post_mean = 5.0
post_var = 2.0
T = 100.0

# Generate data using the pdmp function
xx = pdmp_1d(x_initial, v_initial, post_mean, post_var, T)

# Plotting
plt.plot([entry[0] for entry in xx], [entry[1] for entry in xx])
plt.xlabel('Time')
plt.ylabel('Position')
plt.title('Position vs Time')
plt.grid(True)
plt.show()
