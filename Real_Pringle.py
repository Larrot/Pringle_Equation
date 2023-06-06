import numpy as np
from General_Pringle import General_Pringle


def Diff(u, x, t):
    D = x**(3/8)
    # D = x**(1)
    # D = 0.1
    return D


def u_init(x):
    u_init = 100*x**(-3/8)
    # u_init = 1*x**(-3/8)
    return u_init


def u_left(u, x, t):
    u_left = u[1]
    return u_left


# def u_left(u, x, t):
#     u_left = 0.1
#     return u_left


def u_right(u, x, t):
    u_right = 100*100**(-3/8)
    # u_right = 17.78279410038923
    return u_right

# def u_right(u, x, t):
#     if x >= 200:
#         return 0
#     else:
#         return u[-1]


R_in = -1

# R_in = 0.1
# R_out = 100

R_out = 2

epsilon = 10**(-5)


N = 100

t_end = 1

tau = 10**(-3)

# T = np.arange(0, 0.1, 0.01)
T = [0, 0.01, 0.1, 1]
for i in range(len(T)):
    t_end = T[i]
    u = General_Pringle(Diff, u_init, u_left, u_right,
                        t_end, R_in, R_out, tau, N, epsilon)
    np.savetxt(f'Real_T_{t_end}_N_{N}.txt', u, fmt='%.3e', delimiter=' ')

# dt = np.linspace(0, 10, 40)
# t_end = 0
# tau = 10**(-3)
# b0 = np.zeros(N)
# phi0 = np.zeros(40)
# x = np.logspace(-1, 2, N)


# for i in range(40):
#     t_end = dt[i]
#     u = General_Pringle(Diff, u_init, u_left, u_right, t_end, R_in, R_out, tau, N, epsilon)
#     for j in range(N):
#         b0[j] = 100*u[j]/(10**(3)*0.1**(-3/8))
#     phi0[i] = 2*np.pi*np.trapz(b0, x)

# np.savetxt(f'Phi_T_{t_end}_N_{N}.txt', phi0, fmt='%.3e', delimiter=' ')
