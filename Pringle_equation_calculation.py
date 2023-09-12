import numpy as np
from Pringle_module import General_Pringle


# Степень вязкости
gamma = 1

# Коэффициент диффузии
def Diff(u, x, t):
    D = x**(gamma)
    return D



# Начальные условия
def u_init(x):
    if x <= 100:
        u_init = 100*x**(-3/8)
    else:
        u_init = 0    
    return u_init

# Граничные условия слева
def u_left(u, x, t):
    u_left = u[1]
    return u_left


# Граничные условия справа
def u_right(u, x, t):
    u_right = 0
    return u_right


# Внутренняя граница
R_in = -1

# Внешнаяя граница
R_out = 3

# Предельная ошибка
epsilon = 10**(-5)

# Количество узлов
N = 100

#  Шаг по времени
tau = 10**(-3)

# Расчет для различных t_end
T = [0, 0.01, 0.1, 1]
for i in range(len(T)):
    t_end = T[i]
    u = General_Pringle(Diff, u_init, u_left, u_right,
                        t_end, R_in, R_out, tau, N, epsilon)
    np.savetxt(f'Real_T_{t_end}_N_{N}.txt', u, fmt='%.3e', delimiter=' ')



# u = General_Pringle(Diff, u_init, u_left, u_right,
                    # t_end, R_in, R_out, tau, N, epsilon)
# np.savetxt(f'Real_T_{t_end}_N_{N}.txt', u, fmt='%.3e', delimiter=' ')


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
