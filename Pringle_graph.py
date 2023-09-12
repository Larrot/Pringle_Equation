from matplotlib import pyplot as plt
import numpy as np

# Степень вязкости

gamma = -1
#  Внутреняя граница
R_in = -1
# R_in = 0.1

# Внешняя граница
R_out = 3
# R_out = 1000


# Количество узлов
N = 100

# Пространственная сетка
# x = np.linspace(R_in, R_out, N)
x = np.logspace(R_in, R_out, N)

# Данные расчета
data0 = np.loadtxt(f"Real_T_0_N_{N}.txt")
data1 = np.loadtxt(f"Real_T_0.01_N_{N}.txt")
data2 = np.loadtxt(f"Real_T_0.1_N_{N}.txt")
data3 = np.loadtxt(f"Real_T_1_N_{N}.txt")



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)


ax.set_xlim(0.1, 1000)
# ax.set_ylim(10**(1), 3*10**(2))
# ax.set_ylim(10**(1), 3*10**(2))
ax.set_ylim(10**(-5), 10**(6))
ax.set_xlabel(r'$r/r_0$ ')
# ax.set_xlabel(r'$\tau$ ')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r'$\Sigma,~{г}/{см}^2$ ')
ax.set_title(fr'$\Sigma$ для различных $t. N={N}$. $\nu \sim$'+ '$r^{%s}$' %gamma)




ax.plot(x, data0, '--', color='r', label =r'$t=0$')
ax.plot(x, data1, '--', color='b', label =r'$t=0.01$')
ax.plot(x, data2, '--', color='g', label =r'$t=0.1$')
ax.plot(x, data3, '--', color='k', label =r'$t=1$')





ax.legend(loc='best')
fig.show()
plt.tight_layout()
fig.savefig(f"Pringle_Sigma_R_{gamma}_N_{N}.png", orientation='landscape', dpi=300)
# fig.savefig(f"Pringle_Bz_1_log_{N}.png", dpi=300)
#fig.savefig(f"Pringle_Phi_1_r_log_{N}.png", orientation='landscape', dpi=300)
# fig.savefig(f"Pringle_Phi_1_total_log_40.png", orientation='landscape', dpi=300)
