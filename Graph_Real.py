from matplotlib import pyplot as plt
import numpy as np
import scipy.integrate as integrate


R_in = -1

R_out = 2

# R_in = 0.1

# R_out = 100
N = 1000
# h = (R_out-R_in)/(N-1)

# x = np.linspace(R_in, R_out, N)
x = np.logspace(-1, 2, N)
# HH =np.logspace(R_in, R_out, N)

data0 = np.loadtxt(f"Real_T_0_N_{N}.txt")
data1 = np.loadtxt(f"Real_T_0.01_N_{N}.txt")
data2 = np.loadtxt(f"Real_T_0.1_N_{N}.txt")
data3 = np.loadtxt(f"Real_T_1_N_{N}.txt")




# PHI = np.loadtxt(f"Phi_T_10.0_N_100.txt")
# dt = np.linspace(0, 10, 40)


"""
error = np.zeros(N)
for i in range(N):
    error[i] = abs(data2[i]-y2[i])/y2[i]*100
"""


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.set_xlim(R_in, R_out)

# ax.set_xlim(R_in, R_out)
# ax.set_ylim(10**(1), 10**(2))
# ax.set_ylim(60, 10**(2))
ax.set_xlabel(r'$r/r_0$ ')
# ax.set_xlabel(r'$\tau$ ')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r'$\Sigma,~{г}/{см}^2$ ')
ax.set_title(r'$\Sigma$ для различных $t. N=1000$. $\nu \sim r^{3/8}$')

# # ax.set_title(u'$B_z$ для различнх $t. N=100$')
# ax.set_ylabel(r'$B_z, $ Гс ')

#ax.set_title(u'Полный магнитный поток $\Phi. N=20$')
# ax.set_ylabel(r'$\Phi,~10^{26}~{Мкс} $')

# ax.set_title(u'Относительная масса диска в зависимости от времени. ')

# ax.set_ylabel(r'$M_{disk}/M_{\odot} $')


# b0 = np.zeros(N)
# for i in range(N):
#     b0[i] = 100*data0[i]/(10**(2)*0.1**(-3/8))
# b1 = np.zeros(N)
# for i in range(N):
#     b1[i] = 100*data1[i]/(10**(2)*0.1**(-3/8))
# b2 = np.zeros(N)
# for i in range(N):
#     b2[i] = 100*data2[i]/(10**(2)*0.1**(-3/8))
# b3 = np.zeros(N)
# for i in range(N):
#     b3[i] = 100*data3[i]/(10**(2)*0.1**(-3/8))
    
# phi0 = np.zeros(N)
# for i in range(N):
#     phi0[i] = 2*np.pi*np.trapz(b0[:i], x[:i])

   
# phi1 = np.zeros(N)
# for i in range(N):
#     phi1[i] = 2*np.pi*np.trapz(b1[:i], x[:i])

    
# phi2 = np.zeros(N)
# for i in range(N):
#     phi2[i] = 2*np.pi*np.trapz(b2[:i], x[:i])

   
# phi3 = np.zeros(N)
# for i in range(N):
#     phi3[i] = 2*np.pi*np.trapz(b3[:i], x[:i])


# M_disk = np.zeros(40)

# for i in range(40):
#     M_disk[i] = 0.1-10**(-5)*dt[i]


"""
ax.plot(x, sigma1, '--', color='r', label =r'$t=0.001$')
ax.plot(x, sigma2, '--', color='g', label =r'$t=0.01$')
ax.plot(x, sigma3, '--', color='b', label =r'$t=0.1$')
"""


#ax.plot(x, y1, '-', color='k', label =r'$t=0$')
#ax.plot(x, y2, '-', color='g', label =r'$t=1$')
#ax.plot(x, y3, '-', color='r', label =r'$t=10$')

ax.plot(x, data0, '--', color='k', label =r'$t=0$')
ax.plot(x, data1, '--', color='b', label =r'$t=0.001$')
ax.plot(x, data2, '--', color='g', label =r'$t=0.01$')
ax.plot(x, data3, '--', color='r', label =r'$t=1$')
# # #ax.plot(x, data4, '--', color='y', label =r'$t=100$')



# t = np.arange(0, 2*np.pi, 0.01)
# r = 4
# plt.plot(r*np.sin(t), r*np.cos(t), lw=3)
# plt.axis('equal')
# plt.show()




# ax.plot(x, b0, '--', color='k', label =r'$t=0$')
# ax.plot(x, b1, '--', color='b', label =r'$t=0.1$')
# ax.plot(x, b2, '--', color='g', label =r'$t=1$')
# ax.plot(x, b3, '--', color='r', label =r'$t=10$')

# ax.plot(x[2:], phi0[2:], '--', color='k', label =r'$t=0$')
# ax.plot(x[2:], phi1[2:], '--', color='b', label =r'$t=0.1$')
# ax.plot(x[2:], phi2[2:], '--', color='g', label =r'$t=1$')
# ax.plot(x[2:], phi3[2:], '--', color='r', label =r'$t=10$')


# ax.plot(dt, PHI, '.-', color='k')

# ax.plot(dt, M_disk, '.-', color='k')


#ax.plot(x, error)
#print(max(error))

#ax.plot(x, data1, '--', color='b', label =r'$t=1$')


#ax.plot(x, sigma, '-.', color='g', label = '$\Sigma$')

#ax1=plt.errorbar(x[0:4], data2[0:4], xerr=h, fmt='o',  ecolor='red')
#ax1=plt.errorbar(x[4:], data2[4:], xerr=h,  ecolor='red')


ax.legend(loc='best')
fig.show()
plt.tight_layout()
# fig.savefig(f"Pringle_Sigma_R_m0.2.png", orientation='landscape', dpi=300)
# fig.savefig(f"Pringle_Bz_1_log_{N}.png", dpi=300)
#fig.savefig(f"Pringle_Phi_1_r_log_{N}.png", orientation='landscape', dpi=300)
#fig.savefig(f"Pringle_Phi_1_total_log_40.png", orientation='landscape', dpi=300)
