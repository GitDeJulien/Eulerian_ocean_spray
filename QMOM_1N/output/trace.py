import numpy as np
import matplotlib.pyplot as plt

N_T = 20
N_r = 100

u_bar_1 = []
u_bar_2 = []
u_bar_3 = []
u_bar_4 = []
u_bar_5 = []
LOAD_TIME = np.loadtxt('other/tn.dat')
LOAD_TEMPERATURE = np.loadtxt('other/Temperature.mesh.dat')
LOAD_RADIUS = np.loadtxt('other/radius.mesh.dat')


r1 = LOAD_RADIUS[10]
r2 = LOAD_RADIUS[20]
r3 = LOAD_RADIUS[50]
r4 = LOAD_RADIUS[60]
r5 = LOAD_RADIUS[70]

for i in range(0,2001,1):
    LOAD_DAT = np.loadtxt(f'sol/sol.{i}.dat')
    
    u_bar_1.append(LOAD_DAT[10,2])
    u_bar_2.append(LOAD_DAT[20,2])
    u_bar_3.append(LOAD_DAT[50,2])
    u_bar_4.append(LOAD_DAT[60,2])
    u_bar_5.append(LOAD_DAT[70,2])


u_bar_1 = np.array(u_bar_1)
u_bar_2 = np.array(u_bar_2)
u_bar_3 = np.array(u_bar_3)
u_bar_4 = np.array(u_bar_4)
u_bar_5 = np.array(u_bar_5)

# Affichage
plt.figure()
plt.plot(LOAD_TIME, u_bar_1, ':b', label=f'r={r1:.3E}')
plt.plot(LOAD_TIME, u_bar_2, '--m', label=f'r={r2:.3E}')
plt.plot(LOAD_TIME, u_bar_3, '-.g', label=f'r={r3:.3E}')
plt.plot(LOAD_TIME, u_bar_4, ':k', label=f'r={r4:.3E}')
plt.plot(LOAD_TIME, u_bar_5, '--r', label=f'r={r5:.3E}')
# plt.plot(maillage, mean_value_2 , '-^r', label='HLL (L1 error)')
# plt.plot(maillage, mean_value_3 , '-Dg', label='VFRoe (L1 error)')
# plt.plot(maillage, np.exp(pred_f_1), ':b', label=f'pente {-coefs_1[0]:.3f}')
# plt.plot(maillage, np.exp(pred_f_2), ':r', label=f'pente {-coefs_2[0]:.3f}')
# plt.plot(maillage, np.exp(pred_f_3), ':g', label=f'pente {-coefs_3[0]:.3f}')
plt.xlabel("$t(s)$")
plt.ylabel("$\\bar{u}(t,r)$")
# plt.xscale("log")
# plt.yscale("log")
# plt.xlim([maillage[0], maillage[-1]])
plt.legend()
#plt.show()
plt.savefig("u_bar_time.png")