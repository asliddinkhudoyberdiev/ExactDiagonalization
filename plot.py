import numpy as np
import matplotlib.pyplot as plt

file = np.loadtxt('/home/asliddin/CPP/Exact_DIagonalization/Mit_Florian/Zeitentwicklung_Florian/output_time_evolution/Test_mit_N16_Spins.dat').T

time = file[0, :]
energy = file[1, :]
m_of_spin_0 = file[2, :]
m_sublattice_even = file[3, :]
m_sublattice_odd = file[4, :]

fig, ax = plt.subplots(1, 1)
ax.plot(time, energy, '-r')
ax.set_xlabel('Time t')
ax.set_ylabel('Energy')
ax.grid()
plt.show()

fig, ax = plt.subplots(1, 1)
ax.plot(time, m_of_spin_0, '-b')
ax.set_xlabel('Time t')
ax.set_ylabel('Magnetization of spin 0')
ax.grid()
plt.show()

fig, ax = plt.subplots(1, 1)
ax.plot(time, m_sublattice_even, '-g')
ax.set_xlabel('Time t')
ax.set_ylabel('Magnetization of one sublattice')
ax.grid()
plt.show()

fig, ax = plt.subplots(1, 1)
ax.plot(time, m_sublattice_odd, '-g')
ax.set_xlabel('Time t')
ax.set_ylabel('Magnetization of one sublattice')
ax.grid()
plt.show()
