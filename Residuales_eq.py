import numpy as np

compound = ["CO", "CO2", "CH4", "H2"]
y_i = [0.60175, 0.34575, 0.03375, 0.01875]
W = [1.0, 1.0, 1.0, 1.0]
Tr = [1.0, 1.0, 1.0, 1.0]
Pr = [1.0, 1.0, 1.0, 1.0]

A_i = np.zeros((1, len(compound)))
B_i = np.zeros((1, len(compound)))
alfa_i = np.zeros((1, len(compound)))

print(alfa_i)
for i in range(len(compound)):
    alfa_i[0, i] = ((1.0+(0.480+1.574*W[i]-0.176*W[i]**2)*(1-Tr[i]**(1/2)))**2)
    A_i[0, i] = (0.42748*(Pr[i]/Tr[i]**2)*alfa_i[0, i])
    B_i[0, i] = (0.08664*(Pr[i]/Tr[i]))

A_ij = np.zeros((len(compound), len(compound)))
for i in range(len(compound)):
    for j in range(len(compound)):
        A_ij[i, j] = (1-0)*(A_i[0, i]*A_i[0, j])**(1/2)

A_mix = 0
B_mix = 0

for i in range(len(compound)):
    for j in range(len(compound)):
        A_mix += y_i[i]*y_i[j]*A_ij[i, j]
    B_mix = 1

Z_mix = 2

DH_RT = Z_mix - 1 - (3/2)*(A_mix/B_mix)*np.log(1+(B_mix/Z_mix))
DS_RT = np.log(Z_mix-B_mix) - (1/2)*(A_mix/B_mix)*np.log(1+(B_mix/Z_mix))
