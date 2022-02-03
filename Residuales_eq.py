import numpy as np
import pandas as pd

df = pd.read_excel("Heat_Capacity_Poly_Database.xlsx", sheet_name="Gral")
# print(df.info())
compound = ["CO", "CO2", "CH4", "H2"]
y_i = [0.60175, 0.34575, 0.03375, 0.01875]

R = 8.31447
T = 424.85
P = 5

def extract_data(comp):
    W = np.zeros((1, len(comp)))
    Tc = np.zeros((1, len(comp)))
    Pc = np.zeros((1, len(comp)))
    for i in range(len(comp)):
        x = df.loc[df["Chemical Formula"] == comp[i]]
        W[0, i] = float(x.loc[:, "w"])
        Tc[0, i] = float(x.loc[:, "Tc/K"])
        Pc[0, i] = float(x.loc[:, "Pc/bar"])
    return W, Tc, Pc


W, Tc, Pc = extract_data(compound)

Tr = T/Tc
Pr = P/Pc


A_i = np.zeros((1, len(compound)))
B_i = np.zeros((1, len(compound)))
alfa_i = np.zeros((1, len(compound)))

for i in range(len(compound)):
    alfa_i[0, i] = ((1.0+(0.480+1.574*W[0, i]-0.176*W[0, i]**2)*(1-Tr[0, i]**(1/2)))**2)
    A_i[0, i] = (0.42748*(Pr[0, i]/Tr[0, i]**2)*alfa_i[0, i])
    B_i[0, i] = (0.08664*(Pr[0, i]/Tr[0, i]))

A_ij = np.zeros((len(compound), len(compound)))
for i in range(len(compound)):
    for j in range(len(compound)):
        A_ij[i, j] = (1-0)*(A_i[0, i]*A_i[0, j])**(1/2)

A_mix = 0
B_mix = 0

for i in range(len(compound)):
    for j in range(len(compound)):
        A_mix += y_i[i]*y_i[j]*A_ij[i, j]
    B_mix += y_i[i]*B_i[0, i]

r = -A_mix*B_mix
q = A_mix-B_mix-B_mix**2
p = -1
Z_mix = np.polynomial.polynomial.Polynomial([r, q, p, 1]).roots()
Z_mix = np.real(Z_mix @ np.isclose(Z_mix.imag, 0))
DH_RT = Z_mix - 1 - (3/2)*(A_mix/B_mix)*np.log(1+(B_mix/Z_mix))
DS_RT = np.log(Z_mix-B_mix) - (1/2)*(A_mix/B_mix)*np.log(1+(B_mix/Z_mix))
print(DH_RT*R*T*1E3)
print(DS_RT*R*1E3)
