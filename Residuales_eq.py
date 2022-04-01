import numpy as np
import pandas as pd

# print(df.info())

compound = ["CO", "CO2", "CH4", "H2"]

y_i = [0.60175, 0.34575, 0.03375, 0.01875]

T1 = 424.85  # K
T2 = 273.15  # K
P1 = 5  # bar
P2 = 5  # bar

R = 8.31447  # kJ/kmol*K



def extract_data(comp):
    df = pd.read_excel("Heat_Capacity_Poly_Database.xlsx", sheet_name="Gral")
    df = df.fillna(0)
    W = np.zeros((1, len(comp)))
    Tc = np.zeros((1, len(comp)))
    Pc = np.zeros((1, len(comp)))
    poly = np.zeros((4, len(comp)))
    for i in range(len(comp)):
        x = df.loc[df["Chemical Formula"] == comp[i]]
        W[0, i] = float(x.loc[:, "w"])
        Tc[0, i] = float(x.loc[:, "Tc/K"])
        Pc[0, i] = float(x.loc[:, "Pc/bar"])
        poly[:, i] = np.array([float(x.loc[:, "A"]),
                              float(x.loc[:, "B E3"]),
                              float(x.loc[:, "C E6"]),
                              float(x.loc[:, "D E-5"])])
    return W, Tc, Pc, poly


def hs_residual(comp, y_i, T, P):
    W, Tc, Pc, poly = extract_data(comp)
    Tr = T/Tc
    Pr = P/Pc

    R = 8.31447  # J/mol*K

    A_i = np.zeros((1, len(comp)))
    B_i = np.zeros((1, len(comp)))
    alfa_i = np.zeros((1, len(comp)))

    for i in range(len(comp)):
        alfa_i[0, i] = ((1.0+(0.480+1.574*W[0, i]-0.176*W[0, i]**2)*(1-Tr[0, i]**(1/2)))**2)
        A_i[0, i] = (0.42748*(Pr[0, i]/Tr[0, i]**2)*alfa_i[0, i])
        B_i[0, i] = (0.08664*(Pr[0, i]/Tr[0, i]))

    A_ij = np.zeros((len(comp), len(comp)))
    for i in range(len(comp)):
        for j in range(len(comp)):
            A_ij[i, j] = (1-0)*(A_i[0, i]*A_i[0, j])**(1/2)

    A_mix = 0
    B_mix = 0

    for i in range(len(comp)):
        for j in range(len(comp)):
            A_mix += y_i[i]*y_i[j]*A_ij[i, j]
        B_mix += y_i[i]*B_i[0, i]

    r = -A_mix*B_mix
    q = A_mix-B_mix-B_mix**2
    p = -1
    Z_mix = np.polynomial.polynomial.Polynomial([r, q, p, 1]).roots()
    Z_mix = np.real(Z_mix @ np.isclose(Z_mix.imag, 0))
    DH_RT = Z_mix - 1 - (3/2)*(A_mix/B_mix)*np.log(1+(B_mix/Z_mix))
    DS_RT = np.log(Z_mix-B_mix) - (1/2)*(A_mix/B_mix)*np.log(1+(B_mix/Z_mix))
    return np.array([DH_RT*R*T, DS_RT*R])


def hs_ideal(comp, y_i, T1, T2):
    W, Tc, Pc, poly = extract_data(comp)
    print(poly)
    y = np.array(y_i).reshape(len(y_i), -1)
    print(y.shape)
    print((poly @ y).T)



    # return np.array[DH_RT * R * T, DS_RT * R]


def h_s(T1, T2, P1, P2, comp, y_i):
    hs_residual(comp, y_i, T1, P1)
    hs_ideal(comp, y_i, T1, T2)
    hs_residual(comp, y_i, T2, P2)


# H, S =\

h_s(T1, T2, P1, P2, compound, y_i)

