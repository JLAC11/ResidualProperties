import exergy as ex
import numpy as np


def Cp_SVN_gas(coefs, R=8.31446):
    return lambda T: R * (
        coefs[0]
        + coefs[1] * T * 1e-3
        + coefs[2] * T ** 2 * 1e-6
        + coefs[3] * T ** (-2) * 1e5
    )


def Cp_SVN_liq(coefs, R=8.31446):
    return lambda T: R * (
        coefs[0] + coefs[1] * T * 1e-3 + coefs[2] * T ** 2 * 1e-6
    )


def Cp_Shomate(coefs):  # En J / mol K
    return (
        lambda T: coefs[0]
        + coefs[1] * (T / 1000)
        + coefs[2] * (T / 1000) ** 2
        + coefs[3] * (T / 1000) ** 3
        + coefs[4] / (T / 1000) ** 2
    )


def Cp_constant(coef):
    return lambda T: coef


def read_csv_components(filename):
    pass


def main():
    H2O = ex.Component(  # NIST excepto donde se indique
        name="Water",
        crit_temperature=647.096,
        crit_pressure=220.64,
        acentric_factor=0.3443,
        enthalpy_of_formation=-285_830,
        entropy_of_formation=69.95,
        standard_exergy=0,
        heat_capacity=Cp_Shomate(
            [-203.6060, 1523.290, -3196.413, 2474.455, 3.855326]
        ),
    )
    H2O = ex.StreamComponent(H2O, 0.45)

    AcOH = ex.Component(  # NIST excepto donde se indique
        name="Acetic acid",
        crit_temperature=593.0,
        crit_pressure=57.81,
        acentric_factor=0.467,
        enthalpy_of_formation=-483_520,  # J / mol
        entropy_of_formation=158.0,  # J / mol K
        standard_exergy=0,
        heat_capacity=Cp_constant(123.1),  # J / mol K
    )  # Critical density 5.84 mol/L
    AcOH = ex.StreamComponent(AcOH, 0.015)

    ethanol = ex.Component(  # NIST excepto donde se indique
        name="Ethanol",
        crit_temperature=513.9,
        crit_pressure=61.40,
        acentric_factor=0.645,
        enthalpy_of_formation=-277_000,
        entropy_of_formation=159.86,
        standard_exergy=0,
        heat_capacity=Cp_SVN_liq([33.866, -172.60, 349.17]),  # Smith Van Ness
    )
    ethanol = ex.StreamComponent(ethanol, 0.02)

    acetone = ex.Component(  # NIST excepto donde se indique
        name="Acetone",
        crit_temperature=508.2,
        crit_pressure=47.00,
        acentric_factor=0.304,  # Wikipedia
        enthalpy_of_formation=-249_400,  # J / mol
        entropy_of_formation=200.4,  # J / mol K
        standard_exergy=0,
        heat_capacity=Cp_constant(125.45),  # J / mol K
    )
    acetone = ex.StreamComponent(acetone, 0.512)

    benzene = ex.Component(  # NIST excepto donde se indique
        name="Benzene",
        crit_temperature=562.02,
        crit_pressure=49.07277,
        acentric_factor=0.211,
        enthalpy_of_formation=48_950,  # J / mol
        entropy_of_formation=173.26,  # J / mol K
        standard_exergy=0,
        heat_capacity=Cp_SVN_liq([-0.747, 67.96, -37.78]),
    )  # Critical density: 3.901 mol/L
    benzene = ex.StreamComponent(benzene, 0.003)
    state = ex.ThermalProperties(temperature=310, pressure=1)  # 310 K, 1 bar
    components = [H2O, AcOH, ethanol, acetone, benzene]
    stream = ex.Stream(components, state)
    kij = ex.calculate_interaction_parameters_prausnitz(
        [57.1, 171, 167.1, 209, 259]
    )

    with np.printoptions(precision=4, suppress=True):
        print("Estado:")
        print(state)
        print("Propiedades en exceso en la mezcla (SRK):")
        print(
            ex.SRK.get_residual_properties(
                Tr=stream.get_Trs(),
                Pr=stream.get_Prs(),
                x=stream.get_molar_fractions(),
                w=stream.get_ws(),
                kij=kij,
            )
        )
        print("Valores de 1 - kij:")
        print(1 - kij)
        excess_enthalpy = ex.SRK.get_excess_enthalpy(
            state=stream.state,
            Tc=stream.get_Tcs(),
            Pc=stream.get_Pcs(),
            x=stream.get_molar_fractions(),
            w=stream.get_ws(),
            kij=kij,
        )
        print(f"Excess enthalpy: {excess_enthalpy:,.4f} J/mol")
        print(f"Ideal enthalpy: {stream.ideal_enthalpy():,.4f} J/mol")
        print(
            f"Enthalpy: {stream.enthalpy_with_eos(EOS=ex.SRK, kij=kij,liq=True):,.4f} J/mol"
        )
        print(
            f"Entropy: {stream.entropy_with_eos(EOS=ex.SRK,kij=kij,liq=True):,.4f} J/mol K"
        )

        excess_entropy = ex.SRK.get_excess_entropy(
            state=stream.state,
            Tc=stream.get_Tcs(),
            Pc=stream.get_Pcs(),
            x=stream.get_molar_fractions(),
            w=stream.get_ws(),
            kij=kij,
        )
        print(f"Excess entropy: {excess_entropy:,.4f} (J/mol K) ")
        exergy = ex.liq_exergy_with_eos(stream=stream, EOS=ex.SRK, kij=kij)
        print(f"Exergy: {exergy:,.4f} J/mol")

        enthalpy_srk = (
            ex.SRK.get_residual_properties(
                Tr=stream.get_Trs(),
                Pr=stream.get_Prs(),
                x=stream.get_molar_fractions(),
                w=stream.get_ws(),
                kij=kij,
                liq=True,
            ).enthalpy
            * 8.31446
            * state.temperature
            + stream.ideal_enthalpy()
        )
        print(f"Enthalpy SRK: (Ideal + residual): {enthalpy_srk:,.4f} J/mol")

        print("Fugacity coefficients:")
        print(
            ex.SRK.get_fugacity_coeffs(
                state=stream.state,
                Tr=stream.get_Trs(),
                Pr=stream.get_Prs(),
                x=stream.get_molar_fractions(),
                w=stream.get_ws(),
                kij=kij,
                liq=True,
            )
        )


if __name__ == "__main__":
    main()
