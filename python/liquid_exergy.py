import exergy as ex
import numpy as np


def Cp_SVN_gas(coefs, R=8.31447):
    return lambda T: R * (
        coefs[0]
        + coefs[1] * T * 1e-3
        + coefs[2] * T ** 2 * 1e-6
        + coefs[3] * T ** (-2) * 1e5
    )


def Cp_SVN_liq(coefs, R=8.31447):
    return lambda T: R * (
        coefs[0] + coefs[1] * T * 1e-3 + coefs[2] * T ** 2 * 1e-6
    )


def Cp_constant(coef):
    return lambda T: coef * T


def read_csv_components(filename):
    pass


def main():
    H2O = ex.Component(
        name="Water",
        crit_temperature=647.1,
        crit_pressure=220.55,
        acentric_factor=0.345,
        enthalpy_of_formation=-285_830,
        entropy_of_formation=41,
        standard_exergy=0,
        heat_capacity=Cp_SVN_liq([8.712, 1.25, -0.18]),
    )
    H2O = ex.StreamComponent(H2O, 0.45)
    ethanol = ex.Component(
        name="Ethanol",
        crit_temperature=513.9,
        crit_pressure=61.48,
        acentric_factor=0.645,
        enthalpy_of_formation=-277_690,
        entropy_of_formation=159.9,
        standard_exergy=0,
        heat_capacity=Cp_SVN_liq([33.866, -172.60, 349.17]),
    )
    ethanol = ex.StreamComponent(ethanol, 0.015)
    benzene = ex.Component(
        name="Benzene",
        crit_temperature=562.2,
        crit_pressure=48.98,
        acentric_factor=0.210,
        enthalpy_of_formation=49_080,  # J / mol
        entropy_of_formation=173.26,  # J / mol K
        standard_exergy=0,
        heat_capacity=Cp_SVN_liq([-0.747, 67.96, -37.78]),
    )
    benzene = ex.StreamComponent(benzene, 0.02)
    AcOH = ex.Component(
        name="Acetic acid",
        crit_temperature=592.0,
        crit_pressure=57.86,
        acentric_factor=0.467,
        enthalpy_of_formation=-483_500,  # J / mol
        entropy_of_formation=158.0,  # J / mol K
        standard_exergy=0,
        heat_capacity=Cp_constant(123.1),  # J / mol K
    )
    AcOH = ex.StreamComponent(AcOH, 0.512)
    acetone = ex.Component(
        name="Acetone",
        crit_temperature=513.9,
        crit_pressure=61.48,
        acentric_factor=0.645,
        enthalpy_of_formation=-249_400,  # J / mol
        entropy_of_formation=200.4,  # J / mol K
        standard_exergy=0,
        heat_capacity=Cp_constant(125.5),  # J / mol K
    )
    acetone = ex.StreamComponent(acetone, 0.003)
    state = ex.ThermalProperties(temperature=310, pressure=1)
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
        print("Excess enthalpy: (J/mol)")
        print(
            ex.SRK.get_excess_enthalpy(
                state=stream.state,
                Tc=stream.get_Trs(),
                Pc=stream.get_Prs(),
                x=stream.get_molar_fractions(),
                w=stream.get_ws(),
                kij=kij,
            )
        )
        print(
            f"Enthalpy: {stream.enthalpy_with_eos(EOS=ex.SRK, kij=kij):,.4f} J/mol"
        )
        print(f"Entropy: {stream.entropy_liquid():,.4f} J/mol K")
        print(
            f"Exergy: {ex.liq_exergy_with_eos(stream=stream,EOS=ex.SRK,kij=kij):,.4f} J/mol"
        )


if __name__ == "__main__":
    main()
