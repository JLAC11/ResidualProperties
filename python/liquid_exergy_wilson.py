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
    H2O = ex.StreamComponent(H2O, 0.47)

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
    acetone = ex.StreamComponent(acetone, 0.53)
    state = ex.ThermalProperties(temperature=310, pressure=1)
    stream = ex.Stream([acetone, H2O], state=state)
    params = [
        1840.685,
        5884.506,
    ]  # J / mol, lambda 12-lambda 11, lambda 21- lambda 22
    act_model = ex.WilsonActivityModel(stream=stream, params=params)
    ln_g = act_model.calculate_ln_gamma()
    print(np.exp(ln_g))
    print(f"Excess enthalpy: {act_model.excess_enthalpy():,.4f} J/mol")
    print(f"Excess entropy: {act_model.excess_entropy():,.4f} J/mol K")
    print(f"Ideal enthalpy: {stream.ideal_enthalpy():,.4f} J/mol")
    print(f"Ideal entropy: {act_model.entropy():,.4f} J/mol K")
    print(f"Enthalpy: {act_model.enthalpy():,.4f} J/mol")
    print(f"Entropy: {act_model.entropy():,.4f} J/mol K")
    print(f"Exergy: {act_model.exergy():,.4f} J/mol")
    print(f"Temperatura de referencia: {ex.ReferenceState.temperature}")


if __name__ == "__main__":
    main()
