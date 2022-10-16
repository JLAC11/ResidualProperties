from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, List
import numpy as np
import scipy.integrate as integrate
from scipy.misc import derivative


def main():
    pass


def calculate_interaction_parameters_prausnitz(Vols) -> np.ndarray:
    """Calculates binary interaction parameters according to equation 4.9.6 in
    Prausnitz (insert book name here), on page 89.

    Args:
        Vols (np.array): Critical volumes of species.

    Returns:
        np.ndarray: Matrix of binary interaction parameters for given components.
    """
    n = np.size(Vols)
    kij = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            kij[i][j] = (
                1
                - 8
                * np.sqrt(Vols[i] * Vols[j])
                / (Vols[i] ** (1 / 3) + Vols[j] ** (1 / 3)) ** 3
            )
    kij = kij + kij.T - np.diag(kij)
    return kij


@dataclass
class ThermodynamicalProperties:
    internal_energy: float
    enthalpy: float
    entropy: float
    fugacity_coefficient: float


@dataclass
class ResidualProperties:
    """
    ResidualProperties: object containing residual properties in their adimensional
    form. Only contains internal energy, enthalpy, entropy and fugacity coefficient.
    """

    internal_energy: float
    enthalpy: float
    entropy: float
    fugacity_coefficient: float
    fugacity_coefficients: np.ndarray


@dataclass
class Component:
    name: str
    #    ID: str
    crit_temperature: float
    crit_pressure: float
    acentric_factor: float
    enthalpy_of_formation: float
    entropy_of_formation: float
    standard_exergy: float
    heat_capacity: Callable[[np.number], np.number]


@dataclass
class ThermalProperties:
    temperature: float = 298
    pressure: float = 1


class ReferenceState(ThermalProperties):
    """ReferenceState: state which is referred to. By default: 298 K and 1 bar."""


class RestrictedDeadState(ThermalProperties):
    """
    RestrictedDeadState: state at which system is in equilibrium with the
    environment. By default: 298 K and 1 bar.
    """


def enthalpy_ideal_pure(
    component: Component, T, T_ref=ReferenceState.temperature
):
    Cp = component.heat_capacity
    Href = component.enthalpy_of_formation
    H = Href + integrate.quad(Cp, T_ref, T)[0]
    return H


class ResidualEquationOfState(ABC):
    @staticmethod
    @abstractmethod
    def get_Z_factor(
        Pr: np.ndarray,
        Tr: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray,
        liq: bool = False,
    ):
        """Calculates and returns compressibility factor of equation of state

        Args:
            Pr (ndarray): Reduced pressure of each component
            Tr (ndarray): Reduced temperature of each component

        Returns:
            float: compressibility factor at given reduced temperature and pressure.
        """

    @staticmethod
    @abstractmethod
    def get_residual_properties(
        Pr: np.ndarray,
        Tr: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
        liq: bool = False,
    ) -> ResidualProperties:
        """Calculates and returns residual properties for components..
        Vector values must be supplied as a column vector.

        Args:
            Pr (ndarray): Reduced pressure of each component.
            Tr (ndarray): Reduced temperature of each component.
            x (ndarray): Molar fraction in gas phase of each component.
            w (np.ndarray[float]): Acentric factor of each component.
            kij (np.ndarray[float], optional): Correction factor. Defaults to 0.

        Returns:
            ResidualProperties object.
        """

        pass

    @staticmethod
    @abstractmethod
    def get_fugacity_mixture(Pr: np.ndarray, Tr: np.ndarray, y_i: np.ndarray):
        """Calculates and returns fugacity coefficients for gas mixture

        Args:
            Pr (ndarray): Reduced pressure of each component
            Tr (ndarray): Reduced temperature of each component
            y_i (ndarray): Molar fraction in gas phase of each component.

        Returns:
            float: fugacity coefficient of mixture at given reduced temperature
            and pressure
        """
        pass

    @staticmethod
    @abstractmethod
    def get_excess_enthalpy(
        state: ThermalProperties,
        Pc: np.ndarray,
        Tc: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
        R=8.31447,  # J / mol K
    ):
        """Gets excess enthalpy given equation of state (liquid phase)s

        Args:
            state (ThermalProperties): _description_
            Pr (np.ndarray): _description_
            Tr (np.ndarray): _description_
            x (np.ndarray): _description_
            w (np.ndarray): _description_
            kij (np.ndarray, optional): _description_. Defaults to 0.
            R (float, optional): _description_. Defaults to 8.31447.
        """
        pass

    @staticmethod
    @abstractmethod
    def get_excess_entropy(
        state: ThermalProperties,
        Pc: np.ndarray,
        Tc: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
        R=8.31447,  # J / mol K
    ):
        pass


class SRK(ResidualEquationOfState):
    @staticmethod
    def get_Z_factor(
        Pr: np.ndarray,
        Tr: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
        liq: bool = False,
    ):
        # Make them column vectors
        x = x.reshape(-1, 1)
        Pr = Pr.reshape(-1, 1)
        Tr = Tr.reshape(-1, 1)
        w = w.reshape(-1, 1)
        # Individual coefficients
        alpha = (
            1 + (0.48 + 1.574 * w - 0.176 * w ** 2) * (1 - np.sqrt(Tr))
        ) ** 2
        A = 0.42748 * (Pr / Tr ** 2) * alpha
        # SRK: exponente de 2.
        Bi = 0.08664 * (Pr / Tr)

        # Applying Van der Waals' mixing rule
        Aii = np.sqrt(A @ A.T) * (1 - kij)
        if Aii.size == 1:
            A = Aii * x.T @ x
        else:
            A = x.T @ Aii @ x

        B = x.T @ Bi

        p = -1
        q = (A - B - B ** 2)[0][0]  # Extract value from matrix
        r = (-A * B)[0][0]

        zeta = np.roots([1, p, q, r])
        z = zeta[np.imag(zeta) == 0]
        z = z[np.real(z) > 0]
        z = min(z) if liq is True else max(z)

        lnphi_i = (
            (Bi / B) * (z - 1)
            - np.log(z - B)
            - (A / B)
            * (2 * np.sum(x * Aii, axis=0).reshape(-1, 1) / A - Bi / B)
            * np.log(1 + B / z)
        )
        return z, A, B, lnphi_i

    @staticmethod
    def get_residual_properties(
        Pr: np.ndarray,
        Tr: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = np.array([0]),
        liq: bool = False,
    ) -> ResidualProperties:

        z, A, B, lnphis = SRK.get_Z_factor(Pr, Tr, x, w, kij, liq=liq)
        uu = -3 * A / (2 * B) * np.log(1 + B / z)
        hh = z - 1 + uu
        ss = np.log(z - B) - A / (2 * B) * np.log(1 + B / z)
        theta = (A / B) * np.log(1 + B / z)
        phi = np.exp(z - 1 - np.log(z - B) - theta)  # coeficiente de fugacidad
        props = ResidualProperties(
            uu[0][0], hh[0][0], ss[0][0], phi[0][0], np.exp(lnphis)
        )
        return props

    @staticmethod
    def get_fugacity_coeffs(
        state: ThermalProperties,
        Pr: np.ndarray,
        Tr: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
        liq: bool = False,
    ) -> ResidualProperties:
        z, A, B, lnphis = SRK.get_Z_factor(Pr, Tr, x, w, kij, liq=liq)
        phis = np.exp(lnphis)
        return phis

    @staticmethod
    def get_lnphis(
        state: ThermalProperties,
        Pr: np.ndarray,
        Tr: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
        liq: bool = False,
    ) -> ResidualProperties:
        z, A, B, lnphis = SRK.get_Z_factor(Pr, Tr, x, w, kij, liq=liq)
        return lnphis

    @staticmethod
    def get_excess_enthalpy(
        state: ThermalProperties,
        Pc: np.ndarray,
        Tc: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
        R=8.31446,  # J / mol K
    ):
        """H_ex = R*T^2* dln phi / dT

        Args:
            state (ThermalProperties): _description_
            Pc (np.ndarray): _description_
            Tc (np.ndarray): _description_
            x (np.ndarray): _description_
            w (np.ndarray): _description_
            kij (np.ndarray, optional): _description_. Defaults to 0.
            R (float, optional): _description_. Defaults to 8.31447.

        Returns:
            _type_: _description_
        """
        # Order as column vector
        x = x.reshape(-1, 1)

        der_ln_phi_dT = derivative(
            lambda T: SRK.get_lnphis(
                state=ThermalProperties(temperature=T, pressure=state.pressure),
                Pr=state.pressure / Pc,
                Tr=T / Tc,
                x=x,
                w=w,
                kij=kij,
                liq=True,
            ),
            x0=state.temperature,
            dx=0.05,
        )
        H_ex = -R * state.temperature ** 2 * np.sum(x * der_ln_phi_dT)
        return H_ex

    @staticmethod
    def get_excess_entropy(
        state: ThermalProperties,
        Pc: np.ndarray,
        Tc: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
        R=8.31447,  # J / mol K
    ):
        # Sanitize into column vectors
        x = x.reshape(-1, 1)
        Pc = Pc.reshape(-1, 1)
        Tc = Tc.reshape(-1, 1)
        w = w.reshape(-1, 1)
        der_ln_phi_dT = derivative(
            lambda T: SRK.get_lnphis(
                state=ThermalProperties(temperature=T, pressure=state.pressure),
                Pr=state.pressure / Pc,
                Tr=T / Tc,
                x=x,
                w=w,
                kij=kij,
                liq=True,
            ),
            x0=state.temperature,
            dx=0.1,
        )
        ln_phi = SRK.get_lnphis(
            state=state,
            Pr=state.pressure / Pc,
            Tr=state.temperature / Tc,
            x=x,
            w=w,
            kij=kij,
            liq=True,
        )  # En mezcla

        lnphi_i = np.log(
            SRK.get_residual_properties(
                Pr=state.pressure / Pc,
                Tr=state.temperature / Tc,
                x=x,
                w=w,
                kij=kij,
                liq=True,
            ).fugacity_coefficient
        )  # Especie pura

        S_ex = -R * (sum(x * (der_ln_phi_dT + ln_phi - lnphi_i)))
        return S_ex[0]


class PR(ResidualEquationOfState):
    @staticmethod
    def get_Z_factor(
        Pr: np.ndarray,
        Tr: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
    ) -> np.number:
        # Individual coefficients
        coef_acentric = 0.48 + 1.574 * w - 0.176 * w ** 2
        alpha = (1 + coef_acentric * (1 - np.sqrt(Tr))) ** 2
        A = 0.45724 * (Pr / Tr ** 2) * alpha
        # SRK: exponente de 2.
        B = 0.07780 * (Pr / Tr)

        # Applying Van der Waals' mixing rule
        Aii = np.sqrt(A.T @ A) * (1 - kij)

        A = x.T @ Aii @ x
        B = x.T @ B

        # Solving equations
        p = -1 + B
        q = A - 2 * B - 3 * B ** 2
        r = -A * B + B ** 2 + B ** 3

        gamma = coef_acentric * (np.sqrt(Tr / alpha))

        zeta = np.roots([1, p, q, r])
        z = zeta[np.imag(zeta) == 0]
        z = z[np.real(z) > 0]

        return z, A, B, gamma

    @staticmethod
    def get_residual_properties(
        Pr: np.ndarray,
        Tr: np.ndarray,
        x: np.ndarray,
        w: np.ndarray,
        kij: np.ndarray = 0,
    ) -> ResidualProperties:
        z, A, B, gamma = PR.get_Z_factor(Pr, Tr, x, w, kij)
        uu = (
            -A
            * (1 + gamma)
            / (B * np.sqrt(8))
            * np.log((z + B * (1 + np.sqrt(2))) / (z + B * (1 - np.sqrt(2))))
        )
        hh = z - 1 + uu
        ss = np.log(z - B) - A * gamma / (B * np.sqrt(8)) * np.log(
            (z + B * (1 + np.sqrt(2))) / (z + B * (1 - np.sqrt(2)))
        )
        theta = (
            A
            / (B * np.sqrt(8))
            * np.log((z + B * (1 + np.sqrt(2))) / (z + B * (1 - np.sqrt(2))))
        )
        phi = np.exp(z - 1 - np.log(z - B) - theta)
        return ResidualProperties(uu, hh, ss, phi)


@dataclass
class StreamComponent:
    component: Component
    molar_fraction: float

    def get_Pc(self):
        return self.component.crit_pressure

    def get_Tc(self):
        return self.component.crit_temperature

    def get_w(self):
        return self.component.acentric_factor

    def get_enthalpy(self):
        return self.component.enthalpy_of_formation

    def get_entropy(self):
        return self.component.entropy_of_formation


@dataclass
class Stream:
    components: List[StreamComponent]
    state: ThermalProperties

    def get_Pcs(self):
        return np.array([comp.get_Pc() for comp in self.components])

    def get_Prs(self):
        return self.state.pressure / self.get_Pcs()

    def get_Tcs(self):
        return np.array([comp.get_Tc() for comp in self.components])

    def get_Trs(self):
        return self.state.temperature / self.get_Tcs()

    def get_ws(self):
        return np.array([comp.get_w() for comp in self.components])

    def get_molar_fractions(self):
        return np.array([comp.molar_fraction for comp in self.components])

    def get_enthalpies(self):
        return np.array([comp.get_enthalpy() for comp in self.components])

    def get_entropies(self):
        return np.array([comp.get_entropy() for comp in self.components])

    def ideal_enthalpy(self, T_ref=ReferenceState.temperature):
        Hs = np.array(
            [
                enthalpy_ideal_pure(
                    comp.component, self.state.temperature, T_ref
                )
                for comp in self.components
            ]
        )
        return np.sum(self.get_molar_fractions() * Hs)

    def ideal_entropy(self, R=8.31446):
        x = self.get_molar_fractions()
        S_f = self.get_entropies()
        return np.sum(x * (S_f - R * np.log(x)))

    def enthalpy_with_eos(
        self,
        EOS: ResidualEquationOfState,
        T_ref=ReferenceState.temperature,
        kij=0,
        R=8.31446,
        liq: bool = False,
    ):
        H_ideal = self.ideal_enthalpy(T_ref)
        H_residual = (
            EOS.get_residual_properties(
                Pr=self.get_Prs(),
                Tr=self.get_Trs(),
                x=self.get_molar_fractions(),
                w=self.get_ws(),
                kij=kij,
                liq=liq,
            ).enthalpy
            * R
            * self.state.temperature
        )
        if liq is True:
            H_excess = EOS.get_excess_enthalpy(
                state=self.state,
                Tc=self.get_Tcs(),
                Pc=self.get_Pcs(),
                x=self.get_molar_fractions(),
                w=self.get_ws(),
                kij=kij,
                R=R,
            )
        else:
            H_excess = 0

        return H_ideal + H_residual + H_excess

    def entropy_with_eos(
        self,
        EOS: ResidualEquationOfState,
        T_ref=ReferenceState.temperature,
        kij=0,
        R=8.31446,
        liq: bool = False,
    ):
        S_ideal = self.ideal_entropy(R=R)
        S_residual = (
            EOS.get_residual_properties(
                Pr=self.get_Prs(),
                Tr=self.get_Trs(),
                x=self.get_molar_fractions(),
                w=self.get_ws(),
                kij=kij,
                liq=liq,
            ).entropy
            * R
        )
        if liq is True:
            S_excess = EOS.get_excess_entropy(
                state=self.state,
                Tc=self.get_Tcs(),
                Pc=self.get_Pcs(),
                x=self.get_molar_fractions(),
                w=self.get_ws(),
                kij=kij,
                R=R,
            )
        else:
            S_excess = 0

        return S_ideal + S_residual + S_excess


class IdealGas:
    ref = ReferenceState()

    @staticmethod
    def calculate_entropyPT(
        DeltaS: List[float],
        T: float,
        P: float,
        Cp: Callable,
        y: np.ndarray = np.array([[1]]),
        R=8.31447,
    ) -> float:
        """Calculates pointwise entropy given temperature and pressure

        Args:
            DeltaS (float): Entropy of formation at reference temperature and pressure
            T (float): Temperature
            P (float): Pressure
            Cp (Callable[[float],float]): Heat capacity as function of temperature
            y (np.ndarray): molar fraction in mixture, given as a column vector. Default: 1
            R (float): Ideal gas constant, default: in J/mol K
        Returns:
            Ideal gas entropy at given temperature and pressure.
        """
        Cp_T = lambda x: Cp(x) / x
        s_mix = -R * y.T @ np.log(y)
        s_pressure = -R * np.log(P / IdealGas.ref.pressure)
        s_temp = integrate.quad(Cp_T, IdealGas.ref.temperature, T)
        entropy = y.T @ DeltaS + s_temp + s_mix + s_pressure
        return entropy

    @staticmethod
    def calculate_enthalpy(
        DeltaH: np.ndarray, T: float, Cp: Callable, y
    ) -> float:
        y = y.reshape(-1, 1)
        DeltaH = DeltaH.reshape(-1, 1)
        return y.T @ DeltaH + integrate.quad(Cp, IdealGas.ref.temperature, T)


class Exergy:  # Funciona para gases
    @staticmethod
    def calculate_exergy(
        stream: Stream, method: ResidualEquationOfState, R: float = 8.31447
    ):
        T = stream.state.temperature
        P = stream.state.pressure
        T0 = ReferenceState.temperature
        P0 = ReferenceState.pressure
        H, S, G = Exergy._get_pointwise_enthalpy_and_entropy(
            T, P, stream, method
        )
        H0, S0, G0 = Exergy._get_pointwise_enthalpy_and_entropy(
            T0, P0, stream, method
        )
        physical_exergy = H - H0 - T * (S - S0)

        return physical_exergy

    @staticmethod
    def _get_pointwise_enthalpy_and_entropy(
        T: float,
        P: float,
        stream: Stream,
        method: ResidualEquationOfState,
        R: float = 8.31447,
    ):
        # Get vectorized quantities as column vectors
        y = np.array([x.molar_fraction for x in stream.components]).reshape(
            -1, 1
        )
        Pc = [x.component.crit_pressure for x in stream.components]
        Pr = np.array(Pc).reshape(-1, 1) / P
        Tc = [x.component.crit_temperature for x in stream.components]
        Tr = np.array(Tc).reshape(-1, 1) / T
        w = np.array(
            [x.component.acentric_factor for x in stream.components]
        ).reshape(-1, 1)
        cp_list = [x.component.heat_capacity for x in stream.components]

        DS = [x.component.entropy_of_formation for x in stream.components]
        DS = np.array(DS).reshape(-1, 1)
        DH = [x.component.enthalpy_of_formation for x in stream.components]
        DH = np.array(DH).reshape(-1, 1)

        # Get residual properties
        res_properties = method.get_residual_properties(Pr, Tr, y, w)
        # Calculate ideal gas enthalpy and entropy
        # Get average heat capacity
        cp = lambda T: sum(np.array([Cp(T) for Cp in cp_list]) @ y)
        S = IdealGas.calculate_entropyPT(DS, T, P, cp, y=y, R=R)
        S += res_properties.entropy * R
        H = IdealGas.calculate_enthalpy(DH, T, cp, y=y)
        H += res_properties.enthalpy * R * T
        G = H - T * S
        return H, S, G


def liq_exergy_with_eos(
    stream: Stream,
    EOS: ResidualEquationOfState,
    kij=0,
    R=8.31446,
    s_ref: Stream = None,
):
    if s_ref is None:
        s_ref = Stream(stream.components, ReferenceState())

    H = stream.enthalpy_with_eos(
        EOS=EOS, T_ref=s_ref.state.temperature, kij=kij, R=R, liq=True
    )
    H0 = s_ref.enthalpy_with_eos(
        EOS=EOS, T_ref=s_ref.state.temperature, kij=kij, R=R, liq=True
    )
    S = stream.entropy_with_eos(
        EOS=EOS, T_ref=s_ref.state.temperature, kij=kij, R=R, liq=True
    )
    S0 = s_ref.entropy_with_eos(
        EOS=EOS, T_ref=s_ref.state.temperature, kij=kij, R=R, liq=True
    )
    T = stream.state.temperature
    return H - H0 - T * (S - S0)


if __name__ == "__main__":
    main()
