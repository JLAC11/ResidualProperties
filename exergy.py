from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, List
import numpy as np
import scipy.integrate as integrate


def main():
    pass


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


@dataclass
class StreamComponent:
    component: Component
    molar_fraction: float


@dataclass
class Stream:
    components: List[StreamComponent]
    state: ThermalProperties


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
        DeltaH: List[float], T: float, Cp: Callable, y
    ) -> float:
        return y.T @ DeltaH + integrate.quad(Cp, IdealGas.ref.temperature, T)


class ResidualEquationOfState(ABC):
    @staticmethod
    @abstractmethod
    def get_Z_factor(
        Pr: np.ndarray[float],
        Tr: np.ndarray[float],
        x: np.ndarray[float],
        w: np.ndarray[float],
        kij: np.ndarray[float],
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
        Pr: np.ndarray[float],
        Tr: np.ndarray[float],
        x: np.ndarray[float],
        w: np.ndarray[float],
        kij: np.ndarray[float] = 0,
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


class SRK(ResidualEquationOfState):
    @staticmethod
    def get_Z_factor(
        Pr: np.ndarray[float],
        Tr: np.ndarray[float],
        x: np.ndarray[float],
        w: np.ndarray[float],
        kij: np.ndarray[float] = 0,
    ) -> np.number:
        # Individual coefficients
        alpha = (
            1 + (0.48 + 1.574 * w - 0.176 * w ** 2) * (1 - np.sqrt(Tr))
        ) ** 2
        A = 0.42748 * (Pr / Tr ** 2) * alpha
        # SRK: exponente de 2.
        B = 0.08664 * (Pr / Tr)

        # Applying Van der Waals' mixing rule
        Aii = np.sqrt(A.T @ A) * (1 - kij)

        A = x.T @ Aii @ x
        B = x.T @ B

        p = -1
        q = A - B - B ** 2
        r = -A * B

        zeta = np.roots([1, p, q, r])
        z = zeta[np.imag(zeta) == 0]
        z = z[np.real(z) > 0]

        return z, A, B

    @staticmethod
    def get_residual_properties(
        Pr: np.ndarray[float],
        Tr: np.ndarray[float],
        x: np.ndarray[float],
        w: np.ndarray[float],
        kij: np.ndarray[float] = 0,
    ) -> ResidualProperties:

        z, A, B = SRK.get_Z_factor(Pr, Tr, x, w, kij)
        uu = -3 * A / (2 * B) * np.log(1 + B / z)
        hh = z - 1 + uu
        ss = np.log(z - B) - A / (2 * B) * np.log(1 + B / z)
        theta = (A / B) * np.log(1 + B / z)
        phi = np.exp(z - 1 - np.log(z - B) - theta)  # coeficiente de fugacidad
        props = ResidualProperties(uu, hh, ss, phi)
        return props


class PR(ResidualEquationOfState):
    def get_Z_factor(
        Pr: np.ndarray[float],
        Tr: np.ndarray[float],
        x: np.ndarray[float],
        w: np.ndarray[float],
        kij: np.ndarray[float] = 0,
    ):
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
        Pr: np.ndarray[float],
        Tr: np.ndarray[float],
        x: np.ndarray[float],
        w: np.ndarray[float],
        kij: np.ndarray[float] = 0,
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


class Exergy:
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


if __name__ == "__main__":
    main()
