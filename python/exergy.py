from abc import ABC, abstractmethod
import numpy as np


def main():
    pass


class EquationOfState(ABC):
    @staticmethod
    @abstractmethod
    def get_Z_factor(Pr: np.ndarray, Tr: np.ndarray) -> tuple[float]:
        """Calculates and returns compressibility factor of equation of state

        Args:
            Pr (ndarray): Reduced pressure of each component
            Tr (ndarray): Reduced temperature of each component

        Returns:
            float: compressibility factor at given reduced temperature and pressure.
        """

    @staticmethod
    @abstractmethod
    def get_fugacity_components(
        Pr: np.ndarray, Tr: np.ndarray, y_i: np.ndarray
    ) -> np.ndarray[float]:
        """Calculates and returns compressibility factor of equation of state

        Args:
            Pr (ndarray): Reduced pressure of each component
            Tr (ndarray): Reduced temperature of each component
            y_i (ndarray): Molar fraction in gas phase of each component.

        Returns:
            float: fugacity coefficient of each component at given temperature
            and pressure
        """

        pass


class SRK(EquationOfState):
    @staticmethod
    def get_Z_factor(
        Pr: np.ndarray[float],
        Tr: np.ndarray[float],
        x: np.ndarray[float],
        w: np.ndarray[float],
        kij: np.ndarray[float] = 0,
    ) -> tuple[float, float, float]:
        # Individual coefficients
        alpha = (
            1 + (0.48 + 1.574 * w - 0.176 * w ** 2) * (1 - np.sqrt(Tr))
        ) ** 2
        A = 0.42748 * (Pr / Tr ** 2) * alpha
        # SRK: exponente de 2.
        B = 0.08664 * (Pr / Tr)

        # Applying Van der Waals' mixing rule
        Aii = np.sqrt(A.T @ A) * (1 - kij)
        # ! Verificar
        A = x.T @ Aii @ x
        B = x.T @ B

        p = -1
        q = A - B - B ** 2
        r = -A * B

        zeta = np.roots([1, p, q, r])
        z = zeta[np.imag(zeta) == 0]
        z = z[np.real(z) > 0]

        uu = -3 * A / (2 * B) * np.log(1 + B / z)
        hh = z - 1 + uu
        ss = np.log(z - B) - A / (2 * B) * np.log(1 + B / z)
        theta = (A / B) * np.log(1 + B / z)
        phi = np.exp(z - 1 - np.log(z - B) - theta)  # coeficiente de fugacidad
        return z, A, B

    @staticmethod
    def get_fugacity_components(
        Pr: np.ndarray, Tr: np.ndarray, y_i: np.ndarray
    ) -> np.ndarray[float]:
        pass


if __name__ == "__main__":
    main()
