# -*- coding: utf-8 -*-
from typing import Callable, Self


class SkinFrictionFactor:
    """ Skin friction factors for y+ calculations.

    To-do: implement some from the following (after review);
    https://www.cfd-online.com/Wiki/Skin_friction_coefficient
    """
    @staticmethod
    def laminar(Re):
        """ Laminar limit theoretical value. """
        return 64 / Re

    @staticmethod
    def smooth_wall(Re, check: bool = True):
        """ Blasius smooth wall approximation.

        https://doi.org/10.1007/978-3-662-02239-9_1
        """
        if check and not (4_000 < Re < 100_000):
            print(f"WARNING: out-of-range Re = {Re:.2e} not in [4e3; 1e5]")

        return 0.3164 * Re**(-1/4)


class WallGradingCalculator:
    """ Helper class for estimating first cell thickness given y+. """
    def __init__(self, *,
            L: float,
            U: float,
            rho: float,
            mu: float,
            skin_factor: Callable | None = None
        ) -> None:
        self._rho = rho
        self._U   = U
        self._nu  = mu / rho
        self._Re  = U * L / self._nu

        if skin_factor is not None:
            self.set_skin_factor(skin_factor)

    @classmethod
    def from_solution(cls, obj, L: float, U: float,
                      skin_factor: Callable | None = None) -> Self:
        """ Alternative constructor from dimensionless solution. """
        sol = obj.solution
        return cls(L=L, U=U, rho=sol.density_mass, mu=sol.viscosity,
                   skin_factor=skin_factor)

    def set_skin_factor(self, skin_factor: Callable) -> None:
        self._Cf = skin_factor(self._Re)
        self._tw = self.wall_shear_stress(self._rho, self._U, self._Cf)
        self._ut = self.friction_velocity(self._tw, self._rho)
        # print(f"Skin coefficient .... {self._Cf}")
        # print(f"Wall shear stress ... {self._tw}")
        # print(f"Friction velocity ... {self._ut}")

    @staticmethod
    def wall_shear_stress(rho, U, Cf):
        """ Wall shear stress estimater from friction factor [Pa]. """
        return Cf * rho * U**2

    @staticmethod
    def friction_velocity(tw, rho):
        """ Dimensionless friction velocity. """
        return (tw / rho)**(1/2)

    def first_layer(self, y_plus, skin_factor: Callable | None = None) -> None:
        """ Height of first cell for given y+ value [m]. """
        if skin_factor is not None:
            self.set_skin_factor(skin_factor)

        if not hasattr(self, "_ut"):
            raise ValueError("Please call `set_skin_factor` first")

        return y_plus * self._nu / self._ut
