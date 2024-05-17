"""Hyperstruct."""

from dataclasses import dataclass
from importlib.metadata import version


__version__ = version("hyperstruct")


@dataclass
class Material:
    """The HyperStruct material model.

    The material object contains the basic material allowables.
    It also contains thermomechanical curves for environmental correction factors.
    """

    E: float
    """Young's Modulus, tension."""

    E_c: float
    """Young's Modulus, compression."""

    nu: float
    """Poisson's ratio."""

    F_tu: float
    """Tensile ultimate strength."""

    F_ty: float
    """Tensile yield strength."""

    F_cy: float
    """Compressive yield strength."""

    F_su: float
    """Shear ultimate strength."""

    F_bru: float
    """Bearing ultimate strength."""

    F_bry: float
    """Bearing yield strength."""

    F_en: float
    """Endurance limit."""

    db_r: float
    """Random spectrum decibel level that produces an acoustic fatigue life of 10^9."""


@dataclass
class Component:
    """The basic foundation for sizing.

    An Aircraft is composed of different Assemblies, which contain various Components.
    Components have several base methods and attributes that enable sizing and
    weights estimation.
    """

    material: Material
    """material the cover is made of."""

    def synthesis(self) -> None:
        """The sizing method.

        The sizing method collects all sizing routines and executes them
        in the order of the `routines` list.
        """
        # This doesn't work. It's just a placeholder.
        pass
