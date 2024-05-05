"""Hyperstruct."""

from dataclasses import dataclass


@dataclass
class Material:
    """The HyperStruct material model.

    The material object contains the basic material allowables.
    It also contains thermomechanical curves for environmental correction factors.
    """

    E: float  # Young's Modulus, tension
    E_c: float  # Young's Modulus, compression
    F_tu: float  # Tensile ultimate strength
    F_ty: float  # Tensile yield strength
    F_cy: float  # Compressive yield strength
    F_su: float  # Shear ultimate strength
    F_bru: float  # Bearing ultimate strength
    F_bry: float  # Bearing yield strength
