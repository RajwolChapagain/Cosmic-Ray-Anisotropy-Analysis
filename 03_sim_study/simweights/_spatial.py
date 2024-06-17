# SPDX-FileCopyrightText: © 2022 the SimWeights contributors
#
# SPDX-License-Identifier: BSD-2-Clause
from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray


class CylinderBase:
    """Abstract base class for cylinder pdf classes."""

    def __init__(
        self: CylinderBase,
        length: float,
        radius: float,
        cos_zen_min: float,
        cos_zen_max: float,
        colname: str | None = None,
    ) -> None:
        if cos_zen_min < -1 or cos_zen_max > 1:
            raise ValueError(
                self.__class__.__name__ + ": both cos_zen_min and cos_zen_max must be between -1 and +1",
            )
        if cos_zen_min >= cos_zen_max:
            raise ValueError(self.__class__.__name__ + ": cos_zen_min must be less than cos_zen_max")
        self.length = length
        self.radius = radius
        self.cos_zen_min = cos_zen_min
        self.cos_zen_max = cos_zen_max
        self._side = 2e4 * self.radius * self.length
        self._cap = 1e4 * np.pi * self.radius**2
        self.etendue = float(self._diff_etendue(self.cos_zen_max) - self._diff_etendue(self.cos_zen_min))
        self.columns = (colname,)

    def projected_area(self: CylinderBase, cos_zen: ArrayLike) -> NDArray[np.float64]:
        """Cross sectional area of a cylinder in cm^2.

        As seen from the angle described by cos_zen.
        """
        cosz = np.asarray(cos_zen, dtype=np.float64)
        assert np.all(cosz >= -1)
        assert np.all(cosz <= +1)
        return np.asarray(self._cap * np.abs(cosz) + self._side * np.sqrt(1 - cosz**2), dtype=np.float64)

    def _diff_etendue(self: CylinderBase, cos_zen: ArrayLike) -> NDArray[np.float64]:
        cosz = np.asarray(cos_zen, dtype=np.float64)
        assert np.all(cosz >= -1)
        assert np.all(cosz <= +1)
        return np.asarray(
            np.pi * (self._cap * cosz * np.abs(cosz) + self._side * (cosz * np.sqrt(1 - cosz**2) - np.arccos(cosz))),
            dtype=np.float64,
        )

    def pdf(self: CylinderBase, cos_zen: ArrayLike) -> NDArray[np.float64]:
        """The probability density function for the given zenith angle."""
        raise NotImplementedError

    def __repr__(self: CylinderBase) -> str:
        return f"{self.__class__.__name__}" f"({self.length}, {self.radius}, {self.cos_zen_min}, {self.cos_zen_max})"

    def __eq__(self: CylinderBase, other: object) -> bool:
        return (
            isinstance(other, type(self))
            and self.length == other.length
            and self.radius == other.radius
            and self.cos_zen_min == other.cos_zen_min
            and self.cos_zen_max == other.cos_zen_max
        )


class UniformSolidAngleCylinder(CylinderBase):
    r"""Events are generated uniformly on the surface of a sphere.

    A Cylinder where the the angular distribution was sampled as if it were uniform on
    the surface of a sphere. The area of the location surface is proportional to the
    cross section of the cylinder perpendicular to the direction of the primary.

    The Monte Carlo must have been generated with the following zenith angle intensity:

    .. math::

        I \propto \cos\theta

    """

    def _pdf(self: UniformSolidAngleCylinder, cos_zen: NDArray[np.float64]) -> NDArray[np.float64]:
        return 1 / (2 * np.pi * (self.cos_zen_max - self.cos_zen_min) * self.projected_area(cos_zen))

    def pdf(self: UniformSolidAngleCylinder, cos_zen: ArrayLike) -> NDArray[np.float64]:
        cosz = np.asarray(cos_zen, dtype=np.float64)
        return np.piecewise(cosz, [(cosz >= self.cos_zen_min) & (cosz <= self.cos_zen_max)], [self._pdf])


class NaturalRateCylinder(CylinderBase):
    r"""Angular distribution when zenith distribution matched the natural rate of an isotropic source.

    For a given zenith angle the intensity of particles thrown was proportional to the cross-sectional area
    perpendicular to the direction of the particle. This is the distribution generated by the icetray class
    ``I3Surfaces::Cylinder`` and is what is used for triggered CORSIKA in ``I3PrimaryInjector``.
    It is also what CORSIKA will generate with the ``VOLUMECORR`` option, when the keyword ``DETCFG`` is set
    to :math:`l/(2r)`.

    The Monte Carlo must have been generated with the following zenith angle intensity:

    .. math::

        I \propto \pi\cdot r^2\cdot\sin\theta\cdot(\cos\theta+2/\pi\cdot l/r\cdot\sin\theta)
    """

    def __init__(
        self: NaturalRateCylinder,
        length: float,
        radius: float,
        cos_zen_min: float,
        cos_zen_max: float,
        colname: str | None = None,
    ) -> None:
        super().__init__(length, radius, cos_zen_min, cos_zen_max, colname)
        self._normalization = 1 / self.etendue

    def pdf(self: NaturalRateCylinder, cos_zen: ArrayLike) -> NDArray[np.float64]:
        cosz = np.asarray(cos_zen, dtype=np.float64)
        return np.piecewise(
            cosz,
            [(cosz >= self.cos_zen_min) & (cosz <= self.cos_zen_max)],
            [self._normalization],
        )


class CircleInjector:
    """Particles are generated on a circle perpendicular to the direction of the particle.

    Spatial distribution when particles are injected on the surfused by older neutrino-generator versions
    where the particle is injected in a cylinder that is parallel to momentum vector of the primary.
    The etendue is just the area of the circle times the solid angle.
    """

    def __init__(
        self: CircleInjector, radius: float, cos_zen_min: float, cos_zen_max: float, colname: str | None = None
    ) -> None:
        self.radius = radius
        self.cos_zen_min = cos_zen_min
        self.cos_zen_max = cos_zen_max
        self._cap = 1e4 * np.pi * self.radius**2
        self.etendue = 2 * np.pi * (self.cos_zen_max - self.cos_zen_min) * self._cap
        self._normalization = 1 / self.etendue
        self.columns = (colname,)

    def projected_area(self: CircleInjector, cos_zen: float) -> float:  # noqa: ARG002
        """Returns the cross sectional area of the injection area in cm^2."""
        # pylint: disable=unused-argument
        return self._cap

    def pdf(self: CircleInjector, cos_zen: ArrayLike) -> NDArray[np.float64]:
        """The probability density function for the given zenith angle."""
        cosz = np.asarray(cos_zen, dtype=np.float64)
        return np.piecewise(
            cosz,
            [(cosz >= self.cos_zen_min) & (cosz <= self.cos_zen_max)],
            [self._normalization],
        )

    def __repr__(self: CircleInjector) -> str:
        return f"CircleInjector({self.radius}, {self.cos_zen_min}, {self.cos_zen_max})"

    def __eq__(self: CircleInjector, other: object) -> bool:
        return (
            isinstance(other, CircleInjector)
            and self.radius == other.radius
            and self.cos_zen_min == other.cos_zen_min
            and self.cos_zen_max == other.cos_zen_max
        )


SpatialDist = Union[CylinderBase, CircleInjector]