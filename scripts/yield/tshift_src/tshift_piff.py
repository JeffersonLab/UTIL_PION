#!/usr/bin/env python3

#
# Description:
# ================================================================
# Time-stamp: "2026-03-11 12:12:00 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
##################################################################################
# Created - 10/July/2021, Author - Muhammad Junaid, University of Regina, Canada
##################################################################################

"""Calculate the implied -t shift produced by a missing-mass offset.

Positional command-line arguments
---------------------------------
1. Q2              : positive Q^2 in GeV^2
2. W               : invariant mass W in GeV
3. setting         : right, center, or left
4. beam_energy     : incident-electron TOTAL energy in GeV
5. MMoffset        : MM_data - MM_correct in GeV

Setting convention
------------------
right  -> theta_cm = 10 deg, phi_cm =   0 deg
center -> theta_cm =  0 deg, phi_cm =   0 deg
left   -> theta_cm = 10 deg, phi_cm = 180 deg

Output convention
-----------------
Tshift = (-t)_data - (-t)_correct, in GeV^2.
A positive MMoffset means that the data missing-mass peak is high.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from typing import Final


# Masses used by the original Fortran program, in MeV/c^2.
ELECTRON_MASS_MEV: Final[float] = 0.511
PION_MASS_MEV: Final[float] = 139.57
PROTON_MASS_MEV: Final[float] = 938.27
NEUTRON_MASS_MEV: Final[float] = 939.57

GEV_TO_MEV: Final[float] = 1000.0
GEV2_TO_MEV2: Final[float] = GEV_TO_MEV**2


@dataclass(frozen=True)
class SettingAngles:
    theta_cm_deg: float
    phi_cm_deg: float


SETTINGS: Final[dict[str, SettingAngles]] = {
    "right": SettingAngles(theta_cm_deg=10.0, phi_cm_deg=0.0),
    "center": SettingAngles(theta_cm_deg=0.0, phi_cm_deg=0.0),
    "left": SettingAngles(theta_cm_deg=10.0, phi_cm_deg=180.0),
}


@dataclass(frozen=True)
class TShiftResult:
    setting: str
    theta_cm_deg: float
    phi_cm_deg: float
    beam_energy_mev: float
    mm_offset_mev: float
    electron_scattering_angle_deg: float
    minus_t_correct_gev2: float
    minus_t_data_gev2: float
    tshift_gev2: float


def _finite(name: str, value: float) -> None:
    """Raise a useful error when an input is NaN or infinite."""
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number.")


def _two_body_pion_cm(
    w_mev: float,
    recoil_mass_mev: float,
) -> tuple[float, float]:
    """Return pion CM energy and momentum in MeV for gamma*p -> pi + X."""
    if recoil_mass_mev <= 0.0:
        raise ValueError("The shifted recoil (missing) mass must be positive.")

    if w_mev <= PION_MASS_MEV + recoil_mass_mev:
        raise ValueError(
            "W is below the pi+recoil two-body production threshold for the "
            "requested missing-mass offset."
        )

    pion_energy_cm = (
        w_mev**2 + PION_MASS_MEV**2 - recoil_mass_mev**2
    ) / (2.0 * w_mev)

    momentum_squared = pion_energy_cm**2 - PION_MASS_MEV**2
    # Permit only a tiny negative value from floating-point roundoff.
    scale = max(pion_energy_cm**2, PION_MASS_MEV**2, 1.0)
    if momentum_squared < -1.0e-12 * scale:
        raise ValueError("The requested inputs give an unphysical pion CM momentum.")

    pion_momentum_cm = math.sqrt(max(momentum_squared, 0.0))
    return pion_energy_cm, pion_momentum_cm


def _minus_t_mev2(
    q2_mev2: float,
    w_mev: float,
    theta_cm_rad: float,
    recoil_mass_mev: float,
) -> float:
    """Calculate -t in MeV^2 directly in the gamma*-proton CM frame."""
    pion_energy_cm, pion_momentum_cm = _two_body_pion_cm(
        w_mev=w_mev,
        recoil_mass_mev=recoil_mass_mev,
    )

    # Virtual-photon energy and three-momentum in the gamma*-p CM frame.
    photon_energy_cm = (
        w_mev**2 - PROTON_MASS_MEV**2 - q2_mev2
    ) / (2.0 * w_mev)
    photon_momentum_cm = math.sqrt(photon_energy_cm**2 + q2_mev2)

    # t = (q - p_pi)^2 with metric (+,-,-,-); this expression returns -t.
    return (
        q2_mev2
        - PION_MASS_MEV**2
        + 2.0 * photon_energy_cm * pion_energy_cm
        - 2.0
        * photon_momentum_cm
        * pion_momentum_cm
        * math.cos(theta_cm_rad)
    )


def calculate_tshift(
    q2_gev2: float,
    w_gev: float,
    setting: str,
    beam_energy_gev: float,
    mm_offset_gev: float,
) -> TShiftResult:
    
    """Calculate the missing-mass-implied Tshift.
    Parameters
    ----------
    q2_gev2
        Positive Q^2 in GeV^2.
    w_gev
        W in GeV.
    setting
        One of ``right``, ``center``, or ``left``.
    beam_energy_gev
        Incident-electron total energy in GeV. It is converted to MeV and used
        to verify that the supplied Q^2 and W are accessible at that energy.
    mm_offset_gev
        MM_data - MM_correct in GeV. Positive means the data MM peak is high.
    Returns
    -------
    TShiftResult
        ``tshift_gev2`` is (-t)_data - (-t)_correct in GeV^2.
    """

    for name, value in (
        ("Q2", q2_gev2),
        ("W", w_gev),
        ("beam_energy", beam_energy_gev),
        ("MMoffset", mm_offset_gev),
    ):
        _finite(name, value)

    setting_key = setting.strip().lower()
    if setting_key not in SETTINGS:
        valid = ", ".join(SETTINGS)
        raise ValueError(f"setting must be one of: {valid}.")

    if q2_gev2 < 0.0:
        raise ValueError("Q2 must be non-negative and is expressed in GeV^2.")
    if w_gev <= 0.0:
        raise ValueError("W must be positive.")
    if beam_energy_gev <= 0.0:
        raise ValueError("beam_energy must be positive.")

    angles = SETTINGS[setting_key]
    theta_cm_rad = math.radians(angles.theta_cm_deg)

    # Requested automatic unit conversions.
    q2_mev2 = q2_gev2 * GEV2_TO_MEV2
    w_mev = w_gev * GEV_TO_MEV
    beam_energy_mev = beam_energy_gev * GEV_TO_MEV
    mm_offset_mev = mm_offset_gev * GEV_TO_MEV

    # Electron-side kinematic validation. The beam input is total energy.
    energy_transfer_mev = (
        w_mev**2 + q2_mev2 - PROTON_MASS_MEV**2
    ) / (2.0 * PROTON_MASS_MEV)
    virtual_photon_momentum_mev = math.sqrt(
        energy_transfer_mev**2 + q2_mev2
    )

    scattered_energy_mev = beam_energy_mev - energy_transfer_mev
    if beam_energy_mev <= ELECTRON_MASS_MEV:
        raise ValueError("Beam total energy must exceed the electron rest mass.")
    if scattered_energy_mev <= ELECTRON_MASS_MEV:
        raise ValueError(
            "The scattered-electron energy is below its rest energy; the "
            "requested Q2 and W are not accessible at this beam energy."
        )

    incident_momentum_mev = math.sqrt(
        beam_energy_mev**2 - ELECTRON_MASS_MEV**2
    )
    scattered_momentum_mev = math.sqrt(
        scattered_energy_mev**2 - ELECTRON_MASS_MEV**2
    )

    denominator = 2.0 * incident_momentum_mev * scattered_momentum_mev
    cos_electron_angle = (
        incident_momentum_mev**2
        + scattered_momentum_mev**2
        - virtual_photon_momentum_mev**2
    ) / denominator

    angle_tolerance = 1.0e-10
    if cos_electron_angle < -1.0 - angle_tolerance or cos_electron_angle > 1.0 + angle_tolerance:
        raise ValueError(
            "The requested Q2 and W are not physically accessible at the "
            "specified beam energy."
        )
    cos_electron_angle = min(1.0, max(-1.0, cos_electron_angle))
    electron_scattering_angle_deg = math.degrees(math.acos(cos_electron_angle))

    recoil_mass_correct_mev = NEUTRON_MASS_MEV
    recoil_mass_data_mev = NEUTRON_MASS_MEV + mm_offset_mev

    minus_t_correct_mev2 = _minus_t_mev2(
        q2_mev2=q2_mev2,
        w_mev=w_mev,
        theta_cm_rad=theta_cm_rad,
        recoil_mass_mev=recoil_mass_correct_mev,
    )
    minus_t_data_mev2 = _minus_t_mev2(
        q2_mev2=q2_mev2,
        w_mev=w_mev,
        theta_cm_rad=theta_cm_rad,
        recoil_mass_mev=recoil_mass_data_mev,
    )

    minus_t_correct_gev2 = minus_t_correct_mev2 / GEV2_TO_MEV2
    minus_t_data_gev2 = minus_t_data_mev2 / GEV2_TO_MEV2
    tshift_gev2 = minus_t_data_gev2 - minus_t_correct_gev2

    return TShiftResult(
        setting=setting_key,
        theta_cm_deg=angles.theta_cm_deg,
        phi_cm_deg=angles.phi_cm_deg,
        beam_energy_mev=beam_energy_mev,
        mm_offset_mev=mm_offset_mev,
        electron_scattering_angle_deg=electron_scattering_angle_deg,
        minus_t_correct_gev2=minus_t_correct_gev2,
        minus_t_data_gev2=minus_t_data_gev2,
        tshift_gev2=tshift_gev2,
    )


def _setting_argument(value: str) -> str:
    setting = value.strip().lower()
    if setting not in SETTINGS:
        raise argparse.ArgumentTypeError(
            "setting must be right, center, or left"
        )
    return setting


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate Tshift = (-t)_data - (-t)_correct from a missing-mass "
            "offset. The numerical output is in GeV^2."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("Q2", type=float, help="positive Q^2 in GeV^2")
    parser.add_argument("W", type=float, help="W in GeV")
    parser.add_argument(
        "setting",
        type=_setting_argument,
        help="right, center, or left",
    )
    parser.add_argument(
        "beam_energy",
        type=float,
        help="incident-electron total beam energy in GeV",
    )
    parser.add_argument(
        "MMoffset",
        type=float,
        help="MM_data - MM_correct in GeV; positive means data is high",
    )
    parser.add_argument(
        "--details",
        action="store_true",
        help="print input conversions, angles, -t values, and sign convention",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        result = calculate_tshift(
            q2_gev2=args.Q2,
            w_gev=args.W,
            setting=args.setting,
            beam_energy_gev=args.beam_energy,
            mm_offset_gev=args.MMoffset,
        )
    except ValueError as exc:
        parser.error(str(exc))

    summary = (
        f"The Tshift for MMoffset = {args.MMoffset:.6f} GeV "
        f"is {result.tshift_gev2:+.6f} GeV^2"
    )

    if args.details:
        print(f"setting              = {result.setting}")
        print(f"theta_cm             = {result.theta_cm_deg:.6f} deg")
        print(f"phi_cm               = {result.phi_cm_deg:.6f} deg")
        print(f"beam energy          = {result.beam_energy_mev:.6f} MeV")
        print(f"MMoffset             = {result.mm_offset_mev:+.6f} MeV")
        print(
            "electron angle        = "
            f"{result.electron_scattering_angle_deg:.6f} deg"
        )
        print(
            f"(-t)_correct          = {result.minus_t_correct_gev2:.6f} GeV^2"
        )
        print(
            f"(-t)_data             = {result.minus_t_data_gev2:.6f} GeV^2"
        )
        print(
            "Tshift=data-correct   = "
            f"{result.tshift_gev2:+.6f} GeV^2"
        )
        print(
            "correction            : (-t)_correct = (-t)_data - Tshift"
        )
        print(summary)
    else:
        print(summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())