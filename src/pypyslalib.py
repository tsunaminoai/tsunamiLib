"""
PyPySLALib is licensed under GPLv3; 
PySLALib and SLALib are licensed under GPLv2;
The legal notices are below

---------------------------------------------------------------
PyPySLALib, a full python conversion of the astrometrics library SLALib
Copyright (C) 2022  Gregory Foote

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

---------------------------------------------------------------
f2py-generated wrappers for SLALIB
Copyright (C) 2010 Scott M. Ransom

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

---------------------------------------------------------------
SLALIB is a library of routines intended to make accurate and reliable positional-astronomy applications easier to write
Copyright (C) 1995 P.T.Wallace; Starlink; Rutherford Appleton Laboratory

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

---------------------------------------------------------------
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

import re
import time

import numpy as np

class SLALib:
    # 2 * PI
    D2PI = 6.2831853071795864769252867665590057683943387987502
    # PI
    DPI = 3.141592653589793238462643
    # PI/(12 * 3600) : Seconds to Radians
    DS2R = 7.2722052166430399038487115353692196393452995355905e-5
    # Arcseconds to Radians
    AS2R = 4.848136811095359935899141e-6
    # Epsilon for float comparision
    TINY = 1e-30
    # Light time for 1 AU (sec)
    CR = 499.004782e0
    # Gravitational radius of the Sun x 2 (2*mu/c**2, AU)
    GR2 = 1.974126e-8
    # B1950
    B1950 = 1949.9997904423e0
    # Degrees to radians
    DD2R = 1.745329251994329576923691e-2
    # Arc seconds in a full circle
    TURNAS = 1296000e0
    # Reference epoch (J2000), MJD
    DJM0 = 51544.5e0
    # Days per Julian century
    DJC = 36525e0
    # Mean sidereal rate (at J2000) in radians per (UT1) second
    SR = 7.292115855306589e-5
    # Earth equatorial radius (metres)
    A0 = 6378140e0
    # Reference spheroid flattening factor and useful function
    SPHF = 1e0 / 298.257e0
    SPHB = (1e0 - SPHF) ** 2
    # Astronomical unit in metres
    AU = 1.49597870e11

    @staticmethod
    def dmod(A, B):
        return A % B

    @classmethod
    def dranrm(cls, in_val):
        return cls.dmod(in_val, cls.D2PI)

    @classmethod
    def gmst(cls, in_ut1):
        # Julian centuries from fundamental epoch J2000 to this UT
        tu = (in_ut1 - 51544.5) / 36525.0
        return cls.dranrm(
            cls.dmod(in_ut1, 1.0) * cls.D2PI
            + (24110.54841 + (8640184.812866 + (0.093104 - 6.2e-6 * tu) * tu) * tu)
            * cls.DS2R
        )

    @classmethod
    def hour_angle(cls, in_mjd, in_ra, in_long):
        # Not part of the original library
        return np.rad2deg(cls.dranrm(cls.dranrm(cls.gmst(in_mjd) + in_long) - in_ra))

    @classmethod
    def clyd(cls, input_year, input_month, input_day):
        # +
        # - - - - -
        # C L Y D
        # - - - - -
        # Gregorian calendar to year and day in year (in a Julian calendar
        # aligned to the 20th/21st century Gregorian calendar).
        # Given:
        # IY,IM,ID   i    year, month, day in Gregorian calendar
        # Returned:
        # NY         i    year (re-aligned Julian calendar)
        # ND         i    day in year (1 = January 1st)
        # JSTAT      i    status:
        # 0 = OK
        # 1 = bad year (before -4711)
        # 2 = bad month
        # 3 = bad day (but conversion performed)
        # Notes:
        # 1  This routine exists to support the low-precision routines
        # sla_EARTH, sla_MOON and sla_ECOR.
        # 2  Between 1900 March 1 and 2100 February 28 it returns answers
        # which are consistent with the ordinary Gregorian calendar.
        # Outside this range there will be a discrepancy which increases
        # by one day for every non-leap century year.
        # 3  The essence of the algorithm is first to express the Gregorian
        # date as a Julian Day Number and then to convert this back to
        # a Julian calendar date, with day-in-year instead of month and
        # day.  See 12.92-1 and 12.95-1 in the reference.
        # Reference:  Explanatory Supplement to the Astronomical Almanac,
        # ed P.K.Seidelmann, University Science Books (1992),
        # p604-606.

        return_code = 0
        return_year = 0
        return_day = 0

        # Validate year
        if input_year >= -4711:
            #  Validate month
            if 1 <= input_month <= 12:
                month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

                # Allow for (Gregorian) leap year
                if np.mod(input_year, 4) == 0 and (
                    np.mod(input_year, 100) != 0 or np.mod(input_year, 400) == 0
                ):
                    month_lengths[1] = 29

                # Validate day
                if input_day < 1 or input_day > month_lengths[int(input_month) - 1]:
                    return_code = 3
                # Perform the conversion
                temp_i = (14 - input_month) / 12
                temp_k = input_year - temp_i
                temp_j = (
                    (1461 * (temp_k + 4800)) / 4
                    + (367 * (input_month - 2 + 12 * temp_i)) / 12
                    - (3 * ((temp_k + 4900) / 100)) / 4
                    + input_day
                    - 30660
                )
                temp_k = (temp_j - 1) / 1461
                temp_l = temp_j - 1461 * temp_k
                temp_n = (temp_l - 1) / 365 - temp_l / 1461
                temp_j = ((80 * (temp_l - 365 * temp_n + 30)) / 2447) / 11
                temp_i = temp_n + temp_j
                return_day = (
                    59 + temp_l - 365 * temp_i + ((4 - temp_n) / 4) * (1 - temp_j)
                )
                return_year = 4 * temp_k + temp_i - 4716

            # Bad month
            else:
                return_code = 2
        else:
            # Bad year
            return_code = 1

        return return_year, return_day, return_code

    @classmethod
    def dcc2s(cls, in_coord):
        # +
        #      - - - - - -
        #       D C C 2 S
        #      - - - - - -
        #   Cartesian to spherical coordinates
        #   Given:
        #      V     d(3)   x,y,z vector
        #   Returned:
        #      A,B   d      spherical coordinates in radians
        #   The spherical coordinates are longitude (+ve anticlockwise looking
        #   from the +ve latitude pole) and latitude.  The Cartesian coordinates
        #   are right handed, with the x axis at zero longitude and latitude, and
        #   the z axis at the +ve latitude pole.
        #   If V is null, zero A and B are returned.  At either pole, zero A is
        #   returned.
        #   Last revision:   22 July 2004

        x, y, z = in_coord
        r = np.sqrt(x * x + y * y)

        longitude = np.where(r == 0, 0, np.arctan2(y, x))
        latitude = np.where(z == 0, 0, np.arctan2(z, r))

        return longitude, latitude

    @classmethod
    def dcs2c(cls, ra, dec):
        # +
        #      - - - - - -
        #       D C S 2 C
        #      - - - - - -
        #
        #   Spherical coordinates to direction cosines (double precision)
        #
        #   Given:
        #      A,B       d      spherical coordinates in radians
        #                          (RA,Dec), (long,lat) etc.
        #
        #   Returned:
        #      V         d(3)   x,y,z unit vector
        #
        #   The spherical coordinates are longitude (+ve anticlockwise looking
        #   from the +ve latitude pole) and latitude.  The Cartesian coordinates
        #   are right handed, with the x axis at zero longitude and latitude, and
        #   the z axis at the +ve latitude pole.
        #
        #   Last revision:   26 December 2004

        return np.array(
            [np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)]
        )

    @classmethod
    def deuler(cls, order, phi, theta, psi):
        #
        #      - - - - - - -
        #       D E U L E R
        #      - - - - - - -

        #   Form a rotation matrix from the Euler angles - three successive
        #   rotations about specified Cartesian axes (double precision)

        #   Given:
        #     ORDER   c*(*)   specifies about which axes the rotations occur
        #     PHI     d       1st rotation (radians)
        #     THETA   d       2nd rotation (   "   )
        #     PSI     d       3rd rotation (   "   )

        #   Returned:
        #     RMAT    d(3,3)  rotation matrix

        #   A rotation is positive when the reference frame rotates
        #   anticlockwise as seen looking towards the origin from the
        #   positive region of the specified axis.

        #   The characters of ORDER define which axes the three successive
        #   rotations are about.  A typical value is 'ZXZ', indicating that
        #   RMAT is to become the direction cosine matrix corresponding to
        #   rotations of the reference frame through PHI radians about the
        #   old Z-axis, followed by THETA radians about the resulting X-axis,
        #   then PSI radians about the resulting Z-axis.

        #   The axis names can be any of the following, in any order or
        #   combination:  X, Y, Z, uppercase or lowercase, 1, 2, 3.  Normal
        #   axis labelling/numbering conventions apply;  the xyz (=123)
        #   triad is right-handed.  Thus, the 'ZXZ' example given above
        #   could be written 'zxz' or '313' (or even 'ZxZ' or '3xZ').  ORDER
        #   is terminated by length or by the first unrecognized character.

        #   Fewer than three rotations are acceptable, in which case the later
        #   angle arguments are ignored.  If all rotations are zero, the
        #   identity matrix is produced.

        #   Initialize result matrix
        result_mat = np.identity(3)
        #   Look at each character of axis string until finished
        for i_char, axis in enumerate(order):
            # Initialize rotation matrix for the current rotation
            rot_mat = np.identity(3)
            # Pick up the appropriate Euler angle and take sine & cosine
            if i_char == 1:
                angle = phi
            elif i_char == 2:
                angle = theta
            else:
                angle = psi

            ang_sin = np.sin(angle)
            ang_cos = np.cos(angle)

            # Identify the axis
            if axis in ["X", "x", "1"]:
                # Matrix for x-rotation
                rot_mat[1, 1] = ang_cos
                rot_mat[1, 2] = ang_sin
                rot_mat[2, 1] = -ang_sin
                rot_mat[2, 2] = ang_cos

            elif axis in ["Y", "y", "2"]:
                # Matrix for y-rotation
                rot_mat[0, 0] = ang_cos
                rot_mat[0, 2] = -ang_sin
                rot_mat[2, 0] = ang_sin
                rot_mat[2, 2] = ang_cos

            elif axis in ["Z", "z", "3"]:
                # Matrix for z-rotation
                rot_mat[0, 0] = ang_cos
                rot_mat[0, 1] = ang_sin
                rot_mat[1, 0] = -ang_sin
                rot_mat[1, 1] = ang_cos
            else:
                raise ValueError("Invalid Character received for deuler!")

            result_mat = result_mat @ rot_mat
        return result_mat

    @classmethod
    def epb2d(cls, epb):
        # +
        # - - - - - -
        # E P B 2 D
        #
        # Conversion of Besselian Epoch to Modified Julian Date
        # (double precision)
        # Given:
        # epb      dp       Besselian Epoch
        # The result is the Modified Julian Date (JD - 2400000.5).
        # -

        return 15019.81352e0 + (epb - 1900e0) * 365.242198781e0

    @classmethod
    def djcl(cls, in_mjd):
        # +
        # - - - - -
        # D J C L
        # - - - - -
        # Modified Julian Date to Gregorian year, month, day,
        # and fraction of a day.
        # Given:
        # DJM      dp     modified Julian Date (JD-2400000.5)
        # Returned:
        # IY       int    year
        # IM       int    month
        # ID       int    day
        # FD       dp     fraction of day
        # J        int    status:
        # 0 = OK
        # -1 = unacceptable date (before 4701BC March 1)
        # The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
        # Last revision:   22 July 2004

        # Check if date is acceptable.

        r_year = r_month = r_day = r_frac_day = 0

        if in_mjd <= -2395520e0 or in_mjd >= 1e9:
            return_code = -1
        else:
            # Separate day and fraction.
            frac_day = np.mod(in_mjd, 1e0)
            if frac_day < 0e0:
                frac_day += 1e0
            day_int = np.rint(in_mjd - frac_day)

            # Express day in Gregorian calendar.
            jd = np.rint(day_int) + 2400001

            n4 = 4 * (jd + ((6 * ((4 * jd - 17918) / 146097)) / 4 + 1) / 2 - 37)
            nd10 = 10 * (np.mod(n4 - 237, 1461) / 4) + 5

            r_year = n4 / 1461 - 4712
            r_month = np.mod(nd10 / 306 + 2, 12) + 1
            r_day = np.mod(nd10, 306) / 10 + 1
            r_frac_day = frac_day

            return_code = 0
        return r_year, r_month, r_day, r_frac_day, return_code

    @classmethod
    def djcal(cls, ndp, djm):
        """
        - - - - - -
        D J C A L
        - - - - - -

        Modified Julian Date to Gregorian Calendar, expressed
        in a form convenient for formatting messages (namely
        rounded to a specified precision, and with the fields
        stored in a single array)

        Given:
        NDP      i      number of decimal places of days in fraction
        DJM      d      modified Julian Date (JD-2400000.5)

        Returned:
        IYMDF    i(4)   year, month, day, fraction in Gregorian
        calendar
        J        i      status:  nonzero = out of range

        Any date after 4701BC March 1 is accepted.

        NDP should be 4 or less if internal overflows are to be avoided
        on machines which use 32-bit integers.

        The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).

        Last revision:   22 July 2004
        """
        # Validate.
        iymdf = np.array([])
        if (djm <= -2395520e0) or (djm <= -1e9):
            j = -1
        else:
            j = 0

            # Denominator of fraction.
            NFD = 10 ** np.maximum(ndp, 0)
            FD = NFD

            # Round date and express in units of fraction.
            DF = np.rint(djm * FD)

            # Separate day and fraction.
            F = np.mod(DF, FD)
            if F < 0e0:
                F = F + FD
            D = (DF - F) / FD

            # Express day in Gregorian calendar.
            JD = np.rint(D) + 2400001

            N4 = 4 * (JD + ((2 * ((4 * JD - 17918) / 146097) * 3) / 4 + 1) / 2 - 37)
            ND10 = 10 * (np.mod(N4 - 237, 1461) / 4) + 5

            iymdf = np.array(
                [
                    N4 / 1461 - 4712,
                    np.mod(ND10 / 306 + 2, 12) + 1,
                    np.mod(ND10, 306) / 10 + 1,
                    np.rint(F),
                ]
            )
        return j, iymdf

    @classmethod
    def dmxv(cls, in_mat, in_vec):
        # +
        # - - - - -
        # D M X V
        # - - - - -
        #
        # Performs the 3-D forward unitary transformation:
        #
        # vector VB = matrix DM * vector VA
        #
        # (double precision)
        #
        # Given:
        # DM       dp(3,3)    matrix
        # VA       dp(3)      vector
        #
        # Returned:
        # VB       dp(3)      result vector
        #
        # To comply with the ANSI Fortran 77 standard, VA and VB must be
        # different arrays.  However, the routine is coded so as to work
        # properly on many platforms even if this rule is violated.
        #
        # Last revision:   26 December 2004
        # -

        return np.dot(in_mat, in_vec)

    @classmethod
    def dimxv(cls, in_mat, in_vec):
        """
        - - - - - -
        D I M X V
        - - - - - -

        Performs the 3-D backward unitary transformation:

        vector VB = (inverse of matrix DM) * vector VA

        (double precision)

        (n.b.  the matrix must be unitary, as this routine assumes that
        the inverse and transpose are identical)

        Given:
        DM       dp(3,3)    matrix
        VA       dp(3)      vector

        Returned:
        VB       dp(3)      result vector
        Inverse of matrix DM * vector VA -> vector VW
        """
        return np.dot(in_mat.T, in_vec)

    @classmethod
    def prebn(cls, bess_epoch_0, bess_epoch_1):
        # +
        # - - - - - -
        # P R E B N
        # - - - - - -
        #
        # Generate the matrix of precession between two epochs,
        # using the old, pre-IAU1976, Bessel-Newcomb model, using
        # Kinoshita's formulation (double precision)
        #
        # Given:
        # BEP0    dp         beginning Besselian epoch
        # BEP1    dp         ending Besselian epoch
        #
        # Returned:
        # RMATP  dp(3,3)    precession matrix
        #
        # The matrix is in the sense   V(BEP1)  =  RMATP * V(BEP0)
        #
        # Reference:
        # Kinoshita, H. (1975) 'Formulas for precession', SAO Special
        # Report No. 364, Smithsonian Institution Astrophysical
        # Observatory, Cambridge, Massachusetts.
        #
        # Called:  sla_DEULER
        # -

        # Interval between basic epoch B1850.0 and beginning epoch in TC
        epoch_diff = (bess_epoch_0 - 1850e0) / 100e0

        # Interval over which precession required, in tropical centuries
        interval_precision = (bess_epoch_1 - bess_epoch_0) / 100e0

        # Euler angles
        tas2_r = interval_precision * cls.AS2R
        w = 2303.5548e0 + (1.39720e0 + 0.000059e0 * epoch_diff) * epoch_diff

        zeta = (
            w
            + (0.30242e0 - 0.000269e0 * epoch_diff + 0.017996e0 * interval_precision)
            * interval_precision
        ) * tas2_r
        z = (
            w
            + (1.09478e0 + 0.000387e0 * epoch_diff + 0.018324e0 * interval_precision)
            * interval_precision
        ) * tas2_r
        theta = (
            2005.1125e0
            + (-0.85294e0 - 0.000365e0 * epoch_diff) * epoch_diff
            + (-0.42647e0 - 0.000365e0 * epoch_diff - 0.041802e0 * interval_precision)
            * interval_precision
        ) * tas2_r

        # Rotation matrix
        return cls.deuler("ZYZ", -zeta, theta, -z)

    @classmethod
    def prec(cls, start_epoch, end_epoch):
        # +
        # - - - - -
        # P R E C
        # - - - - -
        #
        # Form the matrix of precession between two epochs (IAU 1976, FK5)
        # (double precision)
        #
        # Given:
        # EP0    dp         beginning epoch
        # EP1    dp         ending epoch
        #
        # Returned:
        # RMATP  dp(3,3)    precession matrix
        #
        # Notes:
        #
        # 1)  The epochs are TDB (loosely ET) Julian epochs.
        #
        # 2)  The matrix is in the sense   V(EP1)  =  RMATP * V(EP0)
        #
        # 3)  Though the matrix method itself is rigorous, the precession
        # angles are expressed through canonical polynomials which are
        # valid only for a limited time span.  There are also known
        # errors in the IAU precession rate.  The absolute accuracy
        # of the present formulation is better than 0.1 arcsec from
        # 1960AD to 2040AD, better than 1 arcsec from 1640AD to 2360AD,
        # and remains below 3 arcsec for the whole of the period
        # 500BC to 3000AD.  The errors exceed 10 arcsec outside the
        # range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
        # 5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.
        # The SLALIB routine sla_PRECL implements a more elaborate
        # model which is suitable for problems spanning several
        # thousand years.
        #
        # References:
        # Lieske,J.H., 1979. Astron.Astrophys.,73,282.
        # equations (6) & (7), p283.
        # Kaplan,G.H., 1981. USNO circular no. 163, pA2.
        #
        # Called:  sla_DEULER
        # -

        # Interval between basic epoch J2000.0 and beginning epoch (JC)
        t0 = (start_epoch - 2000e0) / 100e0

        # Interval over which precession required (JC)
        interval_precision = (end_epoch - start_epoch) / 100e0

        # Euler angles
        tas2_r = interval_precision * cls.AS2R
        w = 2306.2181e0 + (1.39656e0 - 0.000139e0 * t0) * t0

        zeta = (
            w
            + ((0.30188e0 - 0.000344e0 * t0) + 0.017998e0 * interval_precision)
            * interval_precision
        ) * tas2_r
        z = (
            w
            + ((1.09468e0 + 0.000066e0 * t0) + 0.018203e0 * interval_precision)
            * interval_precision
        ) * tas2_r
        theta = (
            (2004.3109e0 + (-0.85330e0 - 0.000217e0 * t0) * t0)
            + ((-0.42665e0 - 0.000217e0 * t0) - 0.041833e0 * interval_precision)
            * interval_precision
        ) * tas2_r

        # Rotation matrix
        return cls.deuler("ZYZ", -zeta, theta, -z)

    @classmethod
    def preces(cls, system_name, start_epoch, end_epoch, input_ra, input_dec):
        # +
        # - - - - - - -
        # P R E C E S
        # - - - - - - -
        #
        # Precession - either FK4 (Bessel-Newcomb, pre IAU 1976) or
        # FK5 (Fricke, post IAU 1976) as required.
        #
        # Given:
        # SYSTEM     char   precession to be applied: 'FK4' or 'FK5'
        # EP0,EP1    dp     starting and ending epoch
        # RA,DC      dp     RA,Dec, mean equator & equinox of epoch EP0
        #
        # Returned:
        # RA,DC      dp     RA,Dec, mean equator & equinox of epoch EP1
        #
        # Called:    sla_DRANRM, sla_PREBN, sla_PREC, sla_DCS2C,
        # sla_DMXV, sla_DCC2S
        #
        # Notes:
        #
        # 1)  Lowercase characters in SYSTEM are acceptable.
        #
        # 2)  The epochs are Besselian if SYSTEM='FK4' and Julian if 'FK5'.
        # For example, to precess coordinates in the old system from
        # equinox 1900.0 to 1950.0 the call would be:
        # CALL sla_PRECES ('FK4', 1900D0, 1950D0, RA, DC)
        #
        # 3)  This routine will NOT correctly convert between the old and
        # the new systems - for example conversion from B1950 to J2000.
        # For these purposes see sla_FK425, sla_FK524, sla_FK45Z and
        # sla_FK54Z.
        #
        # 4)  If an invalid SYSTEM is supplied, values of -99,-99 will
        # be returned for both RA and DC.
        # -

        # Convert to uppercase and validate SYSTEM
        if system_name.upper() == "FK4" or system_name.upper() != "FK4" and system_name.upper() == "FK5":
            # Generate appropriate precession matrix
            if system_name.upper() == "FK4":
                prec_mat = cls.prebn(start_epoch, end_epoch)
            else:
                prec_mat = cls.prec(start_epoch, end_epoch)

            # Convert RA,Dec to x,y,z
            v1 = cls.dcs2c(input_ra, input_dec)

            # Precess
            v2 = cls.dmxv(prec_mat, v1)

            # Back to RA,Dec
            ra, dec = cls.dcc2s(v2)
            ra = cls.dranrm(ra)

        else:
            ra = -99e0
            dec = -99e0
        return ra, dec

    @classmethod
    def precl(cls, epoch_start, epoch_end):
        # +
        # - - - - - -
        # P R E C L
        # - - - - - -
        #
        # Form the matrix of precession between two epochs, using the
        # model of Simon et al (1994), which is suitable for long
        # periods of time.
        #
        # (double precision)
        #
        # Given:
        # EP0    dp         beginning epoch
        # EP1    dp         ending epoch
        #
        # Returned:
        # RMATP  dp(3,3)    precession matrix
        #
        # Notes:
        #
        # 1)  The epochs are TDB Julian epochs.
        #
        # 2)  The matrix is in the sense   V(EP1)  =  RMATP * V(EP0)
        #
        # 3)  The absolute accuracy of the model is limited by the
        # uncertainty in the general precession, about 0.3 arcsec per
        # 1000 years.  The remainder of the formulation provides a
        # precision of 1 mas over the interval from 1000AD to 3000AD,
        # 0.1 arcsec from 1000BC to 5000AD and 1 arcsec from
        # 4000BC to 8000AD.
        #
        # Reference:
        # Simon, J.L. et al., 1994. Astron.Astrophys., 282, 663-683.
        #
        # Called:  sla_DEULER
        # -

        # Interval between basic epoch J2000.0 and beginning epoch (1000JY)
        t0 = (epoch_start - 2000e0) / 1000e0

        # Interval over which precession required (1000JY)
        t = (epoch_end - epoch_start) / 1000e0

        # Euler angles
        tas2_r = t * cls.AS2R
        w = (
            23060.9097e0
            + (
                139.7459e0
                + (-0.0038e0 + (-0.5918e0 + (-0.0037e0 + 0.0007e0 * t0) * t0) * t0) * t0
            )
            * t0
        )

        zeta = (
            w
            + (
                30.2226e0
                + (-0.2523e0 + (-0.3840e0 + (-0.0014e0 + 0.0007e0 * t0) * t0) * t0) * t0
                + (
                    18.0183e0
                    + (-0.1326e0 + (0.0006e0 + 0.0005e0 * t0) * t0) * t0
                    + (
                        -0.0583e0
                        + (-0.0001e0 + 0.0007e0 * t0) * t0
                        + (-0.0285e0 + (-0.0002e0) * t) * t
                    )
                    * t
                )
                * t
            )
            * t
        ) * tas2_r

        z = (
            w
            + (
                109.5270e0
                + (0.2446e0 + (-1.3913e0 + (-0.0134e0 + 0.0026e0 * t0) * t0) * t0) * t0
                + (
                    18.2667e0
                    + (-1.1400e0 + (-0.0173e0 + 0.0044e0 * t0) * t0) * t0
                    + (
                        -0.2821e0
                        + (-0.0093e0 + 0.0032e0 * t0) * t0
                        + (-0.0301e0 + 0.0006e0 * t0 - 0.0001e0 * t) * t
                    )
                    * t
                )
                * t
            )
            * t
        ) * tas2_r

        theta = (
            20042.0207e0
            + (
                -85.3131e0
                + (-0.2111e0 + (0.3642e0 + (0.0008e0 + (-0.0005e0) * t0) * t0) * t0)
                * t0
            )
            * t0
            + (
                -42.6566e0
                + (-0.2111e0 + (0.5463e0 + (0.0017e0 + (-0.0012e0) * t0) * t0) * t0)
                * t0
                + (
                    -41.8238e0
                    + (0.0359e0 + (0.0027e0 + (-0.0001e0) * t0) * t0) * t0
                    + (
                        -0.0731e0
                        + (0.0019e0 + 0.0009e0 * t0) * t0
                        + (-0.0127e0 + 0.0011e0 * t0 + 0.0004e0 * t) * t
                    )
                    * t
                )
                * t
            )
            * t
        ) * tas2_r

        # Rotation matrix
        return cls.deuler("ZYZ", -zeta, theta, -z)

    @classmethod
    def dh2e(cls, az, el, phi):
        # +
        # - - - - -
        # D E 2 H
        # - - - - -
        #
        # Horizon to equatorial coordinates:  Az,El to HA,Dec
        #
        # (double precision)
        #
        # Given:
        # AZ      d     azimuth
        # EL      d     elevation
        # PHI     d     observatory latitude
        #
        # Returned:
        # HA      d     hour angle
        # DEC     d     declination
        #
        # Notes:
        #
        # 1)  All the arguments are angles in radians.
        #
        # 2)  The sign convention for azimuth is north zero, east +pi/2.
        #
        # 3)  HA is returned in the range +/-pi.  Declination is returned
        # in the range +/-pi/2.
        #
        # 4)  The latitude is (in principle) geodetic.  In critical
        # applications, corrections for polar motion should be applied.
        #
        # 5)  In some applications it will be important to specify the
        # correct type of elevation in order to produce the required
        # type of HA,Dec.  In particular, it may be important to
        # distinguish between the elevation as affected by refraction,
        # which will yield the "observed" HA,Dec, and the elevation
        # in vacuo, which will yield the "topocentric" HA,Dec.  If the
        # effects of diurnal aberration can be neglected, the
        # topocentric HA,Dec may be used as an approximation to the
        # "apparent" HA,Dec.
        #
        # 6)  No range checking of arguments is done.
        #
        # 7)  In applications which involve many such calculations, rather
        # than calling the present routine it will be more efficient to
        # use inline code, having previously computed fixed terms such
        # as sine and cosine of latitude.
        # -

        # Useful trig functions
        sa = np.sin(az)
        ca = np.cos(az)
        se = np.sin(el)
        ce = np.cos(el)
        sp = np.sin(phi)
        cp = np.cos(phi)

        # HA,Dec as x,y,z
        x = -ca * ce * sp + se * cp
        y = -sa * ce
        z = ca * ce * cp + se * sp

        # To HA,Dec
        r = np.sqrt(x * x + y * y)

        ha = np.where(r == 0, 0, np.arctan2(y, x))
        dec = np.arctan2(z, r)
        return ha, dec

    @classmethod
    def gmsta(cls, date, ut=None):
        # +
        # - - - - - -
        # G M S T A
        # - - - - - -
        #
        # Conversion from Universal Time to Greenwich mean sidereal time,
        # with rounding errors minimized.
        #
        # double precision
        #
        # Given:
        # DATE    d      UT1 date (MJD: integer part of JD-2400000.5))
        # UT      d      UT1 time (fraction of a day)
        #
        # The result is the Greenwich mean sidereal time (double precision,
        # radians, in the range 0 to 2pi).
        #
        # There is no restriction on how the UT is apportioned between the
        # DATE and UT arguments.  Either of the two arguments could, for
        # example, be zero and the entire date+time supplied in the other.
        # However, the routine is designed to deliver maximum accuracy when
        # the DATE argument is a whole number and the UT lies in the range
        # 0 to 1 (or vice versa).
        #
        # The algorithm is based on the IAU 1982 expression (see page S15 of
        # the 1984 Astronomical Almanac).  This is always described as giving
        # the GMST at 0 hours UT1.  In fact, it gives the difference between
        # the GMST and the UT, the steady 4-minutes-per-day drawing-ahead of
        # ST with respect to UT.  When whole days are ignored, the expression
        # happens to equal the GMST at 0 hours UT1 each day.  Note that the
        # factor 1.0027379... does not appear explicitly but in the form of
        # the coefficient 8640184.812866, which is 86400x36525x0.0027379...
        #
        # In this routine, the entire UT1 (the sum of the two arguments DATE
        # and UT) is used directly as the argument for the standard formula.
        # The UT1 is then added, but omitting whole days to conserve accuracy.
        #
        # See also the routine sla_GMST, which accepts the UT as a single
        # argument.  Compared with sla_GMST, the extra numerical precision
        # delivered by the present routine is unlikely to be important in
        # an absolute sense, but may be useful when critically comparing
        # algorithms and in applications where two sidereal times close
        # together are differenced.
        #
        # Called:  sla_DRANRM
        # -

        if ut is None:
            ut = np.mod(date, 1.0)
            date = np.trunc(date)

        # Seconds of time to radians

        # Julian centuries since J2000.
        date_lt_ut = date < ut
        d1 = np.where(date_lt_ut, date, ut)
        d2 = np.where(date_lt_ut, ut, date)

        T = (d1 + (d2 - 51544.5e0)) / 36525e0

        # GMST at this UT1.
        return cls.dranrm(
            cls.DS2R
            * (
                24110.54841e0
                + (8640184.812866e0 + (0.093104e0 - 6.2e-6 * T) * T) * T
                + 86400e0 * (np.mod(d1, 1e0) + np.mod(d2, 1e0))
            )
        )

    @classmethod
    def eqgal(cls, in_ra, in_dec):
        # +
        # - - - - - -
        # E Q G A L
        # - - - - - -
        #
        # Transformation from J2000.0 equatorial coordinates to
        # IAU 1958 galactic coordinates (double precision)
        #
        # Given:
        # DR,DD       dp       J2000.0 RA,Dec
        #
        # Returned:
        # DL,DB       dp       galactic longitude and latitude L2,B2
        #
        # (all arguments are radians)
        #
        # Called:
        # sla_DCS2C, sla_DMXV, sla_DCC2S, sla_DRANRM, sla_DRANGE
        #
        # Note:
        # The equatorial coordinates are J2000.0.  Use the routine
        # sla_EG50 if conversion from B1950.0 'FK4' coordinates is
        # required.

        # L2,B2 system of galactic coordinates
        #
        # P = 192.25       RA of galactic north pole (mean B1950.0)
        # Q =  62.6        inclination of galactic to mean B1950.0 equator
        # R =  33          longitude of ascending node
        #
        # P,Q,R are degrees
        #
        # Equatorial to galactic rotation matrix (J2000.0), obtained by
        # applying the standard FK4 to FK5 transformation, for zero proper
        # motion in FK5, to the columns of the B1950 equatorial to
        # galactic rotation matrix:

        RMAT = np.array(
            [
                [-0.054875539726e0, -0.873437108010e0, -0.483834985808e0],
                [0.494109453312e0, -0.444829589425e0, 0.746982251810e0],
                [-0.867666135858e0, -0.198076386122e0, 0.455983795705e0],
            ]
        )

        # Spherical to Cartesian
        V1 = cls.dcs2c(in_ra, in_dec)

        # Equatorial to galactic
        V2 = cls.dmxv(RMAT, V1)

        # Cartesian to spherical
        r_gal_l, r_gal_b = cls.dcc2s(V2)

        # Express in conventional ranges
        r_gal_l = cls.dranrm(r_gal_l)
        r_gal_b = cls.drange(r_gal_b)

        return r_gal_l, r_gal_b

    @classmethod
    def drange(cls, in_angle):
        rv = np.mod(in_angle, cls.D2PI)
        rv = np.where(np.abs(rv) >= cls.DPI, rv - np.sign(in_angle) * cls.D2PI, rv)
        return rv

    @classmethod
    def addet(cls, rm, dm, eq):
        # +
        # - - - - - -
        # A D D E T
        # - - - - - -
        #
        # Add the E-terms (elliptic component of annual aberration)
        # to a pre IAU 1976 mean place to conform to the old
        # catalogue convention (double precision)
        #
        # Given:
        # RM,DM     dp     RA,Dec (radians) without E-terms
        # EQ        dp     Besselian epoch of mean equator and equinox
        #
        # Returned:
        # RC,DC     dp     RA,Dec (radians) with E-terms included
        #
        # Note:
        #
        # Most star positions from pre-1984 optical catalogues (or
        # derived from astrometry using such stars) embody the
        # E-terms.  If it is necessary to convert a formal mean
        # place (for example a pulsar timing position) to one
        # consistent with such a star catalogue, then the RA,Dec
        # should be adjusted using this routine.
        #
        # Reference:
        # Explanatory Supplement to the Astronomical Ephemeris,
        # section 2D, page 48.
        #
        # Depends:  sla_ETRMS, sla_DCS2C, sla_DCC2S, sla_DRANRM, sla_DRANGE

        # E-terms vector
        A = cls.etrms(eq)

        # Spherical to Cartesian
        V = cls.dcs2c(rm, dm)

        # Include the E-terms
        V += A

        # Cartesian to spherical
        rc, dc = cls.dcc2s(V)

        # Bring RA into conventional range
        rc = cls.dranrm(rc)

        return rc, dc

    @classmethod
    def etrms(cls, ep):
        # +
        # - - - - - -
        # E T R M S
        # - - - - - -
        #
        # Compute the E-terms (elliptic component of annual aberration)
        # vector (double precision)
        #
        # Given:
        # EP      dp      Besselian epoch
        #
        # Returned:
        # EV      dp(3)   E-terms as (dx,dy,dz)
        #
        # Note the use of the J2000 aberration constant (20.49552 arcsec).
        # This is a reflection of the fact that the E-terms embodied in
        # existing star catalogues were computed from a variety of
        # aberration constants.  Rather than adopting one of the old
        # constants the latest value is used here.
        #
        # References:
        # 1  Smith, C.A. et al., 1989.  Astr.J. 97, 265.
        # 2  Yallop, B.D. et al., 1989.  Astr.J. 97, 274.

        # Julian centuries since B1950
        T = (ep - 1950e0) * 1.00002135903e-2

        # Eccentricity
        E = 0.01673011e0 - (0.00004193e0 + 0.000000126e0 * T) * T

        # Mean obliquity
        E0 = (
            84404.836e0 - (46.8495e0 + (0.00319e0 + 0.00181e0 * T) * T) * T
        ) * cls.AS2R

        # Mean longitude of perihelion
        P = (1015489.951e0 + (6190.67e0 + (1.65e0 + 0.012e0 * T) * T) * T) * cls.AS2R

        # E-terms
        EK = E * 20.49552e0 * cls.AS2R
        CP = np.cos(P)

        return np.array([EK * np.sin(P), -EK * CP * np.cos(E0), -EK * CP * np.sin(E0)])

    @classmethod
    def airmas(cls, zd):
        # +
        # - - - - - - -
        # A I R M A S
        # - - - - - - -
        #
        # Air mass at given zenith distance (double precision)
        #
        # Given:
        # ZD     d     Observed zenith distance (radians)
        #
        # The result is an estimate of the air mass, in units of that
        # at the zenith.
        #
        # Notes:
        #
        # 1)  The "observed" zenith distance referred to above means "as
        # affected by refraction".
        #
        # 2)  Uses Hardie's (1962) polynomial fit to Bemporad's data for
        # the relative air mass, X, in units of thickness at the zenith
        # as tabulated by Schoenberg (1929). This is adequate for all
        # normal needs as it is accurate to better than 0.1% up to X =
        # 6.8 and better than 1% up to X = 10. Bemporad's tabulated
        # values are unlikely to be trustworthy to such accuracy
        # because of variations in density, pressure and other
        # conditions in the atmosphere from those assumed in his work.
        #
        # 3)  The sign of the ZD is ignored.
        #
        # 4)  At zenith distances greater than about ZD = 87 degrees the
        # air mass is held constant to avoid arithmetic overflows.
        #
        # References:
        # Hardie, R.H., 1962, in "Astronomical Techniques"
        # ed. W.A. Hiltner, University of Chicago Press, p180.
        # Schoenberg, E., 1929, Hdb. d. Ap.,
        # Berlin, Julius Springer, 2, 268.

        SECZM1 = 1e0 / (np.cos(np.minimum(1.52e0, np.abs(zd)))) - 1e0
        return 1e0 + SECZM1 * (
            0.9981833e0 - SECZM1 * (0.002875e0 + 0.0008083e0 * SECZM1)
        )

    @classmethod
    def altaz(cls, ha, dec, phi):

        # +
        # - - - - - -
        # A L T A Z
        # - - - - - -
        #
        # Positions, velocities and accelerations for an altazimuth
        # telescope mount.
        #
        # (double precision)
        #
        # Given:
        # HA      d     hour angle
        # DEC     d     declination
        # PHI     d     observatory latitude
        #
        # Returned:
        # A      d     azimuth
        # AD     d        "    velocity
        # ADD    d        "    acceleration
        # E      d     elevation
        # ED     d         "     velocity
        # EDD    d         "     acceleration
        # Q      d     parallactic angle
        # QD     d         "      "   velocity
        # QDD    d         "      "   acceleration
        #
        # Notes:
        #
        # 1)  Natural units are used throughout.  HA, DEC, PHI, AZ, EL
        # and ZD are in radians.  The velocities and accelerations
        # assume constant declination and constant rate of change of
        # hour angle (as for tracking a star);  the units of AZD, ELD
        # and PAD are radians per radian of HA, while the units of AZDD,
        # ELDD and PADD are radians per radian of HA squared.  To
        # convert into practical degree- and second-based units:
        #
        # angles * 360/2pi -> degrees
        # velocities * (2pi/86400)*(360/2pi) -> degree/sec
        # accelerations * ((2pi/86400)**2)*(360/2pi) -> degree/sec/sec
        #
        # Note that the seconds here are sidereal rather than SI.  One
        # sidereal second is about 0.99727 SI seconds.
        #
        # The velocity and acceleration factors assume the sidereal
        # tracking case.  Their respective numerical values are (exactly)
        # 1/240 and (approximately) 1/3300236.9.
        #
        # 2)  Azimuth is returned in the range 0-2pi;  north is zero,
        # and east is +pi/2.  Elevation and parallactic angle are
        # returned in the range +/-pi.  Parallactic angle is +ve for
        # a star west of the meridian and is the angle NP-star-zenith.
        #
        # 3)  The latitude is geodetic as opposed to geocentric.  The
        # hour angle and declination are topocentric.  Refraction and
        # deficiencies in the telescope mounting are ignored.  The
        # purpose of the routine is to give the general form of the
        # quantities.  The details of a real telescope could profoundly
        # change the results, especially close to the zenith.
        #
        # 4)  No range checking of arguments is carried out.
        #
        # 5)  In applications which involve many such calculations, rather
        # than calling the present routine it will be more efficient to
        # use inline code, having previously computed fixed terms such
        # as sine and cosine of latitude, and (for tracking a star)
        # sine and cosine of declination.

        # Useful functions
        SH = np.sin(ha)
        CH = np.cos(ha)
        SD = np.sin(dec)
        CD = np.cos(dec)
        SP = np.sin(phi)
        CP = np.cos(phi)
        CHCD = CH * CD
        SDCP = SD * CP
        X = -CHCD * SP + SDCP
        Y = -SH * CD
        Z = CHCD * CP + SD * SP
        RSQ = X * X + Y * Y
        R = np.sqrt(RSQ)

        # Azimuth and elevation
        A = np.where(RSQ == 0e0, 0e0, np.arctan2(Y, X))
        A = np.where(A < 0e0, A + cls.D2PI, A)
        E = np.arctan2(Z, R)

        # Parallactic angle
        C = CD * SP - CH * SDCP
        S = SH * CP

        Q = np.where(C * C + S * S > 0, np.arctan2(S, C), cls.DPI - ha)

        # Velocities and accelerations (clamped at zenith/nadir)
        R = np.clip(R, cls.TINY, None)
        RSQ = np.clip(RSQ, cls.TINY, None)

        QD = -X * CP / RSQ
        AD = SP + Z * QD
        ED = CP * Y / R
        EDR = ED / R
        ADD = EDR * (Z * SP + (2e0 - RSQ) * QD)
        EDD = -R * QD * AD
        QDD = EDR * (SP + 2e0 * Z * QD)

        return A, AD, ADD, E, ED, EDD, Q, QD, QDD

    @classmethod
    def amp(cls, ra, da, date, eq):
        # +
        # - - - -
        # A M P
        # - - - -
        #
        # Convert star RA,Dec from geocentric apparent to mean place
        #
        # The mean coordinate system is the post IAU 1976 system,
        # loosely called FK5.
        #
        # Given:
        # RA       d      apparent RA (radians)
        # DA       d      apparent Dec (radians)
        # DATE     d      TDB for apparent place (JD-2400000.5)
        # EQ       d      equinox:  Julian epoch of mean place
        #
        # Returned:
        # RM       d      mean RA (radians)
        # DM       d      mean Dec (radians)
        #
        # References:
        # 1984 Astronomical Almanac, pp B39-B41.
        # (also Lederle & Schwan, Astron. Astrophys. 134,
        # 1-6, 1984)
        #
        # Notes:
        #
        # 1)  The distinction between the required TDB and TT is always
        # negligible.  Moreover, for all but the most critical
        # applications UTC is adequate.
        #
        # 2)  Iterative techniques are used for the aberration and light
        # deflection corrections so that the routines sla_AMP (or
        # sla_AMPQK) and sla_MAP (or sla_MAPQK) are accurate inverses;
        # even at the edge of the Sun's disc the discrepancy is only
        # about 1 nanoarcsecond.
        #
        # 3)  Where multiple apparent places are to be converted to mean
        # places, for a fixed date and equinox, it is more efficient to
        # use the sla_MAPPA routine to compute the required parameters
        # once, followed by one call to sla_AMPQK per star.
        #
        # 4)  The accuracy is sub-milliarcsecond, limited by the
        # precession-nutation model (IAU 1976 precession, Shirai &
        # Fukushima 2001 forced nutation and precession corrections).
        #
        # 5)  The accuracy is further limited by the routine sla_EVP, called
        # by sla_MAPPA, which computes the Earth position and velocity
        # using the methods of Stumpff.  The maximum error is about
        # 0.3 mas.
        #
        # Depends:  sla_MAPPA, sla_AMPQK

        amprms = cls.mappa(eq, date)
        rm, dm = cls.ampqk(ra, da, amprms)

        return rm, dm

    @classmethod
    def ampqk(cls, ra, da, amprms):
        # +
        # - - - - - -
        # A M P Q K
        # - - - - - -
        #
        # Convert star RA,Dec from geocentric apparent to mean place
        #
        # The mean coordinate system is the post IAU 1976 system,
        # loosely called FK5.
        #
        # Use of this routine is appropriate when efficiency is important
        # and where many star positions are all to be transformed for
        # one epoch and equinox.  The star-independent parameters can be
        # obtained by calling the sla_MAPPA routine.
        #
        # Given:
        # RA       d      apparent RA (radians)
        # DA       d      apparent Dec (radians)
        #
        # AMPRMS   d(21)  star-independent mean-to-apparent parameters:
        #
        # (0)      time interval for proper motion (Julian years)
        # (1)    barycentric position of the Earth (AU)
        # (2)    heliocentric direction of the Earth (unit vector)
        # (3)      (grav rad Sun)*2/(Sun-Earth distance)
        # (4)   ABV: barycentric Earth velocity in units of c
        # (5)     sqrt(1-v**2) where v=modulus(ABV)
        # (6)  precession/nutation (3,3) matrix
        #
        # Returned:
        # RM       d      mean RA (radians)
        # DM       d      mean Dec (radians)
        #
        # References:
        # 1984 Astronomical Almanac, pp B39-B41.
        # (also Lederle & Schwan, Astron. Astrophys. 134,
        # 1-6, 1984)
        #
        # Note:
        #
        # Iterative techniques are used for the aberration and
        # light deflection corrections so that the routines
        # sla_AMP (or sla_AMPQK) and sla_MAP (or sla_MAPQK) are
        # accurate inverses;  even at the edge of the Sun's disc
        # the discrepancy is only about 1 nanoarcsecond.
        #
        # Depends:  sla_DCS2C, sla_DIMXV, sla_DVDV, sla_DVN, sla_DCC2S,
        # sla_DRANRM

        # Unpack scalar and vector parameters
        GR2E = amprms[3]
        AB1 = amprms[5]

        EHN = amprms[2]
        ABV = amprms[4]

        # Apparent RA,Dec to Cartesian
        P3 = cls.dcs2c(ra, da)

        # Precession and nutation
        P2 = cls.dimxv(amprms[6], P3)

        # Aberration
        AB1P1 = AB1 + 1e0
        P1 = P2.copy()

        for _ in range(2):
            P1DV = cls.dvdv(P1, ABV)
            P1DVP1 = 1e0 + P1DV
            W = 1e0 + P1DV / AB1P1
            P1 = (P1DVP1 * P2 - W * ABV) / AB1

            P3, W = cls.dvn(P1)
            P1 = P3.copy()

        # Light deflection
        P = P1.copy()

        for _ in range(5):
            PDE = cls.dvdv(P, EHN)
            PDEP1 = 1e0 + PDE
            W = PDEP1 - GR2E * PDE
            P = (PDEP1 * P1 - GR2E * EHN) / W

            P2, W = cls.dvn(P)
            P = P2.copy()

        # Mean RA,Dec
        rm, dm = cls.dcc2s(P)
        rm = cls.dranrm(rm)

        return rm, dm

    @classmethod
    def dvdv(cls, va, vb):
        # +
        # - - - - -
        # D V D V
        # - - - - -
        #
        # Scalar product of two 3-vectors
        #
        # Given:
        # VA      dp(3)     first vector
        # VB      dp(3)     second vector

        return np.dot(va, vb)

    @classmethod
    def dvn(cls, v):
        # +
        # - - - -
        # D V N
        # - - - -
        #
        # Normalizes a 3-vector also giving the modulus (double precision)
        #
        # Given:
        # V       d(3)      vector
        #
        # Returned:
        # UV      d(3)      unit vector in direction of V
        # VM      d         modulus of V
        #
        # Notes:
        # If the modulus of V is zero, UV is set to zero as well.

        # Modulus.
        vm = np.sqrt(np.dot(v, v))

        # Normalize the vector.
        vm = np.where(vm <= 0e0, 1e0, vm)

        return v / vm, vm

    @classmethod
    def mappa(cls, eq, date):
        # +
        # - - - - - -
        # M A P P A
        # - - - - - -
        #
        # Compute star-independent parameters in preparation for
        # conversions between mean place and geocentric apparent place.
        #
        # The parameters produced by this routine are required in the
        # parallax, light deflection, aberration, and precession/nutation
        # parts of the mean/apparent transformations.
        #
        # The reference frames and timescales used are post IAU 1976.
        #
        # Given:
        # EQ       d      epoch of mean equinox to be used (Julian)
        # DATE     d      TDB (JD-2400000.5)
        #
        # Returned:
        # AMPRMS   d(7)  star-independent mean-to-apparent parameters:
        #
        # (0)      time interval for proper motion (Julian years)
        # (1)    barycentric position of the Earth (AU)
        # (2)    heliocentric direction of the Earth (unit vector)
        # (3)      (grav rad Sun)*2/(Sun-Earth distance)
        # (4)   ABV: barycentric Earth velocity in units of c
        # (5)     sqrt(1-v**2) where v=modulus(ABV)
        # (6)  precession/nutation (3,3) matrix
        #
        # References:
        # 1984 Astronomical Almanac, pp B39-B41.
        # (also Lederle & Schwan, Astron. Astrophys. 134,
        # 1-6, 1984)
        #
        # Notes:
        #
        # 1)  For DATE, the distinction between the required TDB and TT
        # is always negligible.  Moreover, for all but the most
        # critical applications UTC is adequate.
        #
        # 2)  The vectors AMPRMS(2-4) and AMPRMS(5-7) are referred to
        # the mean equinox and equator of epoch EQ.
        #
        # 3)  The parameters AMPRMS produced by this routine are used by
        # sla_AMPQK, sla_MAPQK and sla_MAPQKZ.
        #
        # 4)  The accuracy is sub-milliarcsecond, limited by the
        # precession-nutation model (IAU 1976 precession, Shirai &
        # Fukushima 2001 forced nutation and precession corrections).
        #
        # 5)  A further limit to the accuracy of routines using the parameter
        # array AMPRMS is imposed by the routine sla_EVP, used here to
        # compute the Earth position and velocity by the methods of
        # Stumpff.  The maximum error in the resulting aberration
        # corrections is about 0.3 milliarcsecond.
        #
        # Called:
        # sla_EPJ         MDJ to Julian epoch
        # sla_EVP         earth position & velocity
        # sla_DVN         normalize vector
        # sla_PRENUT      precession/nutation matrix
        # Time interval for proper motion correction

        amprms = [None for _ in range(7)]
        amprms[0] = cls.epj(date) - eq

        # Get Earth barycentric and heliocentric position and velocity
        ebd, amprms[1], EHD, EH = cls.evp(date, eq)

        # Heliocentric direction of earth (normalized) and modulus
        amprms[2], E = cls.dvn(EH)

        # Light deflection parameter
        amprms[3] = cls.GR2 / E

        # Aberration parameters
        amprms[4] = ebd * cls.CR

        VN, VM = cls.dvn(ebd * cls.CR)
        amprms[5] = np.sqrt(1e0 - VM * VM)

        # Precession/nutation matrix
        amprms[6] = cls.prenut(eq, date)

        return amprms

    @classmethod
    def epj(cls, date):
        # +
        # - - - -
        # E P J
        # - - - -
        #
        # Conversion of Modified Julian Date to Julian Epoch (double precision)
        #
        # Given:
        # DATE     dp       Modified Julian Date (JD - 2400000.5)
        #
        # The result is the Julian Epoch.
        #
        # Reference:
        # Lieske,J.H., 1979. Astron.Astrophys.,73,282.

        return 2000e0 + (date - 51544.5e0) / 365.25

    @classmethod
    def evp(cls, date, deqx):
        # +
        # - - - -
        # E V P
        # - - - -
        #
        # Barycentric and heliocentric velocity and position of the Earth
        #
        # All arguments are double precision
        #
        # Given:
        #
        # DATE          TDB (loosely ET) as a Modified Julian Date
        # (JD-2400000.5)
        #
        # DEQX          Julian Epoch (e.g. 2000.0D0) of mean equator and
        # equinox of the vectors returned.  If DEQX .LE. 0D0,
        # all vectors are referred to the mean equator and
        # equinox (FK5) of epoch DATE.
        #
        # Returned (all 3D Cartesian vectors):
        #
        # DVB,DPB       barycentric velocity, position (AU/s, AU)
        # DVH,DPH       heliocentric velocity, position (AU/s, AU)
        #
        # Depends:  sla_EPJ, sla_PREC
        #
        # Notes:
        #
        # 1  This routine is accurate enough for many purposes but faster and
        # more compact than the sla_EPV routine.  The maximum deviations
        # from the JPL DE96 ephemeris are as follows:
        #
        # barycentric velocity         0.42  m/s
        # barycentric position         6900  km
        #
        # heliocentric velocity        0.42  m/s
        # heliocentric position        1600  km
        #
        # 2  The routine is adapted from the BARVEL and BARCOR subroutines of
        # Stumpff (1980).  Most of the changes are merely cosmetic and do
        # not affect the results at all.  However, some adjustments have
        # been made so as to give results that refer to the IAU 1976 'FK5'
        # equinox and precession, although the differences these changes
        # make relative to the results from Stumpff's original 'FK4' version
        # are smaller than the inherent accuracy of the algorithm.  One
        # minor shortcoming in the original routines that has NOT been
        # corrected is that better numerical accuracy could be achieved if
        # the various polynomial evaluations were nested.
        #
        # Reference:
        #
        # Stumpff, P., Astron.Astrophys.Suppl.Ser. 41, 1-8 (1980).

        # Constants DCFEL(I,K) of fast changing elements
        # I=1                I=2              I=3
        DCFEL = np.array(
            [
                [1.7400353e0, 6.2833195099091e2, 5.2796e-6],
                [6.2565836e0, 6.2830194572674e2, -2.6180e-6],
                [4.7199666e0, 8.3997091449254e3, -1.9780e-5],
                [1.9636505e-1, 8.4334662911720e3, -5.6044e-5],
                [4.1547339e0, 5.2993466764997e1, 5.8845e-6],
                [4.6524223e0, 2.1354275911213e1, 5.6797e-6],
                [4.2620486e0, 7.5025342197656e0, 5.5317e-6],
                [1.4740694e0, 3.8377331909193e0, 5.6093e-6],
            ]
        )

        #
        # Constants DCEPS and CCSEL(I,K) of slowly changing elements
        # I=1           I=2           I=3
        DCEPS = np.array([4.093198e-1, -2.271110e-4, -2.860401e-8])
        CCSEL = np.array(
            [
                [1.675104e-2, -4.179579e-5, -1.260516e-7],
                [2.220221e-1, 2.809917e-2, 1.852532e-5],
                [1.589963, 3.418075e-2, 1.430200e-5],
                [2.994089, 2.590824e-2, 4.155840e-6],
                [8.155457e-1, 2.486352e-2, 6.836840e-6],
                [1.735614, 1.763719e-2, 6.370440e-6],
                [1.968564, 1.524020e-2, -2.517152e-6],
                [1.282417, 8.703393e-3, 2.289292e-5],
                [2.280820, 1.918010e-2, 4.484520e-6],
                [4.833473e-2, 1.641773e-4, -4.654200e-7],
                [5.589232e-2, -3.455092e-4, -7.388560e-7],
                [4.634443e-2, -2.658234e-5, 7.757000e-8],
                [8.997041e-3, 6.329728e-6, -1.939256e-9],
                [2.284178e-2, -9.941590e-5, 6.787400e-8],
                [4.350267e-2, -6.839749e-5, -2.714956e-7],
                [1.348204e-2, 1.091504e-5, 6.903760e-7],
                [3.106570e-2, -1.665665e-4, -1.590188e-7],
            ]
        )

        #
        # Constants of the arguments of the short-period perturbations
        # by the planets:   DCARGS(I,K)
        # I=1               I=2
        DCARGS = np.array(
            [
                [5.0974222e0, -7.8604195454652e2],
                [3.9584962e0, -5.7533848094674e2],
                [1.6338070e0, -1.1506769618935e3],
                [2.5487111e0, -3.9302097727326e2],
                [4.9255514e0, -5.8849265665348e2],
                [1.3363463e0, -5.5076098609303e2],
                [1.6072053e0, -5.2237501616674e2],
                [1.3629480e0, -1.1790629318198e3],
                [5.5657014e0, -1.0977134971135e3],
                [5.0708205e0, -1.5774000881978e2],
                [3.9318944e0, 5.2963464780000e1],
                [4.8989497e0, 3.9809289073258e1],
                [1.3097446e0, 7.7540959633708e1],
                [3.5147141e0, 7.9618578146517e1],
                [3.5413158e0, -5.4868336758022e2],
            ]
        )

        #
        # Amplitudes CCAMPS(N,K) of the short-period perturbations
        # N=1          N=2          N=3          N=4          N=5
        CCAMPS = np.array(
            [
                [-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5, -2.490817e-7],
                [-3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5, -1.823138e-7],
                [6.593466e-7, 1.322572e-5, 9.258695e-6, -4.674248e-7, -3.646275e-7],
                [1.140767e-5, -2.049792e-5, -4.747930e-6, -2.638763e-6, -1.245408e-7],
                [9.516893e-6, -2.748894e-6, -1.319381e-6, -4.549908e-6, -1.864821e-7],
                [7.310990e-6, -1.924710e-6, -8.772849e-7, -3.334143e-6, -1.745256e-7],
                [-2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6, -1.655307e-7],
                [-3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6, -3.736225e-7],
                [3.442177e-7, 2.671323e-6, 1.832858e-6, -2.394688e-7, -3.478444e-7],
                [8.702406e-6, -8.421214e-6, -1.372341e-6, -1.455234e-6, -4.998479e-8],
                [-1.488378e-6, -1.251789e-5, 5.226868e-7, -2.049301e-7, 0.0],
                [-8.043059e-6, -2.991300e-6, 1.473654e-7, -3.154542e-7, 0.0],
                [3.699128e-6, -3.316126e-6, 2.901257e-7, 3.407826e-7, 0.0],
                [2.550120e-6, -1.241123e-6, 9.901116e-8, 2.210482e-7, 0.0],
                [-6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.0],
            ]
        )

        #
        # Constants of the secular perturbations in longitude
        # CCSEC3 and CCSEC(N,K)
        # N=1           N=2           N=3
        CCSEC3 = -7.757020e-8
        CCSEC = np.array(
            [
                [1.289600e-6, 5.550147e-1, 2.076942],
                [3.102810e-5, 4.035027, 3.525565e-1],
                [9.124190e-6, 9.990265e-1, 2.622706],
                [9.793240e-7, 5.508259, 1.559103e1],
            ]
        )

        # Sidereal rate DCSLD in longitude, rate CCSGD in mean anomaly
        DCSLD = 1.990987e-7
        CCSGD = 1.990969e-7

        # Some constants used in the calculation of the lunar contribution
        CCKM = 3.122140e-5
        CCMLD = 2.661699e-6
        CCFDI = 2.399485e-7

        #
        # Constants DCARGM(I,K) of the arguments of the perturbations
        # of the motion of the Moon
        # I=1               I=2
        DCARGM = np.array(
            [
                [5.1679830e0, 8.3286911095275e3],
                [5.4913150e0, -7.2140632838100e3],
                [5.9598530e0, 1.5542754389685e4],
            ]
        )

        #
        # Amplitudes CCAMPM(N,K) of the perturbations of the Moon
        # N=1          N=2           N=3           N=4
        CCAMPM = np.array(
            [
                [1.097594e-1, 2.896773e-7, 5.450474e-2, 1.438491e-7],
                [-2.223581e-2, 5.083103e-8, 1.002548e-2, -2.291823e-8],
                [1.148966e-2, 5.658888e-8, 8.249439e-3, 4.063015e-8],
            ]
        )

        #
        # CCPAMV(K)=A*M*DL/DT (planets), DC1MME=1-MASS(Earth+Moon)
        CCPAMV = np.array([8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12])
        DC1MME = 0.99999696e0

        # CCPAM(K)=A*M(planets), CCIM=INCLINATION(Moon)
        CCPAM = np.array([4.960906e-3, 2.727436e-3, 8.392311e-4, 1.556861e-3])
        CCIM = 8.978749e-2

        #
        # EXECUTION
        # ---------

        # Control parameter IDEQ, and time arguments
        DT = (date - 15019.5e0) / 36525e0

        dt_poly = np.array([1.0, DT, DT * DT])
        # Values of all elements for the instant DATE
        DLOCAL = np.mod(np.dot(DCFEL, dt_poly), cls.D2PI)
        DML = DLOCAL[0]
        FORBEL = DLOCAL[1:]

        DEPS = np.mod(np.dot(DCEPS, dt_poly), cls.D2PI)
        SORBEL = np.mod(np.dot(CCSEL, dt_poly), cls.D2PI)

        # Secular perturbations in longitude
        dt_line = np.array([1.0, DT])
        SN = np.sin(np.mod(np.dot(CCSEC[:, 1:], dt_line), cls.D2PI))

        # Periodic perturbations of the EMB (Earth-Moon barycentre)
        PERTL = (
            CCSEC[0, 0] * SN(1)
            + CCSEC[1, 0] * SN(2)
            + (CCSEC[2, 0] + DT * CCSEC3) * SN(3)
            + CCSEC[3, 0] * SN(4)
        )

        PERTLD = 0.0
        PERTR = 0.0
        PERTRD = 0.0
        for row_ind, (args_row, amps_row) in enumerate(zip(DCARGS, CCAMPS)):
            A = np.mod(np.dot(args_row[:2], dt_line), cls.D2PI)
            COSA = np.cos(A)
            SINA = np.sin(A)
            PERTL += amps_row[0] * COSA + amps_row[1] * SINA
            PERTR += amps_row[2] * COSA + amps_row[3] * SINA
            if row_ind < 11:
                PERTLD += (amps_row[1] * COSA - amps_row[0] * SINA) * amps_row[4]
                PERTRD += (amps_row[3] * COSA - amps_row[2] * SINA) * amps_row[4]

        # Elliptic part of the motion of the EMB
        E = SORBEL[0]
        G = FORBEL[0]
        ESQ = E * E
        DPARAM = 1e0 - ESQ
        PARAM = DPARAM
        TWOE = E + E
        TWOG = G + G
        PHI = TWOE * (
            (1.0 - ESQ * 0.125) * np.sin(G)
            + E * 0.625 * np.sin(TWOG)
            + ESQ * 0.54166667 * np.sin(G + TWOG)
        )

        F = G + PHI
        SINF = np.sin(F)
        COSF = np.cos(F)
        DPSI = DPARAM / (1e0 + (E * COSF))
        PHID = (
            TWOE * CCSGD * ((1.0 + ESQ * 1.5) * COSF + E * (1.25 - SINF * SINF * 0.5))
        )
        PSID = CCSGD * E * SINF / np.sqrt(PARAM)

        # Perturbed heliocentric motion of the EMB
        D1PDRO = 1e0 + PERTR
        DRD = D1PDRO * (PSID + DPSI * PERTRD)
        DRLD = D1PDRO * DPSI * (DCSLD + PHID + PERTLD)
        DTL = np.mod(DML + PHI + PERTL, cls.D2PI)
        DSINLS = np.sin(DTL)
        DCOSLS = np.cos(DTL)
        DXHD = DRD * DCOSLS - DRLD * DSINLS
        DYHD = DRD * DSINLS + DRLD * DCOSLS

        # Influence of eccentricity, evection and variation on the
        # geocentric motion of the Moon
        PERTL = 0.0
        PERTLD = 0.0
        PERTP = 0.0
        PERTPD = 0.0
        for arg_row, amp_row in zip(DCARGM, CCAMPM):
            A = np.mod(np.dot(arg_row, dt_line), cls.D2PI)
            SINA = np.sin(A)
            COSA = np.cos(A)
            PERTL += amp_row[0] * SINA
            PERTLD += amp_row[1] * COSA
            PERTP += amp_row[2] * COSA
            PERTPD -= amp_row[3] * SINA

        # Heliocentric motion of the Earth
        TL = FORBEL[1] + PERTL
        SINLM = np.sin(TL)
        COSLM = np.cos(TL)
        SIGMA = CCKM / (1.0 + PERTP)
        A = SIGMA * (CCMLD + PERTLD)
        B = SIGMA * PERTPD
        DXHD = DXHD + (A * SINLM) + (B * COSLM)
        DYHD = DYHD - (A * COSLM) + (B * SINLM)
        DZHD = -(SIGMA * CCFDI * np.cos(FORBEL[2]))

        # Barycentric motion of the Earth
        DXBD = DXHD * DC1MME
        DYBD = DYHD * DC1MME
        DZBD = DZHD * DC1MME

        SINLP = np.zeros(4)
        COSLP = np.zeros(4)
        for K in range(4):
            PLON = FORBEL[K + 3]
            POMG = SORBEL[K + 1]
            PECC = SORBEL[K + 9]
            TL = np.mod(PLON + 2.0 * PECC * np.sin(PLON - POMG), cls.D2PI)
            SINLP[K] = np.sin(TL)
            COSLP[K] = np.cos(TL)
            DXBD = DXBD + (CCPAMV[K] * (SINLP[K] + PECC * np.sin(POMG)))
            DYBD = DYBD - (CCPAMV[K] * (COSLP[K] + PECC * np.cos(POMG)))
            DZBD = DZBD - (CCPAMV[K] * SORBEL[K + 13] * np.cos(PLON - SORBEL[K + 5]))

        # Transition to mean equator of date
        DCOSEP = np.cos(DEPS)
        DSINEP = np.sin(DEPS)
        DYAHD = DCOSEP * DYHD - DSINEP * DZHD
        DZAHD = DSINEP * DYHD + DCOSEP * DZHD
        DYABD = DCOSEP * DYBD - DSINEP * DZBD
        DZABD = DSINEP * DYBD + DCOSEP * DZBD

        # Heliocentric coordinates of the Earth
        DR = DPSI * D1PDRO
        FLATM = CCIM * np.sin(FORBEL(3))
        A = SIGMA * np.cos(FLATM)
        DXH = DR * DCOSLS - (A * COSLM)
        DYH = DR * DSINLS - (A * SINLM)
        DZH = -(SIGMA * np.sin(FLATM))

        # Barycentric coordinates of the Earth
        DXB = DXH * DC1MME
        DYB = DYH * DC1MME
        DZB = DZH * DC1MME
        for K in range(4):
            FLAT = SORBEL[K + 13] * np.sin(FORBEL[K + 3] - SORBEL[K + 5])
            A = CCPAM[K] * (1.0 - SORBEL[K + 9] * np.cos(FORBEL[K + 3] - SORBEL[K + 1]))
            B = A * np.cos(FLAT)
            DXB = DXB - (B * COSLP[K])
            DYB = DYB - (B * SINLP[K])
            DZB = DZB - (A * np.sin(FLAT))

        # Transition to mean equator of date
        DYAH = DCOSEP * DYH - DSINEP * DZH
        DZAH = DSINEP * DYH + DCOSEP * DZH
        DYAB = DCOSEP * DYB - DSINEP * DZB
        DZAB = DSINEP * DYB + DCOSEP * DZB

        # Copy result components into vectors, correcting for FK4 equinox
        DEPJ = cls.epj(date)
        DEQCOR = cls.AS2R * (0.035e0 + 0.00085e0 * (DEPJ - cls.B1950))
        dvh = np.array([DXHD - DEQCOR * DYAHD, DYAHD + DEQCOR * DXHD, DZAHD])
        dvb = np.array([DXBD - DEQCOR * DYABD, DYABD + DEQCOR * DXBD, DZABD])
        dph = np.array([DXH - DEQCOR * DYAH, DYAH + DEQCOR * DXH, DZAH])
        dpb = np.array([DXB - DEQCOR * DYAB, DYAB + DEQCOR * DXB, DZAB])

        # Was precession to another equinox requested?
        if deqx > 0.0:
            # Yes: compute precession matrix from MJD DATE to Julian epoch DEQX
            DPREMA = cls.prec(DEPJ, deqx)

            # Rotate DVH
            dvh = DPREMA @ dvh

            # Rotate DVB
            dvb = DPREMA @ dvb

            # Rotate DPH
            dph = DPREMA @ dph

            # Rotate DPB
            dpb = DPREMA @ dpb

        return dvb, dpb, dvh, dph

    @classmethod
    def prenut(cls, epoch, date):
        # +
        # - - - - - - -
        # P R E N U T
        # - - - - - - -
        #
        # Form the matrix of precession and nutation (SF2001)
        # (double precision)
        #
        # Given:
        # EPOCH   dp         Julian Epoch for mean coordinates
        # DATE    dp         Modified Julian Date (JD-2400000.5)
        # for true coordinates
        #
        # Returned:
        # RMATPN  dp(3,3)    combined precession/nutation matrix
        #
        # Depends:  sla_PREC, sla_EPJ, sla_NUT, sla_DMXM
        #
        # Notes:
        #
        # 1)  The epoch and date are TDB (loosely ET).  TT will do, or even
        # UTC.
        #
        # 2)  The matrix is in the sense   V(true) = RMATPN * V(mean)
        # -

        # Precession
        RMATP = cls.prec(epoch, cls.epj(date))

        # Nutation
        RMATN = cls.nut(date)

        return cls.dmxm(RMATN, RMATP)

    @classmethod
    def polmo(cls, elongm, phim, xp, yp):
        # +
        # - - - - - -
        # P O L M O
        #
        # Polar motion:  correct site longitude and latitude for polar
        # motion and calculate azimuth difference between celestial and
        # terrestrial poles.
        # Given:
        # elongm   d	  mean longitude of the observer (radians, east +ve)
        # phim	 d	  mean geodetic latitude of the observer (radians)
        # xp	   d	  polar motion x-coordinate (radians)
        # yp	   d	  polar motion y-coordinate (radians)
        # Returned:
        # elong	d	  true longitude of the observer (radians, east +ve)
        # phi	  d	  true geodetic latitude of the observer (radians)
        # daz	  d	  azimuth correction (terrestrial-celestial, radians)
        # Notes:
        # 1)  "Mean" longitude and latitude are the (fixed) values for the
        # site's location with respect to the IERS terrestrial reference
        # frame;  the latitude is geodetic.  TAKE CARE WITH THE LONGITUDE
        # np.sign CONVENTION.  The longitudes used by the present routine
        # are east-positive, in accordance with geographical convention
        # (and right-handed).  In particular, note that the longitudes
        # returned by the sla_OBS routine are west-positive, following
        # astronomical usage, and must be reversed in sign before use in
        # the present routine.
        # 2)  xp and yp are the (changing) coordinates of the Celestial
        # Ephemeris Pole with respect to the IERS Reference Pole.
        # xp is positive along the meridian at longitude 0 degrees,
        # and yp is positive along the meridian at longitude
        # 270 degrees (i.e. 90 degrees west).  Values for xp,yp can
        # be obtained from IERS circulars and equivalent publications;
        # the maximum amplitude observed so far is about 0.3 arcseconds.
        # 3)  "True" longitude and latitude are the (moving) values for
        # the site's location with respect to the celestial ephemeris
        # pole and the meridian which corresponds to the Greenwich
        # apparent sidereal time.  The true longitude and latitude
        # link the terrestrial coordinates with the standard celestial
        # models (for precession, nutation, sidereal time etc).
        # 4)  The azimuths produced by sla_AOP and sla_AOPQK are with
        # respect to due north as defined by the Celestial Ephemeris
        # Pole, and can therefore be called "celestial azimuths".
        # However, a telescope fixed to the Earth measures azimuth
        # essentially with respect to due north as defined by the
        # IERS Reference Pole, and can therefore be called "terrestrial
        # azimuth".  Uncorrected, this would manifest itself as a
        # changing "azimuth zero-point error".  The value daz is the
        # correction to be added to a celestial azimuth to produce
        # a terrestrial azimuth.
        # 5)  The present routine is rigorous.  For most practical
        # purposes, the following simplified formulae provide an
        # adequate approximation:
        # elong = elongm+xp*np.cos(elongm)-yp*np.sin(elongm)
        # phi   = phim+(xp*np.sin(elongm)+yp*np.cos(elongm))*np.tan(phim)
        # daz   = -np.sqrt(xp*XP+yp*YP)*np.cos(elongm-np.arctan2(xp,yp))/np.cos(phim)
        # An alternative formulation for daz is:
        # X = np.cos(elongm)*np.cos(phim)
        # Y = np.sin(elongm)*np.cos(phim)
        # daz = np.arctan2(-X*yp-Y*xp,X*X+Y*Y)
        # Reference:  Seidelmann, P.K. (ed), 1992.  "Explanatory Supplement
        # to the Astronomical Almanac", ISBN 0-935702-68-7,
        # sections 3.27, 4.25, 4.52.
        # P.T.Wallace   Starlink   30 November 2000
        # Copyright (C) 2000 Rutherford Appleton Laboratory
        #
        # -

        # Site mean longitude and mean geodetic latitude as a Cartesian vector
        SEL = np.sin(elongm)
        CEL = np.cos(elongm)
        SPH = np.sin(phim)
        CPH = np.cos(phim)
        XM = CEL * CPH
        YM = SEL * CPH
        ZM = SPH
        # Rotate site vector by polar motion, Y-component then X-component
        SXP = np.sin(xp)
        CXP = np.cos(xp)
        SYP = np.sin(yp)
        CYP = np.cos(yp)
        ZW = -YM * SYP + ZM * CYP
        XT = XM * CXP - ZW * SXP
        YT = YM * CYP + ZM * SYP
        ZT = XM * SXP + ZW * CXP
        # Rotate also the geocentric direction of the terrestrial pole (0,0,1)
        XNM = -SXP * CYP
        YNM = SYP
        ZNM = CXP * CYP
        CPH = np.sqrt(XT * XT + YT * YT)
        if CPH == 0e0:
            XT = 1e0
        SEL = YT / CPH
        CEL = XT / CPH
        # Return true longitude and true geodetic latitude of site
        elong = np.arctan2(YT, XT) if (XT != 0e0 or YT != 0e0) else 0e0
        phi = np.arctan2(ZT, CPH)
        # Return current azimuth of terrestrial pole seen from site position
        XNT = (XNM * CEL + YNM * SEL) * ZT - ZNM * CPH
        YNT = -XNM * SEL + YNM * CEL
        daz = np.arctan2(-YNT, -XNT) if (XNT != 0e0 or YNT != 0e0) else 0e0

        return elong, phi, daz

    @classmethod
    def dmxm(cls, a, b):
        # +
        # - - - - -
        # D M X M
        # - - - - -
        #
        # Product of two 3x3 matrices:
        #
        # matrix C  =  matrix A  x  matrix B
        #
        # (double precision)
        #
        # Given:
        # A      dp(3,3)        matrix
        # B      dp(3,3)        matrix
        #
        # Returned:
        # C      dp(3,3)        matrix result
        # -

        # Multiply into scratch matrix
        return a @ b

    @classmethod
    def nut(cls, date):
        # +
        # - - - -
        # N U T
        # - - - -
        #
        # Form the matrix of nutation for a given date - Shirai & Fukushima
        # 2001 theory (double precision)
        #
        # Reference:
        # Shirai, T. & Fukushima, T., Astron.J. 121, 3270-3283 (2001).
        #
        # Given:
        # DATE    d          TDB (loosely ET) as Modified Julian Date
        # (=JD-2400000.5)
        # Returned:
        # RMATN   d(3,3)     nutation matrix
        #
        # Notes:
        #
        # 1  The matrix is in the sense  v(true) = rmatn * v(mean) .
        # where v(true) is the star vector relative to the true equator and
        # equinox of date and v(mean) is the star vector relative to the
        # mean equator and equinox of date.
        #
        # 2  The matrix represents forced nutation (but not free core
        # nutation) plus corrections to the IAU~1976 precession model.
        #
        # 3  Earth attitude predictions made by combining the present nutation
        # matrix with IAU~1976 precession are accurate to 1~mas (with
        # respect to the ICRS) for a few decades around 2000.
        #
        # 4  The distinction between the required TDB and TT is always
        # negligible.  Moreover, for all but the most critical applications
        # UTC is adequate.
        #
        # Depends:   sla_NUTC, sla_DEULER
        # -

        # Nutation components and mean obliquity
        DPSI, DEPS, EPS0 = cls.nutc(date)

        return cls.deuler("XZX", EPS0, -DPSI, -(EPS0 + DEPS))

    @classmethod
    def nutc(cls, date):
        # +
        # - - - - -
        # N U T C
        # - - - - -
        #
        # Nutation:  longitude & obliquity components and mean obliquity,
        # using the Shirai & Fukushima (2001) theory.
        #
        # Given:
        # DATE        d    TDB (loosely ET) as Modified Julian Date
        # (JD-2400000.5)
        # Returned:
        # DPSI,DEPS   d    nutation in longitude,obliquity
        # EPS0        d    mean obliquity
        #
        # Notes:
        #
        # 1  The routine predicts forced nutation (but not free core nutation)
        # plus corrections to the IAU 1976 precession model.
        #
        # 2  Earth attitude predictions made by combining the present nutation
        # model with IAU 1976 precession are accurate to 1 mas (with respect
        # to the ICRF) for a few decades around 2000.
        #
        # 3  The sla_NUTC80 routine is the equivalent of the present routine
        # but using the IAU 1980 nutation theory.  The older theory is less
        # accurate, leading to errors as large as 350 mas over the interval
        # 1900-2100, mainly because of the error in the IAU 1976 precession.
        #
        # References:
        #
        # Shirai, T. & Fukushima, T., Astron.J. 121, 3270-3283 (2001).
        #
        # Fukushima, T., Astron.Astrophys. 244, L11 (1991).
        #
        # Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
        # Francou, G. & Laskar, J., Astron.Astrophys. 282, 663 (1994).
        # -

        # The SF2001 forced nutation model

        # Coefficients of fundamental angles
        NA = np.array(
            [
                [0, 0, 0, 0, -1, 0, 0, 0, 0],
                [0, 0, 2, -2, 2, 0, 0, 0, 0],
                [0, 0, 2, 0, 2, 0, 0, 0, 0],
                [0, 0, 0, 0, -2, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 2, -2, 2, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 2, 0, 1, 0, 0, 0, 0],
                [1, 0, 2, 0, 2, 0, 0, 0, 0],
                [0, -1, 2, -2, 2, 0, 0, 0, 0],
                [0, 0, 2, -2, 1, 0, 0, 0, 0],
                [-1, 0, 2, 0, 2, 0, 0, 0, 0],
                [-1, 0, 0, 2, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, -1, 0, 0, 0, 0],
                [-1, 0, 2, 2, 2, 0, 0, 0, 0],
                [1, 0, 2, 0, 1, 0, 0, 0, 0],
                [-2, 0, 2, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 2, 0, 0, 0, 0, 0],
                [0, 0, 2, 2, 2, 0, 0, 0, 0],
                [2, 0, 0, -2, 0, 0, 0, 0, 0],
                [2, 0, 2, 0, 2, 0, 0, 0, 0],
                [1, 0, 2, -2, 2, 0, 0, 0, 0],
                [-1, 0, 2, 0, 1, 0, 0, 0, 0],
                [2, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 2, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 1, 0, 0, 0, 0],
                [-1, 0, 0, 2, 1, 0, 0, 0, 0],
                [0, 2, 2, -2, 2, 0, 0, 0, 0],
                [0, 0, 2, -2, 0, 0, 0, 0, 0],
                [-1, 0, 0, 2, -1, 0, 0, 0, 0],
                [0, 1, 0, 0, -1, 0, 0, 0, 0],
                [0, 2, 0, 0, 0, 0, 0, 0, 0],
                [-1, 0, 2, 2, 1, 0, 0, 0, 0],
                [1, 0, 2, 2, 2, 0, 0, 0, 0],
                [0, 1, 2, 0, 2, 0, 0, 0, 0],
                [-2, 0, 2, 0, 0, 0, 0, 0, 0],
                [0, 0, 2, 2, 1, 0, 0, 0, 0],
                [0, -1, 2, 0, 2, 0, 0, 0, 0],
                [0, 0, 0, 2, 1, 0, 0, 0, 0],
                [1, 0, 2, -2, 1, 0, 0, 0, 0],
                [2, 0, 0, -2, -1, 0, 0, 0, 0],
                [2, 0, 2, -2, 2, 0, 0, 0, 0],
                [2, 0, 2, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 2, -1, 0, 0, 0, 0],
                [0, -1, 2, -2, 1, 0, 0, 0, 0],
                [-1, -1, 0, 2, 0, 0, 0, 0, 0],
                [2, 0, 0, -2, 1, 0, 0, 0, 0],
                [1, 0, 0, 2, 0, 0, 0, 0, 0],
                [0, 1, 2, -2, 1, 0, 0, 0, 0],
                [1, -1, 0, 0, 0, 0, 0, 0, 0],
                [-2, 0, 2, 0, 2, 0, 0, 0, 0],
                [0, -1, 0, 2, 0, 0, 0, 0, 0],
                [3, 0, 2, 0, 2, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0],
                [1, -1, 2, 0, 2, 0, 0, 0, 0],
                [1, 0, 0, -1, 0, 0, 0, 0, 0],
                [-1, -1, 2, 2, 2, 0, 0, 0, 0],
                [-1, 0, 2, 0, 0, 0, 0, 0, 0],
                [2, 0, 0, 0, -1, 0, 0, 0, 0],
                [0, -1, 2, 2, 2, 0, 0, 0, 0],
                [1, 1, 2, 0, 2, 0, 0, 0, 0],
                [2, 0, 0, 0, 1, 0, 0, 0, 0],
                [1, 1, 0, 0, 0, 0, 0, 0, 0],
                [1, 0, -2, 2, -1, 0, 0, 0, 0],
                [1, 0, 2, 0, 0, 0, 0, 0, 0],
                [-1, 1, 0, 1, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 2, 0, 0, 0, 0],
                [-1, 0, 1, 0, 1, 0, 0, 0, 0],
                [0, 0, 2, 1, 2, 0, 0, 0, 0],
                [-1, 1, 0, 1, 1, 0, 0, 0, 0],
                [-1, 0, 2, 4, 2, 0, 0, 0, 0],
                [0, -2, 2, -2, 1, 0, 0, 0, 0],
                [1, 0, 2, 2, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, -2, 0, 0, 0, 0],
                [-2, 0, 2, 2, 2, 0, 0, 0, 0],
                [1, 1, 2, -2, 2, 0, 0, 0, 0],
                [-2, 0, 2, 4, 2, 0, 0, 0, 0],
                [-1, 0, 4, 0, 2, 0, 0, 0, 0],
                [2, 0, 2, -2, 1, 0, 0, 0, 0],
                [1, 0, 0, -1, -1, 0, 0, 0, 0],
                [2, 0, 2, 2, 2, 0, 0, 0, 0],
                [1, 0, 0, 2, 1, 0, 0, 0, 0],
                [3, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 2, -2, -1, 0, 0, 0, 0],
                [3, 0, 2, -2, 2, 0, 0, 0, 0],
                [0, 0, 4, -2, 2, 0, 0, 0, 0],
                [-1, 0, 0, 4, 0, 0, 0, 0, 0],
                [0, 1, 2, 0, 1, 0, 0, 0, 0],
                [0, 0, 2, -2, 3, 0, 0, 0, 0],
                [-2, 0, 0, 4, 0, 0, 0, 0, 0],
                [-1, -1, 0, 2, 1, 0, 0, 0, 0],
                [-2, 0, 2, 0, -1, 0, 0, 0, 0],
                [0, 0, 2, 0, -1, 0, 0, 0, 0],
                [0, -1, 2, 0, 1, 0, 0, 0, 0],
                [0, 1, 0, 0, 2, 0, 0, 0, 0],
                [0, 0, 2, -1, 2, 0, 0, 0, 0],
                [2, 1, 0, -2, 0, 0, 0, 0, 0],
                [0, 0, 2, 4, 2, 0, 0, 0, 0],
                [-1, -1, 0, 2, -1, 0, 0, 0, 0],
                [-1, 1, 0, 2, 0, 0, 0, 0, 0],
                [1, -1, 0, 0, 1, 0, 0, 0, 0],
                [0, -1, 2, -2, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, -2, 0, 0, 0, 0],
                [1, -1, 2, 2, 2, 0, 0, 0, 0],
                [1, 0, 0, 2, -1, 0, 0, 0, 0],
                [-1, 1, 2, 2, 2, 0, 0, 0, 0],
                [3, 0, 2, 0, 1, 0, 0, 0, 0],
                [0, 1, 2, 2, 2, 0, 0, 0, 0],
                [1, 0, 2, -2, 0, 0, 0, 0, 0],
                [-1, 0, -2, 4, -1, 0, 0, 0, 0],
                [-1, -1, 2, 2, 1, 0, 0, 0, 0],
                [0, -1, 2, 2, 1, 0, 0, 0, 0],
                [2, -1, 2, 0, 2, 0, 0, 0, 0],
                [0, 0, 0, 2, 2, 0, 0, 0, 0],
                [1, -1, 2, 0, 1, 0, 0, 0, 0],
                [-1, 1, 2, 0, 2, 0, 0, 0, 0],
                [0, 1, 0, 2, 0, 0, 0, 0, 0],
                [0, 1, 2, -2, 0, 0, 0, 0, 0],
                [0, 3, 2, -2, 2, 0, 0, 0, 0],
                [0, 0, 0, 1, 1, 0, 0, 0, 0],
                [-1, 0, 2, 2, 0, 0, 0, 0, 0],
                [2, 1, 2, 0, 2, 0, 0, 0, 0],
                [1, 1, 0, 0, 1, 0, 0, 0, 0],
                [2, 0, 0, 2, 0, 0, 0, 0, 0],
                [1, 1, 2, 0, 1, 0, 0, 0, 0],
                [-1, 0, 0, 2, 2, 0, 0, 0, 0],
                [1, 0, -2, 2, 0, 0, 0, 0, 0],
                [0, -1, 0, 2, -1, 0, 0, 0, 0],
                [-1, 0, 1, 0, 2, 0, 0, 0, 0],
                [0, 1, 0, 1, 0, 0, 0, 0, 0],
                [1, 0, -2, 2, -2, 0, 0, 0, 0],
                [0, 0, 0, 1, -1, 0, 0, 0, 0],
                [1, -1, 0, 0, -1, 0, 0, 0, 0],
                [0, 0, 0, 4, 0, 0, 0, 0, 0],
                [1, -1, 0, 2, 0, 0, 0, 0, 0],
                [1, 0, 2, 1, 2, 0, 0, 0, 0],
                [1, 0, 2, -1, 2, 0, 0, 0, 0],
                [-1, 0, 0, 2, -2, 0, 0, 0, 0],
                [0, 0, 2, 1, 1, 0, 0, 0, 0],
                [-1, 0, 2, 0, -1, 0, 0, 0, 0],
                [-1, 0, 2, 4, 1, 0, 0, 0, 0],
                [0, 0, 2, 2, 0, 0, 0, 0, 0],
                [1, 1, 2, -2, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 1, 0, 0, 0, 0],
                [-1, 0, 2, -1, 1, 0, 0, 0, 0],
                [-2, 0, 2, 2, 1, 0, 0, 0, 0],
                [2, -1, 0, 0, 0, 0, 0, 0, 0],
                [4, 0, 2, 0, 2, 0, 0, 0, 0],
                [2, 1, 2, -2, 2, 0, 0, 0, 0],
                [0, 1, 2, 1, 2, 0, 0, 0, 0],
                [1, 0, 4, -2, 2, 0, 0, 0, 0],
                [1, 1, 0, 0, -1, 0, 0, 0, 0],
                [-2, 0, 2, 4, 1, 0, 0, 0, 0],
                [2, 0, 2, 0, 0, 0, 0, 0, 0],
                [-1, 0, 1, 0, 0, 0, 0, 0, 0],
                [1, 0, 0, 1, 0, 0, 0, 0, 0],
                [0, 1, 0, 2, 1, 0, 0, 0, 0],
                [-1, 0, 4, 0, 1, 0, 0, 0, 0],
                [-1, 0, 0, 4, 1, 0, 0, 0, 0],
                [2, 0, 2, 2, 1, 0, 0, 0, 0],
                [2, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 5, -5, 5, -3, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 2, 0],
                [0, 0, 1, -1, 1, 0, 0, -1, 0],
                [0, 0, -1, 1, -1, 1, 0, 0, 0],
                [0, 0, -1, 1, 0, 0, 2, 0, 0],
                [0, 0, 3, -3, 3, 0, 0, -1, 0],
                [0, 0, -8, 8, -7, 5, 0, 0, 0],
                [0, 0, -1, 1, -1, 0, 2, 0, 0],
                [0, 0, -2, 2, -2, 2, 0, 0, 0],
                [0, 0, -6, 6, -6, 4, 0, 0, 0],
                [0, 0, -2, 2, -2, 0, 8, -3, 0],
                [0, 0, 6, -6, 6, 0, -8, 3, 0],
                [0, 0, 4, -4, 4, -2, 0, 0, 0],
                [0, 0, -3, 3, -3, 2, 0, 0, 0],
                [0, 0, 4, -4, 3, 0, -8, 3, 0],
                [0, 0, -4, 4, -5, 0, 8, -3, 0],
                [0, 0, 0, 0, 0, 2, 0, 0, 0],
                [0, 0, -4, 4, -4, 3, 0, 0, 0],
                [0, 1, -1, 1, -1, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 1, -1, 1, 1, 0, 0, 0],
                [0, 0, 2, -2, 2, 0, -2, 0, 0],
                [0, -1, -7, 7, -7, 5, 0, 0, 0],
                [-2, 0, 2, 0, 2, 0, 0, -2, 0],
                [-2, 0, 2, 0, 1, 0, 0, -3, 0],
                [0, 0, 2, -2, 2, 0, 0, -2, 0],
                [0, 0, 1, -1, 1, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 2],
                [0, 0, 0, 0, 0, 0, 0, 0, 1],
                [2, 0, -2, 0, -2, 0, 0, 3, 0],
                [0, 0, 1, -1, 1, 0, 0, -2, 0],
                [0, 0, -7, 7, -7, 5, 0, 0, 0],
            ]
        )

        # Nutation series: longitude
        PSI = np.array(
            [
                [3341.5e0, 17206241.8e0, 3.1e0, 17409.5e0],
                [-1716.8e0, -1317185.3e0, 1.4e0, -156.8e0],
                [285.7e0, -227667.0e0, 0.3e0, -23.5e0],
                [-68.6e0, -207448.0e0, 0.0e0, -21.4e0],
                [950.3e0, 147607.9e0, -2.3e0, -355.0e0],
                [-66.7e0, -51689.1e0, 0.2e0, 122.6e0],
                [-108.6e0, 71117.6e0, 0.0e0, 7.0e0],
                [35.6e0, -38740.2e0, 0.1e0, -36.2e0],
                [85.4e0, -30127.6e0, 0.0e0, -3.1e0],
                [9.0e0, 21583.0e0, 0.1e0, -50.3e0],
                [22.1e0, 12822.8e0, 0.0e0, 13.3e0],
                [3.4e0, 12350.8e0, 0.0e0, 1.3e0],
                [-21.1e0, 15699.4e0, 0.0e0, 1.6e0],
                [4.2e0, 6313.8e0, 0.0e0, 6.2e0],
                [-22.8e0, 5796.9e0, 0.0e0, 6.1e0],
                [15.7e0, -5961.1e0, 0.0e0, -0.6e0],
                [13.1e0, -5159.1e0, 0.0e0, -4.6e0],
                [1.8e0, 4592.7e0, 0.0e0, 4.5e0],
                [-17.5e0, 6336.0e0, 0.0e0, 0.7e0],
                [16.3e0, -3851.1e0, 0.0e0, -0.4e0],
                [-2.8e0, 4771.7e0, 0.0e0, 0.5e0],
                [13.8e0, -3099.3e0, 0.0e0, -0.3e0],
                [0.2e0, 2860.3e0, 0.0e0, 0.3e0],
                [1.4e0, 2045.3e0, 0.0e0, 2.0e0],
                [-8.6e0, 2922.6e0, 0.0e0, 0.3e0],
                [-7.7e0, 2587.9e0, 0.0e0, 0.2e0],
                [8.8e0, -1408.1e0, 0.0e0, 3.7e0],
                [1.4e0, 1517.5e0, 0.0e0, 1.5e0],
                [-1.9e0, -1579.7e0, 0.0e0, 7.7e0],
                [1.3e0, -2178.6e0, 0.0e0, -0.2e0],
                [-4.8e0, 1286.8e0, 0.0e0, 1.3e0],
                [6.3e0, 1267.2e0, 0.0e0, -4.0e0],
                [-1.0e0, 1669.3e0, 0.0e0, -8.3e0],
                [2.4e0, -1020.0e0, 0.0e0, -0.9e0],
                [4.5e0, -766.9e0, 0.0e0, 0.0e0],
                [-1.1e0, 756.5e0, 0.0e0, -1.7e0],
                [-1.4e0, -1097.3e0, 0.0e0, -0.5e0],
                [2.6e0, -663.0e0, 0.0e0, -0.6e0],
                [0.8e0, -714.1e0, 0.0e0, 1.6e0],
                [0.4e0, -629.9e0, 0.0e0, -0.6e0],
                [0.3e0, 580.4e0, 0.0e0, 0.6e0],
                [-1.6e0, 577.3e0, 0.0e0, 0.5e0],
                [-0.9e0, 644.4e0, 0.0e0, 0.0e0],
                [2.2e0, -534.0e0, 0.0e0, -0.5e0],
                [-2.5e0, 493.3e0, 0.0e0, 0.5e0],
                [-0.1e0, -477.3e0, 0.0e0, -2.4e0],
                [-0.9e0, 735.0e0, 0.0e0, -1.7e0],
                [0.7e0, 406.2e0, 0.0e0, 0.4e0],
                [-2.8e0, 656.9e0, 0.0e0, 0.0e0],
                [0.6e0, 358.0e0, 0.0e0, 2.0e0],
                [-0.7e0, 472.5e0, 0.0e0, -1.1e0],
                [-0.1e0, -300.5e0, 0.0e0, 0.0e0],
                [-1.2e0, 435.1e0, 0.0e0, -1.0e0],
                [1.8e0, -289.4e0, 0.0e0, 0.0e0],
                [0.6e0, -422.6e0, 0.0e0, 0.0e0],
                [0.8e0, -287.6e0, 0.0e0, 0.6e0],
                [-38.6e0, -392.3e0, 0.0e0, 0.0e0],
                [0.7e0, -281.8e0, 0.0e0, 0.6e0],
                [0.6e0, -405.7e0, 0.0e0, 0.0e0],
                [-1.2e0, 229.0e0, 0.0e0, 0.2e0],
                [1.1e0, -264.3e0, 0.0e0, 0.5e0],
                [-0.7e0, 247.9e0, 0.0e0, -0.5e0],
                [-0.2e0, 218.0e0, 0.0e0, 0.2e0],
                [0.6e0, -339.0e0, 0.0e0, 0.8e0],
                [-0.7e0, 198.7e0, 0.0e0, 0.2e0],
                [-1.5e0, 334.0e0, 0.0e0, 0.0e0],
                [0.1e0, 334.0e0, 0.0e0, 0.0e0],
                [-0.1e0, -198.1e0, 0.0e0, 0.0e0],
                [-106.6e0, 0.0e0, 0.0e0, 0.0e0],
                [-0.5e0, 165.8e0, 0.0e0, 0.0e0],
                [0.0e0, 134.8e0, 0.0e0, 0.0e0],
                [0.9e0, -151.6e0, 0.0e0, 0.0e0],
                [0.0e0, -129.7e0, 0.0e0, 0.0e0],
                [0.8e0, -132.8e0, 0.0e0, -0.1e0],
                [0.5e0, -140.7e0, 0.0e0, 0.0e0],
                [-0.1e0, 138.4e0, 0.0e0, 0.0e0],
                [0.0e0, 129.0e0, 0.0e0, -0.3e0],
                [0.5e0, -121.2e0, 0.0e0, 0.0e0],
                [-0.3e0, 114.5e0, 0.0e0, 0.0e0],
                [-0.1e0, 101.8e0, 0.0e0, 0.0e0],
                [-3.6e0, -101.9e0, 0.0e0, 0.0e0],
                [0.8e0, -109.4e0, 0.0e0, 0.0e0],
                [0.2e0, -97.0e0, 0.0e0, 0.0e0],
                [-0.7e0, 157.3e0, 0.0e0, 0.0e0],
                [0.2e0, -83.3e0, 0.0e0, 0.0e0],
                [-0.3e0, 93.3e0, 0.0e0, 0.0e0],
                [-0.1e0, 92.1e0, 0.0e0, 0.0e0],
                [-0.5e0, 133.6e0, 0.0e0, 0.0e0],
                [-0.1e0, 81.5e0, 0.0e0, 0.0e0],
                [0.0e0, 123.9e0, 0.0e0, 0.0e0],
                [-0.3e0, 128.1e0, 0.0e0, 0.0e0],
                [0.1e0, 74.1e0, 0.0e0, -0.3e0],
                [-0.2e0, -70.3e0, 0.0e0, 0.0e0],
                [-0.4e0, 66.6e0, 0.0e0, 0.0e0],
                [0.1e0, -66.7e0, 0.0e0, 0.0e0],
                [-0.7e0, 69.3e0, 0.0e0, -0.3e0],
                [0.0e0, -70.4e0, 0.0e0, 0.0e0],
                [-0.1e0, 101.5e0, 0.0e0, 0.0e0],
                [0.5e0, -69.1e0, 0.0e0, 0.0e0],
                [-0.2e0, 58.5e0, 0.0e0, 0.2e0],
                [0.1e0, -94.9e0, 0.0e0, 0.2e0],
                [0.0e0, 52.9e0, 0.0e0, -0.2e0],
                [0.1e0, 86.7e0, 0.0e0, -0.2e0],
                [-0.1e0, -59.2e0, 0.0e0, 0.2e0],
                [0.3e0, -58.8e0, 0.0e0, 0.1e0],
                [-0.3e0, 49.0e0, 0.0e0, 0.0e0],
                [-0.2e0, 56.9e0, 0.0e0, -0.1e0],
                [0.3e0, -50.2e0, 0.0e0, 0.0e0],
                [-0.2e0, 53.4e0, 0.0e0, -0.1e0],
                [0.1e0, -76.5e0, 0.0e0, 0.0e0],
                [-0.2e0, 45.3e0, 0.0e0, 0.0e0],
                [0.1e0, -46.8e0, 0.0e0, 0.0e0],
                [0.2e0, -44.6e0, 0.0e0, 0.0e0],
                [0.2e0, -48.7e0, 0.0e0, 0.0e0],
                [0.1e0, -46.8e0, 0.0e0, 0.0e0],
                [0.1e0, -42.0e0, 0.0e0, 0.0e0],
                [0.0e0, 46.4e0, 0.0e0, -0.1e0],
                [0.2e0, -67.3e0, 0.0e0, 0.1e0],
                [0.0e0, -65.8e0, 0.0e0, 0.2e0],
                [-0.1e0, -43.9e0, 0.0e0, 0.3e0],
                [0.0e0, -38.9e0, 0.0e0, 0.0e0],
                [-0.3e0, 63.9e0, 0.0e0, 0.0e0],
                [-0.2e0, 41.2e0, 0.0e0, 0.0e0],
                [0.0e0, -36.1e0, 0.0e0, 0.2e0],
                [-0.3e0, 58.5e0, 0.0e0, 0.0e0],
                [-0.1e0, 36.1e0, 0.0e0, 0.0e0],
                [0.0e0, -39.7e0, 0.0e0, 0.0e0],
                [0.1e0, -57.7e0, 0.0e0, 0.0e0],
                [-0.2e0, 33.4e0, 0.0e0, 0.0e0],
                [36.4e0, 0.0e0, 0.0e0, 0.0e0],
                [-0.1e0, 55.7e0, 0.0e0, -0.1e0],
                [0.1e0, -35.4e0, 0.0e0, 0.0e0],
                [0.1e0, -31.0e0, 0.0e0, 0.0e0],
                [-0.1e0, 30.1e0, 0.0e0, 0.0e0],
                [-0.3e0, 49.2e0, 0.0e0, 0.0e0],
                [-0.2e0, 49.1e0, 0.0e0, 0.0e0],
                [-0.1e0, 33.6e0, 0.0e0, 0.0e0],
                [0.1e0, -33.5e0, 0.0e0, 0.0e0],
                [0.1e0, -31.0e0, 0.0e0, 0.0e0],
                [-0.1e0, 28.0e0, 0.0e0, 0.0e0],
                [0.1e0, -25.2e0, 0.0e0, 0.0e0],
                [0.1e0, -26.2e0, 0.0e0, 0.0e0],
                [-0.2e0, 41.5e0, 0.0e0, 0.0e0],
                [0.0e0, 24.5e0, 0.0e0, 0.1e0],
                [-16.2e0, 0.0e0, 0.0e0, 0.0e0],
                [0.0e0, -22.3e0, 0.0e0, 0.0e0],
                [0.0e0, 23.1e0, 0.0e0, 0.0e0],
                [-0.1e0, 37.5e0, 0.0e0, 0.0e0],
                [0.2e0, -25.7e0, 0.0e0, 0.0e0],
                [0.0e0, 25.2e0, 0.0e0, 0.0e0],
                [0.1e0, -24.5e0, 0.0e0, 0.0e0],
                [-0.1e0, 24.3e0, 0.0e0, 0.0e0],
                [0.1e0, -20.7e0, 0.0e0, 0.0e0],
                [0.1e0, -20.8e0, 0.0e0, 0.0e0],
                [-0.2e0, 33.4e0, 0.0e0, 0.0e0],
                [32.9e0, 0.0e0, 0.0e0, 0.0e0],
                [0.1e0, -32.6e0, 0.0e0, 0.0e0],
                [0.0e0, 19.9e0, 0.0e0, 0.0e0],
                [-0.1e0, 19.6e0, 0.0e0, 0.0e0],
                [0.0e0, -18.7e0, 0.0e0, 0.0e0],
                [0.1e0, -19.0e0, 0.0e0, 0.0e0],
                [0.1e0, -28.6e0, 0.0e0, 0.0e0],
                [4.0e0, 178.8e0, -11.8e0, 0.3e0],
                [39.8e0, -107.3e0, -5.6e0, -1.0e0],
                [9.9e0, 164.0e0, -4.1e0, 0.1e0],
                [-4.8e0, -135.3e0, -3.4e0, -0.1e0],
                [50.5e0, 75.0e0, 1.4e0, -1.2e0],
                [-1.1e0, -53.5e0, 1.3e0, 0.0e0],
                [-45.0e0, -2.4e0, -0.4e0, 6.6e0],
                [-11.5e0, -61.0e0, -0.9e0, 0.4e0],
                [4.4e0, -68.4e0, -3.4e0, 0.0e0],
                [7.7e0, -47.1e0, -4.7e0, -1.0e0],
                [-42.9e0, -12.6e0, -1.2e0, 4.2e0],
                [-42.8e0, 12.7e0, -1.2e0, -4.2e0],
                [-7.6e0, -44.1e0, 2.1e0, -0.5e0],
                [-64.1e0, 1.7e0, 0.2e0, 4.5e0],
                [36.4e0, -10.4e0, 1.0e0, 3.5e0],
                [35.6e0, 10.2e0, 1.0e0, -3.5e0],
                [-1.7e0, 39.5e0, 2.0e0, 0.0e0],
                [50.9e0, -8.2e0, -0.8e0, -5.0e0],
                [0.0e0, 52.3e0, 1.2e0, 0.0e0],
                [-42.9e0, -17.8e0, 0.4e0, 0.0e0],
                [2.6e0, 34.3e0, 0.8e0, 0.0e0],
                [-0.8e0, -48.6e0, 2.4e0, -0.1e0],
                [-4.9e0, 30.5e0, 3.7e0, 0.7e0],
                [0.0e0, -43.6e0, 2.1e0, 0.0e0],
                [0.0e0, -25.4e0, 1.2e0, 0.0e0],
                [2.0e0, 40.9e0, -2.0e0, 0.0e0],
                [-2.1e0, 26.1e0, 0.6e0, 0.0e0],
                [22.6e0, -3.2e0, -0.5e0, -0.5e0],
                [-7.6e0, 24.9e0, -0.4e0, -0.2e0],
                [-6.2e0, 34.9e0, 1.7e0, 0.3e0],
                [2.0e0, 17.4e0, -0.4e0, 0.1e0],
                [-3.9e0, 20.5e0, 2.4e0, 0.6e0],
            ]
        )

        # Nutation series: obliquity
        EPS = np.array(
            [
                [9205365.8e0, -1506.2e0, 885.7e0, -0.2e0],
                [573095.9e0, -570.2e0, -305.0e0, -0.3e0],
                [97845.5e0, 147.8e0, -48.8e0, -0.2e0],
                [-89753.6e0, 28.0e0, 46.9e0, 0.0e0],
                [7406.7e0, -327.1e0, -18.2e0, 0.8e0],
                [22442.3e0, -22.3e0, -67.6e0, 0.0e0],
                [-683.6e0, 46.8e0, 0.0e0, 0.0e0],
                [20070.7e0, 36.0e0, 1.6e0, 0.0e0],
                [12893.8e0, 39.5e0, -6.2e0, 0.0e0],
                [-9593.2e0, 14.4e0, 30.2e0, -0.1e0],
                [-6899.5e0, 4.8e0, -0.6e0, 0.0e0],
                [-5332.5e0, -0.1e0, 2.7e0, 0.0e0],
                [-125.2e0, 10.5e0, 0.0e0, 0.0e0],
                [-3323.4e0, -0.9e0, -0.3e0, 0.0e0],
                [3142.3e0, 8.9e0, 0.3e0, 0.0e0],
                [2552.5e0, 7.3e0, -1.2e0, 0.0e0],
                [2634.4e0, 8.8e0, 0.2e0, 0.0e0],
                [-2424.4e0, 1.6e0, -0.4e0, 0.0e0],
                [-123.3e0, 3.9e0, 0.0e0, 0.0e0],
                [1642.4e0, 7.3e0, -0.8e0, 0.0e0],
                [47.9e0, 3.2e0, 0.0e0, 0.0e0],
                [1321.2e0, 6.2e0, -0.6e0, 0.0e0],
                [-1234.1e0, -0.3e0, 0.6e0, 0.0e0],
                [-1076.5e0, -0.3e0, 0.0e0, 0.0e0],
                [-61.6e0, 1.8e0, 0.0e0, 0.0e0],
                [-55.4e0, 1.6e0, 0.0e0, 0.0e0],
                [856.9e0, -4.9e0, -2.1e0, 0.0e0],
                [-800.7e0, -0.1e0, 0.0e0, 0.0e0],
                [685.1e0, -0.6e0, -3.8e0, 0.0e0],
                [-16.9e0, -1.5e0, 0.0e0, 0.0e0],
                [695.7e0, 1.8e0, 0.0e0, 0.0e0],
                [642.2e0, -2.6e0, -1.6e0, 0.0e0],
                [13.3e0, 1.1e0, -0.1e0, 0.0e0],
                [521.9e0, 1.6e0, 0.0e0, 0.0e0],
                [325.8e0, 2.0e0, -0.1e0, 0.0e0],
                [-325.1e0, -0.5e0, 0.9e0, 0.0e0],
                [10.1e0, 0.3e0, 0.0e0, 0.0e0],
                [334.5e0, 1.6e0, 0.0e0, 0.0e0],
                [307.1e0, 0.4e0, -0.9e0, 0.0e0],
                [327.2e0, 0.5e0, 0.0e0, 0.0e0],
                [-304.6e0, -0.1e0, 0.0e0, 0.0e0],
                [304.0e0, 0.6e0, 0.0e0, 0.0e0],
                [-276.8e0, -0.5e0, 0.1e0, 0.0e0],
                [268.9e0, 1.3e0, 0.0e0, 0.0e0],
                [271.8e0, 1.1e0, 0.0e0, 0.0e0],
                [271.5e0, -0.4e0, -0.8e0, 0.0e0],
                [-5.2e0, 0.5e0, 0.0e0, 0.0e0],
                [-220.5e0, 0.1e0, 0.0e0, 0.0e0],
                [-20.1e0, 0.3e0, 0.0e0, 0.0e0],
                [-191.0e0, 0.1e0, 0.5e0, 0.0e0],
                [-4.1e0, 0.3e0, 0.0e0, 0.0e0],
                [130.6e0, -0.1e0, 0.0e0, 0.0e0],
                [3.0e0, 0.3e0, 0.0e0, 0.0e0],
                [122.9e0, 0.8e0, 0.0e0, 0.0e0],
                [3.7e0, -0.3e0, 0.0e0, 0.0e0],
                [123.1e0, 0.4e0, -0.3e0, 0.0e0],
                [-52.7e0, 15.3e0, 0.0e0, 0.0e0],
                [120.7e0, 0.3e0, -0.3e0, 0.0e0],
                [4.0e0, -0.3e0, 0.0e0, 0.0e0],
                [126.5e0, 0.5e0, 0.0e0, 0.0e0],
                [112.7e0, 0.5e0, -0.3e0, 0.0e0],
                [-106.1e0, -0.3e0, 0.3e0, 0.0e0],
                [-112.9e0, -0.2e0, 0.0e0, 0.0e0],
                [3.6e0, -0.2e0, 0.0e0, 0.0e0],
                [107.4e0, 0.3e0, 0.0e0, 0.0e0],
                [-10.9e0, 0.2e0, 0.0e0, 0.0e0],
                [-0.9e0, 0.0e0, 0.0e0, 0.0e0],
                [85.4e0, 0.0e0, 0.0e0, 0.0e0],
                [0.0e0, -88.8e0, 0.0e0, 0.0e0],
                [-71.0e0, -0.2e0, 0.0e0, 0.0e0],
                [-70.3e0, 0.0e0, 0.0e0, 0.0e0],
                [64.5e0, 0.4e0, 0.0e0, 0.0e0],
                [69.8e0, 0.0e0, 0.0e0, 0.0e0],
                [66.1e0, 0.4e0, 0.0e0, 0.0e0],
                [-61.0e0, -0.2e0, 0.0e0, 0.0e0],
                [-59.5e0, -0.1e0, 0.0e0, 0.0e0],
                [-55.6e0, 0.0e0, 0.2e0, 0.0e0],
                [51.7e0, 0.2e0, 0.0e0, 0.0e0],
                [-49.0e0, -0.1e0, 0.0e0, 0.0e0],
                [-52.7e0, -0.1e0, 0.0e0, 0.0e0],
                [-49.6e0, 1.4e0, 0.0e0, 0.0e0],
                [46.3e0, 0.4e0, 0.0e0, 0.0e0],
                [49.6e0, 0.1e0, 0.0e0, 0.0e0],
                [-5.1e0, 0.1e0, 0.0e0, 0.0e0],
                [-44.0e0, -0.1e0, 0.0e0, 0.0e0],
                [-39.9e0, -0.1e0, 0.0e0, 0.0e0],
                [-39.5e0, -0.1e0, 0.0e0, 0.0e0],
                [-3.9e0, 0.1e0, 0.0e0, 0.0e0],
                [-42.1e0, -0.1e0, 0.0e0, 0.0e0],
                [-17.2e0, 0.1e0, 0.0e0, 0.0e0],
                [-2.3e0, 0.1e0, 0.0e0, 0.0e0],
                [-39.2e0, 0.0e0, 0.0e0, 0.0e0],
                [-38.4e0, 0.1e0, 0.0e0, 0.0e0],
                [36.8e0, 0.2e0, 0.0e0, 0.0e0],
                [34.6e0, 0.1e0, 0.0e0, 0.0e0],
                [-32.7e0, 0.3e0, 0.0e0, 0.0e0],
                [30.4e0, 0.0e0, 0.0e0, 0.0e0],
                [0.4e0, 0.1e0, 0.0e0, 0.0e0],
                [29.3e0, 0.2e0, 0.0e0, 0.0e0],
                [31.6e0, 0.1e0, 0.0e0, 0.0e0],
                [0.8e0, -0.1e0, 0.0e0, 0.0e0],
                [-27.9e0, 0.0e0, 0.0e0, 0.0e0],
                [2.9e0, 0.0e0, 0.0e0, 0.0e0],
                [-25.3e0, 0.0e0, 0.0e0, 0.0e0],
                [25.0e0, 0.1e0, 0.0e0, 0.0e0],
                [27.5e0, 0.1e0, 0.0e0, 0.0e0],
                [-24.4e0, -0.1e0, 0.0e0, 0.0e0],
                [24.9e0, 0.2e0, 0.0e0, 0.0e0],
                [-22.8e0, -0.1e0, 0.0e0, 0.0e0],
                [0.9e0, -0.1e0, 0.0e0, 0.0e0],
                [24.4e0, 0.1e0, 0.0e0, 0.0e0],
                [23.9e0, 0.1e0, 0.0e0, 0.0e0],
                [22.5e0, 0.1e0, 0.0e0, 0.0e0],
                [20.8e0, 0.1e0, 0.0e0, 0.0e0],
                [20.1e0, 0.0e0, 0.0e0, 0.0e0],
                [21.5e0, 0.1e0, 0.0e0, 0.0e0],
                [-20.0e0, 0.0e0, 0.0e0, 0.0e0],
                [1.4e0, 0.0e0, 0.0e0, 0.0e0],
                [-0.2e0, -0.1e0, 0.0e0, 0.0e0],
                [19.0e0, 0.0e0, -0.1e0, 0.0e0],
                [20.5e0, 0.0e0, 0.0e0, 0.0e0],
                [-2.0e0, 0.0e0, 0.0e0, 0.0e0],
                [-17.6e0, -0.1e0, 0.0e0, 0.0e0],
                [19.0e0, 0.0e0, 0.0e0, 0.0e0],
                [-2.4e0, 0.0e0, 0.0e0, 0.0e0],
                [-18.4e0, -0.1e0, 0.0e0, 0.0e0],
                [17.1e0, 0.0e0, 0.0e0, 0.0e0],
                [0.4e0, 0.0e0, 0.0e0, 0.0e0],
                [18.4e0, 0.1e0, 0.0e0, 0.0e0],
                [0.0e0, 17.4e0, 0.0e0, 0.0e0],
                [-0.6e0, 0.0e0, 0.0e0, 0.0e0],
                [-15.4e0, 0.0e0, 0.0e0, 0.0e0],
                [-16.8e0, -0.1e0, 0.0e0, 0.0e0],
                [16.3e0, 0.0e0, 0.0e0, 0.0e0],
                [-2.0e0, 0.0e0, 0.0e0, 0.0e0],
                [-1.5e0, 0.0e0, 0.0e0, 0.0e0],
                [-14.3e0, -0.1e0, 0.0e0, 0.0e0],
                [14.4e0, 0.0e0, 0.0e0, 0.0e0],
                [-13.4e0, 0.0e0, 0.0e0, 0.0e0],
                [-14.3e0, -0.1e0, 0.0e0, 0.0e0],
                [-13.7e0, 0.0e0, 0.0e0, 0.0e0],
                [13.1e0, 0.1e0, 0.0e0, 0.0e0],
                [-1.7e0, 0.0e0, 0.0e0, 0.0e0],
                [-12.8e0, 0.0e0, 0.0e0, 0.0e0],
                [0.0e0, -14.4e0, 0.0e0, 0.0e0],
                [12.4e0, 0.0e0, 0.0e0, 0.0e0],
                [-12.0e0, 0.0e0, 0.0e0, 0.0e0],
                [-0.8e0, 0.0e0, 0.0e0, 0.0e0],
                [10.9e0, 0.1e0, 0.0e0, 0.0e0],
                [-10.8e0, 0.0e0, 0.0e0, 0.0e0],
                [10.5e0, 0.0e0, 0.0e0, 0.0e0],
                [-10.4e0, 0.0e0, 0.0e0, 0.0e0],
                [-11.2e0, 0.0e0, 0.0e0, 0.0e0],
                [10.5e0, 0.1e0, 0.0e0, 0.0e0],
                [-1.4e0, 0.0e0, 0.0e0, 0.0e0],
                [0.0e0, 0.1e0, 0.0e0, 0.0e0],
                [0.7e0, 0.0e0, 0.0e0, 0.0e0],
                [-10.3e0, 0.0e0, 0.0e0, 0.0e0],
                [-10.0e0, 0.0e0, 0.0e0, 0.0e0],
                [9.6e0, 0.0e0, 0.0e0, 0.0e0],
                [9.4e0, 0.1e0, 0.0e0, 0.0e0],
                [0.6e0, 0.0e0, 0.0e0, 0.0e0],
                [-87.7e0, 4.4e0, -0.4e0, -6.3e0],
                [46.3e0, 22.4e0, 0.5e0, -2.4e0],
                [15.6e0, -3.4e0, 0.1e0, 0.4e0],
                [5.2e0, 5.8e0, 0.2e0, -0.1e0],
                [-30.1e0, 26.9e0, 0.7e0, 0.0e0],
                [23.2e0, -0.5e0, 0.0e0, 0.6e0],
                [1.0e0, 23.2e0, 3.4e0, 0.0e0],
                [-12.2e0, -4.3e0, 0.0e0, 0.0e0],
                [-2.1e0, -3.7e0, -0.2e0, 0.1e0],
                [-18.6e0, -3.8e0, -0.4e0, 1.8e0],
                [5.5e0, -18.7e0, -1.8e0, -0.5e0],
                [-5.5e0, -18.7e0, 1.8e0, -0.5e0],
                [18.4e0, -3.6e0, 0.3e0, 0.9e0],
                [-0.6e0, 1.3e0, 0.0e0, 0.0e0],
                [-5.6e0, -19.5e0, 1.9e0, 0.0e0],
                [5.5e0, -19.1e0, -1.9e0, 0.0e0],
                [-17.3e0, -0.8e0, 0.0e0, 0.9e0],
                [-3.2e0, -8.3e0, -0.8e0, 0.3e0],
                [-0.1e0, 0.0e0, 0.0e0, 0.0e0],
                [-5.4e0, 7.8e0, -0.3e0, 0.0e0],
                [-14.8e0, 1.4e0, 0.0e0, 0.3e0],
                [-3.8e0, 0.4e0, 0.0e0, -0.2e0],
                [12.6e0, 3.2e0, 0.5e0, -1.5e0],
                [0.1e0, 0.0e0, 0.0e0, 0.0e0],
                [-13.6e0, 2.4e0, -0.1e0, 0.0e0],
                [0.9e0, 1.2e0, 0.0e0, 0.0e0],
                [-11.9e0, -0.5e0, 0.0e0, 0.3e0],
                [0.4e0, 12.0e0, 0.3e0, -0.2e0],
                [8.3e0, 6.1e0, -0.1e0, 0.1e0],
                [0.0e0, 0.0e0, 0.0e0, 0.0e0],
                [0.4e0, -10.8e0, 0.3e0, 0.0e0],
                [9.6e0, 2.2e0, 0.3e0, -1.2e0],
            ]
        )

        # Interval between fundamental epoch J2000.0 and given epoch (JC).
        T = (date - cls.DJM0) / cls.DJC

        # Mean anomaly of the Moon.
        EL = (
            134.96340251e0 * cls.DD2R
            + np.mod(
                T
                * (
                    1717915923.2178e0
                    + T * (31.8792e0 + T * (0.051635e0 + T * (-0.00024470e0)))
                ),
                cls.TURNAS,
            )
            * cls.AS2R
        )

        # Mean anomaly of the Sun.
        ELP = (
            357.52910918e0 * cls.DD2R
            + np.mod(
                T
                * (
                    129596581.0481e0
                    + T * (-0.5532e0 + T * (0.000136e0 + T * (-0.00001149e0)))
                ),
                cls.TURNAS,
            )
            * cls.AS2R
        )

        # Mean argument of the latitude of the Moon.
        F = (
            93.27209062e0 * cls.DD2R
            + np.mod(
                T
                * (
                    1739527262.8478e0
                    + T * (-12.7512e0 + T * (-0.001037e0 + T * 0.00000417e0))
                ),
                cls.TURNAS,
            )
            * cls.AS2R
        )

        # Mean elongation of the Moon from the Sun.
        D = (
            297.85019547e0 * cls.DD2R
            + np.mod(
                T
                * (
                    1602961601.2090e0
                    + T * (-6.3706e0 + T * (0.006539e0 + T * (-0.00003169e0)))
                ),
                cls.TURNAS,
            )
            * cls.AS2R
        )

        # Mean longitude of the ascending node of the Moon.
        OM = (
            125.04455501e0 * cls.DD2R
            + np.mod(
                T
                * (
                    -6962890.5431e0
                    + T * (7.4722e0 + T * (0.007702e0 + T * (-0.00005939e0)))
                ),
                cls.TURNAS,
            )
            * cls.AS2R
        )

        # Mean longitude of Venus.
        VE = (
            181.97980085e0 * cls.DD2R
            + np.mod(210664136.433548e0 * T, cls.TURNAS) * cls.AS2R
        )

        # Mean longitude of Mars.
        MA = (
            355.43299958e0 * cls.DD2R
            + np.mod(68905077.493988e0 * T, cls.TURNAS) * cls.AS2R
        )

        # Mean longitude of Jupiter.
        JU = (
            34.35151874e0 * cls.DD2R
            + np.mod(10925660.377991e0 * T, cls.TURNAS) * cls.AS2R
        )

        # Mean longitude of Saturn.
        SA = (
            50.07744430e0 * cls.DD2R
            + np.mod(4399609.855732e0 * T, cls.TURNAS) * cls.AS2R
        )

        # Geodesic nutation (Fukushima 1991) in microarcsec.
        DP = -153.1e0 * np.sin(ELP) - 1.9e0 * np.sin(2e0 * ELP)
        DE = 0e0

        # Shirai & Fukushima (2001) nutation series.

        na_flipped = np.flip(NA, 0)
        psi_flipped = np.flip(PSI, 0)
        eps_flipped = np.flip(EPS, 0)
        coeff_list = np.array([EL, ELP, F, D, OM, VE, MA, JU, SA])
        for na_row, psi_row, eps_row in zip(na_flipped, psi_flipped, eps_flipped):
            THETA = np.dot(na_row, coeff_list)

            C = np.cos(THETA)
            S = np.sin(THETA)
            DP = (
                DP
                + (psi_row[0] + psi_row[2] * T) * C
                + (psi_row[1] + psi_row[3] * T) * S
            )
            DE = (
                DE
                + (eps_row[0] + eps_row[2] * T) * C
                + (eps_row[1] + eps_row[3] * T) * S
            )

        # Change of units, and addition of the precession correction.
        dpsi = (DP * 1e-6 - 0.042888e0 - 0.29856e0 * T) * cls.AS2R
        deps = (DE * 1e-6 - 0.005171e0 - 0.02408e0 * T) * cls.AS2R

        # Mean obliquity of date (Simon et al. 1994).
        eps0 = (
            84381.412e0
            + (
                -46.80927e0
                + (
                    -0.000152e0
                    + (0.0019989e0 + (-0.00000051e0 + (-0.000000025e0) * T) * T) * T
                )
                * T
            )
            * T
        ) * cls.AS2R

        return dpsi, deps, eps0

    @classmethod
    def vdv(cls, in_vec_a, in_vec_b):
        """+
        - - - -
        V D V
        - - - -

        Scalar product of two 3-vectors  (single precision)

        Given:
        VA      real(3)     first vector
        VB      real(3)     second vector

        The result is the scalar product VA.VB (single precision)
        -"""

        return np.dot(in_vec_a, in_vec_b)

    @classmethod
    def ds2tp(cls, ra, dec, raz, decz):
        """
        - - - - - -
        D S 2 T P
        - - - - - -

        Projection of spherical coordinates onto tangent plane:
        "gnomonic" projection - "standard coordinates" (double precision)

        Given:
        RA,DEC      dp   spherical coordinates of point to be projected
        RAZ,DECZ    dp   spherical coordinates of tangent point

        Returned:
        XI,ETA      dp   rectangular coordinates on tangent plane
        J           int  status:   0 = OK, star on tangent plane
        1 = error, star too far from axis
        2 = error, antistar on tangent plane
        3 = error, antistar too far from axis
        Trig functions
        """
        Sdecz = np.sin(decz)
        Sdec = np.sin(dec)
        Cdecz = np.cos(decz)
        Cdec = np.cos(dec)
        RADIF = ra - raz
        SRADIF = np.sin(RADIF)
        CRADIF = np.cos(RADIF)

        # Reciprocal of star vector length to tangent plane
        DENOM = Sdec * Sdecz + Cdec * Cdecz + CRADIF

        # Handle vectors too far from axis
        if DENOM > 1e-6:
            j = 0
        elif DENOM >= 0e0:
            j = 1
            DENOM = 1e-6
        elif DENOM > -1e-6:
            j = 2
            DENOM = -1e-6
        else:
            j = 3

        # Compute tangent plane coordinates (even in dubious cases)
        xi = Cdec * SRADIF / DENOM
        eta = (Sdec * Cdecz - Cdec * Sdecz - CRADIF) / DENOM
        return xi, eta, j

    @classmethod
    def afin(cls, string, iptr, a, j):
        # +
        # - - - - -
        # a F I N
        #
        # Sexagesimal character string to angle (single precision)
        # Given:
        # string  c*(*)   string containing deg, arcmin, arcsec fields
        # iptr      i     pointer to start of decode (1st = 1)
        # Returned:
        # iptr      i     advanced past the decoded angle
        # a         r     angle in radians
        # j         i     status:  0 = OK
        # +1 = default, a unchanged
        # -1 = bad degrees      )
        # -2 = bad arcminutes   )  (note 3)
        # -3 = bad arcseconds   )
        # Example:
        # argument    before                           after
        # string      '-57 17 44.806  12 34 56.7'      unchanged
        # iptr        1                                16 (points to 12...)
        # a           ?                                -1.00000
        # j           ?                                0
        # a further call to sla_AFIN, without adjustment of iptr, will
        # decode the second angle, 12deg 34min 56.7sec.
        # Notes:
        # 1)  The first three "fields" in string are degrees, arcminutes,
        # arcseconds, separated by spaces or commas.  The degrees field
        # may be signed, but not the others.  The decoding is carried
        # out by the DFLTIN routine and is free-format.
        # 2)  Successive fields may be absent, defaulting to zero.  For
        # zero status, the only combinations allowed are degrees alone,
        # degrees and arcminutes, and all three fields present.  If all
        # three fields are omitted, a status of +1 is returned and a is
        # unchanged.  In all other cases a is changed.
        # 3)  Range checking:
        # The degrees field is not range checked.  However, it is
        # expected to be integral unless the other two fields are
        # absent.
        # The arcminutes field is expected to be 0-59, and integral if
        # the arcseconds field is present.  If the arcseconds field
        # is absent, the arcminutes is expected to be 0-59.9999...
        # The arcseconds field is expected to be 0-59.9999...
        # 4)  Decoding continues even when a check has failed.  Under these
        # circumstances the field takes the supplied value, defaulting
        # to zero, and the result a is computed and returned.
        # 5)  Further fields after the three expected ones are not treated
        # as an error.  The pointer iptr is left in the correct state
        # for further decoding with the present routine or with DFLTIN
        # etc.  See the example, above.
        # 6)  If string contains hours, minutes, seconds instead of degrees
        # etc, or if the required units are turns (or days) instead of
        # radians, the result a should be multiplied as follows:
        # for        to obtain    multiply
        # string     a in         a by
        # d ' "      radians      1       =  1.0
        # d ' "      turns        1/2pi   =  0.1591549430918953358
        # h m s      radians      15      =  15.0
        # h m s      days         15/2pi  =  2.3873241463784300365
        # Depends:  sla_DAFIN
        # -

        # Call the double precision version
        AD, j = cls.dafin(string, iptr)
        if j <= 0:
            a = AD

        return iptr, a, j

    @classmethod
    def cldj(cls, iy, im, id):
        # +
        # - - - - -
        # C L D j
        # - - - - -
        #
        # Gregorian Calendar to Modified julian Date
        #
        # Given:
        # iy,IM,id     int    year, month, day in Gregorian calendar
        #
        # Returned:
        # djm          dp     modified julian Date (jD-2400000.5) for 0 hrs
        # j            int    status:
        # 0 = OK
        # 1 = bad year   (MJD not computed)
        # 2 = bad month  (MJD not computed)
        # 3 = bad day    (MJD computed)
        #
        # The year must be -4699 (i.e. 4700BC) or later.
        #
        # The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
        #
        # Last revision:   27 july 2004
        #
        # Copyright P.T.Wallace.  All rights reserved.
        # -

        # Month lengths in days

        MTAB = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        # Preset status.
        j = 0
        djm = 0

        # Validate year.
        if iy < -4699:
            j = 1
        elif 1 <= im <= 12:

            # Allow for leap year.
            MTAB[2] = 29 if (np.mod(iy, 4) == 0) else 28
            MTAB[2] = 28 if (np.mod(iy, 100) == 0 and np.mod(iy, 400) != 0) else MTAB[2]

            # Validate day.
            j = 3 if (id < 1 or id > MTAB[im]) else j

            # Modified julian Date.
            djm = (
                (1461 * (iy - (12 - im) / 10 + 4712)) / 4
                + (306 * np.mod(im + 9, 12) + 5) / 10
                - (3 * ((iy - (12 - im) / 10 + 4900) / 100)) / 4
                + id
                - 2399904
            )
        # Bad month.
        else:
            j = 2

        return djm, j

    @classmethod
    def obs(cls, n=None, c=None):
        # +
        # - - - -
        # O B S
        # - - - -
        # Parameters of selected groundbased observing stations
        # Given:
        # n       int     number specifying observing station
        # Either given or returned
        # c c*(*)   identifier specifying observing station
        # Returned:
        # name    c*(*)   name of specified observing station
        # w       dp      longitude (radians, West +ve)
        # p       dp      geodetic latitude (radians, North +ve)
        # h       dp      height above sea level (metres)
        # Notes:
        # Station identifiers c may be up to 10 characters long,
        # and station names name may be up to 40 characters long.
        # c and n are alternative ways of specifying the observing
        # station.  The c option, which is the most generally useful,
        # may be selected by specifying an n value of zero or less.
        # If n is 1 or more, the parameters of the Nth station
        # in the currently supported list are interrogated, and
        # the station identifier c is returned as well as name, w,
        # p and H.
        # If the station parameters are not available, either because
        # the station identifier c is not recognized, or because an
        # n value greater than the number of stations supported is
        # given, a name of '?' is returned and c, w, p and h are left
        # in their current states.
        # Programs can obtain a list of all currently supported
        # stations by calling the routine repeatedly, with n=1,2,3...
        # When name='?' is seen, the list of stations has been
        # exhausted.
        # Station numbers, identifiers, names and other details are
        # subject to change and should not be hardwired into
        # application programs.
        # All station identifiers c are uppercase only;  lowercase
        # characters must be converted to uppercase by the calling
        # program.  The station names returned may contain both upper-
        # and lowercase.  All characters up to the first space are
        # checked;  thus an abbreviated ID will return the parameters
        # for the first station in the list which matches the
        # abbreviation supplied, and no station in the list will ever
        # contain embedded spaces.  c must not have leading spaces.
        # IMPORTANT -- BEWARE OF THE LONGITUDE SIGN CONVENTION.  The
        # longitude returned by sla_OBS is west-positive in accordance
        # with astronomical usage.  However, this sign convention is
        # left-handed and is the opposite of the one used by geographers;
        # elsewhere in SLALIB the preferable east-positive convention is
        # used.  In particular, note that for use in sla_AOP, sla_AOPPA
        # and sla_OAP the sign of the longitude must be reversed.
        # Users are urged to inform the author of any improvements
        # they would like to see made.  For example:
        # typographical corrections
        # more accurate parameters
        # better station identifiers or names
        # additional stations
        # P.T.Wallace   Starlink   15 March 2002
        # Copyright (c) 2002 Rutherford Appleton Laboratory
        # -
        AS2R = 0.484813681109535994e-5

        # Table of station identifiers
        CTAB = [
            "AAT",
            "LPO4.2",
            "LPO2.5",
            "LPO1",
            "LICK120",
            "MMT",
            "DAO72",
            "DUPONT",
            "MTHOP1.5",
            "STROMLO74",
            "ANU2.3",
            "GBVA140",
            "TOLOLO4M",
            "TOLOLO1.5M",
            "TIDBINBLA",
            "BLOEMF",
            "BOSQALEGRE",
            "FLAGSTF61",
            "LOWELL72",
            "HARVARD",
            "OKAYAMA",
            "KPNO158",
            "KPNO90",
            "KPNO84",
            "KPNO36FT",
            "KOTTAMIA",
            "ESO3.6",
            "MAUNAK88",
            "UKIRT",
            "QUEBEC1.6",
            "MTEKAR",
            "MTLEMMON60",
            "MCDONLD2.7",
            "MCDONLD2.1",
            "PALOMAR200",
            "PALOMAR60",
            "DUNLAP74",
            "HPROV1.93",
            "HPROV1.52",
            "SANPM83",
            "SAAO74",
            "TAUTNBG",
            "CATALINA61",
            "STEWARD90",
            "USSR6",
            "ARECIBO",
            "CAMB5KM",
            "CAMB1MILE",
            "EFFELSBERG",
            "GBT",
            "JODRELL1",
            "PARKES",
            "VLA",
            "SUGARGROVE",
            "USSR600",
            "NOBEYAMA",
            "JCMT",
            "ESONTT",
            "ST.ANDREWS",
            "APO3.5",
            "KECK1",
            "TAUTSCHM",
            "PALOMAR48",
            "UKST",
            "KISO",
            "ESOSCHM",
            "ATCA",
            "MOPRA",
            "SUBARU",
            "CFHT",
            "KECK2",
            "GEMININ",
            "FCRAO",
            "IRTF",
            "CSO",
            "VLT1",
            "VLT2",
            "VLT3",
            "VLT4",
            "GEMINIS",
            "KOSMA3M",
            "MAGELLAN1",
            "MAGELLAN2",
        ]

        if n is not None:
            c = CTAB[n]

        # Degrees, arcminutes, arcseconds to radians
        WEST = lambda ID, IAM, AS: AS2R * ((60 * (60 * ID + IAM)) + AS)
        NORTH = WEST
        EAST = lambda ID, IAM, AS: -1 * WEST(ID, IAM, AS)
        SOUTH = lambda ID, IAM, AS: -1 * WEST(ID, IAM, AS)

        observatories = {
            # AAT (Observer's Guide)
            "AAT": (
                "Anglo-Australian 3.9m Telescope",
                EAST(149, 3, 57.91),
                SOUTH(31, 16, 37.34),
                1164e0,
            ),
            # WHT (Gemini, April 1987)
            "LPO4.2": (
                "William Herschel 4.2m Telescope",
                WEST(17, 52, 53.9),
                NORTH(28, 45, 38.1),
                2332e0,
            ),
            # INT (Gemini, April 1987)
            "LPO2.5": (
                "Isaac Newton 2.5m Telescope",
                WEST(17, 52, 39.5),
                NORTH(28, 45, 43.2),
                2336e0,
            ),
            # JKT (Gemini, April 1987)
            "LPO1": (
                "Jacobus Kapteyn 1m Telescope",
                WEST(17, 52, 41.2),
                NORTH(28, 45, 39.9),
                2364e0,
            ),
            # Lick 120" (S.L.Allen, private communication, 2002)
            "LICK120": (
                "Lick 120 inch",
                WEST(121, 38, 13.689),
                NORTH(37, 20, 34.931),
                1286e0,
            ),
            # MMT 6.5m conversion (MMT Observatory website)
            "MMT": (
                "MMT 6.5m, Mt Hopkins",
                WEST(110, 53, 4.4),
                NORTH(31, 41, 19.6),
                2608e0,
            ),
            # Victoria B.C. 1.85m (1984 Almanac)
            "DAO72": (
                "DAO Victoria BC 1.85 metre",
                WEST(123, 25, 1.18),
                NORTH(48, 31, 11.9),
                238e0,
            ),
            # Las Campanas (1983 Almanac)
            "DUPONT": (
                "Du Pont 2.5m Telescope, Las Campanas",
                WEST(70, 42, 9.0),
                SOUTH(29, 0, 11.0),
                2280e0,
            ),
            # Mt Hopkins 1.5m (1983 Almanac)
            "MTHOP1.5": (
                "Mt Hopkins 1.5 metre",
                WEST(110, 52, 39.00),
                NORTH(31, 40, 51.4),
                2344e0,
            ),
            # Mt Stromlo 74" (1983 Almanac)
            "STROMLO74": (
                "Mount Stromlo 74 inch",
                EAST(149, 0, 27.59),
                SOUTH(35, 19, 14.3),
                767e0,
            ),
            # ANU 2.3m, SSO (Gary Hovey)
            "ANU2.3": (
                "Siding Spring 2.3 metre",
                EAST(149, 3, 40.3),
                SOUTH(31, 16, 24.1),
                1149e0,
            ),
            # Greenbank 140' (1983 Almanac)
            "GBVA140": (
                "Greenbank 140 foot",
                WEST(79, 50, 9.61),
                NORTH(38, 26, 15.4),
                881e0,
            ),
            # Cerro Tololo 4m (1982 Almanac)
            "TOLOLO4M": (
                "Cerro Tololo 4 metre",
                WEST(70, 48, 53.6),
                SOUTH(30, 9, 57.8),
                2235e0,
            ),
            # Cerro Tololo 1.5m (1982 Almanac)
            "TOLOLO1.5M": (
                "Cerro Tololo 1.5 metre",
                WEST(70, 48, 54.5),
                SOUTH(30, 9, 56.3),
                2225e0,
            ),
            # Tidbinbilla 64m (1982 Almanac)
            "TIDBINBLA": (
                "Tidbinbilla 64 metre",
                EAST(148, 58, 48.20),
                SOUTH(35, 24, 14.3),
                670e0,
            ),
            # Bloemfontein 1.52m (1981 Almanac)
            "BLOEMF": (
                "Bloemfontein 1.52 metre",
                EAST(26, 24, 18.0),
                SOUTH(29, 2, 18.0),
                1387e0,
            ),
            # Bosque Alegre 1.54m (1981 Almanac)
            "BOSQALEGRE": (
                "Bosque Alegre 1.54 metre",
                WEST(64, 32, 48.0),
                SOUTH(31, 35, 53.0),
                1250e0,
            ),
            # USNO 61" astrographic reflector, Flagstaff (1981 Almanac)
            "FLAGSTF61": (
                "USNO 61 inch astrograph, Flagstaff",
                WEST(111, 44, 23.6),
                NORTH(35, 11, 2.5),
                2316e0,
            ),
            # Lowell 72" (1981 Almanac)
            "LOWELL72": (
                "Perkins 72 inch, Lowell",
                WEST(111, 32, 9.3),
                NORTH(35, 5, 48.6),
                2198e0,
            ),
            # Harvard 1.55m (1981 Almanac)
            "HARVARD": (
                "Harvard College Observatory 1.55m",
                WEST(71, 33, 29.32),
                NORTH(42, 30, 19.0),
                185e0,
            ),
            # Okayama 1.88m (1981 Almanac)
            "OKAYAMA": (
                "Okayama 1.88 metre",
                EAST(133, 35, 47.29),
                NORTH(34, 34, 26.1),
                372e0,
            ),
            # Kitt Peak Mayall 4m (1981 Almanac)
            "KPNO158": (
                "Kitt Peak 158 inch",
                WEST(111, 35, 57.61),
                NORTH(31, 57, 50.3),
                2120e0,
            ),
            # Kitt Peak 90 inch (1981 Almanac)
            "KPNO90": (
                "Kitt Peak 90 inch",
                WEST(111, 35, 58.24),
                NORTH(31, 57, 46.9),
                2071e0,
            ),
            # Kitt Peak 84 inch (1981 Almanac)
            "KPNO84": (
                "Kitt Peak 84 inch",
                WEST(111, 35, 51.56),
                NORTH(31, 57, 29.2),
                2096e0,
            ),
            # Kitt Peak 36 foot (1981 Almanac)
            "KPNO36FT": (
                "Kitt Peak 36 foot",
                WEST(111, 36, 51.12),
                NORTH(31, 57, 12.1),
                1939e0,
            ),
            # Kottamia 74" (1981 Almanac)
            "KOTTAMIA": (
                "Kottamia 74 inch",
                EAST(31, 49, 30.0),
                NORTH(29, 55, 54.0),
                476e0,
            ),
            # La Silla 3.6m (1981 Almanac)
            "ESO3.6": (
                "ESO 3.6 metre",
                WEST(70, 43, 36.0),
                SOUTH(29, 15, 36.0),
                2428e0,
            ),
            # Mauna Kea 88 inch        (IfA website, Richard Wainscoat)
            "MAUNAK88": (
                "Mauna Kea 88 inch",
                WEST(155, 28, 9.96),
                NORTH(19, 49, 22.77),
                4213.6e0,
            ),
            # UKIRT (IfA website, Richard Wainscoat)
            "UKIRT": (
                "UK Infra Red Telescope",
                WEST(155, 28, 13.18),
                NORTH(19, 49, 20.75),
                4198.5e0,
            ),
            # Quebec 1.6m (1981 Almanac)
            "QUEBEC1.6": (
                "Quebec 1.6 metre",
                WEST(71, 9, 9.7),
                NORTH(45, 27, 20.6),
                1114e0,
            ),
            # Mt Ekar 1.82m (1981 Almanac)
            "MTEKAR": (
                "Mt Ekar 1.82 metre",
                EAST(11, 34, 15.0),
                NORTH(45, 50, 48.0),
                1365e0,
            ),
            # Mt Lemmon 60" (1981 Almanac)
            "MTLEMMON60": (
                "Mt Lemmon 60 inch",
                WEST(110, 42, 16.9),
                NORTH(32, 26, 33.9),
                2790e0,
            ),
            # Mt Locke 2.7m (1981 Almanac)
            "MCDONLD2.7": (
                "McDonald 2.7 metre",
                WEST(104, 1, 17.60),
                NORTH(30, 40, 17.7),
                2075e0,
            ),
            # Mt Locke 2.1m (1981 Almanac)
            "MCDONLD2.1": (
                "McDonald 2.1 metre",
                WEST(104, 1, 20.10),
                NORTH(30, 40, 17.7),
                2075e0,
            ),
            # Palomar 200" (1981 Almanac)
            "PALOMAR200": (
                "Palomar 200 inch",
                WEST(116, 51, 50.0),
                NORTH(33, 21, 22.0),
                1706e0,
            ),
            # Palomar 60" (1981 Almanac)
            "PALOMAR60": (
                "Palomar 60 inch",
                WEST(116, 51, 31.0),
                NORTH(33, 20, 56.0),
                1706e0,
            ),
            # David Dunlap 74" (1981 Almanac)
            "DUNLAP74": (
                "David Dunlap 74 inch",
                WEST(79, 25, 20.0),
                NORTH(43, 51, 46.0),
                244e0,
            ),
            # Haute Provence 1.93m (1981 Almanac)
            "HPROV1.93": (
                "Haute Provence 1.93 metre",
                EAST(5, 42, 46.75),
                NORTH(43, 55, 53.3),
                665e0,
            ),
            # Haute Provence 1.52m (1981 Almanac)
            "HPROV1.52": (
                "Haute Provence 1.52 metre",
                EAST(5, 42, 43.82),
                NORTH(43, 56, 0.2),
                667e0,
            ),
            # San Pedro Martir 83" (1981 Almanac)
            "SANPM83": (
                "San Pedro Martir 83 inch",
                WEST(115, 27, 47.0),
                NORTH(31, 2, 38.0),
                2830e0,
            ),
            # Sutherland 74" (1981 Almanac)
            "SAAO74": (
                "Sutherland 74 inch",
                EAST(20, 48, 44.3),
                SOUTH(32, 22, 43.4),
                1771e0,
            ),
            # Tautenburg 2m (1981 Almanac)
            "TAUTNBG": (
                "Tautenburg 2 metre",
                EAST(11, 42, 45.0),
                NORTH(50, 58, 51.0),
                331e0,
            ),
            # Catalina 61" (1981 Almanac)
            "CATALINA61": (
                "Catalina 61 inch",
                WEST(110, 43, 55.1),
                NORTH(32, 25, 0.7),
                2510e0,
            ),
            # Steward 90" (1981 Almanac)
            "STEWARD90": (
                "Steward 90 inch",
                WEST(111, 35, 58.24),
                NORTH(31, 57, 46.9),
                2071e0,
            ),
            # Russian 6m (1981 Almanac)
            "USSR6": ("USSR 6 metre", EAST(41, 26, 30.0), NORTH(43, 39, 12.0), 2100e0),
            # Arecibo 1000' (1981 Almanac)
            "ARECIBO": (
                "Arecibo 1000 foot",
                WEST(66, 45, 11.1),
                NORTH(18, 20, 36.6),
                496e0,
            ),
            # Cambridge 5km (1981 Almanac)
            "CAMB5KM": ("Cambridge 5km", EAST(0, 2, 37.23), NORTH(52, 10, 12.2), 17e0),
            # Cambridge 1 mile (1981 Almanac)
            "CAMB1MILE": (
                "Cambridge 1 mile",
                EAST(0, 2, 21.64),
                NORTH(52, 9, 47.3),
                17e0,
            ),
            # Bonn 100m (1981 Almanac)
            "EFFELSBERG": (
                "Effelsberg 100 metre",
                EAST(6, 53, 1.5),
                NORTH(50, 31, 28.6),
                366e0,
            ),
            # Green Bank Telescop
            "100m": (
                "Green Bank Telescope",
                WEST(79, 50, 23.406),
                NORTH(38, 25, 59.236),
                880e0,
            ),
            # Jodrell Bank Mk 1 (1981 Almanac)
            "JODRELL1": (
                "Jodrell Bank 250 foot",
                WEST(2, 18, 25.0),
                NORTH(53, 14, 10.5),
                78e0,
            ),
            # Australia Telescope Parkes Observatory   (Peter te Lintel Hekkert)
            "PARKES": (
                "Parkes 64 metre",
                EAST(148, 15, 44.3591),
                SOUTH(32, 59, 59.8657),
                391.79e0,
            ),
            # VLA (1981 Almanac)
            "VLA": (
                "Very Large Array",
                WEST(107, 37, 3.82),
                NORTH(34, 4, 43.5),
                2124e0,
            ),
            # Sugar Grove 150' (1981 Almanac)
            "SUGARGROVE": (
                "Sugar Grove 150 foot",
                WEST(79, 16, 23.0),
                NORTH(38, 31, 14.0),
                705e0,
            ),
            # Russian 600' (1981 Almanac)
            "USSR600": (
                "USSR 600 foot",
                EAST(41, 35, 25.5),
                NORTH(43, 49, 32.0),
                973e0,
            ),
            # Nobeyama 45 metre mm dish (based on 1981 Almanac entry)
            "NOBEYAMA": (
                "Nobeyama 45 metre",
                EAST(138, 29, 12.0),
                NORTH(35, 56, 19.0),
                1350e0,
            ),
            # James Clerk Maxwell 15 metre mm telescope, Mauna Kea   (IfA website, Richard Wainscoat, height from I.Coulson)
            "JCMT": (
                "JCMT 15 metre",
                WEST(155, 28, 37.20),
                NORTH(19, 49, 22.11),
                4111e0,
            ),
            # ESO 3.5 metre NTT, La Silla (K.Wirenstrand)
            "ESONTT": (
                "ESO 3.5 metre NTT",
                WEST(70, 43, 7.0),
                SOUTH(29, 15, 30.0),
                2377e0,
            ),
            # St Andrews University Observatory (1982 Almanac)
            "ST.ANDREWS": ("St Andrews", WEST(2, 48, 52.5), NORTH(56, 20, 12.0), 30e0),
            # Apache Point 3.5 metre (R.Owen)
            "APO3.5": (
                "Apache Point 3.5m",
                WEST(105, 49, 11.56),
                NORTH(32, 46, 48.96),
                2809e0,
            ),
            # W.M.Keck Observatory, Telescope 1        (William Lupton)
            "KECK1": (
                "Keck 10m Telescope #1",
                WEST(155, 28, 28.99),
                NORTH(19, 49, 33.41),
                4160e0,
            ),
            # Tautenberg Schmidt (1983 Almanac)
            "TAUTSCHM": (
                "Tautenberg 1.34 metre Schmidt",
                EAST(11, 42, 45.0),
                NORTH(50, 58, 51.0),
                331e0,
            ),
            # Palomar Schmidt (1981 Almanac)
            "PALOMAR48": (
                "Palomar 48-inch Schmidt",
                WEST(116, 51, 32.0),
                NORTH(33, 21, 26.0),
                1706e0,
            ),
            # UK Schmidt, Siding Spring (1983 Almanac)
            "UKST": (
                "UK 1.2 metre Schmidt, Siding Spring",
                EAST(149, 4, 12.8),
                SOUTH(31, 16, 27.8),
                1145e0,
            ),
            # Kiso Schmidt, Japan (1981 Almanac)
            "KISO": (
                "Kiso 1.05 metre Schmidt, Japan",
                EAST(137, 37, 42.2),
                NORTH(35, 47, 38.7),
                1130e0,
            ),
            # ESO Schmidt, La Silla (1981 Almanac)
            "ESOSCHM": (
                "ESO 1 metre Schmidt, La Silla",
                WEST(70, 43, 46.5),
                SOUTH(29, 15, 25.8),
                2347e0,
            ),
            # Australia Telescope Compact Array (WGS84 coordinates of Station 35, Mark Calabretta)
            "ATCA": (
                "Australia Telescope Compact Array",
                EAST(149, 33, 0.500),
                SOUTH(30, 18, 46.385),
                236.9e0,
            ),
            # Australia Telescope Mopra Observatory (Peter te Lintel Hekkert)
            "MOPRA": (
                "ATNF Mopra Observatory",
                EAST(149, 5, 58.732),
                SOUTH(31, 16, 4.451),
                850e0,
            ),
            # Subaru telescope, Mauna Kea     (IfA website, Richard Wainscoat)
            "SUBARU": (
                "Subaru 8m telescope",
                WEST(155, 28, 33.67),
                NORTH(19, 49, 31.81),
                4163e0,
            ),
            # Canada-France-Hawaii Telescope, Mauna Kea     (IfA website, Richard Wainscoat)
            "CFHT": (
                "Canada-France-Hawaii 3.6m Telescope",
                WEST(155, 28, 7.95),
                NORTH(19, 49, 30.91),
                4204.1e0,
            ),
            # W.M.Keck Observatory, Telescope 2       (William Lupton)
            "KECK2": (
                "Keck 10m Telescope #2",
                WEST(155, 28, 27.24),
                NORTH(19, 49, 35.62),
                4159.6e0,
            ),
            # Gemini North, Mauna Kea         (IfA website, Richard Wainscoat)
            "GEMININ": (
                "Gemini North 8-m telescope",
                WEST(155, 28, 8.57),
                NORTH(19, 49, 25.69),
                4213.4e0,
            ),
            # Five College Radio Astronomy Observatory     (Tim Jenness)
            "FCRAO": (
                "Five College Radio Astronomy Obs",
                WEST(72, 20, 42.0),
                NORTH(42, 23, 30.0),
                314e0,
            ),
            # NASA Infra Red Telescope Facility  (IfA website, Richard Wainscoat)
            "IRTF": (
                "NASA IR Telescope Facility, Mauna Kea",
                WEST(155, 28, 19.20),
                NORTH(19, 49, 34.39),
                4168.1e0,
            ),
            # Caltech Submillimeter Observatory    (IfA website, Richard Wainscoat; height estimated)
            "CSO": (
                "Caltech Sub-mm Observatory, Mauna Kea",
                WEST(155, 28, 31.79),
                NORTH(19, 49, 20.78),
                4080e0,
            ),
            # ESO VLT, UT1      (ESO website, VLT Whitebook Chapter 2)
            "VLT1": (
                "ESO VLT, Paranal, Chile: UT1",
                WEST(70, 24, 11.642),
                SOUTH(24, 37, 33.117),
                2635.43,
            ),
            # ESO VLT, UT2        (ESO website, VLT Whitebook Chapter 2)
            "VLT2": (
                "ESO VLT, Paranal, Chile: UT2",
                WEST(70, 24, 10.855),
                SOUTH(24, 37, 31.465),
                2635.43,
            ),
            # ESO VLT, UT3  (ESO website, VLT Whitebook Chapter 2)
            "VLT3": (
                "ESO VLT, Paranal, Chile: UT3",
                WEST(70, 24, 9.896),
                SOUTH(24, 37, 30.300),
                2635.43,
            ),
            # ESO VLT, UT4    (ESO website, VLT Whitebook Chapter 2)
            "VLT4": (
                "ESO VLT, Paranal, Chile: UT4",
                WEST(70, 24, 8.000),
                SOUTH(24, 37, 31.000),
                2635.43,
            ),
            # Gemini South, Cerro Pachon       (GPS readings by Patrick Wallace)
            "GEMINIS": (
                "Gemini South 8-m telescope",
                WEST(70, 44, 11.5),
                SOUTH(30, 14, 26.7),
                2738e0,
            ),
            # Cologne Observatory for Submillimeter Astronomy (KOSMA)   (Holger Jakob)
            "KOSMA3M": (
                "KOSMA 3m telescope, Gornergrat",
                EAST(7, 47, 3.48),
                NORTH(45, 58, 59.772),
                3141e0,
            ),
            # Magellan 1, 6.5m telescope at Las Campanas, Chile     (Skip Schaller)
            "MAGELLAN1": (
                "Magellan 1, 6.5m, Las Campanas",
                WEST(70, 41, 31.9),
                SOUTH(29, 0, 51.7),
                2408e0,
            ),
            # Magellan 2, 6.5m telescope at Las Campanas, Chile (Skip Schaller)
            "MAGELLAN2": (
                "Magellan 2, 6.5m, Las Campanas",
                WEST(70, 41, 33.5),
                SOUTH(29, 0, 50.3),
                2408e0,
            ),
        }

        # Exit
        name, w, p, h = observatories.get(c, ("?", None, None, None))

        return c, name, w, p, h

    @classmethod
    def xy2xy(cls, x1, y1, coeffs):
        # +
        # - - - - - -
        # X Y 2 X Y
        # - - - - - -
        #
        # Transform one [X,Y] into another using a linear model of the type
        # produced by the sla_FITXY routine.
        #
        # Given:
        # x1       d        x-coordinate
        # y1       d        y-coordinate
        # coeffs  d(6)      transformation coefficients (see note)
        #
        # Returned:
        # x2       d        x-coordinate
        # y2       d        y-coordinate
        #
        # The model relates two sets of [X,Y] coordinates as follows.
        # Naming the elements of COEFFS:
        #
        # COEFFS[1] = A
        # COEFFS[2] = B
        # COEFFS[3] = C
        # COEFFS[4] = D
        # COEFFS[5] = E
        # COEFFS[6] = F
        #
        # the present routine performs the transformation:
        #
        # x2 = A + B*x1 + C*y1
        # y2 = D + E*x1 + F*y1
        #
        # See also sla_FITXY, sla_PXY, sla_INVF, sla_DCMPF
        #
        # P.T.Wallace   Starlink   5 December 1994
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        # -

        x2 = coeffs[0] + coeffs[1] * x1 + coeffs[2] * y1
        y2 = coeffs[3] + coeffs[4] * x1 + coeffs[5] * y1

        return x2, y2

    @classmethod
    def zd(cls, ha, dec, phi):
        # +
        # - - -
        # Z D
        # - - -
        #
        # ha, Dec to Zenith Distance (double precision)
        #
        # Given:
        # ha     d     Hour Angle in radians
        # dec    d     declination in radians
        # phi    d     observatory latitude in radians
        #
        # The result is in the range 0 to pi.
        #
        # Notes:
        #
        # 1)  The latitude must be geodetic.  In critical applications,
        # corrections for polar motion should be applied.
        #
        # 2)  In some applications it will be important to specify the
        # correct type of hour angle and declination in order to
        # produce the required type of zenith distance.  In particular,
        # it may be important to distinguish between the zenith distance
        # as affected by refraction, which would require the "observed"
        # ha,Dec, and the zenith distance in vacuo, which would require
        # the "topocentric" ha,Dec.  If the effects of diurnal aberration
        # can be neglected, the "apparent" ha,Dec may be used instead of
        # the topocentric ha,Dec.
        #
        # 3)  No range checking of arguments is done.
        #
        # 4)  In applications which involve many zenith distance calculations,
        # rather than calling the present routine it will be more efficient
        # to use inline code, having previously computed fixed terms such
        # as sine and cosine of latitude, and perhaps sine and cosine of
        # declination.
        #
        # P.T.Wallace   Starlink   3 April 1994
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        # -

        SH = np.sin(ha)
        CH = np.cos(ha)
        SD = np.sin(dec)
        CD = np.cos(dec)
        SP = np.sin(phi)
        CP = np.cos(phi)
        X = CH * CD * SP - SD * CP
        Y = SH * CD
        Z = CH * CD * CP + SD * SP
        return np.arctan2(np.sqrt(X * X + Y * Y), Z)

    @classmethod
    def wait(cls, delay):
        # +
        # - - - - -
        # W A I T
        # - - - - -
        #
        # Interval wait
        #
        # !!! Version for: SPARC/SunOS4,
        # SPARC/Solaris2,
        # DEC Mips/Ultrix
        # DEC AXP/Digital Unix
        # Intel/Linux
        # Convex
        #
        # Given:
        # delay     real      delay in seconds
        #
        # Depends:  SLEEP (a Fortran Intrinsic on all above platforms)
        #
        # P.T.Wallace   Starlink   22 January 1998
        #
        # Copyright (C) 1998 Rutherford Appleton Laboratory
        # -
        time.sleep(np.rint(delay))
        return

    @classmethod
    def vxv(cls, va, vb):
        # +
        # - - - -
        # V X V
        # - - - -
        #
        # Vector product of two 3-vectors (single precision)
        #
        # Given:
        # va      real(3)     first vector
        # vb      real(3)     second vector
        #
        # Returned:
        # vc      real(3)     vector result
        #
        # P.T.Wallace   Starlink   March 1986
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        # Form the vector product va cross vb
        return np.cross(va, vb)

    @classmethod
    def vn(cls, v):
        # +
        # - - -
        # v N
        # - - -
        #
        # Normalizes a 3-vector also giving the modulus (single precision)
        #
        # Given:
        # v       real(3)      vector
        #
        # Returned:
        # uv      real(3)      unit vector in direction of v
        # vm      real         modulus of v
        #
        # If the modulus of v is zero, uv is set to zero as well
        #
        # P.T.Wallace   Starlink   23 November 1995
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        # Modulus
        vm = np.sqrt(np.dot(v, v))

        # Normalize the vector
        vm = 1.0 if (vm <= 0.0) else vm
        uv = v / vm

        return uv, vm

    @classmethod
    def vers(cls):
        # +
        #     - - - - -
        #      V E R S
        #     - - - - -
        #
        #  Report the SLALIB version number.
        #
        #  Given:
        #    None
        #
        #  Returned:
        #    VERSION   c*(*)   Version number, in the form 'm.n-r'.
        #                      The major version is m, the minor version n, and
        #                      release r.  The string passed in should be at least
        #                      8 characters in length, to account for the (remote)
        #                      possibility that these numbers will ever go to
        #                      two digits.
        #
        #  Notes:
        #
        #    To obtain the version number in a more easily processed form, see
        #    function sla_veri().
        #
        #    The sla_vers subroutine was introduced in SLALIB version 2.5-1, so
        #    if this function is absent, one can only tell that the release
        #    predates that one.
        return "2.5-4"

    @classmethod
    def veri(cls):
        # +
        #     - - - - -
        #      V E R I
        #     - - - - -
        #
        #  Report the SLALIB version number as an integer.
        #
        #  Given:
        #    None
        #
        #  The result is the SLALIB version number as an integer m*1e6+n*1e3+r,
        #  where m is the major version, n the minor version and r the release
        #  number.
        #
        #  Notes:
        #
        #    To obtain the version number in a printable form, see
        #    subroutine sla_vers(version).
        #
        #    The sla_veri subroutine was introduced in SLALIB version 2.5-1, so
        #    if this function is absent, one can only tell that the release
        #    predates that one.
        #
        #  Norman Gray   Starlink   8 April 2005
        #
        #  Copyright (C) 2005 Council for the Central Laboratory of the
        #  Research Councils
        #
        #  Licence:
        #    This program is free software; you can redistribute it and/or modify
        #    it under the terms of the GNU General Public License as published by
        #    the Free Software Foundation; either version 2 of the License, or
        #    (at your option) any later version.
        #
        #    This program is distributed in the hope that it will be useful,
        #    but WITHOUT ANY WARRANTY; without even the implied warranty of
        #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        #    GNU General Public License for more details.
        #
        #    You should have received a copy of the GNU General Public License
        #    along with this program (see SLA_CONDITIONS); if not, write to the
        #    Free Software Foundation, Inc., 59 Temple Place, Suite 330,
        #    Boston, MA  02111-1307  USA
        #
        # -
        return 2005004

    @classmethod
    def v2tp(cls, v, v0):
        # +
        # - - - - -
        # v 2 T P
        # - - - - -
        #
        # Given the direction cosines of a star and of the tangent point,
        # determine the star's tangent-plane coordinates.
        #
        # (single precision)
        #
        # Given:
        # v         r(3)    direction cosines of star
        # v0        r(3)    direction cosines of tangent point
        #
        # Returned:
        # xi,eta    r       tangent plane coordinates of star
        # j         i       status:   0 = OK
        # 1 = error, star too far from axis
        # 2 = error, antistar on tangent plane
        # 3 = error, antistar too far from axis
        #
        # Notes:
        #
        # 1  If vector v0 is not of unit length, or if vector v is of zero
        # length, the results will be wrong.
        #
        # 2  If v0 points at a pole, the returned xi,eta will be based on the
        # arbitrary assumption that the RA of the tangent point is zero.
        #
        # 3  This routine is the Cartesian equivalent of the routine sla_S2TP.
        #
        # P.T.Wallace   Starlink   27 November 1996
        #
        # Copyright (C) 1996 Rutherford Appleton Laboratory
        #
        # -

        TINY = 1e-6

        X, Y, Z = v
        X0, Y0, Z0 = v0
        R2 = X0 * X0 + Y0 * Y0
        R = np.sqrt(R2)
        if R == 0.0:
            R = 1e-20
            X0 = R

        W = X * X0 + Y * Y0
        D = W + Z * Z0
        if D > TINY:
            j = 0
        elif D >= 0.0:
            j = 1
            D = TINY
        elif D > -TINY:
            j = 2
            D = -TINY
        else:
            j = 3

        D = D * R
        xi = (Y * X0 - X * Y0) / D
        eta = (Z * R2 - Z0 * W) / D
        return xi, eta, j

    @classmethod
    def unpcd(cls, disco, x, y):
        # +
        # - - - - - -
        # U N P C D
        # - - - - - -
        #
        # Remove pincushion/barrel distortion from a distorted [x,y] to give
        # tangent-plane [x,y].
        #
        # Given:
        # disco    d      pincushion/barrel distortion coefficient
        # x,y      d      distorted coordinates
        #
        # Returned:
        # x,y      d      tangent-plane coordinates
        #
        # Notes:
        #
        # 1)  The distortion is of the form RP = R*(1+C*R^2), where R is
        # the radial distance from the tangent point, C is the disco
        # argument, and RP is the radial distance in the presence of
        # the distortion.
        #
        # 2)  For pincushion distortion, C is +ve;  for barrel distortion,
        # C is -ve.
        #
        # 3)  For x,y in "radians" - units of one projection radius,
        # which in the case of a photograph is the focal length of
        # the camera - the following disco values apply:
        #
        # Geometry          disco
        #
        # astrograph         0.0
        # Schmidt           -0.3333
        # AAT PF doublet  +147.069
        # AAT PF triplet  +178.585
        # AAT f/8          +21.20
        # JKT f/8          +13.32
        #
        # 4)  The present routine is a rigorous inverse of the companion
        # routine sla_PCD.  The expression for RP in Note 1 is rewritten
        # in the form x^3+a*x+b=0 and solved by standard techniques.
        #
        # 5)  Cases where the cubic has multiple real roots can sometimes
        # occur, corresponding to extreme instances of barrel distortion
        # where up to three different undistorted [X,Y]s all produce the
        # same distorted [X,Y].  However, only one solution is returned,
        # the one that produces the smallest change in [X,Y].
        #
        # P.T.Wallace   Starlink   3 September 2000
        #
        # Copyright (C) 2000 Rutherford Appleton Laboratory
        # -

        THIRD = 1e0 / 3e0

        # Distance of the point from the origin.
        RP = np.sqrt(x * x + y * y)

        # If zero, or if no distortion, no action is necessary.
        if RP != 0e0 and disco != 0e0:

            # Begin algebraic solution.
            Q = 1e0 / (3e0 * disco)
            R = RP / (2e0 * disco)
            W = Q * Q * Q + R * R

            # Continue if one real root, or three of which only one is
            # positive.
            if W >= 0e0:
                D = np.sqrt(W)
                W = R + D
                S = np.sign(np.abs(W) ** THIRD, W)
                W = R - D
                T = np.sign((np.abs(W)) ** THIRD, W)
                F = S + T
            else:
                # Three different real roots:  use geometrical method instead.
                W = 2e0 / np.sqrt(-3e0 * disco)
                C = 4e0 * RP / (disco * W * W * W)
                S = np.sqrt(1e0 - np.minimum(C * C, 1e0))
                T3 = np.arctan2(S, C)

                # The three solutions.
                F1 = W * np.cos((cls.D2PI - T3) / 3e0)
                F2 = W * np.cos(T3 / 3e0)
                F3 = W * np.cos((cls.D2PI + T3) / 3e0)

                # Pick the one that moves [X,Y] least.
                W1 = np.abs(F1 - RP)
                W2 = np.abs(F2 - RP)
                W3 = np.abs(F3 - RP)
                if W1 < W2 and W1 < W3:
                    F = F1
                elif W1 < W2 or W2 >= W3:
                    F = F3
                else:
                    F = F2
            # Remove the distortion.
            F = F / RP
            x = F * x
            y = F * y

        return x, y

    @classmethod
    def ue2pv(cls, date, u):
        # +
        # - - - - - -
        # u E 2 P V
        # - - - - - -
        #
        # Heliocentric position and velocity of a planet, asteroid or comet,
        # starting from orbital elements in the "universal variables" form.
        #
        # Given:
        # date     d       date, Modified Julian Date (JD-2400000.5)
        #
        # Given and returned:
        # u        d(13)   universal orbital elements (updated; Note 1)
        #
        # given    (1)   combined mass (M+m)
        # "      (2)   total energy of the orbit (alpha)
        # "      (3)   reference (osculating) epoch (t0)
        # "    (4-6)   position at reference epoch (r0)
        # "    (7-9)   velocity at reference epoch (v0)
        # "     (10)   heliocentric distance at reference epoch
        # "     (11)   r0.v0
        # returned  (12)   date (t)
        # "     (13)   universal eccentric anomaly (psi) of date
        #
        # Returned:
        # pv       d(6)    position (AU) and velocity (AU/s)
        # jstat    i       status:  0 = OK
        # -1 = radius vector zero
        # -2 = failed to converge
        #
        # Notes
        #
        # 1  The "universal" elements are those which define the orbit for the
        # purposes of the method of universal variables (see reference).
        # They consist of the combined mass of the two bodies, an epoch,
        # and the position and velocity vectors (arbitrary reference frame)
        # at that epoch.  The parameter set used here includes also various
        # quantities that can, in fact, be derived from the other
        # information.  This approach is taken to avoiding unnecessary
        # computation and loss of accuracy.  The supplementary quantities
        # are (i) alpha, which is proportional to the total energy of the
        # orbit, (ii) the heliocentric distance at epoch, (iii) the
        # outwards component of the velocity at the given epoch, (iv) an
        # estimate of psi, the "universal eccentric anomaly" at a given
        # date and (v) that date.
        #
        # 2  The companion routine is sla_EL2UE.  This takes the conventional
        # orbital elements and transforms them into the set of numbers
        # needed by the present routine.  A single prediction requires one
        # one call to sla_EL2UE followed by one call to the present routine;
        # for convenience, the two calls are packaged as the routine
        # sla_PLANEL.  Multiple predictions may be made by again
        # calling sla_EL2UE once, but then calling the present routine
        # multiple times, which is faster than multiple calls to sla_PLANEL.
        #
        # It is not obligatory to use sla_EL2UE to obtain the parameters.
        # However, it should be noted that because sla_EL2UE performs its
        # own validation, no checks on the contents of the array u are made
        # by the present routine.
        #
        # 3  date is the instant for which the prediction is required.  It is
        # in the TT timescale (formerly Ephemeris Time, ET) and is a
        # Modified Julian Date (JD-2400000.5).
        #
        # 4  The universal elements supplied in the array u are in canonical
        # units (solar masses, AU and canonical days).  The position and
        # velocity are not sensitive to the choice of reference frame.  The
        # sla_EL2UE routine in fact produces coordinates with respect to the
        # J2000 equator and equinox.
        #
        # 5  The algorithm was originally adapted from the EPHSLA program of
        # D.H.P.Jones (private communication, 1996).  The method is based
        # on Stumpff's Universal Variables.
        #
        # Reference:  Everhart, E. & Pitkin, E.T., Am.J.Phys. 51, 712, 1983.
        #
        # P.T.Wallace   Starlink   22 October 2005
        #
        # Copyright (C) 2005 Rutherford Appleton Laboratory
        #
        # -

        # Gaussian gravitational constant (exact)
        GCON = 0.01720209895e0
        # Canonical days to seconds
        CD2S = GCON / 86400e0
        # Test value for solution and maximum number of iterations
        TEST = 1e-13
        NITMAX = 25
        # Unpack the parameters.
        CM = u[0]
        ALPHA = u[1]
        T0 = u[2]
        P0 = u[3:6]
        V0 = u[6:9]

        R0 = u[9]
        SIGMA0 = u[10]
        T = u[11]
        PSI = u[12]

        # Approximately update the universal eccentric anomaly.
        PSI = PSI + (date - T) * GCON / R0

        # Time from reference epoch to date (in Canonical Days: a canonical
        # day is 58.1324409... days, defined as 1/GCON).
        DT = (date - T0) * GCON

        # Refine the universal eccentric anomaly, psi.
        NIT = 1
        W = 1e0
        TOL = 0e0
        FLAST = 0
        jstat = 0
        while np.abs(W) >= TOL:
            # Form half angles until BETA small enough.
            N = 0
            PSJ = PSI
            PSJ2 = PSJ * PSJ
            BETA = ALPHA * PSJ2
            while np.abs(BETA) > 0.7e0:
                N += 1
                BETA = BETA / 4e0
                PSJ = PSJ / 2e0
                PSJ2 = PSJ2 / 4e0

            # Calculate Universal Variables S0,S1,S2,S3 by nested series.
            S3 = (
                PSJ
                * PSJ2
                * (
                    (
                        (
                            (
                                ((BETA / 210e0 + 1e0) * BETA / 156e0 + 1e0)
                                * BETA
                                / 110e0
                                + 1e0
                            )
                            * BETA
                            / 72e0
                            + 1e0
                        )
                        * BETA
                        / 42e0
                        + 1e0
                    )
                    * BETA
                    / 20e0
                    + 1e0
                )
                / 6e0
            )

            S2 = (
                PSJ2
                * (
                    (
                        (
                            (
                                ((BETA / 182e0 + 1e0) * BETA / 132e0 + 1e0)
                                * BETA
                                / 90e0
                                + 1e0
                            )
                            * BETA
                            / 56e0
                            + 1e0
                        )
                        * BETA
                        / 30e0
                        + 1e0
                    )
                    * BETA
                    / 12e0
                    + 1e0
                )
                / 2e0
            )

            S1 = PSJ + ALPHA * S3
            S0 = 1e0 + ALPHA * S2

            # Undo the angle-halving.
            TOL = TEST
            while N > 0:
                S3 = 2e0 * (S0 * S3 + PSJ * S2)
                S2 = 2e0 * S1 * S1
                S1 = 2e0 * S0 * S1
                S0 = 2e0 * S0 * S0 - 1e0
                PSJ = PSJ + PSJ
                TOL += TOL
                N -= 1

            # Values of F and F' corresponding to the current value of psi.
            FF = R0 * S1 + SIGMA0 * S2 + CM * S3 - DT
            R = R0 * S0 + SIGMA0 * S1 + CM * S2

            # If first iteration, create dummy "last F".
            FLAST = FF if (NIT == 1) else FLAST

            # Check for sign change.
            if FF * FLAST < 0e0:

                # Sign change:  get psi adjustment using secant method.
                W = FF * (PLAST - PSI) / (FLAST - FF)
            else:

                # No sign change:  use Newton-Raphson method instead.
                if R == 0e0:
                    jstat = -1
                    break
                W = FF / R

            # Save the last psi and F values.
            PLAST = PSI
            FLAST = FF

            # Apply the Newton-Raphson or secant adjustment to psi.
            PSI = PSI - W

            # Next iteration, unless too many already.
            if NIT > NITMAX:
                jstat = -2
                break
            NIT += 1

        if jstat > -1:
            # Project the position and velocity vectors (scaling velocity to AU/s).
            W = CM * S2
            F = 1e0 - W / R0
            G = DT - CM * S3
            FD = -CM * S1 / (R0 * R)
            GD = 1e0 - W / R
            pv = P0 * F + V0 * G
            pv[3:] = CD2S * (P0 * FD + V0 * GD)

        # Update the parameters to allow speedy prediction of PSI next time.
        u[11] = date
        u[12] = PSI

        # OK exit.
        jstat = 0

        return u, pv, jstat

    @classmethod
    def tpv2c(cls, xi, eta, v):
        # +
        # - - - - - -
        # T P v 2 C
        # - - - - - -
        #
        # Given the tangent-plane coordinates of a star and its direction
        # cosines, determine the direction cosines of the tangent-point.
        #
        # (single precision)
        #
        # Given:
        # xi,eta    r       tangent plane coordinates of star
        # v         r(3)    direction cosines of star
        #
        # Returned:
        # v01       r(3)    direction cosines of tangent point, solution 1
        # v02       r(3)    direction cosines of tangent point, solution 2
        # n         i       number of solutions:
        # 0 = no solutions returned (note 2)
        # 1 = only the first solution is useful (note 3)
        # 2 = both solutions are useful (note 3)
        #
        # Notes:
        #
        # 1  The vector v must be of unit length or the result will be wrong.
        #
        # 2  Cases where there is no solution can only arise near the poles.
        # For example, it is clearly impossible for a star at the pole
        # itself to have a non-zero xi value, and hence it is meaningless
        # to ask where the tangent point would have to be.
        #
        # 3  Also near the poles, cases can arise where there are two useful
        # solutions.  The argument n indicates whether the second of the
        # two solutions returned is useful.  n=1 indicates only one useful
        # solution, the usual case;  under these circumstances, the second
        # solution can be regarded as valid if the vector v02 is interpreted
        # as the "over-the-pole" case.
        #
        # 4  This routine is the Cartesian equivalent of the routine sla_TPS2C.
        #
        # P.T.Wallace   Starlink   5 June 1995
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        X = v[0]
        Y = v[1]
        Z = v[2]
        RXY2 = X * X + Y * Y
        XI2 = xi * xi
        ETA2P1 = eta * eta + 1.0
        SDF = Z * np.sqrt(XI2 + ETA2P1)
        R2 = RXY2 * ETA2P1 - Z * Z * XI2
        v01 = np.zeros(3)
        v02 = np.zeros(3)
        if R2 > 0.0:
            R = np.sqrt(R2)
            C = (SDF * eta + R) / (ETA2P1 * np.sqrt(RXY2 * (R2 + XI2)))
            v01[0] = C * (X * R + Y * xi)
            v01[1] = C * (Y * R - X * xi)
            v01[2] = (SDF - eta * R) / ETA2P1
            R = -R
            C = (SDF * eta + R) / (ETA2P1 * np.sqrt(RXY2 * (R2 + XI2)))
            v02[0] = C * (X * R + Y * xi)
            v02[1] = C * (Y * R - X * xi)
            v02[2] = (SDF - eta * R) / ETA2P1
            n = 1 if np.abs(SDF) < 1.0 else 2
        else:
            n = 0
        return v01, v02, n

    @classmethod
    def ue2el(cls, u):
        # +
        # - - - - - -
        # u e 2 e L
        # - - - - - -
        #
        # Transform universal elements into conventional heliocentric
        # osculating elements.
        #
        # Given:
        # u         d(13)  universal orbital elements (Note 1)
        #
        # (1)  combined mass (M+m)
        # (2)  total energy of the orbit (alpha)
        # (3)  reference (osculating) epoch (t0)
        # (4-6)  position at reference epoch (r0)
        # (7-9)  velocity at reference epoch (v0)
        # (10)  heliocentric distance at reference epoch
        # (11)  r0.v0
        # (12)  date (t)
        # (13)  universal eccentric anomaly (psi) of date, approx
        #
        # jformr    i      requested element set (1-3; Note 3)
        #
        # Returned:
        # jform     d      element set actually returned (1-3; Note 4)
        # epoch     d      epoch of elements (TT MJD)
        # orbinc    d      inclination (radians)
        # anode     d      longitude of the ascending node (radians)
        # perih     d      longitude or argument of perihelion (radians)
        # aorq      d      mean distance or perihelion distance (AU)
        # e         d      eccentricity
        # aorl      d      mean anomaly or longitude (radians, jform=1,2 only)
        # dm        d      daily motion (radians, jform=1 only)
        # jstat     i      status:  0 = OK
        # -1 = illegal combined mass
        # -2 = illegal jformr
        # -3 = position/velocity out of range
        #
        # Notes
        #
        # 1  The "universal" elements are those which define the orbit for the
        # purposes of the method of universal variables (see reference 2).
        # They consist of the combined mass of the two bodies, an epoch,
        # and the position and velocity vectors (arbitrary reference frame)
        # at that epoch.  The parameter set used here includes also various
        # quantities that can, in fact, be derived from the other
        # information.  This approach is taken to avoiding unnecessary
        # computation and loss of accuracy.  The supplementary quantities
        # are (i) alpha, which is proportional to the total energy of the
        # orbit, (ii) the heliocentric distance at epoch, (iii) the
        # outwards component of the velocity at the given epoch, (iv) an
        # estimate of psi, the "universal eccentric anomaly" at a given
        # date and (v) that date.
        #
        # 2  The universal elements are with respect to the mean equator and
        # equinox of epoch J2000.  The orbital elements produced are with
        # respect to the J2000 ecliptic and mean equinox.
        #
        # 3  Three different element-format options are supported:
        #
        # Option jform=1, suitable for the major planets:
        #
        # epoch  = epoch of elements (TT MJD)
        # orbinc = inclination i (radians)
        # anode  = longitude of the ascending node, big omega (radians)
        # perih  = longitude of perihelion, curly pi (radians)
        # aorq   = mean distance, a (AU)
        # e      = eccentricity, e
        # aorl   = mean longitude L (radians)
        # dm     = daily motion (radians)
        #
        # Option jform=2, suitable for minor planets:
        #
        # epoch  = epoch of elements (TT MJD)
        # orbinc = inclination i (radians)
        # anode  = longitude of the ascending node, big omega (radians)
        # perih  = argument of perihelion, little omega (radians)
        # aorq   = mean distance, a (AU)
        # e      = eccentricity, e
        # aorl   = mean anomaly M (radians)
        #
        # Option jform=3, suitable for comets:
        #
        # epoch  = epoch of perihelion (TT MJD)
        # orbinc = inclination i (radians)
        # anode  = longitude of the ascending node, big omega (radians)
        # perih  = argument of perihelion, little omega (radians)
        # aorq   = perihelion distance, q (AU)
        # e      = eccentricity, e
        #
        # 4  It may not be possible to generate elements in the form
        # requested through JFORMR.  The caller is notified of the form
        # of elements actually returned by means of the jform argument:
        #
        # jformr   jform     meaning
        #
        # 1        1       OK - elements are in the requested format
        # 1        2       never happens
        # 1        3       orbit not elliptical
        #
        # 2        1       never happens
        # 2        2       OK - elements are in the requested format
        # 2        3       orbit not elliptical
        #
        # 3        1       never happens
        # 3        2       never happens
        # 3        3       OK - elements are in the requested format
        #
        # 5  The arguments returned for each value of jform (cf Note 6: jform
        # may not be the same as jformr) are as follows:
        #
        # jform         1              2              3
        # epoch         t0             t0             T
        # orbinc        i              i              i
        # anode         Omega          Omega          Omega
        # perih         curly pi       omega          omega
        # aorq          a              a              q
        # e             e              e              e
        # aorl          L              M              -
        # dm            n              -              -
        #
        # where:
        #
        # t0           is the epoch of the elements (MJD, TT)
        # T              "    epoch of perihelion (MJD, TT)
        # i              "    inclination (radians)
        # Omega          "    longitude of the ascending node (radians)
        # curly pi       "    longitude of perihelion (radians)
        # omega          "    argument of perihelion (radians)
        # a              "    mean distance (AU)
        # q              "    perihelion distance (AU)
        # e              "    eccentricity
        # L              "    longitude (radians, 0-2pi)
        # M              "    mean anomaly (radians, 0-2pi)
        # n              "    daily motion (radians)
        # -             means no value is set
        #
        # 6  At very small inclinations, the longitude of the ascending node
        # anode becomes indeterminate and under some circumstances may be
        # set arbitrarily to zero.  Similarly, if the orbit is close to
        # circular, the true anomaly becomes indeterminate and under some
        # circumstances may be set arbitrarily to zero.  In such cases,
        # the other elements are automatically adjusted to compensate,
        # and so the elements remain a valid description of the orbit.
        #
        # References:
        #
        # 1  Sterne, Theodore E., "An Introduction to Celestial Mechanics",
        # Interscience Publishers Inc., 1960.  Section 6.7, p199.
        #
        # 2  Everhart, E. & Pitkin, E.T., Am.J.Phys. 51, 712, 1983.
        #
        # Depends:  sla_PV2EL
        #
        # P.T.Wallace   Starlink   18 March 1999
        #
        # Copyright (C) 1999 Rutherford Appleton Laboratory
        # -

        # Gaussian gravitational constant (exact)

        GCON = 0.01720209895e0

        # Canonical days to seconds

        CD2S = GCON / 86400e0

        # Unpack the universal elements.
        PMASS = u[0] - 1e0
        DATE = u[2]
        PV = np.zeros(6)
        PV[:3] = u
        PV[3:] = u[3:6] * CD2S

        # Convert the position and velocity etc into conventional elements.
        jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm, jstat = cls.pv2el(PV, DATE, PMASS)

        return jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm, jstat

    @classmethod
    def tp2v(cls, xi, eta, v0):
        # +
        # - - - - -
        # T P 2 v
        # - - - - -
        #
        # Given the tangent-plane coordinates of a star and the direction
        # cosines of the tangent point, determine the direction cosines
        # of the star.
        #
        # (single precision)
        #
        # Given:
        # xi,eta    r       tangent plane coordinates of star
        # v0        r(3)    direction cosines of tangent point
        #
        # Returned:
        # v         r(3)    direction cosines of star
        #
        # Notes:
        #
        # 1  If vector v0 is not of unit length, the returned vector v will
        # be wrong.
        #
        # 2  If vector v0 points at a pole, the returned vector v will be
        # based on the arbitrary assumption that the RA of the tangent
        # point is zero.
        #
        # 3  This routine is the Cartesian equivalent of the routine sla_TP2S.
        #
        # P.T.Wallace   Starlink   11 February 1995
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        X = v0[0]
        Y = v0[1]
        Z = v0[2]
        F = np.sqrt(1.0 + xi * xi + eta * eta)
        R = np.sqrt(X * X + Y * Y)
        if R == 0.0:
            R = 1e-20
            X = R

        v = np.zeros(3)
        v[0] = (X - (xi * Y + eta * X * Z) / R) / F
        v[1] = (Y + (xi * X - eta * Y * Z) / R) / F
        v[2] = (Z + eta * R) / F

        return v

    @classmethod
    def tps2c(cls, xi, eta, ra, dec):
        # +
        # - - - - - -
        # T P S 2 C
        # - - - - - -
        #
        # From the tangent plane coordinates of a star of known ra,Dec,
        # determine the ra,Dec of the tangent point.
        #
        # (single precision)
        #
        # Given:
        # xi,eta      r    tangent plane rectangular coordinates
        # ra,dec      r    spherical coordinates
        #
        # Returned:
        # raz1,decz1  r    spherical coordinates of tangent point, solution 1
        # raz2,decz2  r    spherical coordinates of tangent point, solution 2
        # n           i    number of solutions:
        # 0 = no solutions returned (note 2)
        # 1 = only the first solution is useful (note 3)
        # 2 = both solutions are useful (note 3)
        #
        # Notes:
        #
        # 1  The raz1 and raz2 values are returned in the range 0-2pi.
        #
        # 2  Cases where there is no solution can only arise near the poles.
        # For example, it is clearly impossible for a star at the pole
        # itself to have a non-zero xi value, and hence it is
        # meaningless to ask where the tangent point would have to be
        # to bring about this combination of xi and DEC.
        #
        # 3  Also near the poles, cases can arise where there are two useful
        # solutions.  The argument n indicates whether the second of the
        # two solutions returned is useful.  n=1 indicates only one useful
        # solution, the usual case;  under these circumstances, the second
        # solution corresponds to the "over-the-pole" case, and this is
        # reflected in the values of raz2 and decz2 which are returned.
        #
        # 4  The decz1 and decz2 values are returned in the range +/-pi, but
        # in the usual, non-pole-crossing, case, the range is +/-pi/2.
        #
        # 5  This routine is the spherical equivalent of the routine sla_DTPV2C.
        #
        # Depends:  sla_RANORM
        #
        # P.T.Wallace   Starlink   5 June 1995
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        X2 = xi * xi
        Y2 = eta * eta
        SD = np.sin(dec)
        CD = np.cos(dec)
        SDF = SD * np.sqrt(1.0 + X2 + Y2)
        R2 = CD * CD * (1.0 + Y2) - SD * SD * X2
        if R2 >= 0.0:
            R = np.sqrt(R2)
            S = SDF - eta * R
            C = SDF * eta + R
            R = 1.0 if (xi == 0.0 and R == 0.0) else R
            raz1 = cls.ranorm(ra - np.arctan2(xi, R))
            decz1 = np.arctan2(S, C)
            R = -R
            S = SDF - eta * R
            C = SDF * eta + R
            raz2 = cls.ranorm(ra - np.arctan2(xi, R))
            decz2 = np.arctan2(S, C)
            n = 1 if np.abs(SDF) < 1.0 else 2
        else:
            n = 0

        return raz1, decz1, raz2, decz2, n

    @classmethod
    def tp2s(cls, xi, eta, raz, decz):
        # +
        # - - - - -
        # T P 2 S
        # - - - - -
        #
        # Transform tangent plane coordinates into spherical
        # (single precision)
        #
        # Given:
        # xi,eta      real  tangent plane rectangular coordinates
        # raz,decz    real  spherical coordinates of tangent point
        #
        # Returned:
        # ra,dec      real  spherical coordinates (0-2pi,+/-pi/2)
        #
        # Depends:        sla_RANORM
        #
        # P.T.Wallace   Starlink   24 July 1995
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        SDECZ = np.sin(decz)
        CDECZ = np.cos(decz)

        DENOM = CDECZ - eta * SDECZ

        ra = cls.ranorm(np.arctan2(xi, DENOM) + raz)
        dec = np.arctan2(SDECZ + eta * CDECZ, np.sqrt(xi * xi + DENOM * DENOM))

        return ra, dec

    @classmethod
    def svdsol(cls, m, n, mp, nnp, b, u, w, v):
        # +
        # - - - - - - -
        # S v D S O L
        # - - - - - - -
        #
        # From a given vector and the SVD of a matrix (as obtained from
        # the SVD routine), obtain the solution vector (double precision)
        #
        # This routine solves the equation:
        #
        # A . x = b
        #
        # where:
        #
        # A   is a given m (rows) x n (columns) matrix, where m >= n
        # x   is the n-vector we wish to find
        # b   is a given m-vector
        #
        # by means of the Singular Value Decomposition method (SVD).  In
        # this method, the matrix A is first factorised (for example by
        # the routine sla_SVD) into the following components:
        #
        # A = u x w x VT
        #
        # where:
        #
        # A   is the m (rows) x n (columns) matrix
        # u   is an m x n column-orthogonal matrix
        # w   is an n x n diagonal matrix with w(I,I) >= 0
        # VT  is the transpose of an NxN orthogonal matrix
        #
        # Note that m and n, above, are the LOGICAL dimensions of the
        # matrices and vectors concerned, which can be located in
        # arrays of larger PHYSICAL dimensions mp and NP.
        #
        # The solution is found from the expression:
        #
        # x = v . [diag(1/Wj)] . (transpose(u) . b)
        #
        # Notes:
        #
        # 1)  If matrix A is square, and if the diagonal matrix w is not
        # adjusted, the method is equivalent to conventional solution
        # of simultaneous equations.
        #
        # 2)  If M>N, the result is a least-squares fit.
        #
        # 3)  If the solution is poorly determined, this shows up in the
        # SVD factorisation as very small or zero Wj values.  Where
        # a Wj value is small but non-zero it can be set to zero to
        # avoid ill effects.  The present routine detects such zero
        # Wj values and produces a sensible solution, with highly
        # correlated terms kept under control rather than being allowed
        # to elope to infinity, and with meaningful values for the
        # other terms.
        #
        # Given:
        # m,n    i         numbers of rows and columns in matrix A
        # mp,np  i         physical dimensions of array containing matrix A
        # b      d(m)      known vector b
        # u      d(mp,np)  array containing MxN matrix u
        # w      d(n)      NxN diagonal matrix w (diagonal elements only)
        # v      d(np,np)  array containing NxN orthogonal matrix v
        #
        # Returned:
        # work   d(n)      workspace
        # x      d(n)      unknown vector x
        #
        # Reference:
        # Numerical Recipes, section 2.9.
        #
        # P.T.Wallace   Starlink   29 October 1993
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        # Calculate [diag(1/Wj)] . transpose(u) . b (or zero for zero Wj)
        work = np.linalg.multi_dot([np.diag(1 / w), u.T, b])

        # Multiply by matrix v to get result
        x = work @ v
        return work, x

    @classmethod
    def ctf2r(cls, ihour, imin, sec):
        # +
        # - - - - - -
        # C T F 2 R
        # - - - - - -
        #
        # Convert hours, minutes, seconds to radians (single precision)
        #
        # Given:
        # ihour       int       hours
        # imin        int       minutes
        # sec         real      seconds
        #
        # Returned:
        # rad         real      angle in radians
        # j           int       status:  0 = OK
        # 1 = ihour outside range 0-23
        # 2 = imin outside range 0-59
        # 3 = sec outside range 0-59.999...
        #
        # Called:
        # sla_CTF2D
        #
        # Notes:
        #
        # 1)  The result is computed even if any of the range checks
        # fail.
        #
        # 2)  The sign must be dealt with outside this routine.
        #
        # P.T.Wallace   Starlink   November 1984
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        # Turns to radians

        T2R = 6.283185307179586476925287

        # Convert to turns then radians
        days, j = cls.ctf2d(ihour, imin, sec, TURNS)
        rad = T2R * TURNS

        return rad, j

    @classmethod
    def dtf2r(cls, ihour, imin, sec):
        # +
        # - - - - - -
        # D T F 2 R
        # - - - - - -
        #
        # Convert hours, minutes, seconds to radians (double precision)
        #
        # Given:
        # ihour       int       hours
        # imin        int       minutes
        # sec         dp        seconds
        #
        # Returned:
        # rad         dp        angle in radians
        # j           int       status:  0 = OK
        # 1 = ihour outside range 0-23
        # 2 = imin outside range 0-59
        # 3 = sec outside range 0-59.999...
        #
        # Called:
        # sla_DTF2D
        #
        # Notes:
        #
        # 1)  The result is computed even if any of the range checks fail.
        #
        # 2)  The sign must be dealt with outside this routine.
        #
        # P.T.Wallace   Starlink   July 1984
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        # Turns to radians

        T2R = 6.283185307179586476925287e0

        # Convert to turns then radians
        TURNS, j = cls.dtf2d(ihour, imin, sec)
        rad = T2R * TURNS

        return rad, j

    @classmethod
    def dtf2d(cls, IHOUR, IMIN, SEC):
        # +
        # - - - - - -
        # D T F 2 D
        # - - - - - -
        #
        # Convert hours, minutes, seconds to days (double precision)
        #
        # Given:
        # IHOUR       int       hours
        # IMIN        int       minutes
        # SEC         dp        seconds
        #
        # Returned:
        # DAYS        dp        interval in days
        # J           int       status:  0 = OK
        # 1 = IHOUR outside range 0-23
        # 2 = IMIN outside range 0-59
        # 3 = SeC outside range 0-59.999...
        #
        # Notes:
        #
        # 1)  The result is computed even if any of the range checks fail.
        #
        # 2)  The sign must be dealt with outside this routine.
        #
        # P.T.Wallace   Starlink   July 1984
        # -

        # Seconds per day
        D2S = 86400e0

        # Preset status
        J = 0

        # Validate sec, min, hour
        J = 3 if (SEC < 0e0 or SEC >= 60e0) else J
        J = 2 if (IMIN < 0.0 or IMIN > 59) else J
        J = 1 if (IHOUR < 0.0 or IHOUR > 23) else J

        # Compute interval
        DAYS = (60e0 * (60e0 * IHOUR + IMIN) + SEC) / D2S

        return DAYS, J

    @classmethod
    def bear(cls, a1, b1, a2, b2):
        # +
        # - - - - -
        # B E A R
        # - - - - -
        #
        # Bearing (position angle) of one point on a sphere relative to another
        # (single precision)
        #
        # Given:
        # a1,b1    r    spherical coordinates of one point
        # a2,b2    r    spherical coordinates of the other point
        #
        # (The spherical coordinates are RA,Dec, Long,Lat etc, in radians.)
        #
        # The result is the bearing (position angle), in radians, of point
        # a2,b2 as seen from point a1,B1.  It is in the range +/- pi.  If
        # a2,b2 is due east of a1,b1 the bearing is +pi/2.  Zero is returned
        # if the two points are coincident.
        #
        # P.T.Wallace   Starlink   23 March 1991
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        DA = a2 - a1
        Y = np.sin(DA) * np.cos(b2)
        X = np.sin(b2) * np.cos(b1) - np.cos(b2) * np.sin(b1) * np.cos(DA)
        return np.arctan2(Y, X) if X != 0.0 or Y != 0.0 else 0.0

    @classmethod
    def cr2tf(cls, ndp, angle):
        # +
        # - - - - - -
        # C R 2 T F
        # - - - - - -
        #
        # Convert an angle in radians into hours, minutes, seconds
        # (single precision)
        #
        # Given:
        # ndp       int      number of decimal places of seconds
        # angle     real     angle in radians
        #
        # Returned:
        # sign      char     '+' or '-'
        # ihmsf     int(4)   hours, minutes, seconds, fraction
        #
        # Notes:
        #
        # 1)  ndp less than zero is interpreted as zero.
        #
        # 2)  The largest useful value for ndp is determined by the size of
        # angle, the format of
        # machine, and the risk of overflowing ihmsf(4).  On some
        # architectures, for angle up to 2pi, the available floating-point
        # precision corresponds roughly to ndp=3.  This is well below
        # the ultimate limit of ndp=9 set by the capacity of a typical
        # 32-bit ihmsf(4).
        #
        # 3)  The absolute value of angle may exceed 2pi.  In cases where it
        # does not, it is up to the caller to test for and handle the
        # case where angle is very nearly 2pi and rounds up to 24 hours,
        # by testing for IHMSF[1] = 24 and setting ihmsf(1-4) to zero.
        #
        # Depends:  sla_CD2TF
        #
        # Last revision:   26 December 2004
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        #
        #
        # -
        # Turns to radians

        T2R = 6.283185307179586476925287

        # Scale then use days to h,m,s routine
        sign, ihmsf = cls.cd2tf(ndp, angle / T2R)

        return sign, ihmsf

    @classmethod
    def dtps2c(cls, xi, eta, ra, dec):

        # +
        # - - - - - - -
        # D T P S 2 C
        # - - - - - - -
        #
        # From the tangent plane coordinates of a star of known ra,Dec,
        # determine the ra,Dec of the tangent point.
        #
        # (double precision)
        #
        # Given:
        # xi,eta      d    tangent plane rectangular coordinates
        # ra,dec      d    spherical coordinates
        #
        # Returned:
        # raz1,decz1  d    spherical coordinates of tangent point, solution 1
        # raz2,decz2  d    spherical coordinates of tangent point, solution 2
        # n           i    number of solutions:
        # 0 = no solutions returned (note 2)
        # 1 = only the first solution is useful (note 3)
        # 2 = both solutions are useful (note 3)
        #
        # Notes:
        #
        # 1  The raz1 and raz2 values are returned in the range 0-2pi.
        #
        # 2  Cases where there is no solution can only arise near the poles.
        # For example, it is clearly impossible for a star at the pole
        # itself to have a non-zero xi value, and hence it is
        # meaningless to ask where the tangent point would have to be
        # to bring about this combination of xi and DEC.
        #
        # 3  Also near the poles, cases can arise where there are two useful
        # solutions.  The argument n indicates whether the second of the
        # two solutions returned is useful.  n=1 indicates only one useful
        # solution, the usual case;  under these circumstances, the second
        # solution corresponds to the "over-the-pole" case, and this is
        # reflected in the values of raz2 and decz2 which are returned.
        #
        # 4  The decz1 and decz2 values are returned in the range +/-pi, but
        # in the usual, non-pole-crossing, case, the range is +/-pi/2.
        #
        # 5  This routine is the spherical equivalent of the routine sla_DTPV2C.
        #
        # Depends:  sla_DRANRM
        #
        # P.T.Wallace   Starlink   5 June 1995
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -
        raz1 = 0.0
        decz1 = 0.0
        raz2 = 0.0
        decz2 = 0.0

        X2 = xi * xi
        Y2 = eta * eta
        SD = np.sin(dec)
        CD = np.cos(dec)
        SDF = SD * np.sqrt(1e0 + X2 + Y2)
        R2 = CD * CD * (1e0 + Y2) - SD * SD * X2
        if R2 >= 0e0:
            R = np.sqrt(R2)
            S = SDF - eta * R
            C = SDF * eta + R
            R = 1e0 if (xi == 0e0 and R == 0e0) else R
            raz1 = cls.dranrm(ra - np.arctan2(xi, R))
            decz1 = np.arctan2(S, C)
            R = -R
            S = SDF - eta * R
            C = SDF * eta + R
            raz2 = cls.dranrm(ra - np.arctan2(xi, R))
            decz2 = np.arctan2(S, C)
            n = 1 if np.abs(SDF) < 1e0 else 2
        else:
            n = 0

        return raz1, decz1, raz2, decz2, n

    @classmethod
    def ecor(cls, rm, dm, iy, iid, fd):
        # +
        # - - - - -
        # E C O R
        # - - - - -
        #
        # Component of Earth orbit velocity and heliocentric
        # light time in a given direction (single precision)
        #
        # Given:
        # rm,dm    real    mean RA, Dec of date (radians)
        # iy       int     year
        # id       int     day in year (1 = Jan 1st)
        # fd       real    fraction of day
        #
        # Returned:
        # rv       real    component of Earth orbital velocity (km/sec)
        # tl       real    component of heliocentric light time (sec)
        #
        # Notes:
        #
        # 1  The date and time is TDB (loosely ET) in a Julian calendar
        # which has been aligned to the ordinary Gregorian
        # calendar for the interval 1900 March 1 to 2100 February 28.
        # The year and day can be obtained by calling sla_CALYD or
        # sla_CLYD.
        #
        # 2  Sign convention:
        #
        # The velocity component is +ve when the Earth is receding from
        # the given point on the sky.  The light time component is +ve
        # when the Earth lies between the Sun and the given point on
        # the sky.
        #
        # 3  Accuracy:
        #
        # The velocity component is usually within 0.004 km/s of the
        # correct value and is never in error by more than 0.007 km/s.
        # The error in light time correction is about 0.03s at worst,
        # but is usually better than 0.01s. For applications requiring
        # higher accuracy, see the sla_EVP and sla_EPV routines.
        #
        # Depends:  sla_EARTH, sla_CS2C, sla_VDV
        #
        # Last revision:   5 April 2005
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        # License:
        # This program is free software; you can redistribute it and/or modify
        # it under the terms of the GNU General Public License as published by
        # the Free Software Foundation; either version 2 of the License, or
        # (at your option) any later version.
        #
        # This program is distributed in the hope that it will be useful,
        # but WITHOUT ANY WARRANTY; without even the implied warranty of
        # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        # GNU General Public License for more details.
        #
        # You should have received a copy of the GNU General Public License
        # along with this program (see SLA_CONDITIONS); if not, write to the
        # Free Software Foundation, Inc., 59 Temple Place, Suite 330,
        # Boston, MA  02111-1307  USA
        #
        # -
        # AU to km and light sec (1985 Almanac)
        AUKM = 1.4959787066e8
        AUSEC = 499.0047837

        # Sun:Earth position & velocity vector
        pv = cls.earth(iy, iid, fd)

        # Star position vector
        v = cls.cs2c(rm, dm)

        # Velocity component
        rv = -AUKM * cls.vdv(pv[3], v)

        # Light time component
        tl = AUSEC * cls.vdv(pv[0], v)

        return rv, tl

    @classmethod
    def earth(cls, iy, id, fd):
        # +
        # - - - - - -
        # E A R T H
        # - - - - - -
        #
        # Approximate heliocentric position and velocity of the Earth
        #
        # Given:
        # iy       I       year
        # id       I       day in year (1 = Jan 1st)
        # fd       R       fraction of day
        #
        # Returned:
        # pv       R(6)    Earth position & velocity vector
        #
        # Notes:
        #
        # 1  The date and time is TDB (loosely ET) in a Julian calendar
        # which has been aligned to the ordinary Gregorian
        # calendar for the interval 1900 March 1 to 2100 February 28.
        # The year and day can be obtained by calling sla_CALYD or
        # sla_CLYD.
        #
        # 2  The Earth heliocentric 6-vector is mean equator and equinox
        # of date.  Position part, pv(1-3), is in AU;  velocity part,
        # pv(4-6), is in AU/sec.
        #
        # 3  Max/RMS errors 1950-2050:
        # 13/5 E-5 AU = 19200/7600 km in position
        # 47/26 E-10 AU/s = 0.0070/0.0039 km/s in speed
        #
        # 4  More accurate results are obtainable with the routines sla_EVP
        # and sla_EPV.
        #
        # Last revision:   5 April 2005
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        # License:
        # This program is free software; you can redistribute it and/or modify
        # it under the terms of the GNU General Public License as published by
        # the Free Software Foundation; either version 2 of the License, or
        # (at your option) any later version.
        #
        # This program is distributed in the hope that it will be useful,
        # but WITHOUT ANY WARRANTY; without even the implied warranty of
        # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        # GNU General Public License for more details.
        #
        # You should have received a copy of the GNU General Public License
        # along with this program (see SLA_CONDITIONS); if not, write to the
        # Free Software Foundation, Inc., 59 Temple Place, Suite 330,
        # Boston, MA  02111-1307  USA
        #
        # -
        TWOPI = 6.28318530718

        # Mean orbital speed of Earth, AU/s
        SPEED = 1.9913e-7

        # Mean Earth:EMB distance and speed, AU and AU/s
        REMB = 3.12e-5
        SEMB = 8.31e-11

        # Whole years & fraction of year, and years since J1900.0
        YI = float(iy - 1900)
        IY4 = np.mod(np.mod(iy, 4) + 4, 4)
        YF = (float(4 * (id - 1 / (IY4 + 1)) - IY4 - 2) + 4.0 * fd) / 1461.0
        T = YI + YF

        # Geometric mean longitude of Sun
        # (cf 4.881627938+6.283319509911*T MOD 2PI)
        ELM = np.mod(4.881628 + TWOPI * YF + 0.00013420 * T, TWOPI)

        # Mean longitude of perihelion
        GAMMA = 4.908230 + 3.0005e-4 * T

        # Mean anomaly
        EM = ELM - GAMMA

        # Mean obliquity
        EPS0 = 0.40931975 - 2.27e-6 * T

        # Eccentricity
        E = 0.016751 - 4.2e-7 * T
        ESQ = E * E

        # True anomaly
        V = EM + 2.0 * E * np.sin(EM) + 1.25 * ESQ * np.sin(2.0 * EM)

        # True ecliptic longitude
        ELT = V + GAMMA

        # True distance
        R = (1.0 - ESQ) / (1.0 + E * np.cos(V))

        # Moon's mean longitude
        ELMM = np.mod(4.72 + 83.9971 * T, TWOPI)

        # Useful functions
        COSELT = np.cos(ELT)
        SINEPS = np.sin(EPS0)
        COSEPS = np.cos(EPS0)
        W1 = -R * np.sin(ELT)
        W2 = -SPEED * (COSELT + E * np.cos(GAMMA))
        SELMM = np.sin(ELMM)
        CELMM = np.cos(ELMM)

        # Earth position and velocity
        pv = [None for _ in range(6)]
        pv[0] = -R * COSELT - REMB * CELMM
        pv[1] = (W1 - REMB * SELMM) * COSEPS
        pv[2] = W1 * SINEPS
        pv[3] = SPEED * (np.sin(ELT) + E * np.sin(GAMMA)) + SEMB * SELMM
        pv[4] = (W2 - SEMB * CELMM) * COSEPS
        pv[5] = W2 * SINEPS

        return pv

    @classmethod
    def epj2d(cls, epj):
        # +
        # - - - - - -
        # E P J 2 D
        # - - - - - -
        #
        # Conversion of Julian Epoch to Modified Julian Date (double precision)
        #
        # Given:
        # epj      dp       Julian Epoch
        #
        # The result is the Modified Julian Date (JD - 2400000.5).
        #
        # Reference:
        # Lieske,J.H., 1979. Astron.Astrophys.,73,282.
        #
        # P.T.Wallace   Starlink   February 1984
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        return 51544.5e0 + (epj - 2000e0) * 365.25e0

    @classmethod
    def pcd(cls, disco, x, y):
        # +
        # - - - -
        # P C D
        # - - - -
        #
        # Apply pincushion/barrel distortion to a tangent-plane [x,y].
        #
        # Given:
        # disco    d      pincushion/barrel distortion coefficient
        # x,y      d      tangent-plane coordinates
        #
        # Returned:
        # x,y      d      distorted coordinates
        #
        # Notes:
        #
        # 1)  The distortion is of the form RP = R*(1 + C*R**2), where R is
        # the radial distance from the tangent point, C is the disco
        # argument, and RP is the radial distance in the presence of
        # the distortion.
        #
        # 2)  For pincushion distortion, C is +ve;  for barrel distortion,
        # C is -ve.
        #
        # 3)  For x,y in units of one projection radius (in the case of
        # a photographic plate, the focal length), the following
        # disco values apply:
        #
        # Geometry          disco
        #
        # astrograph         0.0
        # Schmidt           -0.3333
        # AAT PF doublet  +147.069
        # AAT PF triplet  +178.585
        # AAT f/8          +21.20
        # JKT f/8          +13.32
        #
        # 4)  There is a companion routine, sla_UNPCD, which performs the
        # inverse operation.
        #
        # P.T.Wallace   Starlink   3 September 2000
        #
        # Copyright (C) 2000 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        F = 1e0 + disco * (x * x + y * y)
        x = x * F
        y = y * F

        return x, y

    @classmethod
    def refro(cls, zobs, hm, tdk, pmb, rh, wl, phi, tlr, eps):

        # +
        # - - - - - -
        # R E F R O
        # - - - - - -
        #
        # Atmospheric refraction for radio and optical/IR wavelengths.
        #
        # Given:
        # zobs    d  observed zenith distance of the source (radian)
        # hm      d  height of the observer above sea level (metre)
        # tdk     d  ambient temperature at the observer (K)
        # pmb     d  pressure at the observer (millibar)
        # rh      d  relative humidity at the observer (range 0-1)
        # wl      d  effective wavelength of the source (micrometre)
        # phi     d  latitude of the observer (radian, astronomical)
        # tlr     d  temperature lapse rate in the troposphere (K/metre)
        # eps     d  precision required to terminate iteration (radian)
        #
        # Returned:
        # ref     d  refraction: in vacuo ZD minus observed ZD (radian)
        #
        # Notes:
        #
        # 1  A suggested value for the tlr argument is 0.0065e0.  The
        # refraction is significantly affected by tlr, and if studies
        # of the local atmosphere have been carried out a better tlr
        # value may be available.  The sign of the supplied tlr value
        # is ignored.
        #
        # 2  A suggested value for the eps argument is 1e-8.  The result is
        # usually at least two orders of magnitude more computationally
        # precise than the supplied eps value.
        #
        # 3  The routine computes the refraction for zenith distances up
        # to and a little beyond 90 deg using the method of Hohenkerk
        # and Sinclair (NAO Technical Notes 59 and 63, subsequently adopted
        # in the Explanatory Supplement, 1992 edition - see section 3.281).
        #
        # 4  The code is a development of the optical/IR refraction subroutine
        # AREF of C.Hohenkerk (HMNAO, September 1984), with extensions to
        # support the radio case.  Apart from merely cosmetic changes, the
        # following modifications to the original HMNAO optical/IR refraction
        # code have been made:
        #
        # .  The angle arguments have been changed to radians.
        #
        # .  Any value of zobs is allowed (see note 6, below).
        #
        # .  Other argument values have been limited to safe values.
        #
        # .  Murray's values for the gas constants have been used
        # (Vectorial Astrometry, Adam Hilger, 1983).
        #
        # .  The numerical integration phase has been rearranged for
        # extra clarity.
        #
        # .  A better model for Ps(T) has been adopted (taken from
        # Gill, Atmosphere-Ocean Dynamics, Academic Press, 1982).
        #
        # .  More accurate expressions for Pwo have been adopted
        # (again from Gill 1982).
        #
        # .  The formula for the water vapour pressure, given the
        # saturation pressure and the relative humidity, is from
        # Crane (1976), expression 2.5.5.
        #
        # .  Provision for radio wavelengths has been added using
        # expressions devised by A.T.Sinclair, RGO (private
        # communication 1989).  The refractivity model currently
        # used is from J.M.Rueger, "Refractive Index Formulae for
        # Electronic Distance Measurement with Radio and Millimetre
        # Waves", in Unisurv Report S-68 (2002), School of Surveying
        # and Spatial Information Systems, University of New South
        # Wales, Sydney, Australia.
        #
        # .  The optical refractivity for dry air is from Resolution 3 of
        # the International Association of Geodesy adopted at the XXIIth
        # General Assembly in Birmingham, UK, 1999.
        #
        # .  Various small changes have been made to gain speed.
        #
        # 5  The radio refraction is chosen by specifying wl > 100 micrometres.
        # Because the algorithm takes no account of the ionosphere, the
        # accuracy deteriorates at low frequencies, below about 30 MHz.
        #
        # 6  Before use, the value of zobs is expressed in the range +/- pi.
        # If this ranged zobs is -ve, the result ref is computed from its
        # absolute value before being made -ve to match.  In addition, if
        # it has an absolute value greater than 93 deg, a fixed ref value
        # equal to the result for zobs = 93 deg is returned, appropriately
        # signed.
        #
        # 7  As in the original Hohenkerk and Sinclair algorithm, fixed values
        # of the water vapour polytrope exponent, the height of the
        # tropopause, and the height at which refraction is negligible are
        # used.
        #
        # 8  The radio refraction has been tested against work done by
        # Iain Coulson, JACH, (private communication 1995) for the
        # James Clerk Maxwell Telescope, Mauna Kea.  For typical conditions,
        # agreement at the 0.1 arcsec level is achieved for moderate ZD,
        # worsening to perhaps 0.5-1.0 arcsec at ZD 80 deg.  At hot and
        # humid sea-level sites the accuracy will not be as good.
        #
        # 9  It should be noted that the relative humidity rh is formally
        # defined in terms of "mixing ratio" rather than pressures or
        # densities as is often stated.  It is the mass of water per unit
        # mass of dry air divided by that for saturated air at the same
        # temperature and pressure (see Gill 1982).
        #
        # 10 The algorithm is designed for observers in the troposphere.  The
        # supplied temperature, pressure and lapse rate are assumed to be
        # for a point in the troposphere and are used to define a model
        # atmosphere with the tropopause at 11km altitude and a constant
        # temperature above that.  However, in practice, the refraction
        # values returned for stratospheric observers, at altitudes up to
        # 25km, are quite usable.
        #
        # Depends:  sla_DRANGE, sla__ATMT, sla__ATMS
        #
        # Last revision:   5 December 2005
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        #
        #
        # -

        #
        # Fixed parameters
        #

        # 93 degrees in radians
        D93 = 1.623156204e0
        # Universal gas constant
        GCR = 8314.32e0
        # Molecular weight of dry air
        DMD = 28.9644e0
        # Molecular weight of water vapour
        DMW = 18.0152e0
        # Mean Earth radius (metre)
        S = 6378120e0
        # Exponent of temperature dependence of water vapour pressure
        DELTA = 18.36e0
        # Height of tropopause (metre)
        HT = 11000e0
        # Upper limit for refractive effects (metre)
        HS = 80000e0
        # Numerical integration: maximum number of strips.
        ISMAX = 16384
        # The refraction integrand

        REFI = lambda DN, RDNDR: RDNDR / (DN + RDNDR)

        # Transform zobs into the normal range.
        ZOBS1 = cls.drange(zobs)
        ZOBS2 = np.minimum(np.abs(ZOBS1), D93)

        # Keep other arguments within safe bounds.
        HMOK = np.minimum(np.maximum(hm, -1e3), HS)
        TDKOK = np.minimum(np.maximum(tdk, 100e0), 500e0)
        PMBOK = np.minimum(np.maximum(pmb, 0e0), 10000e0)
        RHOK = np.minimum(np.maximum(rh, 0e0), 1e0)
        WLOK = np.maximum(wl, 0.1e0)
        ALPHA = np.minimum(np.maximum(np.abs(tlr), 0.001e0), 0.01e0)

        # Tolerance for iteration.
        TOL = np.minimum(np.maximum(np.abs(eps), 1e-12), 0.1e0) / 2e0

        # Decide whether optical/IR or radio case - switch at 100 microns.
        OPTIC = WLOK <= 100e0

        # Set up model atmosphere parameters defined at the observer.
        WLSQ = WLOK * WLOK
        GB = 9.784e0 * (1e0 - 0.0026e0 * np.cos(phi + phi) - 0.00000028e0 * HMOK)
        if OPTIC:
            A = (
                (287.6155e0 + (1.62887e0 + 0.01360e0 / WLSQ) / WLSQ)
                * 273.15e-6
                / 1013.25e0
            )
        else:
            A = 77.6890e-6

        GAMAL = (GB * DMD) / GCR
        GAMMA = GAMAL / ALPHA
        GAMM2 = GAMMA - 2e0
        DELM2 = DELTA - 2e0
        TDC = TDKOK - 273.15e0
        PSAT = 10e0 ** ((0.7859e0 + 0.03477e0 * TDC) / (1e0 + 0.00412e0 * TDC)) * (
            1e0 + PMBOK * (4.5e-6 + 6e-10 * TDC * TDC)
        )

        if PMBOK > 0e0:
            PWO = RHOK * PSAT / (1e0 - (1e0 - RHOK) * PSAT / PMBOK)
        else:
            PWO = 0e0

        W = PWO * (1e0 - DMW / DMD) * GAMMA / (DELTA - GAMMA)
        C1 = A * (PMBOK + W) / TDKOK
        if OPTIC:
            C2 = (A * W + 11.2684e-6 * PWO) / TDKOK
        else:
            C2 = (A * W + 6.3938e-6 * PWO) / TDKOK

        C3 = (GAMMA - 1e0) * ALPHA * C1 / TDKOK
        C4 = (DELTA - 1e0) * ALPHA * C2 / TDKOK
        if OPTIC:
            C5 = 0e0
            C6 = 0e0
        else:
            C5 = 375463e-6 * PWO / TDKOK
            C6 = C5 * DELM2 * ALPHA / (TDKOK * TDKOK)

        # Conditions at the observer.
        R0 = S + HMOK
        TEMPO, DN0, RDNDR0 = cls._atmt(
            R0, TDKOK, ALPHA, GAMM2, DELM2, C1, C2, C3, C4, C5, C6, R0
        )

        SK0 = DN0 * R0 * np.sin(ZOBS2)
        F0 = REFI(DN0, RDNDR0)

        # Conditions in the troposphere at the tropopause.
        RT = S + np.maximum(HT, HMOK)
        TT, DNT, RDNDRT = cls._atmt(
            R0, TDKOK, ALPHA, GAMM2, DELM2, C1, C2, C3, C4, C5, C6, RT
        )

        SINE = SK0 / (RT * DNT)
        ZT = np.arctan2(SINE, np.sqrt(np.maximum(1e0 - SINE * SINE, 0e0)))
        FT = REFI(DNT, RDNDRT)

        # Conditions in the stratosphere at the tropopause.
        DNTS, RDNDRP = cls._atms(RT, TT, DNT, GAMAL, RT)
        SINE = SK0 / (RT * DNTS)
        ZTS = np.arctan2(SINE, np.sqrt(np.maximum(1e0 - SINE * SINE, 0e0)))
        FTS = REFI(DNTS, RDNDRP)

        # Conditions at the stratosphere limit.
        RS = S + HS
        DNS, RDNDRS = cls._atms(RT, TT, DNT, GAMAL, RS)
        SINE = SK0 / (RS * DNS)
        ZS = np.arctan2(SINE, np.sqrt(np.maximum(1e0 - SINE * SINE, 0e0)))
        FS = REFI(DNS, RDNDRS)

        # Variable initialization to avoid compiler warning.
        REFT = 0e0

        # Integrate the refraction integral in two parts;  first in the
        # troposphere (K=1), then in the stratosphere (K=2).

        for K in range(2):
            # Initialize previous refraction to ensure at least two iterations.
            REFOLD = 1e0
            # Start off with 8 strips.
            IS = 8
            # Start Z, Z range, and start and end values.
            if K == 0:
                Z0 = ZOBS2
                ZRANGE = ZT - Z0
                FB = F0
                FF = FT
            else:
                Z0 = ZTS
                ZRANGE = ZS - Z0
                FB = FTS
                FF = FS

            # Sums of odd and even values.
            FO = 0e0
            FE = 0e0

            # First time through the loop we have to do every point.
            N = 1

            # Start of iteration loop (terminates at specified precision).
            LOOP = True
            while LOOP:

                # Strip width.
                H = ZRANGE / float(IS)

                # Initialize distance from Earth centre for quadrature pass.
                if K == 1:
                    R = R0
                else:
                    R = RT

                # One pass (no need to compute evens after first time).
                for I in range(0, IS, N):

                    # Sine of observed zenith distance.
                    SZ = np.sin(Z0 + H * float(I))

                    # Find R (to the nearest metre, maximum four iterations).
                    if SZ > 1e-20:
                        W = SK0 / SZ
                        RG = R
                        DR = 1e6
                        J = 0
                        while np.abs(DR) > 1e0 and J < 4:
                            J = J + 1
                            if K == 1:
                                TG, DN, RDNDR = cls._atmt(
                                    R0,
                                    TDKOK,
                                    ALPHA,
                                    GAMM2,
                                    DELM2,
                                    C1,
                                    C2,
                                    C3,
                                    C4,
                                    C5,
                                    C6,
                                    RG,
                                )

                            else:
                                DN, RDNDR = cls._atms(RT, TT, DNT, GAMAL, RG)

                            DR = (RG * DN - W) / (DN + RDNDR)
                            RG = RG - DR

                        R = RG

                    # Find the refractive index and integrand at R.
                    if K == 1:
                        T, DN, RDNDR = cls._atmt(
                            R0, TDKOK, ALPHA, GAMM2, DELM2, C1, C2, C3, C4, C5, C6, R
                        )
                    else:
                        DN, RDNDR = cls._atms(RT, TT, DNT, GAMAL, R)

                    F = REFI(DN, RDNDR)

                    # Accumulate odd and (first time only) even values.
                    if N == 1 and np.mod(I, 2) == 0:
                        FE = FE + F
                    else:
                        FO = FO + F

                # Evaluate the integrand using Simpson's Rule.
                REFP = H * (FB + 4e0 * FO + 2e0 * FE + FF) / 3e0

                # Has the required precision been achieved (or can't be)?
                if np.abs(REFP - REFOLD) > TOL and IS < ISMAX:

                    # No: prepare for next iteration.

                    # Save current value for convergence test.
                    REFOLD = REFP

                    # Double the number of strips.
                    IS = IS + IS

                    # Sum of all current values = sum of next pass's even values.
                    FE = FE + FO

                    # Prepare for new odd values.
                    FO = 0e0

                    # Skip even values next time.
                    N = 2
                else:

                    # Yes: save troposphere component and terminate the loop.
                    REFT = REFP if (K == 1) else REFT
                    LOOP = False
        # Result.
        ref = REFT + REFP
        ref = -ref if (ZOBS1 < 0e0) else ref

        return ref

    @classmethod
    def pa(cls, ha, dec, phi):
        # +
        # - - -
        # P A
        # - - -
        #
        # ha, Dec to Parallactic Angle (double precision)
        #
        # Given:
        # ha     d     hour angle in radians (geocentric apparent)
        # dec    d     declination in radians (geocentric apparent)
        # phi    d     observatory latitude in radians (geodetic)
        #
        # The result is in the range -pi to +pi
        #
        # Notes:
        #
        # 1)  The parallactic angle at a point in the sky is the position
        # angle of the vertical, i.e. the angle between the direction to
        # the pole and to the zenith.  In precise applications care must
        # be taken only to use geocentric apparent ha,Dec and to consider
        # separately the effects of atmospheric refraction and telescope
        # mount errors.
        #
        # 2)  At the pole a zero result is returned.
        #
        # P.T.Wallace   Starlink   16 August 1994
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        CP = np.cos(phi)
        SQSZ = CP * np.sin(ha)
        CQSZ = np.sin(phi) * np.cos(dec) - CP * np.sin(dec) * np.cos(ha)
        CQSZ = 1e0 if (SQSZ == 0e0 and CQSZ == 0e0) else CQSZ
        return np.arctan2(SQSZ, CQSZ)

    @classmethod
    def daf2r(cls, ideg, iamin, asec):
        # +
        # - - - - - -
        # D A F 2 R
        # - - - - - -
        #
        # Convert degrees, arcminutes, arcseconds to radians
        # (double precision)
        #
        # Given:
        # ideg        int       degrees
        # iamin       int       arcminutes
        # asec        dp        arcseconds
        #
        # Returned:
        # rad         dp        angle in radians
        # j           int       status:  0 = OK
        # 1 = ideg outside range 0-359
        # 2 = iamin outside range 0-59
        # 3 = asec outside range 0-59.999...
        #
        # Notes:
        # 1)  The result is computed even if any of the range checks
        # fail.
        # 2)  The sign must be dealt with outside this routine.
        #
        # P.T.Wallace   Starlink   23 August 1996
        #
        # Copyright (C) 1996 Rutherford Appleton Laboratory
        #
        #
        #
        # -

        # Arc seconds to radians

        AS2R = 0.484813681109535994e-5

        # Preset status
        j = 0

        # Validate arcsec, arcmin, deg
        j = 3 if (asec < 0e0 or asec >= 60e0) else j
        j = 2 if (iamin < 0 or iamin > 59) else j
        j = 1 if (ideg < 0 or ideg > 359) else j

        # Compute angle
        rad = AS2R * (60e0 * (60e0 * ideg + iamin) + asec)

        return rad, j

    @classmethod
    def dmoon(cls, date):
        # +
        # - - - - - -
        # D M O O N
        # - - - - - -
        #
        # Approximate geocentric position and velocity of the Moon
        # (double precision)
        #
        # Given:
        # DATE       D       TDB (loosely ET) as a Modified Julian Date
        # (JD-2400000.5)
        #
        # Returned:
        # PV         D(6)    Moon x,y,z,xdot,ydot,zdot, mean equator and
        # equinox of date (AU, AU/s)
        #
        # Notes:
        #
        # 1  This routine is a full implementation of the algorithm
        # published by Meeus (see reference).
        #
        # 2  Meeus quotes accuracies of 10 arcsec in longitude, 3 arcsec in
        # latitude and 0.2 arcsec in HP (equivalent to about 20 km in
        # distance).  Comparison with JPL DE200 over the interval
        # 1960-2025 gives RMS errors of 3.7 arcsec and 83 mas/hour in
        # longitude, 2.3 arcsec and 48 mas/hour in latitude, 11 km
        # and 81 mm/s in distance.  The maximum errors over the same
        # interval are 18 arcsec and 0.50 arcsec/hour in longitude,
        # 11 arcsec and 0.24 arcsec/hour in latitude, 40 km and 0.29 m/s
        # in distance.
        #
        # 3  The original algorithm is expressed in terms of the obsolete
        # timescale Ephemeris Time.  Either TDB or TT can be used, but
        # not UT without incurring significant errors (30 arcsec at
        # the present time) due to the Moon's 0.5 arcsec/sec movement.
        #
        # 4  The algorithm is based on pre IAU 1976 standards.  However,
        # the result has been moved onto the new (FK5) equinox, an
        # adjustment which is in any case much smaller than the
        # intrinsic accuracy of the procedure.
        #
        # 5  Velocity is obtained by a complete analytical differentiation
        # of the Meeus model.
        #
        # Reference:
        # Meeus, l'Astronomie, June 1984, p348.
        #
        # P.T.Wallace   Starlink   22 January 1998
        #
        # Copyright (C) 1998 Rutherford Appleton Laboratory
        #
        # License:
        # This program is free software; you can redistribute it and/or modify
        # it under the terms of the GNU General Public License as published by
        # the Free Software Foundation; either version 2 of the License, or
        # (at your option) any later version.
        #
        # This program is distributed in the hope that it will be useful,
        # but WITHOUT ANY WARRANTY; without even the implied warranty of
        # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        # GNU General Public License for more details.
        #
        # You should have received a copy of the GNU General Public License
        # along with this program (see SLA_CONDITIONS); if not, write to the
        # Free Software Foundation, Inc., 59 Temple Place, Suite 330,
        # Boston, MA  02111-1307  USA
        #
        # -

        DATE = np.zeros(6)
        PV = np.zeros(6)

        # Degrees, arcseconds and seconds of time to radians
        D2R = 0.0174532925199432957692369e0
        DAS2R = 4.848136811095359935899141e-6
        DS2R = 7.272205216643039903848712e-5

        # Seconds per Julian century (86400*36525)
        CJ = 3155760000e0

        # Julian epoch of B1950
        B1950 = 1949.9997904423e0

        # Earth equatorial radius in AU ( = 6378.137 / 149597870 )
        ERADAU = 4.2635212653763e-5

        #
        # Coefficients for fundamental arguments
        #
        # at J1900:  T**0, T**1, T**2, T**3
        # at epoch:  T**0, T**1
        #
        # Units are degrees for position and Julian centuries for time
        #

        # Moon's mean longitude
        ELP0 = 270.434164e0
        ELP1 = 481267.8831e0
        ELP2 = -0.001133e0
        ELP3 = 0.0000019e0

        # Sun's mean anomaly
        EM0 = 358.475833e0
        EM1 = 35999.0498e0
        EM2 = -0.000150e0
        EM3 = -0.0000033e0

        # Moon's mean anomaly
        EMP0 = 296.104608e0
        EMP1 = 477198.8491e0
        EMP2 = 0.009192e0
        EMP3 = 0.0000144e0

        # Moon's mean elongation
        D0 = 350.737486e0
        D1 = 445267.1142e0
        D2 = -0.001436e0
        D3 = 0.0000019e0

        # Mean distance of the Moon from its ascending node
        F0 = 11.250889e0
        F1 = 483202.0251e0
        F2 = -0.003211e0
        F3 = -0.0000003e0

        # Longitude of the Moon's ascending node
        OM0 = 259.183275e0
        OM1 = -1934.1420e0
        OM2 = 0.002078e0
        OM3 = 0.0000022e0

        # Coefficients for (dimensionless) E factor
        E1 = -0.002495e0
        E2 = -0.00000752e0

        # Coefficients for periodic variations etc
        PAC = 0.000233e0
        PA0 = 51.2e0
        PA1 = 20.2e0
        PBC = -0.001778e0
        PCC = 0.000817e0
        PDC = 0.002011e0
        PEC = 0.003964e0
        PE0 = 346.560e0
        PE1 = 132.870e0
        PE2 = -0.0091731e0

        PFC = 0.001964e0
        PGC = 0.002541e0
        PHC = 0.001964e0
        PIC = -0.024691e0
        PJC = -0.004328e0
        PJ0 = 275.05e0
        PJ1 = -2.30e0
        CW1 = 0.0004664e0
        CW2 = 0.0000754e0

        #
        # Coefficients for Moon position
        #
        # Tx(N)       = coefficient of L, B or P term (deg)
        # ITx(N,1-5)  = coefficients of M, M', D, F, E**n in argument
        #
        NL = 50
        NB = 45
        NP = 31

        TL = np.zeros(NL)
        TB = np.zeros(NB)
        TP = np.zeros(NP)

        ITL = np.zeros((5, NL))
        ITB = np.zeros((5, NB))
        ITP = np.zeros((5, NP))
        #
        # Longitude
        # M   M'  D   F   n
        TL[1] = 6.288750e0
        ITL[1] = [0, 1, 0, 0, 0]

        TL[2] = 1.274018e0
        ITL[2] = [0, -1, 2, 0, 0]

        TL[3] = 0.658309e0
        ITL[3] = [0, 0, 2, 0, 0]

        TL[4] = 0.213616e0
        ITL[4] = [0, 2, 0, 0, 0]

        TL[5] = -0.185596e0
        ITL[5] = [1, 0, 0, 0, 1]

        TL[6] = -0.114336e0
        ITL[6] = [0, 0, 0, 2, 0]

        TL[7] = 0.058793e0
        ITL[7] = [0, -2, 2, 0, 0]

        TL[8] = 0.057212e0
        ITL[8] = [-1, -1, 2, 0, 1]

        TL[9] = 0.053320e0
        ITL[9] = [0, 1, 2, 0, 0]

        TL[10] = 0.045874e0
        ITL[10] = [-1, 0, 2, 0, 1]

        TL[11] = 0.041024e0
        ITL[11] = [-1, 1, 0, 0, 1]

        TL[12] = -0.034718e0
        ITL[12] = [0, 0, 1, 0, 0]

        TL[13] = -0.030465e0
        ITL[13] = [1, 1, 0, 0, 1]

        TL[14] = 0.015326e0
        ITL[14] = [0, 0, 2, -2, 0]

        TL[15] = -0.012528e0
        ITL[15] = [0, 1, 0, 2, 0]

        TL[16] = -0.010980e0
        ITL[16] = [0, -1, 0, 2, 0]

        TL[17] = 0.010674e0
        ITL[17] = [0, -1, 4, 0, 0]

        TL[18] = 0.010034e0
        ITL[18] = [0, 3, 0, 0, 0]

        TL[19] = 0.008548e0
        ITL[19] = [0, -2, 4, 0, 0]

        TL[20] = -0.007910e0
        ITL[20] = [1, -1, 2, 0, 1]

        TL[21] = -0.006783e0
        ITL[21] = [1, 0, 2, 0, 1]

        TL[22] = 0.005162e0
        ITL[22] = [0, 1, -1, 0, 0]

        TL[23] = 0.005000e0
        ITL[23] = [1, 0, 1, 0, 1]

        TL[24] = 0.004049e0
        ITL[24] = [-1, 1, 2, 0, 1]

        TL[25] = 0.003996e0
        ITL[25] = [0, 2, 2, 0, 0]

        TL[26] = 0.003862e0
        ITL[26] = [0, 0, 4, 0, 0]

        TL[27] = 0.003665e0
        ITL[27] = [0, -3, 2, 0, 0]

        TL[28] = 0.002695e0
        ITL[28] = [-1, 2, 0, 0, 1]

        TL[29] = 0.002602e0
        ITL[29] = [0, 1, -2, -2, 0]

        TL[30] = 0.002396e0
        ITL[30] = [-1, -2, 2, 0, 1]

        TL[31] = -0.002349e0
        ITL[31] = [0, 1, 1, 0, 0]

        TL[32] = 0.002249e0
        ITL[32] = [-2, 0, 2, 0, 2]

        TL[33] = -0.002125e0
        ITL[33] = [1, 2, 0, 0, 1]

        TL[34] = -0.002079e0
        ITL[34] = [2, 0, 0, 0, 2]

        TL[35] = 0.002059e0
        ITL[35] = [-2, -1, 2, 0, 2]

        TL[36] = -0.001773e0
        ITL[36] = [0, 1, 2, -2, 0]

        TL[37] = -0.001595e0
        ITL[37] = [0, 0, 2, 2, 0]

        TL[38] = 0.001220e0
        ITL[38] = [-1, -1, 4, 0, 1]

        TL[39] = -0.001110e0
        ITL[39] = [0, 2, 0, 2, 0]

        TL[40] = 0.000892e0
        ITL[40] = [0, 1, -3, 0, 0]

        TL[41] = -0.000811e0
        ITL[41] = [1, 1, 2, 0, 1]

        TL[42] = 0.000761e0
        ITL[42] = [-1, -2, 4, 0, 1]

        TL[43] = 0.000717e0
        ITL[43] = [-2, 1, 0, 0, 2]

        TL[44] = 0.000704e0
        ITL[44] = [-2, 1, -2, 0, 2]

        TL[45] = 0.000693e0
        ITL[45] = [1, -2, 2, 0, 1]

        TL[46] = 0.000598e0
        ITL[46] = [-1, 0, 2, -2, 1]

        TL[47] = 0.000550e0
        ITL[47] = [0, 1, 4, 0, 0]

        TL[48] = 0.000538e0
        ITL[48] = [0, 4, 0, 0, 0]

        TL[49] = 0.000521e0
        ITL[49] = [-1, 0, 4, 0, 1]

        TL[50] = 0.000486e0
        ITL[50] = [0, 2, -1, 0, 0]

        #
        # Latitude
        # M   M'  D   F   n
        TB[1] = 5.128189e0
        ITB[1] = [0, 0, 0, 1, 0]

        TB[2] = 0.280606e0
        ITB[2] = [0, 1, 0, 1, 0]

        TB[3] = 0.277693e0
        ITB[3] = [0, 1, 0, -1, 0]

        TB[4] = 0.173238e0
        ITB[4] = [0, 0, 2, -1, 0]

        TB[5] = 0.055413e0
        ITB[5] = [0, -1, 2, 1, 0]

        TB[6] = 0.046272e0
        ITB[6] = [0, -1, 2, -1, 0]

        TB[7] = 0.032573e0
        ITB[7] = [0, 0, 2, 1, 0]

        TB[8] = 0.017198e0
        ITB[8] = [0, 2, 0, 1, 0]

        TB[9] = 0.009267e0
        ITB[9] = [0, 1, 2, -1, 0]

        TB[10] = 0.008823e0
        ITB[10] = [0, 2, 0, -1, 0]

        TB[11] = 0.008247e0
        ITB[11] = [-1, 0, 2, -1, 1]

        TB[12] = 0.004323e0
        ITB[12] = [0, -2, 2, -1, 0]

        TB[13] = 0.004200e0
        ITB[13] = [0, 1, 2, 1, 0]

        TB[14] = 0.003372e0
        ITB[14] = [-1, 0, -2, 1, 1]

        TB[15] = 0.002472e0
        ITB[15] = [-1, -1, 2, 1, 1]

        TB[16] = 0.002222e0
        ITB[16] = [-1, 0, 2, 1, 1]

        TB[17] = 0.002072e0
        ITB[17] = [-1, -1, 2, -1, 1]

        TB[18] = 0.001877e0
        ITB[18] = [-1, 1, 0, 1, 1]

        TB[19] = 0.001828e0
        ITB[19] = [0, -1, 4, -1, 0]

        TB[20] = -0.001803e0
        ITB[20] = [1, 0, 0, 1, 1]

        TB[21] = -0.001750e0
        ITB[21] = [0, 0, 0, 3, 0]

        TB[22] = 0.001570e0
        ITB[22] = [-1, 1, 0, -1, 1]

        TB[23] = -0.001487e0
        ITB[23] = [0, 0, 1, 1, 0]

        TB[24] = -0.001481e0
        ITB[24] = [1, 1, 0, 1, 1]

        TB[25] = 0.001417e0
        ITB[25] = [-1, -1, 0, 1, 1]

        TB[26] = 0.001350e0
        ITB[26] = [-1, 0, 0, 1, 1]

        TB[27] = 0.001330e0
        ITB[27] = [0, 0, -1, 1, 0]

        TB[28] = 0.001106e0
        ITB[28] = [0, 3, 0, 1, 0]

        TB[29] = 0.001020e0
        ITB[29] = [0, 0, 4, -1, 0]

        TB[30] = 0.000833e0
        ITB[30] = [0, -1, 4, 1, 0]

        TB[31] = 0.000781e0
        ITB[31] = [0, 1, 0, -3, 0]

        TB[32] = 0.000670e0
        ITB[32] = [0, -2, 4, 1, 0]

        TB[33] = 0.000606e0
        ITB[33] = [0, 0, 2, -3, 0]

        TB[34] = 0.000597e0
        ITB[34] = [0, 2, 2, -1, 0]

        TB[35] = 0.000492e0
        ITB[35] = [-1, 1, 2, -1, 1]

        TB[36] = 0.000450e0
        ITB[36] = [0, 2, -2, -1, 0]

        TB[37] = 0.000439e0
        ITB[37] = [0, 3, 0, -1, 0]

        TB[38] = 0.000423e0
        ITB[38] = [0, 2, 2, 1, 0]

        TB[39] = 0.000422e0
        ITB[39] = [0, -3, 2, -1, 0]

        TB[40] = -0.000367e0
        ITB[40] = [1, -1, 2, 1, 1]

        TB[41] = -0.000353e0
        ITB[41] = [1, 0, 2, 1, 1]

        TB[42] = 0.000331e0
        ITB[42] = [0, 0, 4, 1, 0]

        TB[43] = 0.000317e0
        ITB[43] = [-1, 1, 2, 1, 1]

        TB[44] = 0.000306e0
        ITB[44] = [-2, 0, 2, -1, 2]

        TB[45] = -0.000283e0
        ITB[45] = [0, 1, 0, 3, 0]

        #
        # Parallax
        # M   M'  D   F   n
        TP[1] = 0.950724e0
        ITP[1] = [0, 0, 0, 0, 0]

        TP[2] = 0.051818e0
        ITP[2] = [0, 1, 0, 0, 0]

        TP[3] = 0.009531e0
        ITP[3] = [0, -1, 2, 0, 0]

        TP[4] = 0.007843e0
        ITP[4] = [0, 0, 2, 0, 0]

        TP[5] = 0.002824e0
        ITP[5] = [0, 2, 0, 0, 0]

        TP[6] = 0.000857e0
        ITP[6] = [0, 1, 2, 0, 0]

        TP[7] = 0.000533e0
        ITP[7] = [-1, 0, 2, 0, 1]

        TP[8] = 0.000401e0
        ITP[8] = [-1, -1, 2, 0, 1]

        TP[9] = 0.000320e0
        ITP[9] = [-1, 1, 0, 0, 1]

        TP[10] = -0.000271e0
        ITP[10] = [0, 0, 1, 0, 0]

        TP[11] = -0.000264e0
        ITP[11] = [1, 1, 0, 0, 1]

        TP[12] = -0.000198e0
        ITP[12] = [0, -1, 0, 2, 0]

        TP[13] = 0.000173e0
        ITP[13] = [0, 3, 0, 0, 0]

        TP[14] = 0.000167e0
        ITP[14] = [0, -1, 4, 0, 0]

        TP[15] = -0.000111e0
        ITP[15] = [1, 0, 0, 0, 1]

        TP[16] = 0.000103e0
        ITP[16] = [0, -2, 4, 0, 0]

        TP[17] = -0.000084e0
        ITP[17] = [0, 2, -2, 0, 0]

        TP[18] = -0.000083e0
        ITP[18] = [1, 0, 2, 0, 1]

        TP[19] = 0.000079e0
        ITP[19] = [0, 2, 2, 0, 0]

        TP[20] = 0.000072e0
        ITP[20] = [0, 0, 4, 0, 0]

        TP[21] = 0.000064e0
        ITP[21] = [-1, 1, 2, 0, 1]

        TP[22] = -0.000063e0
        ITP[22] = [1, -1, 2, 0, 1]

        TP[23] = 0.000041e0
        ITP[23] = [1, 0, 1, 0, 1]

        TP[24] = 0.000035e0
        ITP[24] = [-1, 2, 0, 0, 1]

        TP[25] = -0.000033e0
        ITP[25] = [0, 3, -2, 0, 0]

        TP[26] = -0.000030e0
        ITP[26] = [0, 1, 1, 0, 0]

        TP[27] = -0.000029e0
        ITP[27] = [0, 0, -2, 2, 0]

        TP[28] = -0.000029e0
        ITP[28] = [1, 2, 0, 0, 1]

        TP[29] = 0.000026e0
        ITP[29] = [-2, 0, 2, 0, 2]

        TP[30] = -0.000023e0
        ITP[30] = [0, 1, -2, 2, 0]

        TP[31] = 0.000019e0
        ITP[31] = [-1, -1, 4, 0, 1]

        # Centuries since J1900
        T = (DATE - 15019.5e0) / 36525e0

        #
        # Fundamental arguments (radians) and derivatives (radians per
        # Julian century) for the current epoch
        #

        # Moon's mean longitude
        ELP = D2R * np.mod(ELP0 + (ELP1 + (ELP2 + ELP3 * T) * T) * T, 360e0)
        DELP = D2R * (ELP1 + (2e0 * ELP2 + 3e0 * ELP3 * T) * T)

        # Sun's mean anomaly
        EM = D2R * np.mod(EM0 + (EM1 + (EM2 + EM3 * T) * T) * T, 360e0)
        DEM = D2R * (EM1 + (2e0 * EM2 + 3e0 * EM3 * T) * T)

        # Moon's mean anomaly
        EMP = D2R * np.mod(EMP0 + (EMP1 + (EMP2 + EMP3 * T) * T) * T, 360e0)
        DEMP = D2R * (EMP1 + (2e0 * EMP2 + 3e0 * EMP3 * T) * T)

        # Moon's mean elongation
        D = D2R * np.mod(D0 + (D1 + (D2 + D3 * T) * T) * T, 360e0)
        DD = D2R * (D1 + (2e0 * D2 + 3e0 * D3 * T) * T)

        # Mean distance of the Moon from its ascending node
        F = D2R * np.mod(F0 + (F1 + (F2 + F3 * T) * T) * T, 360e0)
        DF = D2R * (F1 + (2e0 * F2 + 3e0 * F3 * T) * T)

        # Longitude of the Moon's ascending node
        OM = D2R * np.mod(OM0 + (OM1 + (OM2 + OM3 * T) * T) * T, 360e0)
        DOM = D2R * (OM1 + (2e0 * OM2 + 3e0 * OM3 * T) * T)
        SINOM = np.sin(OM)
        COSOM = np.cos(OM)
        DOMCOM = DOM * COSOM

        # Add the periodic variations
        THETA = D2R * (PA0 + PA1 * T)
        WA = np.sin(THETA)
        DWA = D2R * PA1 * np.cos(THETA)
        THETA = D2R * (PE0 + (PE1 + PE2 * T) * T)
        WB = PEC * np.sin(THETA)
        DWB = D2R * PEC * (PE1 + 2e0 * PE2 * T) * np.cos(THETA)
        ELP = ELP + D2R * (PAC * WA + WB + PFC * SINOM)
        DELP = DELP + D2R * (PAC * DWA + DWB + PFC * DOMCOM)
        EM = EM + D2R * PBC * WA
        DEM = DEM + D2R * PBC * DWA
        EMP = EMP + D2R * (PCC * WA + WB + PGC * SINOM)
        DEMP = DEMP + D2R * (PCC * DWA + DWB + PGC * DOMCOM)
        D = D + D2R * (PDC * WA + WB + PHC * SINOM)
        DD = DD + D2R * (PDC * DWA + DWB + PHC * DOMCOM)
        WOM = OM + D2R * (PJ0 + PJ1 * T)
        DWOM = DOM + D2R * PJ1
        SINWOM = np.sin(WOM)
        COSWOM = np.cos(WOM)
        F = F + D2R * (WB + PIC * SINOM + PJC * SINWOM)
        DF = DF + D2R * (DWB + PIC * DOMCOM + PJC * DWOM * COSWOM)

        # E-factor, and square
        E = 1e0 + (E1 + E2 * T) * T
        DE = E1 + 2e0 * E2 * T
        ESQ = E * E
        DESQ = 2e0 * E * DE

        #
        # Series expansions
        #

        # Longitude
        V = 0e0
        DV = 0e0
        for N in range(NL, 1, -1):
            COEFF = TL[N]
            EMN = ITL[1][N]
            EMPN = ITL[2][N]
            DN = ITL[3][N]
            FN = ITL[4][N]
            I = ITL[5][N]
            if I == 0:
                EN = 1e0
                DEN = 0e0
            elif I == 1:
                EN = E
                DEN = DE
            else:
                EN = ESQ
                DEN = DESQ

            THETA = EMN * EM + EMPN * EMP + DN * D + FN * F
            DTHETA = EMN * DEM + EMPN * DEMP + DN * DD + FN * DF
            FTHETA = np.sin(THETA)
            V = V + COEFF * FTHETA * EN
            DV = DV + COEFF * (np.cos(THETA) * DTHETA * EN + FTHETA * DEN)

        EL = ELP + D2R * V
        DEL = (DELP + D2R * DV) / CJ

        # Latitude
        V = 0e0
        DV = 0e0
        for N in range(NB, 1, -1):
            COEFF = TB[N]
            EMN = ITB[1][N]
            EMPN = ITB[2][N]
            DN = ITB[3][N]
            FN = ITB[4][N]
            I = ITB[5][N]
            if I == 0:
                EN = 1e0
                DEN = 0e0
            elif I == 1:
                EN = E
                DEN = DE
            else:
                EN = ESQ
                DEN = DESQ

            THETA = EMN * EM + EMPN * EMP + DN * D + FN * F
            DTHETA = EMN * DEM + EMPN * DEMP + DN * DD + FN * DF
            FTHETA = np.sin(THETA)
            V = V + COEFF * FTHETA * EN
            DV = DV + COEFF * (np.cos(THETA) * DTHETA * EN + FTHETA * DEN)

        BF = 1e0 - CW1 * COSOM - CW2 * COSWOM
        DBF = CW1 * DOM * SINOM + CW2 * DWOM * SINWOM
        B = D2R * V * BF
        DB = D2R * (DV * BF + V * DBF) / CJ

        # Parallax
        V = 0e0
        DV = 0e0
        for N in range(NP, 1, -1):
            COEFF = TP[N]
            EMN = ITP[1][N]
            EMPN = ITP[2][N]
            DN = ITP[3][N]
            FN = ITP[4][N]
            I = ITP[5][N]
            if I == 0:
                EN = 1e0
                DEN = 0e0
            elif I == 1:
                EN = E
                DEN = DE
            else:
                EN = ESQ
                DEN = DESQ

            THETA = EMN * EM + EMPN * EMP + DN * D + FN * F
            DTHETA = EMN * DEM + EMPN * DEMP + DN * DD + FN * DF
            FTHETA = np.cos(THETA)
            V = V + COEFF * FTHETA * EN
            DV = DV + COEFF * (-np.sin(THETA) * DTHETA * EN + FTHETA * DEN)

        P = D2R * V
        DP = D2R * DV / CJ

        #
        # Transformation into final form
        #

        # Parallax to distance (AU, AU/sec)
        SP = np.sin(P)
        R = ERADAU / SP
        DR = -R * DP * np.cos(P) / SP

        # Longitude, latitude to x,y,z (AU)
        SEL = np.sin(EL)
        CEL = np.cos(EL)
        SB = np.sin(B)
        CB = np.cos(B)
        RCB = R * CB
        RBD = R * DB
        W = RBD * SB - CB * DR
        X = RCB * CEL
        Y = RCB * SEL
        Z = R * SB
        XD = -Y * DEL - W * CEL
        YD = X * DEL - W * SEL
        ZD = RBD * CB + SB * DR

        # Julian centuries since J2000
        T = (DATE - 51544.5e0) / 36525e0

        # Fricke equinox correction
        EPJ = 2000e0 + T * 100e0
        EQCOR = DS2R * (0.035e0 + 0.00085e0 * (EPJ - B1950))

        # Mean obliquity (IAU 1976)
        EPS = DAS2R * (
            84381.448e0 + (-46.8150e0 + (-0.00059e0 + 0.001813e0 * T) * T) * T
        )

        # To the equatorial system, mean of date, FK5 system
        SINEPS = np.sin(EPS)
        COSEPS = np.cos(EPS)
        ES = EQCOR * SINEPS
        EC = EQCOR * COSEPS
        PV[1] = X - EC * Y + ES * Z
        PV[2] = EQCOR * X + Y * COSEPS - Z * SINEPS
        PV[3] = Y * SINEPS + Z * COSEPS
        PV[4] = XD - EC * YD + ES * ZD
        PV[5] = EQCOR * XD + YD * COSEPS - ZD * SINEPS
        PV[6] = YD * SINEPS + ZD * COSEPS

        return PV
    
    @classmethod
    def dpav(cls, V1, V2):
        # +
        # - - - - -
        # D P A V
        # - - - - -
        #
        # Position angle of one celestial direction with respect to another.
        #
        # (double precision)
        #
        # Given:
        # V1    d(3)    direction cosines of one point
        # V2    d(3)    direction cosines of the other point
        #
        # (The coordinate frames correspond to RA,Dec, Long,Lat etc.)
        #
        # The result is the bearing (position angle), in radians, of point
        # V2 with respect to point V1.  It is in the range +/- pi.  The
        # sense is such that if V2 is a small distance east of V1, the
        # bearing is about +pi/2.  Zero is returned if the two points
        # are coincident.
        #
        # V1 and V2 need not be unit vectors.
        #
        # The routine sla_DBEAR performs an equivalent function except
        # that the points are specified in the form of spherical
        # coordinates.
        #
        # Last revision:   16 March 2005
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        # -

        V1 = np.zeros(3)
        V2 = np.zeros(3)

        # The unit vector to point 1.
        X1 = V1[1]
        Y1 = V1[2]
        Z1 = V1[3]
        W = np.sqrt(X1 * X1 + Y1 * Y1 + Z1 * Z1)
        if W != 0e0:
            X1 = X1 / W
            Y1 = Y1 / W
            Z1 = Z1 / W

        # The vector to point 2.
        X2 = V2[1]
        Y2 = V2[2]
        Z2 = V2[3]

        # Position angle.
        SQ = Y2 * X1 - X2 * Y1
        CQ = Z2 * (X1 * X1 + Y1 * Y1) - Z1 * (X2 * X1 + Y2 * Y1)
        CQ = 1e0 if (SQ == 0e0 and CQ == 0e0) else CQ
        return np.arctan2(SQ, CQ)

    @classmethod
    def flotin(cls, STRING, NSTRT, RESLT=0.):
        # +
        # - - - - - - -
        # F L O T I N
        # - - - - - - -
        #
        # Convert free-format input into single precision floating point
        #
        # Given:
        # STRING     c     string containing number to be decoded
        # NSTRT      i     pointer to where decoding is to start
        # RESLT      r     current value of result
        #
        # Returned:
        # NSTRT      i      advanced to next number
        # RESLT      r      result
        # JFLAG      i      status: -1 = -OK, 0 = +OK, 1 = null, 2 = error
        #
        # Called:  sla_DFLTIN
        #
        # Notes:
        #
        # 1     The reason FLOTIN has separate OK status values for +
        # and - is to enable minus zero to be detected.   This is
        # of crucial importance when decoding mixed-radix numbers.
        # For example, an angle expressed as deg, arcmin, arcsec
        # may have a leading minus sign but a zero degrees field.
        #
        # 2     A TAB is interpreted as a space, and lowercase characters
        # are interpreted as uppercase.
        #
        # 3     The basic format is the sequence of fields #^.^@#^, where
        # # is a sign character + or -, ^ means a string of decimal
        # digits, and @, which indicates an exponent, means D or E.
        # Various combinations of these fields can be omitted, and
        # embedded blanks are permissible in certain places.
        #
        # 4     Spaces:
        #
        # .  Leading spaces are ignored.
        #
        # .  Embedded spaces are allowed only after +, -, D or E,
        # and after the decomal point if the first sequence of
        # digits is absent.
        #
        # .  Trailing spaces are ignored;  the first signifies
        # end of decoding and subsequent ones are skipped.
        #
        # 5     eelimiters:
        #
        # .  Any character other than +,-,0-9,.,D,E or space may be
        # used to signal the end of the number and terminate
        # decoding.
        #
        # .  Comma is recognized by FLOTIN as a special case;  it
        # is skipped, leaving the pointer on the next character.
        # See 13, below.
        #
        # 6     Both signs are optional.  The default is +.
        #
        # 7     The mantissa ^.^ defaults to 1.
        #
        # 8     The exponent @#^ defaults to e0.
        #
        # 9     The strings of decimal digits may be of any length.
        #
        # 10    The decimal point is optional for whole numbers.
        #
        # 11    A "null result" occurs when the string of characters being
        # decoded does not begin with +,-,0-9,.,D or E, or consists
        # entirely of spaces.  When this condition is detected, JFLAG
        # is set to 1 and RESLT is left untouched.
        #
        # 12    NSTRT = 1 for the first character in the string.
        #
        # 13    On return from FLOTIN, NSTRT is set ready for the next
        # decode - following trailing blanks and any comma.  If a
        # delimiter other than comma is being used, NSTRT must be
        # incremented before the next call to FLOTIN, otherwise
        # all subsequent calls will return a null result.
        #
        # 14    errors (JFLAG=2) occur when:
        #
        # .  a +, -, D or E is left unsatisfied;  or
        #
        # .  the decimal point is present without at least
        # one decimal digit before or after it;  or
        #
        # .  an exponent more than 100 has been presented.
        #
        # 15    When an error has been detected, NSTRT is left
        # pointing to the character following the last
        # one used before the error came to light.  This
        # may be after the point at which a more sophisticated
        # program could have detected the error.  For example,
        # FLOTIN does not detect that '1E999' is unacceptable
        # (on a computer where this is so) until the entire number
        # has been decoded.
        #
        # 16    Certain highly unlikely combinations of mantissa &
        # exponent can cause arithmetic faults during the
        # decode, in some cases despite the fact that they
        # together could be construed as a valid number.
        #
        # 17    eecoding is left to right, one pass.
        #
        # 18    See also eFLTIN and INTIN
        #
        # P.T.Wallace   Starlink   23 November 1995
        #
        # Copyright (C) 1995 Rutherford Appleton Laboratory
        # -

        # Call the double precision version
        STRING = re.sub(r"\s+", r" ", STRING).replace("D","e").replace("E","e")
        STRARR = [xx.strip() for xx in re.split(r"[^+-0-9.e]", STRING)]
        JFLAG = 0
        try:
            RESLT = float(STRARR[NSTRT])
        except ValueError:
            JFLAG = 2
        NSTRT += 1
        return NSTRT, RESLT, JFLAG
    
    @classmethod
    def dfltin(cls, STRING, NSTRT, DRESLT=0.):
        # +
        # - - - - - - -
        # D F L T I N
        # - - - - - - -
        #
        # Convert free-format input into double precision floating point
        #
        # Given:
        # STRING     c     string containing number to be decoded
        # NSTRT      i     pointer to where decoding is to start
        # DRESLT     d     current value of result
        #
        # Returned:
        # NSTRT      i      advanced to next number
        # DRESLT     d      result
        # JFLAG      i      status: -1 = -OK, 0 = +OK, 1 = null, 2 = error
        #
        # Notes:
        #
        # 1     The reason eFLTIN has separate OK status values for +
        # and - is to enable minus zero to be detected.   This is
        # of crucial importance when decoding mixed-radix numbers.
        # For example, an angle expressed as deg, arcmin, arcsec
        # may have a leading minus sign but a zero degrees field.
        #
        # 2     A TAB is interpreted as a space, and lowercase characters
        # are interpreted as uppercase.
        #
        # 3     The basic format is the sequence of fields #^.^@#^, where
        # # is a sign character + or -, ^ means a string of decimal
        # digits, and @, which indicates an exponent, means D or E.
        # Various combinations of these fields can be omitted, and
        # embedded blanks are permissible in certain places.
        #
        # 4     Spaces:
        #
        # .  Leading spaces are ignored.
        #
        # .  Embedded spaces are allowed only after +, -, D or E,
        # and after the decomal point if the first sequence of
        # digits is absent.
        #
        # .  Trailing spaces are ignored;  the first signifies
        # end of decoding and subsequent ones are skipped.
        #
        # 5     eelimiters:
        #
        # .  Any character other than +,-,0-9,.,D,E or space may be
        # used to signal the end of the number and terminate
        # decoding.
        #
        # .  Comma is recognized by DFLTIN as a special case;  it
        # is skipped, leaving the pointer on the next character.
        # See 13, below.
        #
        # 6     Both signs are optional.  The default is +.
        #
        # 7     The mantissa ^.^ defaults to 1.
        #
        # 8     The exponent @#^ defaults to e0.
        #
        # 9     The strings of decimal digits may be of any length.
        #
        # 10    The decimal point is optional for whole numbers.
        #
        # 11    A "null result" occurs when the string of characters being
        # decoded does not begin with +,-,0-9,.,D or E, or consists
        # entirely of spaces.  When this condition is detected, JFLAG
        # is set to 1 and DRESLT is left untouched.
        #
        # 12    NSTRT = 1 for the first character in the string.
        #
        # 13    On return from eFLTIN, NSTRT is set ready for the next
        # decode - following trailing blanks and any comma.  If a
        # delimiter other than comma is being used, NSTRT must be
        # incremented before the next call to DFLTIN, otherwise
        # all subsequent calls will return a null result.
        #
        # 14    errors (JFLAG=2) occur when:
        #
        # .  a +, -, D or E is left unsatisfied;  or
        #
        # .  the decimal point is present without at least
        # one decimal digit before or after it;  or
        #
        # .  an exponent more than 100 has been presented.
        #
        # 15    When an error has been detected, NSTRT is left
        # pointing to the character following the last
        # one used before the error came to light.  This
        # may be after the point at which a more sophisticated
        # program could have detected the error.  For example,
        # DFLTIN does not detect that '1D999' is unacceptable
        # (on a computer where this is so) until the entire number
        # has been decoded.
        #
        # 16    Certain highly unlikely combinations of mantissa &
        # exponent can cause arithmetic faults during the
        # decode, in some cases despite the fact that they
        # together could be construed as a valid number.
        #
        # 17    eecoding is left to right, one pass.
        #
        # 18    See also FLOTIN and INTIN
        #
        # Called:  sla__IDCHF
        #
        # P.T.Wallace   Starlink   18 March 1999
        #
        # Copyright (C) 1999 Rutherford Appleton Laboratory
        #
        # License:
        # This program is free software; you can redistribute it and/or modify
        # it under the terms of the GNU General Public License as published by
        # the Free Software Foundation; either version 2 of the License, or
        # (at your option) any later version.
        #
        # This program is distributed in the hope that it will be useful,
        # but WITHOUT ANY WARRANTY; without even the implied warranty of
        # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        # GNU General Public License for more details.
        #
        # You should have received a copy of the GNU General Public License
        # along with this program (see SLA_CONDITIONS); if not, write to the
        # Free Software Foundation, Inc., 59 Temple Place, Suite 330,
        # Boston, MA  02111-1307  USA
        #
        # -
        
        STRING = re.sub(r"\s+", r" ", STRING).replace("D", "e").replace("E", "e")
        STRARR = [xx.strip() for xx in re.split(r"[^+-0-9.e]", STRING)]
        JFLAG = 0
        try:
            DRESLT = float(STRARR[NSTRT])
        except ValueError:
            JFLAG = 2
        NSTRT += 1
        return NSTRT, DRESLT, JFLAG

    @classmethod
    def cs2c(cls, A, B):
        # +
        # - - - - -
        # C S 2 C
        # - - - - -
        #
        # Spherical coordinates to direction cosines (single precision)
        #
        # Given:
        # A,B      real      spherical coordinates in radians
        # (RA,Dec), (long,lat) etc.
        #
        # Returned:
        # V        real(3)   x,y,z unit vector
        #
        # The spherical coordinates are longitude (+ve anticlockwise looking
        # from the +ve latitude pole) and latitude.  The Cartesian coordinates
        # are right handed, with the x axis at zero longitude and latitude, and
        # the z axis at the +ve latitude pole.
        #
        # Last revision:   22 July 2004
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        # License:
        # This program is free software; you can redistribute it and/or modify
        # it under the terms of the GNU General Public License as published by
        # the Free Software Foundation; either version 2 of the License, or
        # (at your option) any later version.
        #
        # This program is distributed in the hope that it will be useful,
        # but WITHOUT ANY WARRANTY; without even the implied warranty of
        # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        # GNU General Public License for more details.
        #
        # You should have received a copy of the GNU General Public License
        # along with this program (see SLA_CONDITIONS); if not, write to the
        # Free Software Foundation, Inc., 59 Temple Place, Suite 330,
        # Boston, MA  02111-1307  USA
        #
        # -
    
        V = np.zeros(3)
    
        COSB = np.cos(B)
    
        V[1] = np.cos(A) * COSB
        V[2] = np.sin(A) * COSB
        V[3] = np.sin(B)
    
        return V

    @classmethod
    def cr2af(cls, NDP, ANGLE):
        # +
        # - - - - - -
        # C R 2 A F
        # - - - - - -
        #
        # Convert an angle in radians into degrees, arcminutes, arcseconds
        # (single precision)
        #
        # Given:
        # NDP       int      number of decimal places of arcseconds
        # ANGLE     real     angle in radians
        #
        # Returned:
        # SIGN      char     '+' or '-'
        # IDMSF     int(4)   degrees, arcminutes, arcseconds, fraction
        #
        # Notes:
        #
        # 1)  NeP less than zero is interpreted as zero.
        #
        # 2)  The largest useful value for NeP is determined by the size of
        # ANGLE, the format of REAL floating-point numbers on the target
        # machine, and the risk of overflowing IDMSF(4).  On some
        # architectures, for ANGLE up to 2pi, the available floating-
        # point precision corresponds roughly to NDP=3.  This is well
        # below the ultimate limit of NDP=9 set by the capacity of a
        # typical 32-bit IDMSF(4).
        #
        # 3)  The absolute value of ANGLe may exceed 2pi.  In cases where it
        # does not, it is up to the caller to test for and handle the
        # case where ANGLE is very nearly 2pi and rounds up to 360 deg,
        # by testing for IDMSF(1)=360 and setting IDMSF(1-4) to zero.
        #
        # Called:  sla_CD2TF
        #
        # Last revision:   26 December 2004
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        # License:
        # This program is free software; you can redistribute it and/or modify
        # it under the terms of the GNU General Public License as published by
        # the Free Software Foundation; either version 2 of the License, or
        # (at your option) any later version.
        #
        # This program is distributed in the hope that it will be useful,
        # but WITHOUT ANY WARRANTY; without even the implied warranty of
        # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        # GNU General Public License for more details.
        #
        # You should have received a copy of the GNU General Public License
        # along with this program (see SLA_CONDITIONS); if not, write to the
        # Free Software Foundation, Inc., 59 Temple Place, Suite 330,
        # Boston, MA  02111-1307  USA
        #
        # -

        # Hours to degrees * radians to turns
        F = 15.0 / 6.283185307179586476925287

        # Scale then use days to h,m,s routine
        SIGN, IDMSF = cls.cd2tf(NDP, ANGLE * F)

        return SIGN, IDMSF

    @classmethod
    def dd2tf(cls, NDP, DAYS):
        # +
        # - - - - - -
        # D D 2 T F
        # - - - - - -
        #
        # Convert an interval in days into hours, minutes, seconds
        # (double precision)
        #
        # Given:
        # NDP      i      number of decimal places of seconds
        # DAYS     d      interval in days
        #
        # Returned:
        # SIGN     c      '+' or '-'
        # IHMSF    i(4)   hours, minutes, seconds, fraction
        #
        # Notes:
        #
        # 1)  NeP less than zero is interpreted as zero.
        #
        # 2)  The largest useful value for NeP is determined by the size
        # of DAYS, the format of DOUBLE PRECISION floating-point numbers
        # on the target machine, and the risk of overflowing IHMSF(4).
        # On some architectures, for DAYS up to 1D0, the available
        # floating-point precision corresponds roughly to NDP=12.
        # However, the practical limit is NDP=9, set by the capacity of
        # a typical 32-bit IHMSF(4).
        #
        # 3)  The absolute value of eAYS may exceed 1e0.  In cases where it
        # does not, it is up to the caller to test for and handle the
        # case where DAYS is very nearly 1D0 and rounds up to 24 hours,
        # by testing for IHMSF(1)=24 and setting IHMSF(1-4) to zero.
        #
        # Last revision:   26 December 2004
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        # -
    
        IHMSF = np.zeros(4, dtype=int)
    
        # Days to seconds
        D2S = 86400e0
    
        # Handle sign
        SIGN = "+" if DAYS >= 0e0 else "-"
        # Field units in terms of least significant figure
        NRS = 1
        for _ in range(NDP):
            NRS = NRS * 10
    
        RS = NRS
        RM = RS * 60e0
        RH = RM * 60e0
    
        # Round interval and express in smallest units required
        A = np.rint(RS * D2S * np.abs(DAYS))
    
        # Separate into fields
        AH = np.trunc(A / RH)
        A = A - AH * RH
        AM = np.trunc(A / RM)
        A = A - AM * RM
        AS = np.trunc(A / RS)
        AF = A - AS * RS
    
        # Return results
        IHMSF[1] = np.maximum(np.rint(AH), 0)
        IHMSF[2] = np.maximum(np.minimum(np.rint(AM), 59), 0)
        IHMSF[3] = np.maximum(np.minimum(np.rint(AS), 59), 0)
        IHMSF[4] = np.maximum(np.rint(np.minimum(AF, RS - 1e0)), 0)
    
        return SIGN, IHMSF

    @classmethod
    def cd2tf(cls, NDP, DAYS):
        # +
        # - - - - - -
        # C D 2 T F
        # - - - - - -
        #
        # Convert an interval in days into hours, minutes, seconds
        #
        # (single precision)
        #
        # Given:
        # NDP       int      number of decimal places of seconds
        # DAYS      real     interval in days
        #
        # Returned:
        # SIGN      char     '+' or '-'
        # IHMSF     int(4)   hours, minutes, seconds, fraction
        #
        # Notes:
        #
        # 1)  NeP less than zero is interpreted as zero.
        #
        # 2)  The largest useful value for NeP is determined by the size of
        # DAYS, the format of REAL floating-point numbers on the target
        # machine, and the risk of overflowing IHMSF(4).  On some
        # architectures, for DAYS up to 1.0, the available floating-
        # point precision corresponds roughly to NDP=3.  This is well
        # below the ultimate limit of NDP=9 set by the capacity of a
        # typical 32-bit IHMSF(4).
        #
        # 3)  The absolute value of eAYS may exceed 1.0.  In cases where it
        # does not, it is up to the caller to test for and handle the
        # case where DAYS is very nearly 1.0 and rounds up to 24 hours,
        # by testing for IHMSF(1)=24 and setting IHMSF(1-4) to zero.
        #
        # Called:  sla_DD2TF
        #
        # Last revision:   26 December 2004
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        # License:
        # This program is free software; you can redistribute it and/or modify
        # it under the terms of the GNU General Public License as published by
        # the Free Software Foundation; either version 2 of the License, or
        # (at your option) any later version.
        #
        # This program is distributed in the hope that it will be useful,
        # but WITHOUT ANY WARRANTY; without even the implied warranty of
        # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        # GNU General Public License for more details.
        #
        # You should have received a copy of the GNU General Public License
        # along with this program (see SLA_CONDITIONS); if not, write to the
        # Free Software Foundation, Inc., 59 Temple Place, Suite 330,
        # Boston, MA  02111-1307  USA
        #
        # -

        # Call double precision version
        SIGN, IHMSF = cls.dd2tf(NDP, DAYS)

        return SIGN, IHMSF
    
    @classmethod
    def intin(cls, STRING, NSTRT, IRESLT):
        # +
        # - - - - - -
        # I N T I N
        # - - - - - -
        #
        # Convert free-format input into an integer
        #
        # Given:
        # STRING     c     string containing number to be decoded
        # NSTRT      i     pointer to where decoding is to start
        # IRESLT     i     current value of result
        #
        # Returned:
        # NSTRT      i      advanced to next number
        # IRESLT     i      result
        # JFLAG      i      status: -1 = -OK, 0 = +OK, 1 = null, 2 = error
        #
        # Called:  sla__IDCHI
        #
        # Notes:
        #
        # 1     The reason INTIN has separate OK status values for +
        # and - is to enable minus zero to be detected.   This is
        # of crucial importance when decoding mixed-radix numbers.
        # For example, an angle expressed as deg, arcmin, arcsec
        # may have a leading minus sign but a zero degrees field.
        #
        # 2     A TAB is interpreted as a space.
        #
        # 3     The basic format is the sequence of fields #^, where
        # # is a sign character + or -, and ^ means a string of
        # decimal digits.
        #
        # 4     Spaces:
        #
        # .  Leading spaces are ignored.
        #
        # .  Spaces between the sign and the number are allowed.
        #
        # .  Trailing spaces are ignored;  the first signifies
        # end of decoding and subsequent ones are skipped.
        #
        # 5     eelimiters:
        #
        # .  Any character other than +,-,0-9 or space may be
        # used to signal the end of the number and terminate
        # decoding.
        #
        # .  Comma is recognized by INTIN as a special case;  it
        # is skipped, leaving the pointer on the next character.
        # See 9, below.
        #
        # 6     The sign is optional.  The default is +.
        #
        # 7     A "null result" occurs when the string of characters being
        # decoded does not begin with +,- or 0-9, or consists
        # entirely of spaces.  When this condition is detected, JFLAG
        # is set to 1 and IRESLT is left untouched.
        #
        # 8     NSTRT = 1 for the first character in the string.
        #
        # 9     On return from INTIN, NSTRT is set ready for the next
        # decode - following trailing blanks and any comma.  If a
        # delimiter other than comma is being used, NSTRT must be
        # incremented before the next call to INTIN, otherwise
        # all subsequent calls will return a null result.
        #
        # 10    errors (JFLAG=2) occur when:
        #
        # .  there is a + or - but no number;  or
        #
        # .  the number is greater than BIG (defined below).
        #
        # 11    When an error has been detected, NSTRT is left
        # pointing to the character following the last
        # one used before the error came to light.
        #
        # 12    See also FLOTIN and eFLTIN.
        #
        # P.T.Wallace   Starlink   27 April 1998
        #
        # Copyright (C) 1998 Rutherford Appleton Laboratory
        #
        # -
        STRING = re.sub(r"\s+", r" ", STRING).replace("D", "e").replace("E", "e")
        STRARR = [xx.strip() for xx in re.split(r"[^+-0-9.e]", STRING)]
        JFLAG = 0
        try:
            IRESLT = int(STRARR[NSTRT])
        except ValueError:
            JFLAG = 2
        NSTRT += 1
        
        return IRESLT, JFLAG

    @classmethod
    def geoc(cls, P, H):
        # +
        # - - - - -
        # G E O C
        # - - - - -
        #
        # Convert geodetic position to geocentric (double precision)
        #
        # Given:
        # P     dp     latitude (geodetic, radians)
        # H     dp     height above reference spheroid (geodetic, metres)
        #
        # Returned:
        # R     dp     distance from Earth axis (AU)
        # Z     dp     distance from plane of Earth equator (AU)
        #
        # Notes:
        #
        # 1  Geocentric latitude can be obtained by evaluating ATAN2(Z,R).
        #
        # 2  IAU 1976 constants are used.
        #
        # Reference:
        #
        # Green,R.M., Spherical Astronomy, CUP 1985, p98.
        #
        # Last revision:   22 July 2004
        #
        # Copyright P.T.Wallace.  All rights reserved.
        #
        # -

        # Earth equatorial radius (metres)
        A0 = 6378140e0

        # Reference spheroid flattening factor and useful function
        F = 1e0 / 298.257e0
        B = (1e0 - F) ** 2

        # Astronomical unit in metres
        AU = 1.49597870e11

        # Geodetic to geocentric conversion
        SP = np.sin(P)
        CP = np.cos(P)
        C = 1e0 / np.sqrt(CP * CP + B * SP * SP)
        S = B * C
        R = (A0 * C + H) * CP / AU
        Z = (A0 * S + H) * SP / AU

        return R, Z


    @classmethod
    def dmat(cls, N, A, Y):
        # +
        # - - - - -
        # D M A T
        # - - - - -
        #
        # Matrix inversion & solution of simultaneous equations
        # (double precision)
        #
        # For the set of n simultaneous equations in n unknowns:
        # A.Y = X
        #
        # where:
        # A is a non-singular N x N matrix
        # Y is the vector of N unknowns
        # X is the known vector
        #
        # DMATRX computes:
        # the inverse of matrix A
        # the determinant of matrix A
        # the vector of N unknowns
        #
        # Arguments:
        #
        # symbol  type   dimension           before              after
        #
        # N      i                    no. of unknowns       unchanged
        # A      d      (N,N)             matrix             inverse
        # Y      d       (N)            known vector      solution vector
        # D      d                           -             determinant
        # * JF     i                           -           singularity flag
        # IW     i       (N)                 -              workspace
        #
        # * JF is the singularity flag.  If the matrix is non-singular, JF=0
        # is returned.  If the matrix is singular, JF=-1 & D=0D0 are
        # returned.  In the latter case, the contents of array A on return
        # are undefined.
        #
        # Algorithm:
        # Gaussian elimination with partial pivoting.
        #
        # Speed:
        # Very fast.
        #
        # Accuracy:
        # Fairly accurate - errors 1 to 4 times those of routines optimized
        # for accuracy.
        #
        # P.T.Wallace   Starlink   4 December 2001
        #
        # Copyright (C) 2001 Rutherford Appleton Laboratory
        #
        #
        #
        # -
        IW = np.zeros(N, dtype=int)
        SFA = 1e-20

        JF = 0
        D = 1e0
        for K in range(N):
            AMX = np.abs(A[K, K])
            IMX = K
            if K != N:
                for I in range(K + 0, N):
                    T = np.abs(A[I, K])
                    if T > AMX:
                        AMX = T
                        IMX = I

            if AMX < SFA:
                JF = -1
            else:
                if IMX != K:
                    for J in range(N):
                        T = A[K, J]
                        A[K, J] = A[IMX, J]
                        A[IMX, J] = T

                    T = Y[K]
                    Y[K] = Y[IMX]
                    Y[IMX] = T
                    D = -D

                IW[K] = IMX
                AKK = A[K, K]
                D = D * AKK
                if np.abs(D) < SFA:
                    JF = -1
                else:
                    AKK = 1e0 / AKK
                    A[K, K] = AKK
                    for J in range(N):
                        A[K, J] = A[K, J] * AKK if (J != K) else A[K, J]

                    YK = Y[K] * AKK
                    Y[K] = YK
                    for I in range(N):
                        AIK = A[I, K]
                        if I != K:
                            for J in range(N):
                                A[I, J] = (
                                    A[I, J] - AIK * A[K, J] if (J != K) else A[I, J]
                                )

                            Y[I] = Y[I] - AIK * YK

                    for I in range(N):
                        A[I, K] = -A[I, K] * AKK if (I != K) else A[I, K]

        if JF != 0:
            D = 0e0
        else:
            for K in range(N):
                NP1MK = N + 1 - K
                KI = IW[NP1MK]
                if NP1MK != KI:
                    for I in range(N):
                        T = A[I, NP1MK]
                        A[I, NP1MK] = A[I, KI]
                        A[I, KI] = T

        return D, JF, IW