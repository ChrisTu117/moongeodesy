"""High-precision ellipsoid formulas in mpmath, used as references that do not
come from GeographicLib.

The rhumb-line, GeodesicExact, polygon and intersection code in
src/geodesic is ported from GeographicLib, so tables made with GeographicLib's
own tools check a port against its source. This module computes the same
quantities from textbook definitions instead, with the integrals done by
mpmath.quad (tanh-sinh quadrature) at 40+ significant digits:

* meridian arc m(phi) = integral of the meridian radius M(phi);
* isometric latitude psi(phi) = asinh(tan phi) - e atanh(e sin phi);
* area between the equator and the parallel phi, per radian of longitude,
  A(phi) = integral of M(phi) N(phi) cos(phi) (closed form, checked against
  the quadrature in self_check());
* geodesics by Bessel's reduction to the auxiliary sphere (Helmert 1880,
  ch. 5): s = b * integral sqrt(1 + k^2 sin^2 sigma), and
  lambda = omega - f (2 - f) sin(alpha0) * integral d sigma /
  (1 + (1 - f) sqrt(1 + k^2 sin^2 sigma)), both evaluated by quadrature;
* the inverse geodesic by solving the same two integrals for the auxiliary
  longitude (Vincenty's 1975 scheme, with the series replaced by quadrature);
* the area under an edge by Green's theorem, S12 = integral of A(phi) dlambda
  along the curve.

None of GeographicLib's series, elliptic-integral or Fourier code is used.
Inputs are taken as the exact binary value of the double the MoonBit test
passes in (mp.mpf(float(x))). Angles in the public functions are degrees.

Requires: pip install mpmath
"""
import mpmath as mp

mp.mp.dps = 45

PI = mp.pi


def D(x):
    """Exact value of the double nearest x; an mpf is kept as it is."""
    if isinstance(x, mp.mpf):
        return x
    return mp.mpf(float(x))


def deg(x):
    return x * 180 / PI


def rad(x):
    return x * PI / 180


def wrap180(x):
    """Reduce degrees to (-180, 180]."""
    y = mp.fmod(x, 360)
    if y > 180:
        y -= 360
    elif y <= -180:
        y += 360
    return y


def ang_diff(lon1, lon2):
    """lon2 - lon1 reduced to (-180, 180] (east-going on opposite meridians)."""
    return wrap180(lon2 - lon1)


def quad(fn, pts):
    v, err = mp.quad(fn, pts, error=True, maxdegree=10)
    scale = max(abs(v), mp.mpf(1))
    if err > mp.mpf(10) ** (-36) * scale * max(1, len(pts)):
        raise ArithmeticError("quadrature error estimate %s for value %s" % (err, v))
    return v


class Ellipsoid:
    """Ellipsoid of revolution with equatorial radius a and flattening f.

    a and f are doubles (the same doubles the MoonBit code receives).
    f > 0 is oblate, f < 0 prolate.
    """

    def __init__(self, a, f):
        self.a = D(a)
        self.f = D(f)
        self.b = self.a * (1 - self.f)
        self.e2 = self.f * (2 - self.f)
        self.ep2 = self.e2 / (1 - self.f) ** 2
        self._Q = None

    # --- basic radii and latitudes ---

    def M(self, phi):
        """Meridian radius of curvature."""
        return self.a * (1 - self.e2) / (1 - self.e2 * mp.sin(phi) ** 2) ** mp.mpf(1.5)

    def N(self, phi):
        """Prime-vertical radius of curvature."""
        return self.a / mp.sqrt(1 - self.e2 * mp.sin(phi) ** 2)

    def atanhe_over_e(self, x):
        """atanh(e x) / e, real for both oblate and prolate ellipsoids."""
        if self.e2 > 0:
            e = mp.sqrt(self.e2)
            return mp.atanh(e * x) / e
        if self.e2 < 0:
            e = mp.sqrt(-self.e2)
            return mp.atan(e * x) / e
        return x

    def meridian_arc(self, phi):
        """Distance along the meridian from the equator to latitude phi (rad)."""
        if phi == 0:
            return mp.mpf(0)
        return quad(self.M, [0, phi])

    def quarter_meridian(self):
        if self._Q is None:
            self._Q = self.meridian_arc(PI / 2)
        return self._Q

    def lat_from_arc(self, m):
        """Inverse of meridian_arc for |m| <= quarter meridian (Newton)."""
        Q = self.quarter_meridian()
        if abs(m) >= Q:
            return mp.sign(m) * PI / 2
        phi = m / Q * PI / 2
        for _ in range(100):
            d = (self.meridian_arc(phi) - m) / self.M(phi)
            phi -= d
            if abs(d) < mp.mpf(10) ** (-mp.mp.dps + 3):
                break
        else:
            raise ArithmeticError("lat_from_arc did not converge")
        return phi

    def isometric(self, phi):
        """Isometric latitude psi; +-inf at the poles."""
        if abs(phi) >= PI / 2:
            return mp.inf * mp.sign(phi)
        s = mp.sin(phi)
        return mp.asinh(mp.tan(phi)) - self.e2 * self.atanhe_over_e(s)

    def area_fn(self, phi):
        """A(phi): area between the equator and parallel phi per radian of
        longitude, the integral of M N cos(phi) from 0 to phi."""
        s = mp.sin(phi)
        return self.b ** 2 / 2 * (s / (1 - self.e2 * s * s) + self.atanhe_over_e(s))

    def total_area(self):
        return 4 * PI * self.area_fn(PI / 2)

    def self_check(self):
        """The closed form A(phi) against quadrature of its definition, and the
        total area against the textbook surface-of-revolution integral."""
        for p in (mp.mpf("0.3"), mp.mpf("-1.1"), PI / 2):
            q = quad(lambda t: self.M(t) * self.N(t) * mp.cos(t), [0, p])
            assert abs(q - self.area_fn(p)) < mp.mpf(10) ** (-30) * self.a ** 2
        # surface of revolution: 2 pi * integral of x ds over the meridian
        tot = 2 * PI * quad(lambda t: self.N(t) * mp.cos(t) * self.M(t), [-PI / 2, 0, PI / 2])
        assert abs(tot - self.total_area()) < mp.mpf(10) ** (-30) * self.a ** 2

    # --- rhumb lines ---

    def rhumb_inverse(self, lat1, lon1, lat2, lon2):
        """Rhumb line from (lat1, lon1) to (lat2, lon2), shortest way round.

        Returns (azi12 deg, s12, area12). azi12 is None when undefined (both
        points at the same pole).
        """
        lat1, lon1, lat2, lon2 = D(lat1), D(lon1), D(lat2), D(lon2)
        lam12 = rad(ang_diff(lon1, lon2))
        phi1, phi2 = rad(lat1), rad(lat2)
        m12 = self.meridian_arc(phi2) - self.meridian_arc(phi1)
        pole1 = abs(lat1) == 90
        pole2 = abs(lat2) == 90
        if pole1 or pole2:
            # The line is a meridian; all of the longitude change happens at the
            # pole, so the area is lam12 * A(pole) (the limit of the formula below).
            if pole1 and pole2 and lat1 == lat2:
                return None, mp.mpf(0), lam12 * self.area_fn(phi1)
            azi = mp.mpf(0) if m12 > 0 else mp.mpf(180)
            if pole1 and pole2:
                raise ValueError("pole to opposite pole is not used")
            area = lam12 * self.area_fn(phi1 if pole1 else phi2)
            return azi, abs(m12), area
        psi12 = self.isometric(phi2) - self.isometric(phi1)
        azi = deg(mp.atan2(lam12, psi12))
        if psi12 == 0:
            s12 = abs(lam12) * self.N(phi1) * mp.cos(phi1)
            area = lam12 * self.area_fn(phi1)
        else:
            s12 = abs(m12) * mp.sqrt(1 + (lam12 / psi12) ** 2)
            area = self.rhumb_area(phi1, phi2, lam12)
        return azi, s12, area

    def rhumb_area(self, phi1, phi2, lam12):
        """Green's theorem along the rhumb: dlambda = (lam12 / psi12) dpsi, so
        S12 = lam12 / psi12 * integral A(phi(psi)) dpsi. Integrated over phi
        with dpsi/dphi = M / (N cos phi)."""
        if phi1 == phi2:
            return lam12 * self.area_fn(phi1)
        psi12 = self.isometric(phi2) - self.isometric(phi1)
        integ = quad(lambda p: self.area_fn(p) * self.M(p) / (self.N(p) * mp.cos(p)), [phi1, phi2])
        return lam12 / psi12 * integ

    def rhumb_direct(self, lat1, lon1, azi12, s12, unroll=False):
        """Travel s12 along constant azimuth azi12.

        Returns (lat2, lon2, area12, lon12) in degrees / m^2; lon2, area12 and
        lon12 are None when the line passes a pole. lon2 is reduced to
        (-180, 180] unless unroll; lon12 is the unrolled change in longitude.
        """
        lat1, lon1, azi12, s12 = D(lat1), D(lon1), D(azi12), D(s12)
        assert abs(lat1) < 90
        phi1 = rad(lat1)
        alp = rad(azi12)
        calp, salp = mp.cos(alp), mp.sin(alp)
        if azi12 in (90, -90, 270, -270):
            calp = mp.mpf(0)
        if azi12 in (0, 180, -180, 360):
            salp = mp.mpf(0)
        Q = self.quarter_meridian()
        m2 = self.meridian_arc(phi1) + s12 * calp
        if abs(m2) > Q:
            # Past a pole: reflect the rectifying arc; longitude undefined.
            m2 = mp.fmod(m2, 4 * Q)
            if m2 > 2 * Q:
                m2 -= 4 * Q
            elif m2 < -2 * Q:
                m2 += 4 * Q
            if m2 > Q:
                m2 = 2 * Q - m2
            elif m2 < -Q:
                m2 = -2 * Q - m2
            return deg(self.lat_from_arc(m2)), None, None, None
        phi2 = self.lat_from_arc(m2) if calp != 0 else phi1
        if calp == 0:
            lam12 = s12 * salp / (self.N(phi1) * mp.cos(phi1))
        else:
            lam12 = salp / calp * (self.isometric(phi2) - self.isometric(phi1))
        area = self.rhumb_area(phi1, phi2, lam12)
        lon2 = lon1 + deg(lam12) if unroll else wrap180(wrap180(lon1) + deg(lam12))
        return deg(phi2), lon2, area, deg(lam12)

    # --- geodesics on the auxiliary sphere (Bessel / Helmert) ---

    def _w(self, k2, sig):
        return mp.sqrt(1 + k2 * mp.sin(sig) ** 2)

    @staticmethod
    def _splits(s1, s2):
        """Interval [s1, s2] split at the vertices sigma = pi/2 + n pi."""
        lo, hi = (s1, s2) if s1 <= s2 else (s2, s1)
        pts = [lo]
        n = mp.ceil((lo - PI / 2) / PI)
        while True:
            v = PI / 2 + n * PI
            if v >= hi:
                break
            if v > lo:
                pts.append(v)
            n += 1
        pts.append(hi)
        return pts if s1 <= s2 else pts[::-1]

    def _int(self, fn, s1, s2):
        if s1 == s2:
            return mp.mpf(0)
        return quad(fn, self._splits(s1, s2))

    def _dist_int(self, k2, s1, s2):
        """b * integral_{s1}^{s2} sqrt(1 + k2 sin^2)."""
        return self.b * self._int(lambda t: self._w(k2, t), s1, s2)

    def _lam_int(self, k2, s1, s2):
        """integral_{s1}^{s2} (2 - f) / (1 + (1 - f) sqrt(1 + k2 sin^2))."""
        f = self.f
        return self._int(lambda t: (2 - f) / (1 + (1 - f) * self._w(k2, t)), s1, s2)

    def _omega(self, sa0, sig):
        """Auxiliary-sphere longitude from the node, tan(omega) = sa0 tan(sigma),
        continuous in sigma."""
        if sa0 == 0:
            # A meridian: omega jumps by pi at each pole.
            return PI * mp.floor((sig + PI / 2) / PI)
        s = abs(sa0)
        d = mp.atan2(s * mp.sin(sig), mp.cos(sig)) - sig
        d = d - 2 * PI * mp.nint(d / (2 * PI))
        om = sig + d
        return om if sa0 > 0 else -om

    def _lat_azi(self, sa0, ca0, sig):
        sbet = ca0 * mp.sin(sig)
        cbet = mp.hypot(sa0, ca0 * mp.cos(sig))
        phi = mp.atan2(sbet, (1 - self.f) * cbet)
        alp = mp.atan2(sa0, ca0 * mp.cos(sig))
        return phi, alp

    def _area_int(self, sa0, ca0, k2, s1, s2):
        """S12 = integral of A(phi) dlambda along the geodesic, with
        dlambda/dsigma = sa0 [1 / cos^2 beta - f (2 - f) / (1 + (1 - f) w)]."""
        if sa0 == 0:
            # A meridian: no longitude change unless it passes a pole.
            return mp.mpf(0) if len(self._splits(s1, s2)) == 2 else None
        f = self.f

        def integrand(t):
            st = mp.sin(t)
            sbet = ca0 * st
            c2bet = 1 - sbet * sbet
            sphi = sbet / mp.sqrt(sbet * sbet + (1 - f) ** 2 * c2bet)
            w = mp.sqrt(1 + k2 * st * st)
            dlam = sa0 * (1 / c2bet - f * (2 - f) / (1 + (1 - f) * w))
            A = self.b ** 2 / 2 * (sphi / (1 - self.e2 * sphi * sphi) + self.atanhe_over_e(sphi))
            return A * dlam

        return self._int(integrand, s1, s2)

    def _start(self, lat1, azi1):
        phi1, alp1 = rad(lat1), rad(azi1)
        sbet1 = (1 - self.f) * mp.sin(phi1)
        cbet1 = mp.cos(phi1)
        r = mp.hypot(sbet1, cbet1)
        sbet1, cbet1 = sbet1 / r, cbet1 / r
        salp1, calp1 = mp.sin(alp1), mp.cos(alp1)
        if azi1 in (0, 180, -180):
            salp1 = mp.mpf(0)
        if azi1 in (90, -90):
            calp1 = mp.mpf(0)
        sa0 = salp1 * cbet1
        ca0 = mp.hypot(calp1, salp1 * sbet1)
        sig1 = mp.atan2(sbet1, calp1 * cbet1)
        k2 = self.ep2 * ca0 * ca0
        return sa0, ca0, sig1, k2

    def geod_direct(self, lat1, lon1, azi1, s12, area=True):
        """Direct geodesic problem by quadrature.

        Returns dict with lat2, lon2 (reduced to (-180, 180]), azi2 (deg),
        and S12 (m^2, None for a meridian).
        """
        lat1, lon1, azi1, s12 = D(lat1), D(lon1), D(azi1), D(s12)
        assert abs(lat1) < 90
        sa0, ca0, sig1, k2 = self._start(lat1, azi1)
        # Solve b * int_{sig1}^{sig2} w = s12 for sig2 (Newton).
        wbar = mp.sqrt(1 + k2 / 2)
        sig2 = sig1 + s12 / (self.b * wbar)
        for _ in range(200):
            d = (self._dist_int(k2, sig1, sig2) - s12) / (self.b * self._w(k2, sig2))
            sig2 -= d
            if abs(d) < mp.mpf(10) ** (-mp.mp.dps + 4):
                break
        else:
            raise ArithmeticError("geod_direct: sigma did not converge")
        om12 = self._omega(sa0, sig2) - self._omega(sa0, sig1)
        lam12 = om12 - self.f * sa0 * self._lam_int(k2, sig1, sig2)
        phi2, alp2 = self._lat_azi(sa0, ca0, sig2)
        S12 = self._area_int(sa0, ca0, k2, sig1, sig2) if area else None
        return {
            "lat2": deg(phi2),
            "lon2": wrap180(lon1 + deg(lam12)),
            "lon2_unrolled": lon1 + deg(lam12),
            "azi2": deg(alp2),
            "S12": S12,
            "sig12": sig2 - sig1,
        }

    def geod_inverse(self, lat1, lon1, lat2, lon2):
        """Inverse geodesic problem by quadrature: find the auxiliary-sphere
        longitude omega12 for which the longitude integral gives lon2 - lon1
        (secant iteration started from omega12 = lambda12), as in Vincenty's
        1975 scheme but with exact integrals. Valid away from nearly
        antipodal points (|lambda12| well below 180 deg).

        Returns dict with s12, azi1, azi2 (deg), S12, lam12 (rad).
        """
        lat1, lon1, lat2, lon2 = D(lat1), D(lon1), D(lat2), D(lon2)
        lam12 = rad(ang_diff(lon1, lon2))
        f = self.f
        phi1, phi2 = rad(lat1), rad(lat2)

        def red(phi):
            sb = (1 - f) * mp.sin(phi)
            cb = mp.cos(phi)
            r = mp.hypot(sb, cb)
            return sb / r, cb / r

        sbet1, cbet1 = red(phi1)
        sbet2, cbet2 = red(phi2)

        def geom(om):
            som, com = mp.sin(om), mp.cos(om)
            y = cbet2 * som
            x = cbet1 * sbet2 - sbet1 * cbet2 * com
            alp1 = mp.atan2(y, x)
            salp1, calp1 = mp.sin(alp1), mp.cos(alp1)
            sa0 = salp1 * cbet1
            ca0 = mp.hypot(calp1, salp1 * sbet1)
            sig1 = mp.atan2(sbet1, calp1 * cbet1)
            sig12 = mp.atan2(mp.hypot(x, y), sbet1 * sbet2 + cbet1 * cbet2 * com)
            k2 = self.ep2 * ca0 * ca0
            return alp1, sa0, ca0, sig1, sig1 + sig12, k2

        def F(om):
            _, sa0, ca0, s1, s2, k2 = geom(om)
            return om - f * sa0 * self._lam_int(k2, s1, s2) - lam12

        om = mp.findroot(F, (lam12, lam12 * (1 + f / 2) + mp.mpf(10) ** -20), solver="secant",
                         tol=mp.mpf(10) ** (-2 * mp.mp.dps + 10), maxsteps=200)
        assert abs(F(om)) < mp.mpf(10) ** (-mp.mp.dps + 6)
        alp1, sa0, ca0, s1, s2, k2 = geom(om)
        _, alp2 = self._lat_azi(sa0, ca0, s2)
        return {
            "s12": self._dist_int(k2, s1, s2),
            "azi1": deg(alp1),
            "azi2": deg(alp2),
            "S12": self._area_int(sa0, ca0, k2, s1, s2),
            "sig12": s2 - s1,
            "lam12": lam12,
        }


def lit(v):
    """MoonBit literal for the double nearest v (None -> NaN)."""
    if v is None:
        return "@double.not_a_number"
    x = float(v)
    if x != x:
        return "@double.not_a_number"
    r = repr(x)
    if r == "-0.0":
        return "-0.0"
    if "e" in r:
        m, e = r.split("e")
        if "." not in m:
            m += ".0"
        return m + "e" + e
    if "." not in r:
        return r + ".0"
    return r
