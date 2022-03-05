// Package offers two main functionalities: solving the "direct" and the "inverse" problems.
//
// Direct:
//
// Place a second point, given the first point, an azimuth, and a distance.
//
// # Arguments
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
//
// # Returns
//
// There are a variety of outputs associated with this calculation. We save computation by
// only calculating the outputs you need. You can get any and all of the following
//
//  - lat2 latitude of point 2 [degrees].
//  - lon2 longitude of point 2 [degrees].
//  - azi2 (forward) azimuth at point 2 [degrees].
//  - m12 reduced length of geodesic (meters).
//  - M12 geodesic scale of point 2 relative to point 1 [dimensionless].
//  - M21 geodesic scale of point 1 relative to point 2 [dimensionless].
//  - S12 area under the geodesic [meters^2]
//  - a12 arc length of between point 1 and point 2 [degrees].
//
// Call the appropriate function to get the output you need. The following functions are
// available for solving the Direct problem. Each describes what it returns
//
// - DirectCalcLatLon -> calculate latitude and longitude
//
// - DirectCalcLatLonAzi -> calculate latitude, longitude, and azimuth
//
// - DirectCalcLatLonAziReducedLength -> calculate latitude, longitude, azimuth, and
// reduced length of the geodesic
//
// - DirectCalcLatLonAziGeodesicScales -> calculate latitude, longitude, azimuth, and
// the geodesic scales
//
// - DirectCalcLatLonAziReducedLengthGeodesicScales -> calculate latitude, longitude,
// azimuth, reduced length, and the geodesic scales
//
// - DirectCalcAll -> calculates all of the above plus Area under the geodesic and the
// arc length between point 1 and point 2
//
// =====================================================================================
// =====================================================================================
// =====================================================================================
// Indirect:
//
// Measure the distance (and other values) between two points.
//
// # Arguments
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
//
// # Returns
//
// There are a variety of outputs associated with this calculation. We save computation by
// only calculating the outputs you need. You can get any and all of the following
//
// - s12 distance between point 1 and point 2 (meters).
// - azi1 azimuth at point 1 [degrees].
// - azi2 (forward) azimuth at point 2 [degrees].
// - m12 reduced length of geodesic (meters).
// - M12 geodesic scale of point 2 relative to point 1 [dimensionless].
// - M21 geodesic scale of point 1 relative to point 2 [dimensionless].
// - S12 area under the geodesic [meters^2]
// - a12 arc length of between point 1 and point 2 [degrees].
//
// Call the appropriate function to get the output you need. The following functions are
// available for solving the Direct problem. Each describes what it returns
//
//
//
//  `lat1` and `lat2` should be in the range [&minus;90&deg;, 90&deg;].
//  The values of `azi1` and `azi2` returned are in the range
//  [&minus;180&deg;, 180&deg;].
//
// If either point is at a pole, the azimuth is defined by keeping the
// longitude fixed, writing `lat` = &plusmn;(90&deg; &minus; &epsilon;),
// and taking the limit &epsilon; &rarr; 0+.
//
// The solution to the inverse problem is found using Newton's method.  If
// this fails to converge (this is very unlikely in geodetic applications
// but does occur for very eccentric ellipsoids), then the bisection method
// is used to refine the solution.
package geographiclibgo

import "math"

type DirectAndInverse interface {
	DirectCalcAll(lat1_deg, lon1_deg, azi1_deg, s12_m float64) AllDirectResults
	DirectCalcLatLon(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLon
	DirectCalcLatLonAzi(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAzi
	DirectCalcLatLonAziGeodesicScales(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAziGeodesicScales
	DirectCalcLatLonAziReducedLength(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAziReducedLength
	DirectCalcLatLonAziReducedLengthGeodesicScales(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAziReducedLengthGeodesicScales
	DirectCalcWithCapabilities(lat1_deg, lon1_deg, azi1_deg, s12_m float64, capabilities uint64) AllDirectResults
	DirectLineWithCapabilities(lat1_deg, lon1_deg, azi1_deg, s12_m float64, capabilities uint64) GeodesicLine
	EqualtorialRadius() float64
	Flattening() float64
	InverseCalcAll(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) AllInverseResults
	InverseCalcAzimuthsArcLength(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) AzimuthsArcLength
	InverseCalcDistance(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) float64
	InverseCalcDistanceArcLength(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) DistanceArcLength
	InverseCalcDistanceAzimuths(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) DistanceAzimuths
	InverseCalcDistanceAzimuthsArcLength(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) DistanceAzimuthsArcLength
	InverseCalcDistanceAzimuthsArcLengthReducedLength(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) DistanceAzimuthsArcLengthReducedLength
	InverseCalcDistanceAzimuthsArcLengthReducedLengthScales(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) DistanceAzimuthsArcLengthReducedLengthScales
	InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64, capabilities uint64) AllInverseResults
	InverseLineWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64, capabilities uint64) GeodesicLine
	LineWithCapabilities(lat1_deg, lon1_deg, azi1_deg float64, capabilities uint64) GeodesicLine
}

const WGS84_A float64 = 6378137.0

// Evaluating this as 1000000000.0 / (298257223563f_64) reduces the
// round-off error by about 10%.  However, expressing the flattening as
// 1/298.257223563 is well ingrained.
const WGS84_F float64 = 1.0 / ((298257223563.0) / 1000000000.0)

const _GEODESIC_ORDER int64 = 6
const nC3x_ int64 = 15
const nC4x_ int64 = 21

func coeff_A3() [18]float64 {
	return [18]float64{
		-3.0, 128.0, -2.0, -3.0, 64.0, -1.0, -3.0, -1.0, 16.0, 3.0, -1.0, -2.0, 8.0, 1.0, -1.0, 2.0,
		1.0, 1.0,
	}
}

func coeff_C3() [45]float64 {
	return [45]float64{
		3.0, 128.0, 2.0, 5.0, 128.0, -1.0, 3.0, 3.0, 64.0, -1.0, 0.0, 1.0, 8.0, -1.0, 1.0, 4.0, 5.0,
		256.0, 1.0, 3.0, 128.0, -3.0, -2.0, 3.0, 64.0, 1.0, -3.0, 2.0, 32.0, 7.0, 512.0, -10.0, 9.0,
		384.0, 5.0, -9.0, 5.0, 192.0, 7.0, 512.0, -14.0, 7.0, 512.0, 21.0, 2560.0,
	}
}

func coeff_C4() [77]float64 {
	return [77]float64{
		97.0, 15015.0, 1088.0, 156.0, 45045.0, -224.0, -4784.0, 1573.0, 45045.0, -10656.0, 14144.0,
		-4576.0, -858.0, 45045.0, 64.0, 624.0, -4576.0, 6864.0, -3003.0, 15015.0, 100.0, 208.0, 572.0,
		3432.0, -12012.0, 30030.0, 45045.0, 1.0, 9009.0, -2944.0, 468.0, 135135.0, 5792.0, 1040.0,
		-1287.0, 135135.0, 5952.0, -11648.0, 9152.0, -2574.0, 135135.0, -64.0, -624.0, 4576.0, -6864.0,
		3003.0, 135135.0, 8.0, 10725.0, 1856.0, -936.0, 225225.0, -8448.0, 4992.0, -1144.0, 225225.0,
		-1440.0, 4160.0, -4576.0, 1716.0, 225225.0, -136.0, 63063.0, 1024.0, -208.0, 105105.0, 3584.0,
		-3328.0, 1144.0, 315315.0, -128.0, 135135.0, -2560.0, 832.0, 405405.0, 128.0, 99099.0,
	}
}

type Geodesic struct {
	a     float64
	f     float64
	f1    float64
	e2    float64
	ep2   float64
	n     float64
	b     float64
	c2    float64
	etol2 float64

	GEODESIC_ORDER int64
	nC3x_          int64
	nC4x_          int64
	maxit1_        uint64
	maxit2_        uint64

	_A3x [_GEODESIC_ORDER]float64
	_C3x [nC3x_]float64
	_C4x [nC4x_]float64

	tiny_    float64
	tol0_    float64
	tol1_    float64
	tol2_    float64
	tolb_    float64
	xthresh_ float64
}

func NewGeodesic(a, f float64) Geodesic {
	var maxit1_ uint64 = 20
	maxit2_ := maxit1_ + _DIGITS + 10
	tiny_ := math.Sqrt(get_min_val())
	tol0_ := get_epsilon()
	tol1_ := 200.0 * tol0_
	tol2_ := math.Sqrt(tol0_)
	tolb_ := tol0_ * tol2_
	xthresh_ := 1000.0 * tol2_

	_f1 := 1.0 - f
	_e2 := f * (2.0 - f)
	_ep2 := _e2 / sq(_f1)
	_n := f / (2.0 - f)
	_b := a * _f1

	var is_f_neg float64
	if f < 0.0 {
		is_f_neg = -1.0
	} else {
		is_f_neg = 1.0
	}

	to_mul := eatanhe(1.0, is_f_neg*math.Sqrt(math.Abs(_e2))) / _e2
	if _e2 == 0.0 {
		to_mul = 1.0
	}
	_c2 := (sq(a) + sq(_b)*to_mul) / 2.0
	_etol2 := 0.1 * tol2_ / math.Sqrt(math.Max(math.Abs(f), 0.001)*math.Min((1.0-f/2.0), 1.0)/2.0)

	_A3x := [_GEODESIC_ORDER]float64{}
	_C3x := [nC3x_]float64{}
	_C4x := [nC4x_]float64{}

	// Call a3coeff
	var o int64 = 0
	k := 0

	coefa3 := coeff_A3()
	for j := _GEODESIC_ORDER - 1; j >= 0; j-- {
		m := int64(math.Min(float64(j), float64(_GEODESIC_ORDER-j-1)))
		_A3x[k] = polyval(m, coefa3[o:], _n) / coefa3[o+m+1]
		k += 1
		o += m + 2
	}

	// c3coeff
	o = 0
	k = 0

	coefc3 := coeff_C3()
	for l := 1; l < int(_GEODESIC_ORDER); l++ {
		for j := int(_GEODESIC_ORDER) - 1; j >= l; j-- {
			m := int64(math.Min(float64(j), float64(int(_GEODESIC_ORDER)-j-1)))
			_C3x[k] = polyval(m, coefc3[o:], _n) / coefc3[o+m+1]
			k += 1
			o += m + 2
		}

	}

	// c4coeff
	o = 0
	k = 0

	coefc4 := coeff_C4()
	for l := 0; l < int(_GEODESIC_ORDER); l++ {
		for j := int(_GEODESIC_ORDER) - 1; j >= l; j-- {
			m := int64(int(_GEODESIC_ORDER) - j - 1)
			_C4x[k] = polyval(m, coefc4[o:], _n) / coefc4[(o+m+1)]
			k += 1
			o += m + 2
		}
	}

	return Geodesic{
		a,
		f,
		_f1,
		_e2,
		_ep2,
		_n,
		_b,
		_c2,
		_etol2,

		_GEODESIC_ORDER,
		nC3x_,
		nC4x_,
		maxit1_,
		maxit2_,

		_A3x,
		_C3x,
		_C4x,

		tiny_,
		tol0_,
		tol1_,
		tol2_,
		tolb_,
		xthresh_,
	}
}

func Wgs84() Geodesic {
	return NewGeodesic(WGS84_A, WGS84_F)
}

func (g Geodesic) EqualtorialRadius() float64 {
	return g.a
}

func (g Geodesic) Flattening() float64 {
	return g.f
}

func (g Geodesic) _A3f(eps float64) float64 {
	return polyval(int64(_GEODESIC_ORDER-1), g._A3x[:], eps)
}

func (g Geodesic) _C3f(eps float64, c []float64) {
	mult := 1.0
	o := 0
	for l := 1; l < int(_GEODESIC_ORDER); l++ {
		m := int(_GEODESIC_ORDER) - l - 1
		mult *= eps
		c[l] = mult * polyval(int64(m), g._C3x[o:], eps)
		o += m + 1
	}
}

func (g Geodesic) _C4f(eps float64, c []float64) {
	mult := 1.0
	o := 0
	for l := 0; l < int(_GEODESIC_ORDER); l++ {
		m := int(_GEODESIC_ORDER) - l - 1
		c[l] = mult * polyval(int64(m), g._C4x[o:], eps)
		o += m + 1
		mult *= eps
	}
}

func (g Geodesic) _Lengths(
	eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2 float64,
	outmask uint64,
	c1a []float64,
	c2a []float64,
) (float64, float64, float64, float64, float64) {
	outmask &= OUT_MASK
	s12b := math.NaN()
	m12b := math.NaN()
	m0 := math.NaN()
	M12 := math.NaN()
	M21 := math.NaN()

	A1 := 0.0
	A2 := 0.0
	m0x := 0.0
	J12 := 0.0

	if outmask&(DISTANCE|REDUCEDLENGTH|GEODESICSCALE) != 0 {
		A1 = a1m1f(eps, _GEODESIC_ORDER)
		c1f(eps, c1a, int(_GEODESIC_ORDER))
		if outmask&(REDUCEDLENGTH|GEODESICSCALE) != 0 {
			A2 = a2m1f(eps, _GEODESIC_ORDER)
			c2f(eps, c2a, int(_GEODESIC_ORDER))
			m0x = A1 - A2
			A2 = 1.0 + A2
		}
		A1 = 1.0 + A1
	}

	if outmask&DISTANCE != 0 {
		B1 := sin_cos_series(true, ssig2, csig2, c1a) - sin_cos_series(true, ssig1, csig1, c1a)
		s12b = A1 * (sig12 + B1)
		if outmask&(REDUCEDLENGTH|GEODESICSCALE) != 0 {
			B2 := sin_cos_series(true, ssig2, csig2, c2a) - sin_cos_series(true, ssig1, csig1, c2a)
			J12 = m0x*sig12 + (A1*B1 - A2*B2)
		}
	} else if outmask&(REDUCEDLENGTH|GEODESICSCALE) != 0 {
		for l := 1; l <= int(_GEODESIC_ORDER); l++ {
			c2a[l] = A1*c1a[l] - A2*c2a[l]
		}
		J12 = m0x*sig12 + (sin_cos_series(true, ssig2, csig2, c2a) - sin_cos_series(true, ssig1, csig1, c2a))
	}

	if outmask&REDUCEDLENGTH != 0 {
		m0 = m0x
		// J12 is wrong
		m12b = dn2*(csig1*ssig2) - dn1*(ssig1*csig2) - csig1*csig2*J12
	}

	if outmask&GEODESICSCALE != 0 {
		csig12 := csig1*csig2 + ssig1*ssig2
		t := g.ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
		M12 = csig12 + (t*ssig2-csig2*J12)*ssig1/dn1
		M21 = csig12 - (t*ssig1-csig1*J12)*ssig2/dn2
	}
	return s12b, m12b, m0, M12, M21
}

func (g Geodesic) _InverseStart(
	sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12 float64,
	c1a []float64,
	c2a []float64,
) (float64, float64, float64, float64, float64, float64) {
	sig12 := -1.0
	salp2 := math.NaN()
	calp2 := math.NaN()
	dnm := math.NaN()

	var somg12 float64
	var comg12 float64

	sbet12 := sbet2*cbet1 - cbet2*sbet1
	cbet12 := cbet2*cbet1 + sbet2*sbet1

	sbet12a := sbet2 * cbet1
	sbet12a += cbet2 * sbet1

	shortline := cbet12 >= 0.0 && sbet12 < 0.5 && cbet2*lam12 < 0.5
	if shortline {
		sbetm2 := sq(sbet1 + sbet2)
		sbetm2 /= sbetm2 + sq(cbet1+cbet2)
		dnm = math.Sqrt(1.0 + g.ep2*sbetm2)
		omg12 := lam12 / (g.f1 * dnm)
		somg12 = math.Sin(omg12)
		comg12 = math.Cos(omg12)
	} else {
		somg12 = slam12
		comg12 = clam12
	}

	salp1 := cbet2 * somg12

	calp1 := sbet12a - cbet2*sbet1*sq(somg12)/(1.0-comg12)
	if comg12 >= 0.0 {
		calp1 = sbet12 + cbet2*sbet1*sq(somg12)/(1.0+comg12)
	}

	ssig12 := math.Hypot(salp1, calp1)
	csig12 := sbet1*sbet2 + cbet1*cbet2*comg12

	if shortline && (ssig12 < g.etol2) {
		salp2 = cbet1 * somg12
		var to_mul float64
		if comg12 >= 0.0 {
			to_mul = sq(somg12) / (1.0 + comg12)
		} else {
			to_mul = 1.0 - comg12
		}
		calp2 = sbet12 - cbet1*sbet2*to_mul

		salp2, calp2 = norm(salp2, calp2)
		sig12 = math.Atan2(ssig12, csig12)
	} else if math.Abs(g.n) > 0.1 || csig12 >= 0.0 || ssig12 >= 6.0*math.Abs(g.n)*math.Pi*sq(cbet1) {
	} else {
		var x float64
		var y float64
		var betscale float64
		var lamscale float64
		lam12x := math.Atan2(-slam12, -clam12)
		if g.f >= 0.0 {
			k2 := sq(sbet1) * g.ep2
			eps := k2 / (2.0*(1.0+math.Sqrt(1.0+k2)) + k2)
			lamscale = g.f * cbet1 * g._A3f(eps) * math.Pi
			betscale = lamscale * cbet1
			x = lam12x / lamscale
			y = sbet12a / betscale
		} else {
			cbet12a := cbet2*cbet1 - sbet2*sbet1
			bet12a := math.Atan2(sbet12, cbet12a)
			_, m12b, m0, _, _ := g._Lengths(
				g.n,
				math.Pi+bet12a,
				sbet1,
				-cbet1,
				dn1,
				sbet2,
				cbet2,
				dn2,
				cbet1,
				cbet2,
				REDUCEDLENGTH,
				c1a,
				c2a,
			)
			x = -1.0 + m12b/(cbet1*cbet2*m0*math.Pi)
			var betscale float64
			if x < -0.01 {
				betscale = sbet12a / x
			} else {
				betscale = -g.f * sq(cbet1) * math.Pi
			}
			lamscale = betscale / cbet1
			y = lam12x / lamscale
		}
		if y > -g.tol1_ && x > -1.0-g.xthresh_ {
			if g.f >= 0.0 {
				salp1 = math.Min(-x, 1.0)
				calp1 = -math.Sqrt(1.0 - sq(salp1))
			} else {
				var to_compare float64
				if x > -g.tol1_ {
					to_compare = 0.0
				} else {
					to_compare = -1.0
				}
				calp1 = math.Max(x, to_compare)
				salp1 = math.Sqrt(1.0 - sq(calp1))
			}
		} else {
			k := astroid(x, y)
			var to_mul float64
			if g.f >= 0.0 {
				to_mul = -x * k / (1.0 + k)
			} else {
				to_mul = -y * (1.0 + k) / k
			}

			omg12a := lamscale * to_mul
			somg12 = math.Sin(omg12a)
			comg12 = -math.Cos(omg12a)
			salp1 = cbet2 * somg12
			calp1 = sbet12a - cbet2*sbet1*sq(somg12)/(1.0-comg12)
		}
	}

	if !(salp1 <= 0.0) {
		salp1, calp1 = norm(salp1, calp1)
	} else {
		salp1 = 1.0
		calp1 = 0.0
	}
	return sig12, salp1, calp1, salp2, calp2, dnm
}

func (g Geodesic) _Lambda12(
	sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1, slam120, clam120 float64,
	diffp bool,
	c1a []float64,
	c2a []float64,
	c3a []float64,
) (float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64) {

	if sbet1 == 0.0 && calp1 == 0.0 {
		calp1 = -g.tiny_
	}
	salp0 := salp1 * cbet1
	calp0 := math.Hypot(calp1, salp1*sbet1)

	ssig1 := sbet1
	somg1 := salp0 * sbet1
	csig1 := calp1 * cbet1
	comg1 := calp1 * cbet1
	ssig1, csig1 = norm(ssig1, csig1)

	var salp2 float64
	if cbet2 != cbet1 {
		salp2 = salp0 / cbet2
	} else {
		salp2 = salp1
	}

	var to_add float64
	if cbet1 < -sbet1 {
		to_add = (cbet2 - cbet1) * (cbet1 + cbet2)
	} else {
		to_add = (sbet1 - sbet2) * (sbet1 + sbet2)
	}

	calp2 := math.Abs(calp1)
	if cbet2 != cbet1 || math.Abs(sbet2) != -sbet1 {
		calp2 = math.Sqrt(sq(calp1*cbet1)+to_add) / cbet2

	}

	ssig2 := sbet2
	somg2 := salp0 * sbet2
	csig2 := calp2 * cbet2
	comg2 := calp2 * cbet2
	ssig2, csig2 = norm(ssig2, csig2)

	sig12 := math.Atan2(math.Max(csig1*ssig2-ssig1*csig2, 0.0), csig1*csig2+ssig1*ssig2)
	somg12 := math.Max((comg1*somg2 - somg1*comg2), 0.0)
	comg12 := comg1*comg2 + somg1*somg2
	eta := math.Atan2(somg12*clam120-comg12*slam120, comg12*clam120+somg12*slam120)

	k2 := sq(calp0) * g.ep2
	eps := k2 / (2.0*(1.0+math.Sqrt(1.0+k2)) + k2)
	g._C3f(eps, c3a)
	B312 := sin_cos_series(true, ssig2, csig2, c3a) - sin_cos_series(true, ssig1, csig1, c3a)
	domg12 := -g.f * g._A3f(eps) * salp0 * (sig12 + B312)
	lam12 := eta + domg12

	var dlam12 float64
	if diffp {
		if calp2 == 0.0 {
			dlam12 = -2.0 * g.f1 * dn1 / sbet1
		} else {
			_, res, _, _, _ := g._Lengths(
				eps,
				sig12,
				ssig1,
				csig1,
				dn1,
				ssig2,
				csig2,
				dn2,
				cbet1,
				cbet2,
				REDUCEDLENGTH,
				c1a,
				c2a,
			)
			dlam12 = res
			dlam12 *= g.f1 / (calp2 * cbet2)
		}
	} else {
		dlam12 = math.NaN()
	}
	return lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12

}

func (g Geodesic) _gen_inverse_azi(
	lat1, lon1, lat2, lon2 float64,
	outmask uint64,
) (
	a12 float64,
	s12 float64,
	azi1 float64,
	azi2 float64,
	m12 float64,
	M12 float64,
	M21 float64,
	S12 float64,
) {

	azi1 = math.NaN()
	azi2 = math.NaN()
	outmask &= OUT_MASK

	a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12 := g._gen_inverse(
		lat1, lon1, lat2, lon2, outmask,
	)

	if outmask&AZIMUTH != 0 {
		azi1 = atan2_deg(salp1, calp1)
		azi2 = atan2_deg(salp2, calp2)
	}
	return a12, s12, azi1, azi2, m12, M12, M21, S12
}

func (g Geodesic) _gen_inverse(lat1, lon1, lat2, lon2 float64, outmask uint64) (
	a12 float64,
	s12 float64,
	salp1 float64,
	calp1 float64,
	salp2 float64,
	calp2 float64,
	m12 float64,
	M12 float64,
	M21 float64,
	S12 float64,
) {
	a12 = math.NaN()
	s12 = math.NaN()
	m12 = math.NaN()
	M12 = math.NaN()
	M21 = math.NaN()
	S12 = math.NaN()
	outmask &= OUT_MASK

	lon12, lon12s := ang_diff(lon1, lon2)
	var lonsign float64
	if lon12 >= 0.0 {
		lonsign = 1.0
	} else {
		lonsign = -1.0
	}

	lon12 = lonsign * ang_round(lon12)
	lon12s = ang_round((180.0 - lon12) - lonsign*lon12s)
	lam12 := lon12 * DEG2RAD
	var slam12 float64
	var clam12 float64
	if lon12 > 90.0 {
		slam12, clam12 = sincosd(lon12s)
		clam12 = -clam12
	} else {
		slam12, clam12 = sincosd(lon12)
	}
	lat1 = ang_round(lat_fix(lat1))
	lat2 = ang_round(lat_fix(lat2))

	var swapp float64
	if math.Abs(lat1) < math.Abs(lat2) {
		swapp = -1.0
	} else {
		swapp = 1.0
	}

	if swapp < 0.0 {
		lonsign *= -1.0
		lat2, lat1 = lat1, lat2
	}

	var latsign float64
	if lat1 < 0.0 {
		latsign = 1.0
	} else {
		latsign = -1.0
	}
	lat1 *= latsign
	lat2 *= latsign

	sbet1, cbet1 := sincosd(lat1)
	sbet1 *= g.f1

	sbet1, cbet1 = norm(sbet1, cbet1)
	cbet1 = math.Max(cbet1, g.tiny_)

	sbet2, cbet2 := sincosd(lat2)
	sbet2 *= g.f1

	sbet2, cbet2 = norm(sbet2, cbet2)
	cbet2 = math.Max(cbet2, g.tiny_)

	if cbet1 < -sbet1 {
		if cbet2 == cbet1 {
			if sbet2 < 0.0 {
				sbet2 = sbet1
			} else {
				sbet2 = -sbet1
			}
		}
	} else if math.Abs(sbet2) == -sbet1 {
		cbet2 = cbet1
	}

	dn1 := math.Sqrt(1.0 + g.ep2*sq(sbet1))
	dn2 := math.Sqrt(1.0 + g.ep2*sq(sbet2))

	const CARR_SIZE uint64 = uint64(_GEODESIC_ORDER) + 1
	C1a := [CARR_SIZE]float64{}
	C2a := [CARR_SIZE]float64{}
	C3a := [_GEODESIC_ORDER]float64{}

	meridian := lat1 == -90.0 || slam12 == 0.0
	calp1 = 0.0
	salp1 = 0.0
	calp2 = 0.0
	salp2 = 0.0
	ssig1 := 0.0
	csig1 := 0.0
	ssig2 := 0.0
	csig2 := 0.0
	var sig12 float64
	s12x := 0.0
	m12x := 0.0

	if meridian {
		calp1 = clam12
		salp1 = slam12
		calp2 = 1.0
		salp2 = 0.0

		ssig1 = sbet1
		csig1 = calp1 * cbet1
		ssig2 = sbet2
		csig2 = calp2 * cbet2

		sig12 = math.Atan2(math.Max((csig1*ssig2-ssig1*csig2), 0.0), csig1*csig2+ssig1*ssig2)
		res1, res2, _, res4, res5 := g._Lengths(
			g.n,
			sig12,
			ssig1,
			csig1,
			dn1,
			ssig2,
			csig2,
			dn2,
			cbet1,
			cbet2,
			outmask|DISTANCE|REDUCEDLENGTH,
			C1a[:],
			C2a[:],
		)
		s12x = res1
		m12x = res2
		M12 = res4
		M21 = res5

		if sig12 < 1.0 || m12x >= 0.0 {
			if sig12 < 3.0*g.tiny_ {
				sig12 = 0.0
				m12x = 0.0
				s12x = 0.0
			}
			m12x *= g.b
			s12x *= g.b
			a12 = sig12 * RAD2DEG
		} else {
			meridian = false
		}
	}

	somg12 := 2.0
	comg12 := 0.0
	omg12 := 0.0
	var dnm float64
	eps := 0.0

	if !meridian && sbet1 == 0.0 && (g.f <= 0.0 || lon12s >= g.f*180.0) {
		calp1 = 0.0
		calp2 = 0.0
		salp1 = 1.0
		salp2 = 1.0

		s12x = g.a * lam12
		sig12 = lam12 / g.f1
		omg12 = lam12 / g.f1
		m12x = g.b * math.Sin(sig12)
		if outmask&GEODESICSCALE != 0 {
			M12 = math.Cos(sig12)
			M21 = math.Cos(sig12)
		}
		a12 = lon12 / g.f1
	} else if !meridian {
		res1, res2, res3, res4, res5, res6 := g._InverseStart(
			sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, C1a[:], C2a[:],
		)
		sig12 = res1
		salp1 = res2
		calp1 = res3
		salp2 = res4
		calp2 = res5
		dnm = res6

		if sig12 >= 0.0 {
			s12x = sig12 * g.b * dnm
			m12x = sq(dnm) * g.b * math.Sin(sig12/dnm)
			if outmask&GEODESICSCALE != 0 {
				M12 = math.Cos(sig12 / dnm)
				M21 = math.Cos(sig12 / dnm)
			}
			a12 = sig12 * RAD2DEG
			omg12 = lam12 / (g.f1 * dnm)
		} else {
			tripn := false
			tripb := false
			salp1a := g.tiny_
			calp1a := 1.0
			salp1b := g.tiny_
			calp1b := -1.0
			domg12 := 0.0

			for numit := uint64(0); numit < g.maxit2_; numit++ {
				res1, res2, res3, res4, res5, res6, res7, res8, res9, res10, res11 := g._Lambda12(
					sbet1,
					cbet1,
					dn1,
					sbet2,
					cbet2,
					dn2,
					salp1,
					calp1,
					slam12,
					clam12,
					numit < g.maxit1_,
					C1a[:],
					C2a[:],
					C3a[:],
				)
				v := res1
				salp2 = res2
				calp2 = res3
				sig12 = res4
				ssig1 = res5
				csig1 = res6
				ssig2 = res7
				csig2 = res8
				eps = res9
				domg12 = res10
				dv := res11

				var to_mul float64
				if tripn {
					to_mul = 8.0
				} else {
					to_mul = 1.0
				}
				if tripb || !(math.Abs(v) >= to_mul*g.tol0_) {
					break
				}
				if v > 0.0 && (numit > g.maxit1_ || calp1/salp1 > calp1b/salp1b) {
					salp1b = salp1
					calp1b = calp1
				} else if v < 0.0 && (numit > g.maxit1_ || calp1/salp1 < calp1a/salp1a) {
					salp1a = salp1
					calp1a = calp1
				}
				if numit < g.maxit1_ && dv > 0.0 {
					dalp1 := -v / dv
					sdalp1 := math.Sin(dalp1)
					cdalp1 := math.Cos(dalp1)
					nsalp1 := salp1*cdalp1 + calp1*sdalp1
					if nsalp1 > 0.0 && math.Abs(dalp1) < math.Pi {
						calp1 = calp1*cdalp1 - salp1*sdalp1
						salp1 = nsalp1
						salp1, calp1 = norm(salp1, calp1)
						tripn = math.Abs(v) <= 16.0*g.tol0_
						continue
					}
				}

				salp1 = (salp1a + salp1b) / 2.0
				calp1 = (calp1a + calp1b) / 2.0
				salp1, calp1 = norm(salp1, calp1)
				tripn = false
				tripb = math.Abs(salp1a-salp1)+(calp1a-calp1) < g.tolb_ || math.Abs(salp1-salp1b)+(calp1-calp1b) < g.tolb_
			}
			var to_cmp uint64
			if outmask&(REDUCEDLENGTH|GEODESICSCALE) != 0 {
				to_cmp = DISTANCE
			} else {
				to_cmp = EMPTY
			}

			lengthmask := outmask | to_cmp
			res1, res2, _, res4, res5 = g._Lengths(
				eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, lengthmask,
				C1a[:], C2a[:],
			)
			s12x = res1
			m12x = res2
			M12 = res4
			M21 = res5

			m12x *= g.b
			s12x *= g.b
			a12 = sig12 * RAD2DEG
			if outmask&AREA != 0 {
				sdomg12 := math.Sin(domg12)
				cdomg12 := math.Cos(domg12)
				somg12 = slam12*cdomg12 - clam12*sdomg12
				comg12 = clam12*cdomg12 + slam12*sdomg12
			}
		}
	}

	if outmask&DISTANCE != 0 {
		s12 = 0.0 + s12x
	}
	if outmask&REDUCEDLENGTH != 0 {
		m12 = 0.0 + m12x
	}

	if outmask&AREA != 0 {
		salp0 := salp1 * cbet1
		calp0 := math.Hypot(calp1, salp1*sbet1)
		if calp0 != 0.0 && salp0 != 0.0 {
			ssig1 = sbet1
			csig1 = calp1 * cbet1
			ssig2 = sbet2
			csig2 = calp2 * cbet2
			k2 := sq(calp0) * g.ep2
			eps = k2 / (2.0*(1.0+math.Sqrt(1.0+k2)) + k2)
			A4 := sq(g.a) * calp0 * salp0 * g.e2
			ssig1, csig1 = norm(ssig1, csig1)
			ssig2, csig2 = norm(ssig2, csig2)
			C4a := [_GEODESIC_ORDER]float64{}
			g._C4f(eps, C4a[:])
			B41 := sin_cos_series(false, ssig1, csig1, C4a[:])
			B42 := sin_cos_series(false, ssig2, csig2, C4a[:])
			S12 = A4 * (B42 - B41)
		} else {
			S12 = 0.0
		}

		if !meridian && somg12 > 1.0 {
			somg12, comg12 = math.Sincos(omg12)
		}

		var alp12 float64
		if !meridian && comg12 > -0.7071 && sbet2-sbet1 < 1.75 {
			domg12 := 1.0 + comg12
			dbet1 := 1.0 + cbet1
			dbet2 := 1.0 + cbet2
			alp12 = 2.0 * math.Atan2(somg12*(sbet1*dbet2+sbet2*dbet1), domg12*(sbet1*sbet2+dbet1*dbet2))
		} else {
			salp12 := salp2*calp1 - calp2*salp1
			calp12 := calp2*calp1 + salp2*salp1

			if salp12 == 0.0 && calp12 < 0.0 {
				salp12 = g.tiny_ * calp1
				calp12 = -1.0
			}
			alp12 = math.Atan2(salp12, calp12)
		}
		S12 += g.c2 * alp12
		S12 *= swapp * lonsign * latsign
		S12 += 0.0
	}

	if swapp < 0.0 {
		salp2, salp1 = salp1, salp2
		calp2, calp1 = calp1, calp2
		if outmask&GEODESICSCALE != 0 {
			M21, M12 = M12, M21
		}
	}
	salp1 *= swapp * lonsign
	calp1 *= swapp * latsign
	salp2 *= swapp * lonsign
	calp2 *= swapp * latsign

	return a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12
}

// _gen_direct returns (a12, lat2, lon2, azi2, s12, m12, M12, M21, S12, outmask)
func (g Geodesic) _gen_direct(
	lat1 float64,
	lon1 float64,
	azi1 float64,
	arcmode bool,
	s12_a12 float64,
	outmask uint64,
) (float64, float64, float64, float64, float64, float64, float64, float64, float64, uint64) {
	if !arcmode {
		outmask |= DISTANCE_IN
	}

	line := NewGeodesicLineWithCapability(g, lat1, lon1, azi1, outmask)
	a12, lat2, lon2, azi2, s12, m12, M12, M21, S12 := line._gen_position(arcmode, s12_a12, outmask)

	return a12, lat2, lon2, azi2, s12, m12, M12, M21, S12, outmask
}

// LatLon represents latitude and longitude of a point. All units in degrees
type LatLon struct {
	LatDeg, LonDeg float64
}

// DirectCalcLatLon gets the lat and lon of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcLatLon(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLon {
	capabilities := LATITUDE | LONGITUDE
	_, lat2, lon2, _, _, _, _, _, _, _ := g._gen_direct(
		lat1_deg, lon1_deg, azi1_deg, false, s12_m, capabilities,
	)
	return LatLon{LatDeg: lat2, LonDeg: lon2}
}

// LatLonAzi represents latitude, longitude, and azimuth of a point. All units in degrees
type LatLonAzi struct {
	LatDeg, LonDeg, AziDeg float64
}

// DirectCalcLatLonAzi gets the lat, lon, and azimuth of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcLatLonAzi(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAzi {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH
	_, lat2, lon2, azi2, _, _, _, _, _, _ := g._gen_direct(
		lat1_deg, lon1_deg, azi1_deg, false, s12_m, capabilities,
	)
	return LatLonAzi{LatDeg: lat2, LonDeg: lon2, AziDeg: azi2}
}

// LatLonAziReducedLength: All units in degrees and meters
type LatLonAziReducedLength struct {
	LatDeg         float64 // Latitude [degrees]
	LonDeg         float64 // Longitude [degrees]
	AziDeg         float64 // Azimuth [degrees]
	ReducedLengthM float64 // Reduced length of the geodesic [meters]
}

// DirectCalcLatLonAziReducedLength gets the lat, lon, azimuth, and reduced length of geodesic
// of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcLatLonAziReducedLength(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAziReducedLength {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH | REDUCEDLENGTH
	_, lat2, lon2, azi2, _, m12, _, _, _, _ := g._gen_direct(
		lat1_deg, lon1_deg, azi1_deg, false, s12_m, capabilities,
	)
	return LatLonAziReducedLength{LatDeg: lat2, LonDeg: lon2, AziDeg: azi2, ReducedLengthM: m12}
}

// LatLonAziGeodesicScales: All units in degrees. Scales are dimensionless
type LatLonAziGeodesicScales struct {
	LatDeg float64 // Latitude [degrees]
	LonDeg float64 // Longitude [degrees]
	AziDeg float64 // Azimuth [degrees]
	M12    float64 // Geodesic scale of point 2 relative to point 1 [dimensionless]
	M21    float64 // Geodesic scale of point 1 relative to point 2 [dimensionless]
}

// DirectCalcLatLonAziGeodesicScales gets the lat, lon, azimuth, and geodesic scales,
// based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcLatLonAziGeodesicScales(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAziGeodesicScales {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH | GEODESICSCALE
	_, lat2, lon2, azi2, _, _, M12, M21, _, _ := g._gen_direct(
		lat1_deg, lon1_deg, azi1_deg, false, s12_m, capabilities,
	)
	return LatLonAziGeodesicScales{LatDeg: lat2, LonDeg: lon2, AziDeg: azi2, M12: M12, M21: M21}
}

// LatLonAziReducedLengthGeodesicScales: All units in degrees and meters. Scales are dimensionless
type LatLonAziReducedLengthGeodesicScales struct {
	LatDeg         float64 // Latitude [degrees]
	LonDeg         float64 // Longitude [degrees]
	AziDeg         float64 // Azimuth [degrees]
	ReducedLengthM float64 // Reduced length of the geodesic [meters]
	M12            float64 // Geodesic scale of point 2 relative to point 1 [dimensionless]
	M21            float64 // Geodesic scale of point 1 relative to point 2 [dimensionless]
}

// DirectCalcLatLonAziReducedLengthGeodesicScales gets the lat, lon, azimuth, reduced length,
// and geodesic scales based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcLatLonAziReducedLengthGeodesicScales(
	lat1_deg, lon1_deg, azi1_deg, s12_m float64,
) LatLonAziReducedLengthGeodesicScales {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH | REDUCEDLENGTH | GEODESICSCALE
	_, lat2, lon2, azi2, _, m12, M12, M21, _, _ := g._gen_direct(
		lat1_deg, lon1_deg, azi1_deg, false, s12_m, capabilities,
	)
	return LatLonAziReducedLengthGeodesicScales{
		LatDeg:         lat2,
		LonDeg:         lon2,
		AziDeg:         azi2,
		ReducedLengthM: m12,
		M12:            M12,
		M21:            M21,
	}
}

// AllDirectResults contains all information that can be computed from the direct method
// latitude, longitude, azimuth, reduced length, geodesic scales, area under the geodesic,
// and arc length between point 1 and point 2
type AllDirectResults struct {
	LatDeg         float64 // Latitude [degrees]
	LonDeg         float64 // Longitude [degrees]
	AziDeg         float64 // Azimuth [degrees]
	ReducedLengthM float64 // Reduced length of the geodesic [meters]
	M12            float64 // Geodesic scale of point 2 relative to point 1 [dimensionless]
	M21            float64 // Geodesic scale of point 1 relative to point 2 [dimensionless]
	S12M2          float64 // Area under the geodesic [meters^2]
	A12Deg         float64 // Arc length between point 1 and point 2 [degrees]
}

// DirectCalcAll calculates everything possible for the direct method. Takes inputs
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcAll(lat1_deg, lon1_deg, azi1_deg, s12_m float64) AllDirectResults {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH | REDUCEDLENGTH | GEODESICSCALE | AREA
	a12, lat2, lon2, azi2, _, m12, M12, M21, S12, _ := g._gen_direct(
		lat1_deg, lon1_deg, azi1_deg, false, s12_m, capabilities,
	)
	return AllDirectResults{
		LatDeg:         lat2,
		LonDeg:         lon2,
		AziDeg:         azi2,
		ReducedLengthM: m12,
		M12:            M12,
		M21:            M21,
		S12M2:          S12,
		A12Deg:         a12,
	}
}

// DirectCalcWithCapabilities allows the user to specify which capabilites they wish to use.
// This function is useful if you want some other subset of capabilities than those offered
// by the other DirectCalc...() methods.
// Takes inputs
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - capabilities - One or more of the capabilities constant as defined in the file
//     geodesiccapability.go. Usually, they are OR'd together, e.g. LATITUDE | LONGITUDE
func (g Geodesic) DirectCalcWithCapabilities(
	lat1_deg, lon1_deg, azi1_deg, s12_m float64,
	capabilities uint64,
) AllDirectResults {
	a12, lat2, lon2, azi2, _, m12, M12, M21, S12, _ := g._gen_direct(
		lat1_deg, lon1_deg, azi1_deg, false, s12_m, capabilities,
	)
	return AllDirectResults{
		LatDeg:         lat2,
		LonDeg:         lon2,
		AziDeg:         azi2,
		ReducedLengthM: m12,
		M12:            M12,
		M21:            M21,
		S12M2:          S12,
		A12Deg:         a12,
	}
}

// InverseCalcDistance returns the distance from point 1 to point 2 in meters. Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g Geodesic) InverseCalcDistance(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) float64 {
	capabilities := DISTANCE
	_, s12, _, _, _, _, _, _ := g._gen_inverse_azi(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return s12
}

type DistanceArcLength struct {
	DistanceM    float64 // distance between point 1 and point 2 [meters]
	ArcLengthDeg float64 // arc length between point 1 and point 2 [degrees]
}

// InverseCalcDistanceArcLength returns the distance from one point to the next, and the
// arc length between the points. Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g Geodesic) InverseCalcDistanceArcLength(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) DistanceArcLength {
	capabilities := DISTANCE
	a12, s12, _, _, _, _, _, _ := g._gen_inverse_azi(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return DistanceArcLength{DistanceM: s12, ArcLengthDeg: a12}
}

type DistanceAzimuths struct {
	DistanceM   float64 // distance between point 1 and point 2 [meters]
	Azimuth1Deg float64 // azimuth at point 1 [degrees]
	Azimuth2Deg float64 // (forward) azimuth at point 2 [degrees]
}

// InverseCalcDistanceAzimuths returns the distance from one point to the next, and the
// azimuths. Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g Geodesic) InverseCalcDistanceAzimuths(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) DistanceAzimuths {
	capabilities := DISTANCE | AZIMUTH
	_, s12, azi1, azi2, _, _, _, _ := g._gen_inverse_azi(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return DistanceAzimuths{DistanceM: s12, Azimuth1Deg: azi1, Azimuth2Deg: azi2}
}

type AzimuthsArcLength struct {
	Azimuth1Deg  float64 // azimuth at point 1 [degrees]
	Azimuth2Deg  float64 // (forward) azimuth at point 2 [degrees]
	ArcLengthDeg float64 // arc length between point 1 and point 2 [degrees]
}

// InverseCalcAzimuthsArcLength returns the azimuth at point 1, the azimuth at point 2,
// and the arc length between the points. Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g Geodesic) InverseCalcAzimuthsArcLength(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) AzimuthsArcLength {
	capabilities := AZIMUTH
	a12, _, azi1, azi2, _, _, _, _ := g._gen_inverse_azi(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return AzimuthsArcLength{Azimuth1Deg: azi1, Azimuth2Deg: azi2, ArcLengthDeg: a12}
}

type DistanceAzimuthsArcLength struct {
	DistanceM    float64 // distance between point 1 and point 2 [meters]
	Azimuth1Deg  float64 // azimuth at point 1 [degrees]
	Azimuth2Deg  float64 // (forward) azimuth at point 2 [degrees]
	ArcLengthDeg float64 // arc length between point 1 and point 2 [degrees]
}

// InverseCalcDistanceAzimuthsArcLength returns the distance from one point to the next,
// the azimuth at point 1, the azimuth at point 2, and the arc length between the points.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g Geodesic) InverseCalcDistanceAzimuthsArcLength(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) DistanceAzimuthsArcLength {
	capabilities := DISTANCE | AZIMUTH
	a12, s12, azi1, azi2, _, _, _, _ := g._gen_inverse_azi(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return DistanceAzimuthsArcLength{
		DistanceM: s12, Azimuth1Deg: azi1, Azimuth2Deg: azi2, ArcLengthDeg: a12,
	}
}

type DistanceAzimuthsArcLengthReducedLength struct {
	DistanceM      float64 // distance between point 1 and point 2 [meters]
	Azimuth1Deg    float64 // azimuth at point 1 [degrees]
	Azimuth2Deg    float64 // (forward) azimuth at point 2 [degrees]
	ArcLengthDeg   float64 // arc length between point 1 and point 2 [degrees]
	ReducedLengthM float64 // reduced length of geodesic [meters]
}

// InverseCalcDistanceAzimuthsArcLengthReducedLength returns the distance from one point
// to the next, the azimuth at point 1, the azimuth at point 2, the arc length
// between the points, and the reduceed length of the geodesic.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g Geodesic) InverseCalcDistanceAzimuthsArcLengthReducedLength(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) DistanceAzimuthsArcLengthReducedLength {
	capabilities := DISTANCE | AZIMUTH | REDUCEDLENGTH
	a12, s12, azi1, azi2, m12, _, _, _ := g._gen_inverse_azi(
		lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities,
	)

	return DistanceAzimuthsArcLengthReducedLength{
		DistanceM:      s12,
		Azimuth1Deg:    azi1,
		Azimuth2Deg:    azi2,
		ArcLengthDeg:   a12,
		ReducedLengthM: m12,
	}
}

type DistanceAzimuthsArcLengthReducedLengthScales struct {
	DistanceM      float64 // distance between point 1 and point 2 [meters]
	Azimuth1Deg    float64 // azimuth at point 1 [degrees]
	Azimuth2Deg    float64 // (forward) azimuth at point 2 [degrees]
	ArcLengthDeg   float64 // arc length between point 1 and point 2 [degrees]
	ReducedLengthM float64 // reduced length of geodesic [meters]
	M12            float64 // geodesic scale of point 2 relative to point 1 [dimensionless]
	M21            float64 // geodesic scale of point 1 relative to point 2 [dimensionless]
}

// InverseCalcDistanceAzimuthsArcLengthReducedLengthScales returns everything described
// by the `DistanceAzimuthsArcLengthReducedLengthScales` type.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g Geodesic) InverseCalcDistanceAzimuthsArcLengthReducedLengthScales(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) DistanceAzimuthsArcLengthReducedLengthScales {
	capabilities := DISTANCE | AZIMUTH | REDUCEDLENGTH | GEODESICSCALE
	a12, s12, azi1, azi2, m12, M12, M21, _ := g._gen_inverse_azi(
		lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities,
	)

	return DistanceAzimuthsArcLengthReducedLengthScales{
		DistanceM:      s12,
		Azimuth1Deg:    azi1,
		Azimuth2Deg:    azi2,
		ArcLengthDeg:   a12,
		ReducedLengthM: m12,
		M12:            M12,
		M21:            M21,
	}
}

type AllInverseResults struct {
	DistanceM      float64 // distance between point 1 and point 2 [meters]
	Azimuth1Deg    float64 // azimuth at point 1 [degrees]
	Azimuth2Deg    float64 // (forward) azimuth at point 2 [degrees]
	ArcLengthDeg   float64 // arc length between point 1 and point 2 [degrees]
	ReducedLengthM float64 // reduced length of geodesic [meters]
	M12            float64 // geodesic scale of point 2 relative to point 1 [dimensionless]
	M21            float64 // geodesic scale of point 1 relative to point 2 [dimensionless]
	S12M2          float64 // area under the geodesic [meters^2]
}

// InverseCalcAll returns everything described in the `AllInverseResults` results type.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g Geodesic) InverseCalcAll(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) AllInverseResults {
	capabilities := DISTANCE | AZIMUTH | REDUCEDLENGTH | GEODESICSCALE | AREA
	a12, s12, azi1, azi2, m12, M12, M21, S12 := g._gen_inverse_azi(
		lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities,
	)

	return AllInverseResults{
		DistanceM:      s12,
		Azimuth1Deg:    azi1,
		Azimuth2Deg:    azi2,
		ArcLengthDeg:   a12,
		ReducedLengthM: m12,
		M12:            M12,
		M21:            M21,
		S12M2:          S12,
	}
}

// InverseCalcWithCapabilities allows the user to specify which capabilites they wish to use.
// This function is useful if you want some other subset of capabilities than those offered
// by the other InverseCalc...() methods.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
// - capabilities - One or more of the capabilities constant as defined in the file
//     geodesiccapability.go. Usually, they are OR'd together, e.g. LATITUDE | LONGITUDE
func (g Geodesic) InverseCalcWithCapabilities(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
	capabilities uint64,
) AllInverseResults {
	a12, s12, azi1, azi2, m12, M12, M21, S12 := g._gen_inverse_azi(
		lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities,
	)

	return AllInverseResults{
		DistanceM:      s12,
		Azimuth1Deg:    azi1,
		Azimuth2Deg:    azi2,
		ArcLengthDeg:   a12,
		ReducedLengthM: m12,
		M12:            M12,
		M21:            M21,
		S12M2:          S12,
	}
}

// InverseLineWithCapabilities: define a GeodesicLine struct in terms of the inverse geodesic
// problem.
// This function sets point 3 of the GeodesicLine to correspond to point 2 of the
// inverse geodesic problem.
func (g Geodesic) InverseLineWithCapabilities(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
	capabilities uint64,
) GeodesicLine {
	a12, _, salp1, calp1, _, _, _, _, _, _ := g._gen_inverse(lat1_deg, lon1_deg, lat2_deg, lon2_deg, 0)
	azi1 := atan2_deg(salp1, calp1)

	if capabilities&(OUT_MASK&DISTANCE_IN) != 0 {
		capabilities |= DISTANCE
	}

	line := new_geodesic_line_all_options(g, lat1_deg, lon1_deg, azi1, capabilities, salp1, calp1)
	line.set_arc(a12)
	return line
}

func (g Geodesic) _gen_direct_line(
	lat1_deg, lon1_deg, azi1_deg float64,
	arcmode bool,
	s12_a12_m float64,
	capabilities uint64,
) GeodesicLine {
	// Automatically supply DISTANCE_IN if necessary
	if !arcmode {
		capabilities |= DISTANCE_IN
	}
	line := NewGeodesicLineWithCapability(
		g,
		lat1_deg,
		lon1_deg,
		azi1_deg,
		capabilities,
	)
	if arcmode {
		line.set_arc(s12_a12_m)
	} else {
		line.set_distance(s12_a12_m)
	}
	return line
}

// DirectLineWithCapabilities defines a GeodesicLine struct in terms of the the direct
// geodesic problem specified in terms of spherical arc length.
// This function sets point 3 of the GeodesicLine to correspond to point 2 of the
// direct geodesic problem
func (g Geodesic) DirectLineWithCapabilities(
	lat1_deg, lon1_deg, azi1_deg, s12_m float64,
	capabilities uint64,
) GeodesicLine {
	return g._gen_direct_line(lat1_deg, lon1_deg, azi1_deg, false, s12_m, capabilities)
}

// Line returns a GeodesicLine. This allows points along a geodesic starting at
// lat1_deg, lon1_deg with azimuth azi1_deg to be found.
func (g Geodesic) LineWithCapabilities(
	lat1_deg, lon1_deg, azi1_deg float64,
	capabilities uint64,
) GeodesicLine {
	return NewGeodesicLineWithCapability(
		g,
		lat1_deg,
		lon1_deg,
		azi1_deg,
		capabilities,
	)
}
