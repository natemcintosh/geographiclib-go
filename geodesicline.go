package geographiclibgo

import "math"

type GeodesicLine struct {
	tiny_  float64
	_A1m1  float64
	_A2m1  float64
	_A3c   float64
	_A4    float64
	_B11   float64
	_B21   float64
	_B31   float64
	_B41   float64
	_C1a   [GEODESIC_ORDER + 1]float64
	_C1pa  [GEODESIC_ORDER + 1]float64
	_C2a   [GEODESIC_ORDER + 1]float64
	_C3a   [GEODESIC_ORDER]float64
	_C4a   [GEODESIC_ORDER]float64
	_b     float64
	_c2    float64
	_calp0 float64
	_csig1 float64
	_comg1 float64
	_ctau1 float64
	_dn1   float64
	_f1    float64
	_k2    float64
	_salp0 float64
	_somg1 float64
	_ssig1 float64
	_stau1 float64
	a13    float64
	a      float64
	azi1   float64
	calp1  float64
	caps   uint64
	f      float64
	lat1   float64
	lon1   float64
	s13    float64
	salp1  float64
}

// NewGeodesicLine creates a GeodesicLine, with `caps` of STANDARD | DISTANCE_IN
// If you do not wish to specify `salp1` and/or `calp1`, set them as math.NaN()
func NewGeodesicLine(
	geod Geodesic,
	lat1, lon1, azi1 float64,
	salp1, calp1 float64,
) GeodesicLine {
	// Specify default `caps`
	caps := STANDARD | DISTANCE_IN

	return NewGeodesicLineWithCaps(geod, lat1, lon1, azi1, caps, salp1, calp1)
}

// NewGeodesicLineWithCaps is the same as NewGeodesicLine but the user specifies a
// `caps` field
func NewGeodesicLineWithCaps(
	geod Geodesic,
	lat1, lon1, azi1 float64,
	caps uint64,
	salp1, calp1 float64,
) GeodesicLine {
	// This was taken from geodesic, putting it here for convenience
	tiny_ := math.Sqrt(Get_min_val())

	a := geod.a
	f := geod.f
	_b := geod._b
	_c2 := geod._c2
	_f1 := geod._f1
	caps |= LATITUDE | AZIMUTH | LONG_UNROLL

	if math.IsNaN(salp1) || math.IsNaN(calp1) {
		azi1 = Ang_normalize(azi1)
		salp1, calp1 = Sincosd(Ang_round(azi1))
	}

	lat1 = Lat_fix(lat1)

	sbet1, cbet1 := Sincosd(Ang_round(lat1))
	sbet1 *= _f1
	sbet1, cbet1 = Norm(sbet1, cbet1)
	cbet1 = math.Max(tiny_, cbet1)
	_dn1 := math.Sqrt(1.0 + geod._ep2*Sq(sbet1))
	_salp0 := salp1 * cbet1
	_calp0 := math.Hypot(calp1, salp1*sbet1)
	_ssig1 := sbet1
	_somg1 := _salp0 * sbet1

	var _csig1 float64
	if sbet1 != 0.0 || calp1 != 0.0 {
		_csig1 = cbet1 * calp1
	} else {
		_csig1 = 1.0
	}

	_comg1 := _csig1
	_ssig1, _csig1 = Norm(_ssig1, _csig1)
	_k2 := Sq(_calp0) * geod._ep2
	eps := _k2 / (2.0*(1.0+math.Sqrt(1.0+_k2)) + _k2)

	_A1m1 := 0.0
	var _C1a [GEODESIC_ORDER + 1]float64
	_B11 := 0.0
	_stau1 := 0.0
	_ctau1 := 0.0

	if caps&CAP_C1 != 0 {
		_A1m1 = _A1m1f(eps, geod.GEODESIC_ORDER)
		_C1f(eps, _C1a[:], int(geod.GEODESIC_ORDER))
		_B11 = Sin_cos_series(true, _ssig1, _csig1, _C1a[:])
		s := math.Sin(_B11)
		c := math.Cos(_B11)
		_stau1 = _ssig1*c + _csig1*s
		_ctau1 = _csig1*c - _ssig1*s
	}

	var _C1pa [GEODESIC_ORDER + 1]float64
	if caps&CAP_C1p != 0 {
		_C1pf(eps, _C1pa[:], int(geod.GEODESIC_ORDER))
	}

	_A2m1 := 0.0
	var _C2a [GEODESIC_ORDER + 1]float64
	_B21 := 0.0
	if caps&CAP_C2 != 0 {
		_A2m1 = _A2m1f(eps, geod.GEODESIC_ORDER)
		_C2f(eps, _C2a[:], int(geod.GEODESIC_ORDER))
		_B21 = Sin_cos_series(true, _ssig1, _csig1, _C2a[:])
	}

	var _C3a [GEODESIC_ORDER]float64
	_A3c := 0.0
	_B31 := 0.0
	if caps&CAP_C3 != 0 {
		geod._C3f(eps, _C3a[:])
		_A3c = -f * _salp0 * geod._A3f(eps)
		_B31 = Sin_cos_series(true, _ssig1, _csig1, _C3a[:])
	}

	var _C4a [GEODESIC_ORDER]float64
	_A4 := 0.0
	_B41 := 0.0
	if caps&CAP_C4 != 0 {
		geod._C4f(eps, _C4a[:])
		_A4 = Sq(a) * _calp0 * _salp0 * geod._e2
		_B41 = Sin_cos_series(false, _ssig1, _csig1, _C4a[:])
	}

	s13 := math.NaN()
	a13 := math.NaN()

	return GeodesicLine{
		tiny_:  tiny_,
		_A1m1:  _A1m1,
		_A2m1:  _A2m1,
		_A3c:   _A3c,
		_A4:    _A4,
		_B11:   _B11,
		_B21:   _B21,
		_B31:   _B31,
		_B41:   _B41,
		_C1a:   _C1a,
		_C1pa:  _C1pa,
		_comg1: _comg1,
		_C2a:   _C2a,
		_C3a:   _C3a,
		_C4a:   _C4a,
		_b:     _b,
		_c2:    _c2,
		_calp0: _calp0,
		_csig1: _csig1,
		_ctau1: _ctau1,
		_dn1:   _dn1,
		_f1:    _f1,
		_k2:    _k2,
		_salp0: _salp0,
		_somg1: _somg1,
		_ssig1: _ssig1,
		_stau1: _stau1,
		a:      a,
		a13:    a13,
		azi1:   azi1,
		calp1:  calp1,
		caps:   caps,
		f:      f,
		lat1:   lat1,
		lon1:   lon1,
		s13:    s13,
		salp1:  salp1,
	}
}

// _gen_position returns (a12, lat2, lon2, azi2, s12, m12, M12, M21, S12)
func (g GeodesicLine) _gen_position(arcmode bool, s12_a12 float64, outmask uint64) (
	float64, float64, float64, float64, float64, float64, float64, float64, float64,
) {
	a12 := math.NaN()
	lat2 := math.NaN()
	lon2 := math.NaN()
	azi2 := math.NaN()
	s12 := math.NaN()
	m12 := math.NaN()
	M12 := math.NaN()
	M21 := math.NaN()
	S12 := math.NaN()

	outmask &= g.caps & OUT_MASK
	if !(arcmode || (g.caps&(OUT_MASK&DISTANCE_IN) != 0)) {
		return a12, lat2, lon2, azi2, s12, m12, M12, M21, S12
	}

	B12 := 0.0
	AB1 := 0.0
	var sig12 float64
	var ssig12 float64
	var csig12 float64
	var ssig2 float64
	var csig2 float64

	if arcmode {
		sig12 = s12_a12 * DEG2RAD
		ssig12, csig12 = Sincosd(s12_a12)

	} else {
		// tau12 = s12_a12 / (g._b * (1 + g._A1m1))
		tau12 := s12_a12 / (g._b * (1.0 + g._A1m1))

		s := math.Sin(tau12)
		c := math.Cos(tau12)

		B12 = -Sin_cos_series(
			true,
			g._stau1*c+g._ctau1*s,
			g._ctau1*c-g._stau1*s,
			g._C1pa[:],
		)
		sig12 = tau12 - (B12 - g._B11)
		ssig12 = math.Sin(sig12)
		csig12 = math.Cos(sig12)
		if math.Abs(g.f) > 0.01 {
			ssig2 = g._ssig1*csig12 + g._csig1*ssig12
			csig2 = g._csig1*csig12 - g._ssig1*ssig12
			B12 = Sin_cos_series(true, ssig2, csig2, g._C1a[:])
			serr := (1.0+g._A1m1)*(sig12+(B12-g._B11)) - s12_a12/g._b
			sig12 = sig12 - serr/math.Sqrt(1.0+g._k2*Sq(ssig2))
			ssig12 = math.Sin(sig12)
			csig12 = math.Cos(sig12)
		}
	}

	ssig2 = g._ssig1*csig12 + g._csig1*ssig12
	csig2 = g._csig1*csig12 - g._ssig1*ssig12
	dn2 := math.Sqrt(1.0 + g._k2*Sq(ssig2))
	if outmask&(DISTANCE|REDUCEDLENGTH|GEODESICSCALE) != 0 {
		if arcmode || math.Abs(g.f) > 0.01 {
			B12 = Sin_cos_series(true, ssig2, csig2, g._C1a[:])
		}
		AB1 = (1.0 + g._A1m1) * (B12 - g._B11)
	}

	sbet2 := g._calp0 * ssig2
	cbet2 := math.Hypot(g._salp0, g._calp0*csig2)
	if cbet2 == 0.0 {
		cbet2 = g.tiny_
		csig2 = g.tiny_
	}
	salp2 := g._salp0
	calp2 := g._calp0 * csig2

	if outmask&DISTANCE != 0 {
		if arcmode {
			s12 = g._b * ((1.0+g._A1m1)*sig12 + AB1)
		} else {
			s12 = s12_a12
		}
	}

	if outmask&LONGITUDE != 0 {
		somg2 := g._salp0 * ssig2
		comg2 := csig2

		var E float64
		if g._salp0 < 0.0 {
			E = -1.0
		} else {
			E = 1.0
		}

		var omg12 float64
		if outmask&LONG_UNROLL != 0 {
			omg12 = E * (sig12 - (math.Atan2(ssig2, csig2) - math.Atan2(g._ssig1, g._csig1)) + (math.Atan2((E*somg2), comg2) - math.Atan2((E*g._somg1), g._comg1)))
		} else {
			omg12 = math.Atan2((somg2*g._comg1 - comg2*g._somg1), (comg2*g._comg1 + somg2*g._somg1))
		}
		lam12 := omg12 + g._A3c*(sig12+(Sin_cos_series(true, ssig2, csig2, g._C3a[:])-g._B31))
		lon12 := lam12 * RAD2DEG

		if outmask&LONG_UNROLL != 0 {
			lon2 = g.lon1 + lon12
		} else {
			lon2 = Ang_normalize(
				Ang_normalize(g.lon1) + Ang_normalize(lon12),
			)
		}
	}

	if outmask&LATITUDE != 0 {
		lat2 = Atan2_deg(sbet2, g._f1*cbet2)
	}

	if outmask&AZIMUTH != 0 {
		azi2 = Atan2_deg(salp2, calp2)
	}

	if outmask&(REDUCEDLENGTH|GEODESICSCALE) != 0 {
		B22 := Sin_cos_series(true, ssig2, csig2, g._C2a[:])
		AB2 := (1.0 + g._A2m1) * (B22 - g._B21)
		J12 := (g._A1m1-g._A2m1)*sig12 + (AB1 - AB2)
		if outmask&REDUCEDLENGTH != 0 {
			m12 = g._b * ((dn2*(g._csig1*ssig2) - g._dn1*(g._ssig1*csig2)) - g._csig1*csig2*J12)
		}
		if outmask&GEODESICSCALE != 0 {
			t := g._k2 * (ssig2 - g._ssig1) * (ssig2 + g._ssig1) / (g._dn1 + dn2)
			M12 = csig12 + (t*ssig2-csig2*J12)*g._ssig1/g._dn1
			M21 = csig12 - (t*g._ssig1-g._csig1*J12)*ssig2/dn2
		}
	}

	if outmask&AREA != 0 {
		B42 := Sin_cos_series(false, ssig2, csig2, g._C4a[:])
		var salp12 float64
		var calp12 float64
		if g._calp0 == 0.0 || g._salp0 == 0.0 {
			salp12 = salp2*g.calp1 - calp2*g.salp1
			calp12 = calp2*g.calp1 + salp2*g.salp1
		} else {
			var to_mul float64
			if csig12 <= 0.0 {
				to_mul = g._csig1*(1.0-csig12) + ssig12*g._ssig1
			} else {
				to_mul = ssig12 * (g._csig1*ssig12/(1.0+csig12) + g._ssig1)
			}
			salp12 = g._calp0 * g._salp0 * to_mul

			calp12 = Sq(g._salp0) + Sq(g._calp0)*g._csig1*csig2
		}
		S12 = g._c2*math.Atan2(salp12, calp12) + g._A4*(B42-g._B41)
	}

	if arcmode {
		a12 = s12_a12
	} else {
		a12 = sig12 * RAD2DEG
	}
	return a12, lat2, lon2, azi2, s12, m12, M12, M21, S12
}

type PositionResultStandard struct {
	Lat1Deg   float64 // Latitude of point 1 [degrees]
	Lon1Deg   float64 // Longitude of point 1 [degrees]
	Azi1Deg   float64 // Azimuth of point 1 [degrees]
	Lat2Deg   float64 // Latitude of point 2 [degrees]
	Lon2Deg   float64 // Longitude of point 2 [degrees]
	Azi2Deg   float64 // Azimuth of point 2 [degrees]
	DistanceM float64 // Distance from point 1 to point 2 [meters]
}

// PositionStandard finds the position on the line given s12_m [meters]. It uses the
// STANDARD capabilities, and returns a PositionResultStandard struct
func (g GeodesicLine) PositionStandard(s12_m float64) PositionResultStandard {
	outmask := STANDARD
	_, lat2, lon2, azi2, _, _, _, _, _ := g._gen_position(false, s12_m, outmask)

	return PositionResultStandard{
		Lat1Deg:   g.lat1,
		Lon1Deg:   g.lon1,
		Azi1Deg:   g.azi1,
		Lat2Deg:   lat2,
		Lon2Deg:   lon2,
		Azi2Deg:   azi2,
		DistanceM: s12_m,
	}
}

type PositionResult struct {
	Lat1Deg        float64 // Latitude of point 1 [degrees]
	Lon1Deg        float64 // Longitude of point 1 [degrees]
	Azi1Deg        float64 // Azimuth of point 1 [degrees]
	Lat2Deg        float64 // Latitude of point 2 [degrees]
	Lon2Deg        float64 // Longitude of point 2 [degrees]
	Azi2Deg        float64 // Azimuth of point 2 [degrees]
	DistanceM      float64 // Distance from point 1 to point 2 [meters]
	ArcLengthDeg   float64 // Arc length between point 1 and point 2 [degrees]
	ReducedLengthM float64 // Reduced length of the geodesic [meters]
	M12            float64 // Geodesic scale of point 2 relative to point 1 [dimensionless]
	M21            float64 // Geodesic scale of point 1 relative to point 2 [dimensionless]
	S12M2          float64 // Area under the geodesic [meters^2]
}

// PositionWithCapabilities finds the position on the line given s12_m [meters]. It uses
// whatever capabilities are handed in. Any results not asked for with the capabilities
// will be math.NaN()
func (g GeodesicLine) PositionWithCapabilities(s12_m float64, outmask uint64) PositionResult {
	a12, lat2, lon2, azi2, s12, m12, M12, M21, S12 := g._gen_position(false, s12_m, outmask)

	var outlon1 float64
	if outmask&LONG_UNROLL != 0 {
		outlon1 = g.lon1
	} else {
		outlon1 = Ang_normalize(g.lon1)
	}

	return PositionResult{
		Lat1Deg:        g.lat1,
		Lon1Deg:        outlon1,
		Azi1Deg:        g.azi1,
		Lat2Deg:        lat2,
		Lon2Deg:        lon2,
		Azi2Deg:        azi2,
		DistanceM:      s12,
		ArcLengthDeg:   a12,
		ReducedLengthM: m12,
		M12:            M12,
		M21:            M21,
		S12M2:          S12,
	}
}
