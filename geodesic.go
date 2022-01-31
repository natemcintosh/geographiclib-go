package geographiclibgo

import "math"

const WGS84_A float64 = 6378137.0

// Evaluating this as 1000000000.0 / (298257223563f_64) reduces the
// round-off error by about 10%.  However, expressing the flattening as
// 1/298.257223563 is well ingrained.
const WGS84_F float64 = 1.0 / ((298257223563.0) / 1000000000.0)

const GEODESIC_ORDER int64 = 6
const nC3x_ int64 = 15
const nC4x_ int64 = 21

func COEFF_A3() [18]float64 {
	return [18]float64{
		-3.0, 128.0, -2.0, -3.0, 64.0, -1.0, -3.0, -1.0, 16.0, 3.0, -1.0, -2.0, 8.0, 1.0, -1.0, 2.0,
		1.0, 1.0,
	}
}

func COEFF_C3() [45]float64 {
	return [45]float64{
		3.0, 128.0, 2.0, 5.0, 128.0, -1.0, 3.0, 3.0, 64.0, -1.0, 0.0, 1.0, 8.0, -1.0, 1.0, 4.0, 5.0,
		256.0, 1.0, 3.0, 128.0, -3.0, -2.0, 3.0, 64.0, 1.0, -3.0, 2.0, 32.0, 7.0, 512.0, -10.0, 9.0,
		384.0, 5.0, -9.0, 5.0, 192.0, 7.0, 512.0, -14.0, 7.0, 512.0, 21.0, 2560.0,
	}
}

func COEFF_C4() [77]float64 {
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
	a      float64
	f      float64
	_f1    float64
	_e2    float64
	_ep2   float64
	_n     float64
	_b     float64
	_c2    float64
	_etol2 float64

	GEODESIC_ORDER int64
	nC3x_          int64
	nC4x_          int64
	maxit1_        uint64
	maxit2_        uint64

	_A3x [GEODESIC_ORDER]float64
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
	maxit2_ := maxit1_ + DIGITS + 10
	tiny_ := math.Sqrt(Get_min_val())
	tol0_ := Get_epsilon()
	tol1_ := 200.0 * tol0_
	tol2_ := math.Sqrt(tol0_)
	tolb_ := tol0_ * tol2_
	xthresh_ := 1000.0 * tol2_

	_f1 := 1.0 - f
	_e2 := f * (2.0 - f)
	_ep2 := _e2 / Sq(_f1)
	_n := f / (2.0 - f)
	_b := a * _f1

	var is_f_neg float64
	if f < 0.0 {
		is_f_neg = -1.0
	} else {
		is_f_neg = 1.0
	}

	to_mul := Eatanhe(1.0, is_f_neg*math.Sqrt(math.Abs(_e2))) / _e2
	if _e2 == 0.0 {
		to_mul = 1.0
	}
	_c2 := (Sq(a) + Sq(_b)*to_mul) / 2.0
	_etol2 := 0.1 * tol2_ / math.Sqrt(math.Max(math.Abs(f), 0.001)*math.Min((1.0-f/2.0), 1.0)/2.0)

	_A3x := [GEODESIC_ORDER]float64{}
	_C3x := [nC3x_]float64{}
	_C4x := [nC4x_]float64{}

	// Call a3coeff
	var o int64 = 0
	k := 0

	coefa3 := COEFF_A3()
	for j := GEODESIC_ORDER - 1; j >= 0; j-- {
		m := int64(math.Min(float64(j), float64(GEODESIC_ORDER-j-1)))
		_A3x[k] = Polyval(m, coefa3[o:], _n) / coefa3[o+m+1]
		k += 1
		o += m + 2
	}

	// c3coeff
	o = 0
	k = 0

	coefc3 := COEFF_C3()
	for l := 1; l < int(GEODESIC_ORDER); l++ {
		for j := int(GEODESIC_ORDER) - 1; j >= l; j-- {
			m := int64(math.Min(float64(j), float64(int(GEODESIC_ORDER)-j-1)))
			_C3x[k] = Polyval(m, coefc3[o:], _n) / coefc3[o+m+1]
			k += 1
			o += m + 2
		}

	}

	// c4coeff
	o = 0
	k = 0

	coefc4 := COEFF_C4()
	for l := 0; l < int(GEODESIC_ORDER); l++ {
		for j := int(GEODESIC_ORDER) - 1; j >= l; j-- {
			m := int64(int(GEODESIC_ORDER) - j - 1)
			_C4x[k] = Polyval(m, coefc4[o:], _n) / coefc4[(o+m+1)]
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

		GEODESIC_ORDER,
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
	return Polyval(int64(GEODESIC_ORDER-1), g._A3x[:], eps)
}

func (g Geodesic) _C3f(eps float64, c []float64) {
	mult := 1.0
	o := 0
	for l := 1; l < int(GEODESIC_ORDER); l++ {
		m := int(GEODESIC_ORDER) - l - 1
		mult *= eps
		c[l] = mult * Polyval(int64(m), g._C3x[o:], eps)
		o += m + 1
	}
}

func (g Geodesic) _C4f(eps float64, c []float64) {
	mult := 1.0
	o := 0
	for l := 0; l < int(GEODESIC_ORDER); l++ {
		m := int(GEODESIC_ORDER) - l - 1
		c[l] = mult * Polyval(int64(m), g._C4x[o:], eps)
		o += m + 1
		mult *= eps
	}
}

func (g Geodesic) _Lengths(
	eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2 float64,
	outmask uint64,
	C1a []float64,
	C2a []float64,
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
		A1 = _A1m1f(eps, GEODESIC_ORDER)
		_C1f(eps, C1a, int(GEODESIC_ORDER))
		if outmask&(REDUCEDLENGTH|GEODESICSCALE) != 0 {
			A2 = _A2m1f(eps, GEODESIC_ORDER)
			_C2f(eps, C2a, int(GEODESIC_ORDER))
			m0x = A1 - A2
			A2 = 1.0 + A2
		}
		A1 = 1.0 + A1
	}

	if outmask&DISTANCE != 0 {
		B1 := Sin_cos_series(true, ssig2, csig2, C1a) - Sin_cos_series(true, ssig1, csig1, C1a)
		s12b = A1 * (sig12 + B1)
		if outmask&(REDUCEDLENGTH|GEODESICSCALE) != 0 {
			B2 := Sin_cos_series(true, ssig2, csig2, C2a) - Sin_cos_series(true, ssig1, csig1, C2a)
			J12 = m0x*sig12 + (A1*B1 - A2*B2)
		}
	} else if outmask&(REDUCEDLENGTH|GEODESICSCALE) != 0 {
		for l := 1; l <= int(GEODESIC_ORDER); l++ {
			C2a[l] = A1*C1a[l] - A2*C2a[l]
		}
		J12 = m0x*sig12 + (Sin_cos_series(true, ssig2, csig2, C2a) - Sin_cos_series(true, ssig1, csig1, C2a))
	}

	if outmask&REDUCEDLENGTH != 0 {
		m0 = m0x
		// J12 is wrong
		m12b = dn2*(csig1*ssig2) - dn1*(ssig1*csig2) - csig1*csig2*J12
	}

	if outmask&GEODESICSCALE != 0 {
		csig12 := csig1*csig2 + ssig1*ssig2
		t := g._ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
		M12 = csig12 + (t*ssig2-csig2*J12)*ssig1/dn1
		M21 = csig12 - (t*ssig1-csig1*J12)*ssig2/dn2
	}
	return s12b, m12b, m0, M12, M21
}

func (g Geodesic) _InverseStart(
	sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12 float64,
	C1a []float64,
	C2a []float64,
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
		sbetm2 := Sq(sbet1 + sbet2)
		sbetm2 /= sbetm2 + Sq(cbet1+cbet2)
		dnm = math.Sqrt(1.0 + g._ep2*sbetm2)
		omg12 := lam12 / (g._f1 * dnm)
		somg12 = math.Sin(omg12)
		comg12 = math.Cos(omg12)
	} else {
		somg12 = slam12
		comg12 = clam12
	}

	salp1 := cbet2 * somg12

	calp1 := sbet12a - cbet2*sbet1*Sq(somg12)/(1.0-comg12)
	if comg12 >= 0.0 {
		calp1 = sbet12 + cbet2*sbet1*Sq(somg12)/(1.0+comg12)
	}

	ssig12 := math.Hypot(salp1, calp1)
	csig12 := sbet1*sbet2 + cbet1*cbet2*comg12

	if shortline && (ssig12 < g._etol2) {
		salp2 = cbet1 * somg12
		var to_mul float64
		if comg12 >= 0.0 {
			to_mul = Sq(somg12) / (1.0 + comg12)
		} else {
			to_mul = 1.0 - comg12
		}
		calp2 = sbet12 - cbet1*sbet2*to_mul

		salp2, calp2 = Norm(salp2, calp2)
		sig12 = math.Atan2(ssig12, csig12)
	} else if math.Abs(g._n) > 0.1 || csig12 >= 0.0 || ssig12 >= 6.0*math.Abs(g._n)*math.Pi*Sq(cbet1) {
	} else {
		var x float64
		var y float64
		var betscale float64
		var lamscale float64
		lam12x := math.Atan2(-slam12, -clam12)
		if g.f >= 0.0 {
			k2 := Sq(sbet1) * g._ep2
			eps := k2 / (2.0*(1.0+math.Sqrt(1.0+k2)) + k2)
			lamscale = g.f * cbet1 * g._A3f(eps) * math.Pi
			betscale = lamscale * cbet1
			x = lam12x / lamscale
			y = sbet12a / betscale
		} else {
			cbet12a := cbet2*cbet1 - sbet2*sbet1
			bet12a := math.Atan2(sbet12, cbet12a)
			_, m12b, m0, _, _ := g._Lengths(
				g._n,
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
				C1a,
				C2a,
			)
			x = -1.0 + m12b/(cbet1*cbet2*m0*math.Pi)
			var betscale float64
			if x < -0.01 {
				betscale = sbet12a / x
			} else {
				betscale = -g.f * Sq(cbet1) * math.Pi
			}
			lamscale = betscale / cbet1
			y = lam12x / lamscale
		}
		if y > -g.tol1_ && x > -1.0-g.xthresh_ {
			if g.f >= 0.0 {
				salp1 = math.Min(-x, 1.0)
				calp1 = -math.Sqrt(1.0 - Sq(salp1))
			} else {
				var to_compare float64
				if x > -g.tol1_ {
					to_compare = 0.0
				} else {
					to_compare = -1.0
				}
				calp1 = math.Max(x, to_compare)
				salp1 = math.Sqrt(1.0 - Sq(calp1))
			}
		} else {
			k := Astroid(x, y)
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
			calp1 = sbet12a - cbet2*sbet1*Sq(somg12)/(1.0-comg12)
		}
	}

	if !(salp1 <= 0.0) {
		salp1, calp1 = Norm(salp1, calp1)
	} else {
		salp1 = 1.0
		calp1 = 0.0
	}
	return sig12, salp1, calp1, salp2, calp2, dnm
}

func (g Geodesic) _Lambda12(
	sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1, slam120, clam120 float64,
	diffp bool,
	C1a []float64,
	C2a []float64,
	C3a []float64,
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
	ssig1, csig1 = Norm(ssig1, csig1)

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
		calp2 = math.Sqrt(Sq(calp1*cbet1)+to_add) / cbet2

	}

	ssig2 := sbet2
	somg2 := salp0 * sbet2
	csig2 := calp2 * cbet2
	comg2 := calp2 * cbet2
	ssig2, csig2 = Norm(ssig2, csig2)

	sig12 := math.Atan2(math.Max(csig1*ssig2-ssig1*csig2, 0.0), csig1*csig2+ssig1*ssig2)
	somg12 := math.Max((comg1*somg2 - somg1*comg2), 0.0)
	comg12 := comg1*comg2 + somg1*somg2
	eta := math.Atan2(somg12*clam120-comg12*slam120, comg12*clam120+somg12*slam120)

	k2 := Sq(calp0) * g._ep2
	eps := k2 / (2.0*(1.0+math.Sqrt(1.0+k2)) + k2)
	g._C3f(eps, C3a)
	B312 := Sin_cos_series(true, ssig2, csig2, C3a) - Sin_cos_series(true, ssig1, csig1, C3a)
	domg12 := -g.f * g._A3f(eps) * salp0 * (sig12 + B312)
	lam12 := eta + domg12

	var dlam12 float64
	if diffp {
		if calp2 == 0.0 {
			dlam12 = -2.0 * g._f1 * dn1 / sbet1
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
				C1a,
				C2a,
			)
			dlam12 = res
			dlam12 *= g._f1 / (calp2 * cbet2)
		}
	} else {
		dlam12 = math.NaN()
	}
	return lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12

}

// _gen_inverse_azi returns (a12, s12, azi1, azi2, m12, M12, M21, S12)
func (g Geodesic) _gen_inverse_azi(
	lat1, lon1, lat2, lon2 float64,
	outmask uint64,
) (float64, float64, float64, float64, float64, float64, float64, float64) {

	azi1 := math.NaN()
	azi2 := math.NaN()
	outmask &= OUT_MASK

	a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12 := g._gen_inverse(
		lat1, lon1, lat2, lon2, outmask,
	)

	if outmask&AZIMUTH != 0 {
		azi1 = Atan2deg(salp1, calp1)
		azi2 = Atan2deg(salp2, calp2)
	}
	return a12, s12, azi1, azi2, m12, M12, M21, S12
}

// _gen_inverse returns (a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12)
func (g Geodesic) _gen_inverse(lat1, lon1, lat2, lon2 float64, outmask uint64) (
	float64, float64, float64, float64, float64, float64, float64, float64, float64, float64,
) {
	a12 := math.NaN()
	s12 := math.NaN()
	m12 := math.NaN()
	M12 := math.NaN()
	M21 := math.NaN()
	S12 := math.NaN()
	outmask &= OUT_MASK

	lon12, lon12s := Ang_diff(lon1, lon2)
	var lonsign float64
	if lon12 >= 0.0 {
		lonsign = 1.0
	} else {
		lonsign = -1.0
	}

	lon12 = lonsign * Ang_round(lon12)
	lon12s = Ang_round((180.0 - lon12) - lonsign*lon12s)
	lam12 := lon12 * DEG2RAD
	var slam12 float64
	var clam12 float64
	if lon12 > 90.0 {
		s, c := Sincosd(lon12s)
		slam12 = s
		clam12 = c
		clam12 = -clam12
	} else {
		s, c := Sincosd(lon12)
		slam12 = s
		clam12 = c
	}
	lat1 = Ang_round(Lat_fix(lat1))
	lat2 = Ang_round(Lat_fix(lat2))

	var swapp float64
	if math.Abs(lat1) < math.Abs(lat2) {
		swapp = -1.0
	} else {
		swapp = 1.0
	}

	var latsign float64
	if lat1 < 0.0 {
		latsign = 1.0
	} else {
		latsign = -1.0
	}
	lat1 *= latsign
	lat2 *= latsign

	sbet1, cbet1 := Sincosd(lat1)
	sbet1 *= g._f1

	sbet1, cbet1 = Norm(sbet1, cbet1)
	cbet1 = math.Max(cbet1, g.tiny_)

	sbet2, cbet2 := Sincosd(lat2)
	sbet2 *= g._f1

	sbet2, cbet2 = Norm(sbet2, cbet2)
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

	dn1 := math.Sqrt(1.0 + g._ep2*Sq(sbet1))
	dn2 := math.Sqrt(1.0 + g._ep2*Sq(sbet2))

	const CARR_SIZE uint64 = uint64(GEODESIC_ORDER) + 1
	C1a := [CARR_SIZE]float64{}
	C2a := [CARR_SIZE]float64{}
	C3a := [GEODESIC_ORDER]float64{}

	meridian := lat1 == -90.0 || slam12 == 0.0
	calp1 := 0.0
	salp1 := 0.0
	calp2 := 0.0
	salp2 := 0.0
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
			g._n,
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
			m12x *= g._b
			s12x *= g._b
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
		sig12 = lam12 / g._f1
		omg12 = lam12 / g._f1
		m12x = g._b * math.Sin(sig12)
		if outmask&GEODESICSCALE != 0 {
			M12 = math.Cos(sig12)
			M21 = math.Cos(sig12)
		}
		a12 = lon12 / g._f1
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
			s12x = sig12 * g._b * dnm
			m12x = Sq(dnm) * g._b * math.Sin(sig12/dnm)
			if outmask&GEODESICSCALE != 0 {
				M12 = math.Cos(sig12 / dnm)
				M21 = math.Cos(sig12 / dnm)
			}
			a12 = sig12 * RAD2DEG
			omg12 = lam12 / (g._f1 * dnm)
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
						salp1, calp1 = Norm(salp1, calp1)
						tripn = math.Abs(v) <= 16.0*g.tol0_
						continue
					}
				}

				salp1 = (salp1a + salp1b) / 2.0
				calp1 = (calp1a + calp1b) / 2.0
				salp1, calp1 = Norm(salp1, calp1)
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

			m12x *= g._b
			s12x *= g._b
			a12 = sig12 * DEG2RAD
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
			k2 := Sq(calp0) * g._ep2
			eps = k2 / (2.0*(1.0+math.Sqrt(1.0+k2)) + k2)
			A4 := Sq(g.a) * calp0 * salp0 * g._e2
			ssig1, csig1 = Norm(ssig1, csig1)
			ssig2, csig2 = Norm(ssig2, csig2)
			C4a := [GEODESIC_ORDER]float64{}
			g._C4f(eps, C4a[:])
			B41 := Sin_cos_series(false, ssig1, csig1, C4a[:])
			B42 := Sin_cos_series(false, ssig2, csig2, C4a[:])
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
		S12 += g._c2 * alp12
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
