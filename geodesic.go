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

}
