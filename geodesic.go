package geographiclibgo

import "math"

const WGS84_A float64 = 6378137.0

// Evaluating this as 1000000000.0 / (298257223563f64) reduces the
// round-off error by about 10%.  However, expressing the flattening as
// 1/298.257223563 is well ingrained.
const WGS84_F float64 = 1.0 / ((298257223563.0) / 1000000000.0)

const GEODESIC_ORDER int64 = 6
const nC3x_ int64 = 15
const nC4x_ int64 = 21

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
