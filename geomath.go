package geographiclibgo

import "math"

const DIGITS uint64 = 53
const TWO float64 = 2.0

const RAD2DEG float64 = 180.0 / math.Pi
const DEG2RAD float64 = math.Pi / 180.0

func Get_epsilon() float64 {
	return math.Pow(TWO, 1.0-float64(DIGITS))
}

func Get_min_val() float64 {
	return math.SmallestNonzeroFloat64
}

// Sq: square a number
func Sq(x float64) float64 {
	return math.Pow(x, 2.0)
}

func Cbrt(x float64) float64 {
	y := math.Pow(math.Abs(x), 1.0/3.0)

	// return y if x > 0 else (-y if x < 0 else x)
	if x > 0 {
		return y
	} else if x < 0 {
		return -y
	} else {
		return x
	}
}

// Norm: normalize a two-vector
func Norm(x, y float64) (float64, float64) {
	r := math.Sqrt(x*x + y*y)
	return (x / r), (y / r)
}

// Sum: error free transformation of a sum
func Sum(u, v float64) (float64, float64) {
	s := u + v
	up := s - v
	vpp := s - up
	up -= u
	vpp -= v
	t := -(up + vpp)
	return s, t
}

// Polyval: evaluate a polynomial
func Polyval(n int64, p []float64, x float64) float64 {
	if n < 0 {
		return 0.0
	} else {
		y := p[0]
		for _, val := range p[1 : n+1] {
			y = y*x + val
		}
		return y
	}
}

// Ang_round: round an angle so that small values underflow to 0
func Ang_round(x float64) float64 {
	// The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
	// for reals = 0.7 pm on the earth if x is an angle in degrees.  (This
	// is about 1000 times more resolution than we get with angles around 90
	// degrees.)  We use this to avoid having to deal with near singular
	// cases when x is non-zero but tiny (e.g., 1.0e-200).
	z := 1.0 / 16.0
	y := math.Abs(x)
	// The compiler mustn't "simplify" z - (z - y) to y
	if y < z {
		y = z - (z - y)
	}
	if x == 0.0 {
		return 0.0
	} else {
		if x < 0.0 {
			return -y
		} else {
			return y
		}
	}
}

// Ang_normalize: reduce angle to (-180,180]
func Ang_normalize(x float64) float64 {
	// y = Math.remainder(x, 360)
	// return 180 if y == -180 else y
	y := math.Remainder(x, 360.0)
	if y == -180.0 {
		return 180.0
	} else {
		return y
	}
}

// Lat_fix: replace angles outside [-90,90] with NaN
func Lat_fix(x float64) float64 {
	if math.Abs(x) > 90.0 {
		return math.NaN()
	} else {
		return x
	}
}

// Ang_diff: compute y - x and reduce to [-180,180] accurately
func Ang_diff(x, y float64) (float64, float64) {
	d, t := Sum(Ang_normalize(-x), Ang_normalize(y))
	d = Ang_normalize(d)
	if d == 180.0 && t > 0.0 {
		return Sum(-180.0, t)
	} else {
		return Sum(d, t)
	}
}

// Sincosd: compute sine and cosine of x in degrees
func Sincosd(x float64) (float64, float64) {
	// r = math.fmod(x, 360) if Math.isfinite(x) else Math.nan
	r := math.NaN()
	if !math.IsInf(x, 0) {
		r = math.Mod(x, 360.0)
	}

	// q = 0 if Math.isnan(r) else int(round(r / 90))
	q := math.Round(r / 90.0)
	if math.IsNaN(r) {
		q = 0
	}

	// r -= 90 * q; r = math.radians(r)
	r -= 90.0 * q
	r *= DEG2RAD

	// s = math.sin(r); c = math.cos(r)
	s := math.Sin(r)
	c := math.Cos(r)

	// q = q % 4
	q = math.Mod(q, 4.0)

	// if q == 1:
	//     s, c =  c, -s
	// elif q == 2:
	//     s, c = -s, -c
	// elif q == 3:
	//     s, c = -c,  s

	if q < 0 {
		q += 4
	}

	if q == 1 {
		s, c = c, -s
	} else if q == 2 {
		s, c = -s, -c
	} else if q == 3 {
		s, c = -c, s
	}

	// # Remove the minus sign on -0.0 except for sin(-0.0).
	// # On Windows 32-bit with python 2.7, math.fmod(-0.0, 360) = +0.0
	// # (x, c) here fixes this bug.  See also Math::sincosd in the C++ library.
	// # AngNormalize has a similar fix.
	//     s, c = (x, c) if x == 0 else (0.0+s, 0.0+c)
	// return s, c
	if x == 0.0 {
		s = x
	} else {
		s, c = 0.0+s, 0.0+c
	}

	return s, c
}

// Atan2_deg: compute the arc tangent of y/x in degrees
func Atan2_deg(y_deg, x_deg float64) float64 {
	// First convert to radians.
	y_rad := y_deg * DEG2RAD
	x_rad := x_deg * DEG2RAD

	// Then use the atan2 function.
	return math.Atan2(y_rad, x_rad) * RAD2DEG
}

func Eatanhe(x float64, es float64) float64 {
	if es > 0.0 {
		return es * math.Atanh(es*x)
	} else {
		return -es * math.Atan(es*x)
	}
}

// Sin_cos_series: functions that used to be inside Geodesic
func Sin_cos_series(sinp bool, sinx float64, cosx float64, c []float64) float64 {
	k := len(c)

	to_sub := 0
	if sinp {
		to_sub = 1
	}
	n := k - to_sub

	var ar float64 = 2.0 * (cosx - sinx) * (cosx + sinx)
	y1 := 0.0
	y0 := 0.0
	if n&1 != 0 {
		k -= 1
		y0 = c[k]
	}

	n /= 2

	for n > 0 {
		n -= 1
		k -= 1
		y1 = ar*y0 - y1 + c[k]
		k -= 1
		y0 = ar*y1 - y0 + c[k]
	}
	if sinp {
		return 2.0 * sinx * cosx * y0
	} else {
		return cosx * (y0 - y1)
	}
}

// Astroid: solve astroid equation
func Astroid(x, y float64) float64 {
	p := Sq(x)
	q := Sq(y)
	r := (p + q - 1.0) / 6.0
	if !(q == 0.0 && r <= 0.0) {
		s := p * q / 4.0
		r2 := Sq(r)
		r3 := r * r2
		disc := s * (s + 2.0*r3)
		u := r
		if disc >= 0.0 {
			t3 := s + r3

			if t3 < 0.0 {
				t3 += -math.Sqrt(disc)
			} else {
				t3 += math.Sqrt(disc)
			}

			t := Cbrt(t3) // we could use built-in math.Cbrt

			to_add := 0.0
			if t != 0.0 {
				to_add = r2 / t
			}
			u += t + to_add
		} else {
			ang := math.Atan2(math.Sqrt(-disc), -(s + r3))
			u += 2.0 * r * math.Cos(ang/3.0)
		}
		v := math.Sqrt((Sq(u) + q))
		uv := u + v
		if u < 0.0 {
			uv = q / (v + u)
		}
		w := (uv - q) / (2.0 * v)
		return uv / (math.Sqrt(uv+Sq(w)) + w)
	} else {
		return 0.0
	}
}

func _A1m1f(eps float64, geodesic_order int64) float64 {
	COEFF := [5]float64{1.0, 4.0, 64.0, 0.0, 256.0}
	m := geodesic_order / 2
	t := Polyval(m, COEFF[:], Sq(eps)) / COEFF[(m+1)]
	return (t + eps) / (1.0 - eps)
}

func _C1f(eps float64, c []float64, geodesic_order int) {

	COEFF := [18]float64{
		-1.0, 6.0, -16.0, 32.0, -9.0, 64.0, -128.0, 2048.0, 9.0, -16.0, 768.0, 3.0, -5.0, 512.0,
		-7.0, 1280.0, -7.0, 2048.0,
	}
	eps2 := Sq(eps)
	d := eps
	var o int64 = 0

	for l := 1; l <= geodesic_order; l++ {
		m := int64((geodesic_order - l) / 2)
		c[l] = d * Polyval(m, COEFF[o:], eps2) / COEFF[(o+m+1)]
		o += m + 2
		d *= eps
	}
}

func _C1pf(eps float64, c []float64, geodesic_order int) {

	COEFF := [18]float64{
		205.0, -432.0, 768.0, 1536.0, 4005.0, -4736.0, 3840.0, 12288.0, -225.0, 116.0, 384.0,
		-7173.0, 2695.0, 7680.0, 3467.0, 7680.0, 38081.0, 61440.0,
	}

	eps2 := Sq(eps)
	d := eps
	var o int64 = 0

	for l := 1; l <= geodesic_order; l++ {
		m := int64((geodesic_order - l) / 2)
		c[l] = d * Polyval(m, COEFF[o:], eps2) / COEFF[(o+m+1)]
		o += m + 2
		d *= eps
	}
}

func _A2m1f(eps float64, geodesic_order int64) float64 {

	COEFF := []float64{-11.0, -28.0, -192.0, 0.0, 256.0}
	var m int64 = geodesic_order / 2
	t := Polyval(m, COEFF, Sq(eps)) / COEFF[(m+1)]
	return (t - eps) / (1.0 + eps)
}

func _C2f(eps float64, c []float64, geodesic_order int) {

	COEFF := [18]float64{
		1.0, 2.0, 16.0, 32.0, 35.0, 64.0, 384.0, 2048.0, 15.0, 80.0, 768.0, 7.0, 35.0, 512.0, 63.0,
		1280.0, 77.0, 2048.0,
	}
	eps2 := Sq(eps)
	d := eps
	var o int64 = 0

	for l := 1; l <= geodesic_order; l++ {
		m := int64((geodesic_order - l) / 2)
		c[l] = d * Polyval(m, COEFF[o:], eps2) / COEFF[(o+m+1)]
		o += m + 2
		d *= eps
	}
}
