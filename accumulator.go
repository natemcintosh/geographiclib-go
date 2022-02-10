package geographiclibgo

import "math"

// Accumulator allows a running sum of float64s
type Accumulator struct {
	s, t float64
}

func NewAccumulatorFromS(s float64) Accumulator {
	p := Accumulator{}
	p.SetS(s)
	return p
}

func NewAccumulatorFromST(s, t float64) Accumulator {
	p := Accumulator{}
	p.SetST(s, t)
	return p
}

// SetS sets s to the argument, and t to 0
func (p *Accumulator) SetS(s float64) {
	p.s = s
	p.t = 0.0
}

// SetST sets both s and t
func (p *Accumulator) SetST(s, t float64) {
	p.s = s
	p.t = t
}

// Add a value to the sum
func (p *Accumulator) Add(y float64) {
	// Here's Shewchuk's solution...
	// hold exact sum as [s, t, u]
	y, u := Sum(y, p.t)    // Accumulate starting at
	p.s, p.t = Sum(y, p.s) // least significant end

	// Start is _s, _t decreasing and non-adjacent.  Sum is now (s + t + u)
	// exactly with s, t, u non-adjacent and in decreasing order (except
	// for possible zeros).  The following code tries to normalize the
	// result.  Ideally, we want _s = round(s+t+u) and _u = round(s+t+u -
	// _s).  The follow does an approximate job (and maintains the
	// decreasing non-adjacent property).  Here are two "failures" using
	// 3-bit floats:
	//
	// Case 1: _s is not equal to round(s+t+u) -- off by 1 ulp
	// [12, -1] - 8 -> [4, 0, -1] -> [4, -1] = 3 should be [3, 0] = 3
	//
	// Case 2: _s+_t is not as close to s+t+u as it shold be
	// [64, 5] + 4 -> [64, 8, 1] -> [64,  8] = 72 (off by 1)
	//                    should be [80, -7] = 73 (exact)
	//
	// "Fixing" these problems is probably not worth the expense.  The
	// representation inevitably leads to small errors in the accumulated
	// values.  The additional errors illustrated here amount to 1 ulp of
	// the less significant word during each addition to the Accumulator
	// and an additional possible error of 1 ulp in the reported sum.
	//
	// Incidentally, the "ideal" representation described above is not
	// canonical, because _s = round(_s + _t) may not be true.  For
	// example, with 3-bit floats:
	//
	// [128, 16] + 1 -> [160, -16] -- 160 = round(145).
	// But [160, 0] - 16 -> [128, 16] -- 128 = round(144).

	if p.s == 0.0 { // This implies t == 0,
		p.s = u // so result is u
	} else {
		p.t += u // otherwise just accumulate u to t.
	}
}

// Sum returns sum + y
func (p *Accumulator) Sum(y float64) float64 {
	if y == 0.0 {
		return p.s
	} else {
		b := NewAccumulatorFromST(p.s, p.t)
		b.Add(y)
		return b.s
	}
}

// Negate negates the sum
func (p *Accumulator) Negate() {
	p.s *= -1
	p.t *= -1
}

// Remainder on division by y
func (p *Accumulator) Remainder(y float64) {
	p.s = math.Remainder(p.s, y)
	p.Add(0.0)
}
