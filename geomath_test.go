package geographiclibgo

import (
	"math"
	"testing"
)

const float64EqualityThreshold = 1e-9

func almost_equal(a, b float64) bool {
	return math.Abs(a-b) < float64EqualityThreshold
}

func TestSincosd(t *testing.T) {
	testCases := []struct {
		desc string
		in   float64
		out1 float64
		out2 float64
	}{
		{
			desc: "a negative number",
			in:   -77.03196,
			out1: -0.9744953925159129,
			out2: 0.22440750870961693,
		},
		{
			desc: "a positive number",
			in:   69.48894,
			out1: 0.9366045700708676,
			out2: 0.3503881837653281,
		},
		{
			desc: "-1",
			in:   -1.0,
			out1: -0.01745240643728351,
			out2: 0.9998476951563913,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			sin, cos := Sincosd(tC.in)
			if !almost_equal(sin, tC.out1) || !almost_equal(cos, tC.out2) {
				t.Errorf("Sincosd(%v) = %v, %v; want %v, %v", tC.in, sin, cos, tC.out1, tC.out2)
			}
		})
	}
}
