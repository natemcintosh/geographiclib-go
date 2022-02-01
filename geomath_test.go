package geographiclibgo

import (
	"math"
	"testing"
)

const float64EqualityThreshold = 4e-6

func almost_equal(a, b, threshold float64) bool {
	return math.Abs(a-b) < threshold
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
			if !almost_equal(sin, tC.out1, float64EqualityThreshold) || !almost_equal(cos, tC.out2, float64EqualityThreshold) {
				t.Errorf("Sincosd(%v) = %v, %v; want %v, %v", tC.in, sin, cos, tC.out1, tC.out2)
			}
		})
	}
}

func BenchmarkSincosd(b *testing.B) {
	benchmarks := []struct {
		desc string
		in   float64
	}{
		{
			desc: "a negative number",
			in:   -77.03196,
		},
		{
			desc: "a positive number",
			in:   69.48894,
		},
		{
			desc: "-1",
			in:   -1.0,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				Sincosd(bm.in)
			}

		})
	}
}

func Test_A1m1f(t *testing.T) {
	exptected := 0.1404582405272727
	got := _A1m1f(0.12, 6)
	if !almost_equal(exptected, got, float64EqualityThreshold) {
		t.Errorf("_A1m1f(0.12, 6) = %v; want %v", got, exptected)
	}
}

func Benchmark_A1m1f(b *testing.B) {
	benchmarks := []struct {
		desc string
		in1  float64
		in2  int64
	}{
		{
			desc: "1",
			in1:  0.12,
			in2:  6,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_A1m1f(bm.in1, bm.in2)
			}

		})
	}
}

func TestAstroid(t *testing.T) {
	expected := 23.44475767500982
	got := Astroid(21.0, 12.0)
	if !almost_equal(expected, got, float64EqualityThreshold) {
		t.Errorf("Astroid(21.0, 12.0) = %v; want %v", got, expected)
	}
}

func BenchmarkAstroid(b *testing.B) {
	benchmarks := []struct {
		desc string
		in1  float64
		in2  float64
	}{
		{
			desc: "1",
			in1:  21.0,
			in2:  12.0,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				Astroid(bm.in1, bm.in2)
			}

		})
	}
}

func TestSin_cos_series(t *testing.T) {
	testCases := []struct {
		desc string
		sinp bool
		sinx float64
		cosx float64
		c    []float64
		want float64
	}{
		{
			desc: "first",
			sinp: false,
			sinx: -0.8928657853278468,
			cosx: 0.45032287238256896,
			c: []float64{
				0.6660771734724675,
				1.5757752625233906e-05,
				3.8461688963148916e-09,
				1.3040960748120204e-12,
				5.252912023008548e-16,
				2.367770858285795e-19,
			},
			want: 0.29993425660538664,
		},
		{
			desc: "second",
			sinp: false,
			sinx: -0.8928657853278468,
			cosx: 0.45032287238256896,
			c: []float64{
				0.0, 1.0, 2.0, 3.0, 4.0, 5.0,
			},
			want: 1.8998562852254026,
		},
		{
			desc: "third",
			sinp: true,
			sinx: -0.8928657853278468,
			cosx: 0.45032287238256896,
			c: []float64{0.0,
				-0.0003561309485314716,
				-3.170731714689771e-08,
				-7.527972480734327e-12,
				-2.5133854116682488e-15,
				-1.0025061462383107e-18,
				-4.462794158625518e-22,
			},
			want: 0.00028635444718997857,
		},
		{
			desc: "fourth",
			sinp: true,
			sinx: 0.12,
			cosx: 0.21,
			c:    []float64{1.0, 2.0},
			want: 0.10079999999999999,
		},
		{
			desc: "fifth",
			sinp: true,
			sinx: -0.024679833885152578,
			cosx: 0.9996954065111039,
			c: []float64{
				0.0,
				-0.0008355098973052918,
				-1.7444619952659748e-07,
				-7.286557795511902e-11,
				-3.80472772706481e-14,
				-2.2251271876594078e-17,
				1.2789961247944744e-20,
			},
			want: 4.124513511893872e-05,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := Sin_cos_series(tC.sinp, tC.sinx, tC.cosx, tC.c)
			if !almost_equal(tC.want, got, float64EqualityThreshold) {
				t.Errorf("Sin_cos_series(%v, %v, %v, %v) = %v; want %v", tC.sinp, tC.sinx, tC.cosx, tC.c, got, tC.want)
			}
		})
	}
}

func BenchmarkSin_cos_series(b *testing.B) {
	benchmarks := []struct {
		desc string
		sinp bool
		sinx float64
		cosx float64
		c    []float64
	}{
		{
			desc: "first",
			sinp: false,
			sinx: -0.8928657853278468,
			cosx: 0.45032287238256896,
			c: []float64{
				0.6660771734724675,
				1.5757752625233906e-05,
				3.8461688963148916e-09,
				1.3040960748120204e-12,
				5.252912023008548e-16,
				2.367770858285795e-19,
			},
		},
		{
			desc: "second",
			sinp: false,
			sinx: -0.8928657853278468,
			cosx: 0.45032287238256896,
			c: []float64{
				0.0, 1.0, 2.0, 3.0, 4.0, 5.0,
			},
		},
		{
			desc: "third",
			sinp: true,
			sinx: -0.8928657853278468,
			cosx: 0.45032287238256896,
			c: []float64{0.0,
				-0.0003561309485314716,
				-3.170731714689771e-08,
				-7.527972480734327e-12,
				-2.5133854116682488e-15,
				-1.0025061462383107e-18,
				-4.462794158625518e-22,
			},
		},
		{
			desc: "fourth",
			sinp: true,
			sinx: 0.12,
			cosx: 0.21,
			c:    []float64{1.0, 2.0},
		},
		{
			desc: "fifth",
			sinp: true,
			sinx: -0.024679833885152578,
			cosx: 0.9996954065111039,
			c: []float64{
				0.0,
				-0.0008355098973052918,
				-1.7444619952659748e-07,
				-7.286557795511902e-11,
				-3.80472772706481e-14,
				-2.2251271876594078e-17,
				1.2789961247944744e-20,
			},
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				Sin_cos_series(bm.sinp, bm.sinx, bm.cosx, bm.c)
			}

		})
	}
}

func Test_C1f(t *testing.T) {
	c := []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
	_C1f(0.12, c, 6)
	want_c := []float64{
		1.0,
		-0.059676777599999994,
		-0.000893533122,
		-3.57084e-05,
		-2.007504e-06,
		-1.3607999999999999e-07,
		-1.0205999999999999e-08,
	}
	// Compare each number in c and want_c
	for i := 0; i < len(c); i++ {
		if !almost_equal(c[i], want_c[i], float64EqualityThreshold) {
			t.Errorf("c[%v] = %v; want %v", i, c[i], want_c[i])
		}
	}
}

func Benchmark_C1f(b *testing.B) {
	benchmarks := []struct {
		desc           string
		eps            float64
		c              []float64
		geodesic_order int
	}{
		{
			desc:           "1",
			eps:            0.12,
			c:              []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
			geodesic_order: 6,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_C1f(bm.eps, bm.c, bm.geodesic_order)
			}

		})
	}
}

func Test_C1pf(t *testing.T) {
	c := []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
	_C1pf(0.12, c, 6)
	want_c := []float64{
		1.0,
		0.059517321000000005,
		0.004421053215,
		0.0005074200000000001,
		6.997613759999999e-05,
		1.1233080000000001e-05,
		1.8507366e-06,
	}
	// Compare each number in c and want_c
	for i := 0; i < len(c); i++ {
		if !almost_equal(c[i], want_c[i], float64EqualityThreshold) {
			t.Errorf("c[%v] = %v; want %v", i, c[i], want_c[i])
		}
	}
}

func Benchmark_C1pf(b *testing.B) {
	benchmarks := []struct {
		desc           string
		eps            float64
		c              []float64
		geodesic_order int
	}{
		{
			desc:           "1",
			eps:            0.12,
			c:              []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
			geodesic_order: 6,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_C1pf(bm.eps, bm.c, bm.geodesic_order)
			}

		})
	}
}

func Test_C2f(t *testing.T) {
	c := []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
	_C2f(0.12, c, 6)
	want_c := []float64{
		1.0,
		0.0601087776,
		0.00270653103,
		0.000180486,
		1.4215824e-05,
		1.22472e-06,
		1.12266e-07,
	}
	// Compare each number in c and want_c
	for i := 0; i < len(c); i++ {
		if !almost_equal(c[i], want_c[i], float64EqualityThreshold) {
			t.Errorf("c[%v] = %v; want %v", i, c[i], want_c[i])
		}
	}
}

func Benchmark_C2f(b *testing.B) {
	benchmarks := []struct {
		desc           string
		eps            float64
		c              []float64
		geodesic_order int
	}{
		{
			desc:           "1",
			eps:            0.12,
			c:              []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
			geodesic_order: 6,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_C2f(bm.eps, bm.c, bm.geodesic_order)
			}

		})
	}
}

func Test_A2m1f(t *testing.T) {
	got := _A2m1f(0.12, 6)
	want := -0.11680607884285714
	if !almost_equal(got, want, float64EqualityThreshold) {
		t.Errorf("_A2m1f(%v) = %v; want %v", 0.12, got, want)
	}
}

func Benchmark_A2m1f(b *testing.B) {
	benchmarks := []struct {
		desc           string
		eps            float64
		geodesic_order int64
	}{
		{
			desc:           "1",
			eps:            0.12,
			geodesic_order: 6,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_A2m1f(bm.eps, bm.geodesic_order)
			}

		})
	}
}

func TestAng_diff(t *testing.T) {
	testCases := []struct {
		desc string
		x    float64
		y    float64
		d    float64
		t    float64
	}{
		{
			desc: "1",
			x:    0.0,
			y:    1.0,
			d:    1.0,
			t:    0.0,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			d, tvar := Ang_diff(tC.x, tC.y)

			if !f64_equals(d, tC.d) {
				t.Errorf("d = %v; want %v", d, tC.d)
			}

			if !f64_equals(tvar, tC.t) {
				t.Errorf("t = %v; want %v", tvar, tC.t)
			}
		})
	}
}

func TestAng_normalize(t *testing.T) {
	testCases := []struct {
		desc string
		in   float64
		want float64
	}{
		{
			desc: "0",
			in:   0,
			want: 0,
		},
		{
			desc: "1",
			in:   1,
			want: 1,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := Ang_normalize(tC.in)
			if !f64_equals(got, tC.want) {
				t.Errorf("Ang_normalize(%v) = %v; want %v", tC.in, got, tC.want)
			}
		})
	}
}

func TestAtan2_deg(t *testing.T) {
	testCases := []struct {
		desc       string
		x, y, want float64
	}{
		{
			desc: "1",
			x:    0.82632552653105184,
			y:    0.56319279487860996,
			want: 55.723110355324408,
		},
		{
			desc: "2",
			x:    0.82632951065581195,
			y:    0.56318694926225565,
			want: 55.72351567783663,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := Atan2_deg(tC.x, tC.y)

			if !f64_equals(tC.want, got) {
				t.Errorf("Atan2_deg(%v, %v) = %v; want %v", tC.x, tC.y, got, tC.want)
			}
		})
	}
}
