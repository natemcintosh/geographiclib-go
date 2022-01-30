package geographiclibgo

import (
	"math"
	"testing"
)

// f64_equals tests if two float64 values are equal to within a small epsilon.
// If `want` is NaN, then `got` must also be NaN.
// Otherwise just a regular float comparison is performed.
func f64_equals(want, got float64) bool {
	// If want is NaN, then got should be NaN
	if math.IsNaN(want) {
		if !math.IsNaN(got) {
			return false
		}
	} else {
		if !almost_equal(got, want) {
			return false
		}
	}
	return true
}

func TestGeodesicInit(t *testing.T) {
	geod := Wgs84()
	if !almost_equal(geod.a, 6378137.0) {
		t.Errorf("Wgs84() failed: a = %v, expected %v", geod.a, 6378137.0)
	}
	if !almost_equal(geod.f, 0.0033528106647474805) {
		t.Errorf("Wgs84() failed: f = %v, exptected %v", geod.f, 0.0033528106647474805)
	}

	if !almost_equal(geod._f1, 0.9966471893352525) {
		t.Errorf("Wgs84() failed: _f1 = %v, exptected %v", geod._f1, 0.9966471893352525)
	}
	if !almost_equal(geod._e2, 0.0066943799901413165) {
		t.Errorf("Wgs84() failed: _e2 = %v, exptected %v", geod._e2, 0.0066943799901413165)
	}

	if !almost_equal(geod._ep2, 0.006739496742276434) {
		t.Errorf("Wgs84() failed: _ep2 = %v, exptected %v", geod._ep2, 0.006739496742276434)
	}

	if !almost_equal(geod._n, 0.0016792203863837047) {
		t.Errorf("Wgs84() failed: _n = %v, exptected %v", geod._n, 0.0016792203863837047)
	}

	if !almost_equal(geod._b, 6356752.314245179) {
		t.Errorf("Wgs84() failed: _b = %v, exptected %v", geod._b, 6356752.314245179)
	}

	if !almost_equal(geod._c2, 40589732499314.76) {
		t.Errorf("Wgs84() failed: _c2 = %v, exptected %v", geod._c2, 40589732499314.76)
	}

	if !almost_equal(geod._etol2, 3.6424611488788524e-08) {
		t.Errorf("Wgs84() failed: _etol2 = %v, exptected %v", geod._etol2, 3.6424611488788524e-08)
	}

	// Compare all the elements of geod._A3x
	want_A3x := []float64{-0.0234375,
		-0.046927475637074494,
		-0.06281503005876607,
		-0.2502088451303832,
		-0.49916038980680816,
		1.0,
	}
	for i := 0; i < len(geod._A3x); i++ {
		if !almost_equal(geod._A3x[i], want_A3x[i]) {
			t.Errorf("Wgs84() failed: _A3[%v] = %v, exptected %v", i, geod._A3x[i], want_A3x[i])
		}
	}

	want_C3x := []float64{0.0234375,
		0.03908873781853724,
		0.04695366939653196,
		0.12499964752736174,
		0.24958019490340408,
		0.01953125,
		0.02345061890926862,
		0.046822392185686165,
		0.062342661206936094,
		0.013671875,
		0.023393770302437927,
		0.025963026642854565,
		0.013671875,
		0.01362595881755982,
		0.008203125,
	}
	for i := 0; i < len(geod._C3x); i++ {
		if !almost_equal(geod._C3x[i], want_C3x[i]) {
			t.Errorf("Wgs84() failed: _C3x[%v] = %v, exptected %v", i, geod._C3x[i], want_C3x[i])
		}
	}

	want_C4x := []float64{
		0.00646020646020646,
		0.0035037627212872787,
		0.034742279454780166,
		-0.01921732223244865,
		-0.19923321555984239,
		0.6662190894642603,
		0.000111000111000111,
		0.003426620602971002,
		-0.009510765372597735,
		-0.01893413691235592,
		0.0221370239510936,
		0.0007459207459207459,
		-0.004142006291321442,
		-0.00504225176309005,
		0.007584982177746079,
		-0.0021565735851450138,
		-0.001962613370670692,
		0.0036104265913438913,
		-0.0009472009472009472,
		0.0020416649913317735,
		0.0012916376552740189,
	}
	for i := 0; i < len(geod._C4x); i++ {
		if !almost_equal(geod._C4x[i], want_C4x[i]) {
			t.Errorf("Wgs84() failed: _C4x[%v] = %v, exptected %v", i, geod._C4x[i], want_C4x[i])
		}
	}
}

func BenchmarkWgs84(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Wgs84()
	}
}

func Test_A3f(t *testing.T) {
	geod := Wgs84()
	if !almost_equal(geod._A3f(0.12), 0.9363788874000158) {
		t.Errorf("A3f() failed: %v, exptected %v", geod._A3f(0.12), 0.9363788874000158)
	}
}

func Benchmark_A3f(b *testing.B) {
	geod := Wgs84()
	for i := 0; i < b.N; i++ {
		geod._A3f(0.12)
	}
}

func Test_C3f(t *testing.T) {
	geod := Wgs84()
	c := []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
	geod._C3f(0.12, c)
	want := []float64{
		1.0,
		0.031839442894193756,
		0.0009839921354137713,
		5.0055242248766214e-05,
		3.1656788204092044e-06,
		2.0412e-07,
		7.0,
	}
	// Compare each element of `c` with each corresponding element of `want`
	for i := 0; i < len(c); i++ {
		if !almost_equal(c[i], want[i]) {
			t.Errorf("C3f() at index %d failed: %v, exptected %v", i, c[i], want[i])
		}
	}
}

func Benchmark_C3f(b *testing.B) {
	geod := Wgs84()
	c := []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
	for i := 0; i < b.N; i++ {
		geod._C3f(0.12, c)
	}
}

func Test_C4f(t *testing.T) {
	geod := Wgs84()
	c := []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
	geod._C4f(0.12, c)
	want := []float64{
		0.6420952961066771,
		0.0023680700061156517,
		9.96704067834604e-05,
		5.778187189466089e-06,
		3.9979026199316593e-07,
		3.2140078103714466e-08,
		7.0,
	}

	// Compare each element of `c` with each corresponding element of `want`
	for i := 0; i < len(c); i++ {
		if !almost_equal(c[i], want[i]) {
			t.Errorf("C3f() at index %d failed: %v, exptected %v", i, c[i], want[i])
		}
	}
}

func Benchmark_C4f(b *testing.B) {
	geod := Wgs84()
	c := []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
	for i := 0; i < b.N; i++ {
		geod._C4f(0.12, c)
	}
}

func Test_Lengths(t *testing.T) {
	geod := Wgs84()

	testCases := []struct {
		desc    string
		eps     float64
		sig12   float64
		ssig1   float64
		csig1   float64
		dn1     float64
		ssig2   float64
		csig2   float64
		dn2     float64
		cbet1   float64
		cbet2   float64
		outmask uint64
		c1a     []float64
		c2a     []float64
		want1   float64
		want2   float64
		want3   float64
		want4   float64
		want5   float64
	}{
		{
			desc:    "1",
			eps:     0.0008355095326524276,
			sig12:   0.024682339962725352,
			ssig1:   -0.024679833885152578,
			csig1:   0.9996954065111039,
			dn1:     1.0000010195104125,
			ssig2:   0.0,
			csig2:   1.0,
			dn2:     1.0,
			cbet1:   0.9998487145115275,
			cbet2:   1.0,
			outmask: 4101,
			c1a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			c2a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			want1:   math.NaN(),
			want2:   0.024679842274314294,
			want3:   0.0016717180169067588,
			want4:   math.NaN(),
			want5:   math.NaN(),
		},
		{
			desc:    "2",
			eps:     0.0008355096040059597,
			sig12:   0.024682338906797385,
			ssig1:   -0.02467983282954624,
			csig1:   0.9996954065371639,
			dn1:     1.0000010195104125,
			ssig2:   0.0,
			csig2:   1.0,
			dn2:     1.0,
			cbet1:   0.9998487145115275,
			cbet2:   1.0,
			outmask: 4101,
			c1a: []float64{
				0.0,
				-0.00041775465696698233,
				-4.362974596862037e-08,
				-1.2151022357848552e-11,
				-4.7588881620421004e-15,
				-2.226614930167366e-18,
				-1.1627237498131586e-21,
			},
			c2a: []float64{
				0.0,
				-0.0008355098973052918,
				-1.7444619952659748e-07,
				-7.286557795511902e-11,
				-3.80472772706481e-14,
				-2.2251271876594078e-17,
				1.2789961247944744e-20,
			},
			want1: math.NaN(),
			want2: 0.02467984121870759,
			want3: 0.0016717181597332804,
			want4: math.NaN(),
			want5: math.NaN(),
		},
		{
			desc:    "3",
			eps:     0.0008355096040059597,
			sig12:   0.024682338906797385,
			ssig1:   -0.02467983282954624,
			csig1:   0.9996954065371639,
			dn1:     1.0000010195104125,
			ssig2:   0.0,
			csig2:   1.0,
			dn2:     1.0,
			cbet1:   0.9998487145115275,
			cbet2:   1.0,
			outmask: 1920,
			c1a: []float64{
				0.0,
				-0.00041775469264372037,
				-4.362975342068502e-08,
				-1.215102547098435e-11,
				-4.758889787701359e-15,
				-2.2266158809456692e-18,
				-1.1627243456014359e-21,
			},
			c2a: []float64{
				0.0,
				-0.0008355099686589174,
				-1.744462293162189e-07,
				-7.286559662008413e-11,
				-3.804729026574989e-14,
				-2.2251281376754273e-17,
				1.2789967801615795e-20,
			},
			want1: 0.024682347295447677,
			want2: math.NaN(),
			want3: math.NaN(),
			want4: math.NaN(),
			want5: math.NaN(),
		},
		{
			desc:    "4",
			eps:     0.0007122620325664751,
			sig12:   1.405117407023628,
			ssig1:   -0.8928657853278468,
			csig1:   0.45032287238256896,
			dn1:     1.0011366173804046,
			ssig2:   0.2969032234925426,
			csig2:   0.9549075745221299,
			dn2:     1.0001257451360057,
			cbet1:   0.8139459053827204,
			cbet2:   0.9811634781422108,
			outmask: 1920,
			c1a: []float64{
				0.0,
				-0.0003561309485314716,
				-3.170731714689771e-08,
				-7.527972480734327e-12,
				-2.5133854116682488e-15,
				-1.0025061462383107e-18,
				-4.462794158625518e-22,
			},
			c2a: []float64{
				0.0,
				-0.0007122622584701569,
				-1.2678416507678478e-07,
				-4.514641118748122e-11,
				-2.0096353119518367e-14,
				-1.0019350865558619e-17,
				4.90907357448807e-21,
			},
			want1: 1.4056304412645388,
			want2: math.NaN(),
			want3: math.NaN(),
			want4: math.NaN(),
			want5: math.NaN(),
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			c1a := make([]float64, len(tC.c1a))
			copy(c1a, tC.c1a)
			c2a := make([]float64, len(tC.c2a))
			copy(c2a, tC.c2a)
			got1, got2, got3, got4, got5 := geod._Lengths(
				tC.eps, tC.sig12, tC.ssig1, tC.csig1, tC.dn1, tC.ssig2, tC.csig2, tC.dn2,
				tC.cbet1, tC.cbet2, tC.outmask, c1a, c2a)

			// Compare all return values
			if !f64_equals(tC.want1, got1) {
				t.Errorf("Lengths() got1 = %v, want %v", got1, tC.want1)
			}

			if !f64_equals(tC.want2, got2) {
				t.Errorf("Lengths() got2 = %v, want %v", got2, tC.want2)
			}

			if !f64_equals(tC.want3, got3) {
				t.Errorf("Lengths() got3 = %v, want %v", got3, tC.want3)
			}

			if !f64_equals(tC.want4, got4) {
				t.Errorf("Lengths() got4 = %v, want %v", got4, tC.want4)
			}

			if !f64_equals(tC.want5, got5) {
				t.Errorf("Lengths() got5 = %v, want %v", got5, tC.want5)
			}

		})
	}
}

func Benchmark_Lengths(b *testing.B) {
	geod := Wgs84()

	benchmarks := []struct {
		desc    string
		eps     float64
		sig12   float64
		ssig1   float64
		csig1   float64
		dn1     float64
		ssig2   float64
		csig2   float64
		dn2     float64
		cbet1   float64
		cbet2   float64
		outmask uint64
		c1a     []float64
		c2a     []float64
	}{
		{
			desc:    "1",
			eps:     0.0008355095326524276,
			sig12:   0.024682339962725352,
			ssig1:   -0.024679833885152578,
			csig1:   0.9996954065111039,
			dn1:     1.0000010195104125,
			ssig2:   0.0,
			csig2:   1.0,
			dn2:     1.0,
			cbet1:   0.9998487145115275,
			cbet2:   1.0,
			outmask: 4101,
			c1a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			c2a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
		},
		{
			desc:    "2",
			eps:     0.0008355096040059597,
			sig12:   0.024682338906797385,
			ssig1:   -0.02467983282954624,
			csig1:   0.9996954065371639,
			dn1:     1.0000010195104125,
			ssig2:   0.0,
			csig2:   1.0,
			dn2:     1.0,
			cbet1:   0.9998487145115275,
			cbet2:   1.0,
			outmask: 4101,
			c1a: []float64{
				0.0,
				-0.00041775465696698233,
				-4.362974596862037e-08,
				-1.2151022357848552e-11,
				-4.7588881620421004e-15,
				-2.226614930167366e-18,
				-1.1627237498131586e-21,
			},
			c2a: []float64{
				0.0,
				-0.0008355098973052918,
				-1.7444619952659748e-07,
				-7.286557795511902e-11,
				-3.80472772706481e-14,
				-2.2251271876594078e-17,
				1.2789961247944744e-20,
			},
		},
		{
			desc:    "3",
			eps:     0.0008355096040059597,
			sig12:   0.024682338906797385,
			ssig1:   -0.02467983282954624,
			csig1:   0.9996954065371639,
			dn1:     1.0000010195104125,
			ssig2:   0.0,
			csig2:   1.0,
			dn2:     1.0,
			cbet1:   0.9998487145115275,
			cbet2:   1.0,
			outmask: 1920,
			c1a: []float64{
				0.0,
				-0.00041775469264372037,
				-4.362975342068502e-08,
				-1.215102547098435e-11,
				-4.758889787701359e-15,
				-2.2266158809456692e-18,
				-1.1627243456014359e-21,
			},
			c2a: []float64{
				0.0,
				-0.0008355099686589174,
				-1.744462293162189e-07,
				-7.286559662008413e-11,
				-3.804729026574989e-14,
				-2.2251281376754273e-17,
				1.2789967801615795e-20,
			},
		},
		{
			desc:    "4",
			eps:     0.0007122620325664751,
			sig12:   1.405117407023628,
			ssig1:   -0.8928657853278468,
			csig1:   0.45032287238256896,
			dn1:     1.0011366173804046,
			ssig2:   0.2969032234925426,
			csig2:   0.9549075745221299,
			dn2:     1.0001257451360057,
			cbet1:   0.8139459053827204,
			cbet2:   0.9811634781422108,
			outmask: 1920,
			c1a: []float64{
				0.0,
				-0.0003561309485314716,
				-3.170731714689771e-08,
				-7.527972480734327e-12,
				-2.5133854116682488e-15,
				-1.0025061462383107e-18,
				-4.462794158625518e-22,
			},
			c2a: []float64{
				0.0,
				-0.0007122622584701569,
				-1.2678416507678478e-07,
				-4.514641118748122e-11,
				-2.0096353119518367e-14,
				-1.0019350865558619e-17,
				4.90907357448807e-21,
			},
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				geod._Lengths(
					bm.eps, bm.sig12, bm.ssig1, bm.csig1, bm.dn1, bm.ssig2, bm.csig2, bm.dn2,
					bm.cbet1, bm.cbet2, bm.outmask, bm.c1a, bm.c2a)
			}

		})
	}
}

func Test_InverseStart(t *testing.T) {
	geod := Wgs84()

	testCases := []struct {
		desc   string
		sbet1  float64
		cbet1  float64
		dn1    float64
		sbet2  float64
		cbet2  float64
		dn2    float64
		lam12  float64
		slam12 float64
		clam12 float64
		C1a    []float64
		C2a    []float64
		want1  float64
		want2  float64
		want3  float64
		want4  float64
		want5  float64
		want6  float64
	}{
		{
			desc:   "1",
			sbet1:  -0.017393909556108908,
			cbet1:  0.9998487145115275,
			dn1:    1.0000010195104125,
			sbet2:  0.0,
			cbet2:  1.0,
			dn2:    1.0,
			lam12:  0.017453292519943295,
			slam12: 0.01745240643728351,
			clam12: 0.9998476951563913,
			C1a:    []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			C2a:    []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			want1:  -1.0,
			want2:  0.7095310092765433,
			want3:  0.7046742132893822,
			want4:  math.NaN(),
			want5:  math.NaN(),
			want6:  1.0000002548969817,
		},
		{
			desc:   "2",
			sbet1:  -0.017393909556108908,
			cbet1:  0.9998487145115275,
			dn1:    1.0000010195104125,
			sbet2:  0.0,
			cbet2:  1.0,
			dn2:    1.0,
			lam12:  0.017453292519943295,
			slam12: 0.01745240643728351,
			clam12: 0.9998476951563913,
			C1a:    []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			C2a:    []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			want1:  -1.0,
			want2:  0.7095310092765433,
			want3:  0.7046742132893822,
			want4:  math.NaN(),
			want5:  math.NaN(),
			want6:  1.0000002548969817,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got1, got2, got3, got4, got5, got6 := geod._InverseStart(
				tC.sbet1, tC.cbet1, tC.dn1, tC.sbet2, tC.cbet2, tC.dn2,
				tC.lam12, tC.slam12, tC.clam12, tC.C1a, tC.C2a,
			)

			// Compare all return values
			if !f64_equals(tC.want1, got1) {
				t.Errorf("_InverseStart() got1 = %v, want %v", got1, tC.want1)
			}

			if !f64_equals(tC.want2, got2) {
				t.Errorf("_InverseStart() got2 = %v, want %v", got2, tC.want2)
			}

			if !f64_equals(tC.want3, got3) {
				t.Errorf("_InverseStart() got3 = %v, want %v", got3, tC.want3)
			}

			if !f64_equals(tC.want4, got4) {
				t.Errorf("_InverseStart() got4 = %v, want %v", got4, tC.want4)
			}

			if !f64_equals(tC.want5, got5) {
				t.Errorf("_InverseStart() got5 = %v, want %v", got5, tC.want5)
			}

			if !f64_equals(tC.want6, got6) {
				t.Errorf("_InverseStart() got6 = %v, want %v", got6, tC.want6)
			}
		})
	}
}

func Benchmark_InverseStart(b *testing.B) {
	geod := Wgs84()

	benchmarks := []struct {
		desc   string
		sbet1  float64
		cbet1  float64
		dn1    float64
		sbet2  float64
		cbet2  float64
		dn2    float64
		lam12  float64
		slam12 float64
		clam12 float64
		C1a    []float64
		C2a    []float64
	}{
		{
			desc:   "1",
			sbet1:  -0.017393909556108908,
			cbet1:  0.9998487145115275,
			dn1:    1.0000010195104125,
			sbet2:  0.0,
			cbet2:  1.0,
			dn2:    1.0,
			lam12:  0.017453292519943295,
			slam12: 0.01745240643728351,
			clam12: 0.9998476951563913,
			C1a:    []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			C2a:    []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				geod._InverseStart(
					bm.sbet1, bm.cbet1, bm.dn1, bm.sbet2, bm.cbet2, bm.dn2,
					bm.lam12, bm.slam12, bm.clam12, bm.C1a, bm.C2a,
				)
			}

		})
	}
}
