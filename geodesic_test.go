package geographiclibgo

import (
	"encoding/csv"
	"errors"
	"io"
	"log"
	"math"
	"os"
	"strconv"
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
		if !almost_equal(got, want, float64EqualityThreshold) {
			return false
		}
	}
	return true
}

func TestGeodesicInit(t *testing.T) {
	geod := Wgs84()
	if !almost_equal(geod.a, 6378137.0, float64EqualityThreshold) {
		t.Errorf("Wgs84() failed: a = %v, expected %v", geod.a, 6378137.0)
	}
	if !almost_equal(geod.f, 0.0033528106647474805, float64EqualityThreshold) {
		t.Errorf("Wgs84() failed: f = %v, exptected %v", geod.f, 0.0033528106647474805)
	}

	if !almost_equal(geod._f1, 0.9966471893352525, float64EqualityThreshold) {
		t.Errorf("Wgs84() failed: _f1 = %v, exptected %v", geod._f1, 0.9966471893352525)
	}
	if !almost_equal(geod._e2, 0.0066943799901413165, float64EqualityThreshold) {
		t.Errorf("Wgs84() failed: _e2 = %v, exptected %v", geod._e2, 0.0066943799901413165)
	}

	if !almost_equal(geod._ep2, 0.006739496742276434, float64EqualityThreshold) {
		t.Errorf("Wgs84() failed: _ep2 = %v, exptected %v", geod._ep2, 0.006739496742276434)
	}

	if !almost_equal(geod._n, 0.0016792203863837047, float64EqualityThreshold) {
		t.Errorf("Wgs84() failed: _n = %v, exptected %v", geod._n, 0.0016792203863837047)
	}

	if !almost_equal(geod._b, 6356752.314245179, float64EqualityThreshold) {
		t.Errorf("Wgs84() failed: _b = %v, exptected %v", geod._b, 6356752.314245179)
	}

	if !almost_equal(geod._c2, 40589732499314.76, float64EqualityThreshold) {
		t.Errorf("Wgs84() failed: _c2 = %v, exptected %v", geod._c2, 40589732499314.76)
	}

	if !almost_equal(geod._etol2, 3.6424611488788524e-08, float64EqualityThreshold) {
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
		if !almost_equal(geod._A3x[i], want_A3x[i], float64EqualityThreshold) {
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
		if !almost_equal(geod._C3x[i], want_C3x[i], float64EqualityThreshold) {
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
		if !almost_equal(geod._C4x[i], want_C4x[i], float64EqualityThreshold) {
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
	if !almost_equal(geod._A3f(0.12), 0.9363788874000158, float64EqualityThreshold) {
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
		if !almost_equal(c[i], want[i], float64EqualityThreshold) {
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
		if !almost_equal(c[i], want[i], float64EqualityThreshold) {
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

func Test_Lambda12(t *testing.T) {
	geod := Wgs84()

	testCases := []struct {
		desc    string
		sbet1   float64
		cbet1   float64
		dn1     float64
		sbet2   float64
		cbet2   float64
		dn2     float64
		salp1   float64
		calp1   float64
		slam120 float64
		clam120 float64
		diffp   bool
		C1a     []float64
		C2a     []float64
		C3a     []float64
		want1   float64
		want2   float64
		want3   float64
		want4   float64
		want5   float64
		want6   float64
		want7   float64
		want8   float64
		want9   float64
		want10  float64
		want11  float64
	}{
		{
			desc:    "1",
			sbet1:   -0.017393909556108908,
			cbet1:   0.9998487145115275,
			dn1:     1.0000010195104125,
			sbet2:   0.0,
			cbet2:   1.0,
			dn2:     1.0,
			salp1:   0.7095310092765433,
			calp1:   0.7046742132893822,
			slam120: 0.01745240643728351,
			clam120: 0.9998476951563913,
			diffp:   true,
			C1a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			C2a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			C3a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
			want1:   1.4834408705897495e-09,
			want2:   0.7094236675312185,
			want3:   0.7047822783999007,
			want4:   0.024682339962725352,
			want5:   -0.024679833885152578,
			want6:   0.9996954065111039,
			want7:   0.0,
			want8:   1.0,
			want9:   0.0008355095326524276,
			want10:  -5.8708496511415445e-05,
			want11:  0.034900275148485,
		},
		{
			desc:    "2",
			sbet1:   -0.017393909556108908,
			cbet1:   0.9998487145115275,
			dn1:     1.0000010195104125,
			sbet2:   0.0,
			cbet2:   1.0,
			dn2:     1.0,
			salp1:   0.7095309793242709,
			calp1:   0.7046742434480923,
			slam120: 0.01745240643728351,
			clam120: 0.9998476951563913,
			diffp:   true,
			C1a: []float64{
				0.0,
				-0.00041775465696698233,
				-4.362974596862037e-08,
				-1.2151022357848552e-11,
				-4.7588881620421004e-15,
				-2.226614930167366e-18,
				-1.1627237498131586e-21,
			},
			C2a: []float64{
				0.0,
				-0.0008355098973052918,
				-1.7444619952659748e-07,
				-7.286557795511902e-11,
				-3.80472772706481e-14,
				-2.2251271876594078e-17,
				1.2789961247944744e-20,
			},
			C3a: []float64{
				0.0,
				0.00020861391868413911,
				4.3547247296823945e-08,
				1.515432276542012e-11,
				6.645637323698485e-15,
				3.3399223952510497e-18,
			},
			want1:  6.046459990680098e-17,
			want2:  0.7094236375834774,
			want3:  0.7047823085448635,
			want4:  0.024682338906797385,
			want5:  -0.02467983282954624,
			want6:  0.9996954065371639,
			want7:  0.0,
			want8:  1.0,
			want9:  0.0008355096040059597,
			want10: -5.870849152149326e-05,
			want11: 0.03490027216297455,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got1, got2, got3, got4, got5, got6, got7, got8, got9, got10, got11 := geod._Lambda12(
				tC.sbet1, tC.cbet1, tC.dn1, tC.sbet2, tC.cbet2, tC.dn2,
				tC.salp1, tC.calp1, tC.slam120, tC.clam120, tC.diffp,
				tC.C1a, tC.C2a, tC.C3a,
			)

			if !f64_equals(tC.want1, got1) {
				t.Errorf("_Lambda12() got1 = %v, want %v", got1, tC.want1)
			}

			if !f64_equals(tC.want2, got2) {
				t.Errorf("_Lambda12() got2 = %v, want %v", got2, tC.want2)
			}

			if !f64_equals(tC.want3, got3) {
				t.Errorf("_Lambda12() got3 = %v, want %v", got3, tC.want3)
			}

			if !f64_equals(tC.want4, got4) {
				t.Errorf("_Lambda12() got4 = %v, want %v", got4, tC.want4)
			}

			if !f64_equals(tC.want5, got5) {
				t.Errorf("_Lambda12() got5 = %v, want %v", got5, tC.want5)
			}

			if !f64_equals(tC.want6, got6) {
				t.Errorf("_Lambda12() got6 = %v, want %v", got6, tC.want6)
			}

			if !f64_equals(tC.want7, got7) {
				t.Errorf("_Lambda12() got7 = %v, want %v", got7, tC.want7)
			}

			if !f64_equals(tC.want8, got8) {
				t.Errorf("_Lambda12() got8 = %v, want %v", got8, tC.want8)
			}

			if !f64_equals(tC.want9, got9) {
				t.Errorf("_Lambda12() got9 = %v, want %v", got9, tC.want9)
			}

			if !f64_equals(tC.want10, got10) {
				t.Errorf("_Lambda12() got10 = %v, want %v", got10, tC.want10)
			}

			if !f64_equals(tC.want11, got11) {
				t.Errorf("_Lambda12() got11 = %v, want %v", got11, tC.want11)
			}
		})
	}
}

func Benchmark_Lambda12(b *testing.B) {
	geod := Wgs84()

	benchmarks := []struct {
		desc    string
		sbet1   float64
		cbet1   float64
		dn1     float64
		sbet2   float64
		cbet2   float64
		dn2     float64
		salp1   float64
		calp1   float64
		slam120 float64
		clam120 float64
		diffp   bool
		C1a     []float64
		C2a     []float64
		C3a     []float64
	}{
		{
			desc:    "1",
			sbet1:   -0.017393909556108908,
			cbet1:   0.9998487145115275,
			dn1:     1.0000010195104125,
			sbet2:   0.0,
			cbet2:   1.0,
			dn2:     1.0,
			salp1:   0.7095310092765433,
			calp1:   0.7046742132893822,
			slam120: 0.01745240643728351,
			clam120: 0.9998476951563913,
			diffp:   true,
			C1a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			C2a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
			C3a:     []float64{0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
		},
		{
			desc:    "2",
			sbet1:   -0.017393909556108908,
			cbet1:   0.9998487145115275,
			dn1:     1.0000010195104125,
			sbet2:   0.0,
			cbet2:   1.0,
			dn2:     1.0,
			salp1:   0.7095309793242709,
			calp1:   0.7046742434480923,
			slam120: 0.01745240643728351,
			clam120: 0.9998476951563913,
			diffp:   true,
			C1a: []float64{
				0.0,
				-0.00041775465696698233,
				-4.362974596862037e-08,
				-1.2151022357848552e-11,
				-4.7588881620421004e-15,
				-2.226614930167366e-18,
				-1.1627237498131586e-21,
			},
			C2a: []float64{
				0.0,
				-0.0008355098973052918,
				-1.7444619952659748e-07,
				-7.286557795511902e-11,
				-3.80472772706481e-14,
				-2.2251271876594078e-17,
				1.2789961247944744e-20,
			},
			C3a: []float64{
				0.0,
				0.00020861391868413911,
				4.3547247296823945e-08,
				1.515432276542012e-11,
				6.645637323698485e-15,
				3.3399223952510497e-18,
			},
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				geod._Lambda12(
					bm.sbet1, bm.cbet1, bm.dn1, bm.sbet2, bm.cbet2, bm.dn2,
					bm.salp1, bm.calp1, bm.slam120, bm.clam120, bm.diffp,
					bm.C1a, bm.C2a, bm.C3a,
				)
			}

		})
	}
}

func Test_gen_inverse(t *testing.T) {
	testCases := []struct {
		desc                                                     string
		geodesic                                                 Geodesic
		lat1, lon1, lat2, lon2                                   float64
		outmask                                                  uint64
		a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12 float64
	}{
		{
			desc:     "1",
			geodesic: Wgs84(),
			lat1:     0.0,
			lon1:     0.0,
			lat2:     1.0,
			lon2:     1.0,
			outmask:  STANDARD,
			a12:      1.4141938478710363,
			s12:      156899.56829134026,
			salp1:    0.7094236375834774,
			calp1:    0.7047823085448635,
			salp2:    0.7095309793242709,
			calp2:    0.7046742434480923,
			m12:      math.NaN(),
			M12:      math.NaN(),
			M21:      math.NaN(),
			S12:      math.NaN(),
		},
		{
			desc:     "2",
			geodesic: Wgs84(),
			lat1:     0.0,
			lon1:     539.0,
			lat2:     0.0,
			lon2:     181.0,
			outmask:  STANDARD,
			a12:      2.0067281796419527,
			s12:      222638.98158654713,
			salp1:    1,
			calp1:    0,
			salp2:    1,
			calp2:    0,
			m12:      math.NaN(),
			M12:      math.NaN(),
			M21:      math.NaN(),
			S12:      math.NaN(),
		},
		{
			desc:     "3",
			geodesic: Wgs84(),
			lat1:     0.0,
			lon1:     539.0,
			lat2:     0.0,
			lon2:     181.0,
			outmask:  STANDARD | LONG_UNROLL,
			a12:      2.0067281796419527,
			s12:      222638.98158654713,
			salp1:    1,
			calp1:    0,
			salp2:    1,
			calp2:    0,
			m12:      math.NaN(),
			M12:      math.NaN(),
			M21:      math.NaN(),
			S12:      math.NaN(),
		},
		{
			desc:     "4",
			geodesic: NewGeodesic(6.4e6, 0.0),
			lat1:     1,
			lon1:     2,
			lat2:     3,
			lon2:     4,
			outmask:  AREA,
			a12:      2.8274936682015439,
			s12:      math.NaN(),
			salp1:    0.70651412839435823,
			calp1:    0.70769893766993908,
			salp2:    0.70737595703509237,
			calp2:    0.70683750283122859,
			m12:      math.NaN(),
			M12:      math.NaN(),
			M21:      math.NaN(),
			S12:      49911046114.952187,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12 := tC.geodesic._gen_inverse(
				tC.lat1, tC.lon1, tC.lat2, tC.lon2, tC.outmask,
			)

			if !f64_equals(tC.a12, a12) {
				t.Errorf("_gen_inverse() a12 = %v, want %v", a12, tC.a12)
			}

			if !f64_equals(tC.s12, s12) {
				t.Errorf("_gen_inverse() s12 = %v, want %v", s12, tC.s12)
			}

			if !f64_equals(tC.salp1, salp1) {
				t.Errorf("_gen_inverse() salp1 = %v, want %v", salp1, tC.salp1)
			}

			if !f64_equals(tC.calp1, calp1) {
				t.Errorf("_gen_inverse() calp1 = %v, want %v", calp1, tC.calp1)
			}

			if !f64_equals(tC.salp2, salp2) {
				t.Errorf("_gen_inverse() salp2 = %v, want %v", salp2, tC.salp2)
			}

			if !f64_equals(tC.calp2, calp2) {
				t.Errorf("_gen_inverse() calp2 = %v, want %v", calp2, tC.calp2)
			}

			if !f64_equals(tC.m12, m12) {
				t.Errorf("_gen_inverse() m12 = %v, want %v", m12, tC.m12)
			}

			if !f64_equals(tC.M12, M12) {
				t.Errorf("_gen_inverse() M12 = %v, want %v", M12, tC.M12)
			}

			if !f64_equals(tC.M21, M21) {
				t.Errorf("_gen_inverse() M21 = %v, want %v", M21, tC.M21)
			}

			if !f64_equals(tC.S12, S12) {
				t.Errorf("_gen_inverse() S12 = %v, want %v", S12, tC.S12)
			}
		})
	}
}

func Benchmark_gen_inverse(b *testing.B) {
	benchmarks := []struct {
		desc                   string
		geodesic               Geodesic
		lat1, lon1, lat2, lon2 float64
		outmask                uint64
	}{
		{
			desc:     "1",
			geodesic: Wgs84(),
			lat1:     0.0,
			lon1:     0.0,
			lat2:     1.0,
			lon2:     1.0,
			outmask:  STANDARD,
		},
		{
			desc:     "2",
			geodesic: Wgs84(),
			lat1:     0.0,
			lon1:     539.0,
			lat2:     0.0,
			lon2:     181.0,
			outmask:  STANDARD,
		},
		{
			desc:     "3",
			geodesic: Wgs84(),
			lat1:     0.0,
			lon1:     539.0,
			lat2:     0.0,
			lon2:     181.0,
			outmask:  STANDARD | LONG_UNROLL,
		},
		{
			desc:     "4",
			geodesic: NewGeodesic(6.4e6, 0.0),
			lat1:     1,
			lon1:     2,
			lat2:     3,
			lon2:     4,
			outmask:  AREA,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				bm.geodesic._gen_inverse(
					bm.lat1, bm.lon1, bm.lat2, bm.lon2, bm.outmask,
				)
			}

		})
	}
}

func Test_gen_inverse_azi(t *testing.T) {
	testCases := []struct {
		desc                                     string
		geodesic                                 Geodesic
		lat1, lon1, lat2, lon2                   float64
		outmask                                  uint64
		a12, s12, azi1, azi2, m12, M12, M21, S12 float64
	}{
		{
			desc:     "1",
			geodesic: Wgs84(),
			lat1:     54.1589,
			lon1:     15.3872,
			lat2:     54.1591,
			lon2:     15.3877,
			outmask:  ALL,
			a12:      0.00035549325126265914,
			s12:      39.527686386021514,
			azi1:     55.723110355324408,
			azi2:     55.72351567783663,
			m12:      39.527686385839779,
			M12:      0.99999999998083666,
			M21:      0.99999999998083666,
			S12:      286698586.30231881,
		},
		{
			desc:     "2",
			geodesic: Wgs84(),
			lat1:     20.001,
			lon1:     0.0,
			lat2:     20.001,
			lon2:     0.0,
			outmask:  ALL,
			a12:      0,
			s12:      0,
			azi1:     180.0,
			azi2:     180.0,
			m12:      0,
			M12:      0.99999999999999988,
			M21:      0.99999999999999988,
			S12:      0,
		},
		{
			desc:     "3",
			geodesic: Wgs84(),
			lat1:     90.0,
			lon1:     0.0,
			lat2:     90.0,
			lon2:     180.0,
			outmask:  ALL,
			a12:      0,
			s12:      0,
			azi1:     0.0,
			azi2:     180.0,
			m12:      0,
			M12:      1.0,
			M21:      1.0,
			S12:      127516405431022.11,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			a12, s12, azi1, azi2, m12, M12, M21, S12 := tC.geodesic._gen_inverse_azi(
				tC.lat1, tC.lon1, tC.lat2, tC.lon2, tC.outmask,
			)

			if !f64_equals(tC.a12, a12) {
				t.Errorf("_gen_inverse() a12 = %v, want %v", a12, tC.a12)
			}

			if !f64_equals(tC.s12, s12) {
				t.Errorf("_gen_inverse() s12 = %v, want %v", s12, tC.s12)
			}

			if !f64_equals(tC.azi1, azi1) {
				t.Errorf("_gen_inverse() azi1 = %v, want %v", azi1, tC.azi1)
			}

			if !f64_equals(tC.azi2, azi2) {
				t.Errorf("_gen_inverse() azi2 = %v, want %v", azi2, tC.azi2)
			}

			if !f64_equals(tC.m12, m12) {
				t.Errorf("_gen_inverse() m12 = %v, want %v", m12, tC.m12)
			}

			if !f64_equals(tC.M12, M12) {
				t.Errorf("_gen_inverse() M12 = %v, want %v", M12, tC.M12)
			}

			if !f64_equals(tC.M21, M21) {
				t.Errorf("_gen_inverse() M21 = %v, want %v", M21, tC.M21)
			}

			if !f64_equals(tC.S12, S12) {
				t.Errorf("_gen_inverse() S12 = %v, want %v", S12, tC.S12)
			}

		})
	}
}

func Benchmark_gen_inverse_azi(b *testing.B) {
	benchmarks := []struct {
		desc                   string
		geodesic               Geodesic
		lat1, lon1, lat2, lon2 float64
		outmask                uint64
	}{
		{
			desc:     "1",
			geodesic: Wgs84(),
			lat1:     54.1589,
			lon1:     15.3872,
			lat2:     54.1591,
			lon2:     15.3877,
			outmask:  ALL,
		},
		{
			desc:     "2",
			geodesic: Wgs84(),
			lat1:     20.001,
			lon1:     0.0,
			lat2:     20.001,
			lon2:     0.0,
			outmask:  ALL,
		},
		{
			desc:     "3",
			geodesic: Wgs84(),
			lat1:     90.0,
			lon1:     0.0,
			lat2:     90.0,
			lon2:     180.0,
			outmask:  ALL,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				bm.geodesic._gen_inverse_azi(
					bm.lat1, bm.lon1, bm.lat2, bm.lon2, bm.outmask,
				)
			}

		})
	}
}

func Test_gen_direct(t *testing.T) {
	testCases := []struct {
		desc                  string
		geod                  Geodesic
		lat1, lon1, azi1      float64
		arcmode               bool
		s12_a12               float64
		outmask               uint64
		a12, lat2, lon2, azi2 float64
		s12, m12              float64
		M12, M21, S12         float64
	}{
		{
			desc:    "geodsolve28",
			geod:    NewGeodesic(6.4e6, 0.1),
			lat1:    1.0,
			lon1:    2.0,
			azi1:    10.0,
			arcmode: false,
			s12_a12: 5e6,
			outmask: STANDARD,
			a12:     48.555706903681603,
			lat2:    51.434124142195124,
			lon2:    12.486884576686398,
			azi2:    15.178966424439388,
			s12:     5000000,
			m12:     math.NaN(),
			M12:     math.NaN(),
			M21:     math.NaN(),
			S12:     math.NaN(),
		},
		{
			desc:    "geodsolve61",
			geod:    Wgs84(),
			lat1:    45.0,
			lon1:    0.0,
			azi1:    -0.000000000000000003,
			arcmode: false,
			s12_a12: 1e7,
			outmask: STANDARD | LONG_UNROLL,
			a12:     89.88610143063056,
			lat2:    45.306319097990396,
			lon2:    -180,
			azi2:    -180,
			s12:     1.0e7,
			m12:     math.NaN(),
			M12:     math.NaN(),
			M21:     math.NaN(),
			S12:     math.NaN(),
		},
		{
			desc:    "geodsolve17",
			geod:    NewGeodesic(6.4e6, -1/150.0),
			lat1:    40.0,
			lon1:    -75.0,
			azi1:    -10.0,
			arcmode: false,
			s12_a12: 2e7,
			outmask: STANDARD | LONG_UNROLL,
			a12:     178.44443367991155,
			lat2:    -38.469547002608067,
			lon2:    -254.81220904664144,
			azi2:    -170.21964804643204,
			s12:     2.0e7,
			m12:     math.NaN(),
			M12:     math.NaN(),
			M21:     math.NaN(),
			S12:     math.NaN(),
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			a12, lat2, lon2, azi2, s12, m12, M12, M21, S12, _ := tC.geod._gen_direct(
				tC.lat1, tC.lon1, tC.azi1, tC.arcmode, tC.s12_a12, tC.outmask,
			)

			if !f64_equals(tC.a12, a12) {
				t.Errorf("_gen_direct() a12 = %v, want %v", a12, tC.a12)
			}

			if !f64_equals(tC.lat2, lat2) {
				t.Errorf("_gen_direct() lat2 = %v, want %v", lat2, tC.lat2)
			}

			if !f64_equals(tC.lon2, lon2) {
				t.Errorf("_gen_direct() lon2 = %v, want %v", lon2, tC.lon2)
			}

			if !f64_equals(tC.azi2, azi2) {
				t.Errorf("_gen_direct() azi2 = %v, want %v", azi2, tC.azi2)
			}

			if !f64_equals(tC.s12, s12) {
				t.Errorf("_gen_direct() s12 = %v, want %v", s12, tC.s12)
			}

			if !f64_equals(tC.m12, m12) {
				t.Errorf("_gen_direct() m12 = %v, want %v", m12, tC.m12)
			}

			if !f64_equals(tC.M12, M12) {
				t.Errorf("_gen_direct() M12 = %v, want %v", M12, tC.M12)
			}

			if !f64_equals(tC.M21, M21) {
				t.Errorf("_gen_direct() M21 = %v, want %v", M21, tC.M21)
			}

			if !f64_equals(tC.S12, S12) {
				t.Errorf("_gen_direct() S12 = %v, want %v", S12, tC.S12)
			}

		})
	}
}

func Benchmark_gen_direct(b *testing.B) {
	benchmarks := []struct {
		desc             string
		geod             Geodesic
		lat1, lon1, azi1 float64
		arcmode          bool
		s12_a12          float64
		outmask          uint64
	}{
		{
			desc:    "geodsolve28",
			geod:    NewGeodesic(6.4e6, 0.1),
			lat1:    1.0,
			lon1:    2.0,
			azi1:    10.0,
			arcmode: false,
			s12_a12: 5e6,
			outmask: STANDARD,
		},
		{
			desc:    "geodsolve61",
			geod:    Wgs84(),
			lat1:    45.0,
			lon1:    0.0,
			azi1:    -0.000000000000000003,
			arcmode: false,
			s12_a12: 1e7,
			outmask: STANDARD | LONG_UNROLL,
		},
		{
			desc:    "geodsolve17",
			geod:    NewGeodesic(6.4e6, -1/150.0),
			lat1:    40.0,
			lon1:    -75.0,
			azi1:    -10.0,
			arcmode: false,
			s12_a12: 2e7,
			outmask: STANDARD | LONG_UNROLL,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				bm.geod._gen_direct(
					bm.lat1, bm.lon1, bm.azi1, bm.arcmode, bm.s12_a12, bm.outmask,
				)
			}

		})
	}
}

func TestDirectBadInputs(t *testing.T) {
	geod := Wgs84()
	testCases := []struct {
		desc    string
		lat     float64
		lon     float64
		azi     float64
		s12     float64
		wantlat float64
		wantlon float64
		wantazi float64
	}{
		{
			desc:    "inf s12",
			lat:     0.0,
			lon:     0.0,
			azi:     90.0,
			s12:     math.Inf(0),
			wantlat: math.NaN(),
			wantlon: math.NaN(),
			wantazi: math.NaN(),
		},
		{
			desc:    "nan s12",
			lat:     0.0,
			lon:     0.0,
			azi:     90.0,
			s12:     math.NaN(),
			wantlat: math.NaN(),
			wantlon: math.NaN(),
			wantazi: math.NaN(),
		},
		{
			desc:    "inf azi",
			lat:     0.0,
			lon:     0.0,
			azi:     math.Inf(0),
			s12:     1000.0,
			wantlat: math.NaN(),
			wantlon: math.NaN(),
			wantazi: math.NaN(),
		},
		{
			desc:    "nan azi",
			lat:     0.0,
			lon:     0.0,
			azi:     math.NaN(),
			s12:     1000.0,
			wantlat: math.NaN(),
			wantlon: math.NaN(),
			wantazi: math.NaN(),
		},
		{
			desc:    "inf lon",
			lat:     0.0,
			lon:     math.Inf(0),
			azi:     90.0,
			s12:     1000.0,
			wantlat: 0.0,
			wantlon: math.NaN(),
			wantazi: 90,
		},
		{
			desc:    "nan lon",
			lat:     0.0,
			lon:     math.NaN(),
			azi:     90.0,
			s12:     1000.0,
			wantlat: 0.0,
			wantlon: math.NaN(),
			wantazi: 90,
		},
		{
			desc:    "inf lat",
			lat:     math.Inf(0),
			lon:     0.0,
			azi:     90.0,
			s12:     1000.0,
			wantlat: math.NaN(),
			wantlon: math.NaN(),
			wantazi: math.NaN(),
		},
		{
			desc:    "nan lat",
			lat:     math.NaN(),
			lon:     0.0,
			azi:     90.0,
			s12:     1000.0,
			wantlat: math.NaN(),
			wantlon: math.NaN(),
			wantazi: math.NaN(),
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			latlonazi := geod.DirectCalcLatLonAzi(tC.lat, tC.lon, tC.azi, tC.s12)

			if !f64_equals(tC.wantlat, latlonazi.LatDeg) {
				t.Errorf("Direct() lat = %v; want %v", latlonazi.LatDeg, tC.wantlat)
			}

			if !f64_equals(tC.wantlon, latlonazi.LonDeg) {
				t.Errorf("Direct() lon = %v; want %v", latlonazi.LonDeg, tC.wantlon)
			}

			if !f64_equals(tC.wantazi, latlonazi.AziDeg) {
				t.Errorf("Direct() lon = %v; want %v", latlonazi.AziDeg, tC.wantazi)
			}
		})
	}
}

func BenchmarkDirectBadInputs(b *testing.B) {
	geod := Wgs84()
	benchmarks := []struct {
		desc string
		lat  float64
		lon  float64
		azi  float64
		s12  float64
	}{
		{
			desc: "inf s12",
			lat:  0.0,
			lon:  0.0,
			azi:  90.0,
			s12:  math.Inf(0),
		},
		{
			desc: "nan s12",
			lat:  0.0,
			lon:  0.0,
			azi:  90.0,
			s12:  math.NaN(),
		},
		{
			desc: "inf azi",
			lat:  0.0,
			lon:  0.0,
			azi:  math.Inf(0),
			s12:  1000.0,
		},
		{
			desc: "nan azi",
			lat:  0.0,
			lon:  0.0,
			azi:  math.NaN(),
			s12:  1000.0,
		},
		{
			desc: "inf lon",
			lat:  0.0,
			lon:  math.Inf(0),
			azi:  90.0,
			s12:  1000.0,
		},
		{
			desc: "nan lon",
			lat:  0.0,
			lon:  math.NaN(),
			azi:  90.0,
			s12:  1000.0,
		},
		{
			desc: "inf lat",
			lat:  math.Inf(0),
			lon:  0.0,
			azi:  90.0,
			s12:  1000.0,
		},
		{
			desc: "nan lat",
			lat:  math.NaN(),
			lon:  0.0,
			azi:  90.0,
			s12:  1000.0,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				geod.DirectCalcLatLonAzi(bm.lat, bm.lon, bm.azi, bm.s12)
			}
		})
	}
}

type dat_struct struct {
	lat1     float64 // latitude at point 1, lat1 (degrees, exact)
	lon1     float64 // longitude at point 1, lon1 (degrees, always 0)
	azi1     float64 // azimuth at point 1, azi1 (clockwise from north in degrees, exact)
	lat2     float64 // latitude at point 2, lat2 (degrees, accurate to 10−18 deg)
	lon2     float64 // longitude at point 2, lon2 (degrees, accurate to 10−18 deg)
	azi2     float64 // azimuth at point 2, azi2 (degrees, accurate to 10−18 deg)
	s12_dist float64 // geodesic distance from point 1 to point 2, s12 (meters, exact)
	a12      float64 // arc distance on the auxiliary sphere, a12 (degrees, accurate to 10−18 deg)
	m12      float64 // reduced length of the geodesic, m12 (meters, accurate to 0.1 pm)
	s12_area float64 // the area under the geodesic, S12 (m2, accurate to 1 mm2)
}

func parse_dat_struct(record []string) (dat_struct, error) {
	var result dat_struct
	var err error

	result.lat1, err = strconv.ParseFloat(record[0], 64)
	if err != nil {
		return result, err
	}

	result.lon1, err = strconv.ParseFloat(record[1], 64)
	if err != nil {
		return result, err
	}

	result.azi1, err = strconv.ParseFloat(record[2], 64)
	if err != nil {
		return result, err
	}

	result.lat2, err = strconv.ParseFloat(record[3], 64)
	if err != nil {
		return result, err
	}

	result.lon2, err = strconv.ParseFloat(record[4], 64)
	if err != nil {
		return result, err
	}

	result.azi2, err = strconv.ParseFloat(record[5], 64)
	if err != nil {
		return result, err
	}

	result.s12_dist, err = strconv.ParseFloat(record[6], 64)
	if err != nil {
		return result, err
	}

	result.a12, err = strconv.ParseFloat(record[7], 64)
	if err != nil {
		return result, err
	}

	result.m12, err = strconv.ParseFloat(record[8], 64)
	if err != nil {
		return result, err
	}

	result.s12_area, err = strconv.ParseFloat(record[9], 64)
	if err != nil {
		return result, err
	}

	return result, err

}

func read_dot_dat(filename string) []dat_struct {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		log.Fatalf("Could not open file %s", filename)
	}
	defer file.Close()

	// Create the CSV reader and the slice of data to fill
	reader := csv.NewReader(file)
	reader.Comma = ' '
	reader.FieldsPerRecord = 10
	data := make([]dat_struct, 0)

	// For each row
	for {
		// Read the record
		record, err := reader.Read()

		// Check for end of file
		if err == io.EOF {
			break
		}

		// Check for error
		if err != nil {
			log.Fatal(err)
		}

		// Parse everything to floats
		test_case, err := parse_dat_struct(record)
		if err != nil {
			log.Fatal(err)
		}

		data = append(data, test_case)

	}

	return data
}

func TestDirect100(t *testing.T) {
	// Create a geod with which to do the computations
	geod := Wgs84()

	// Get the set of 100 test cases from the test_fixtures folder
	filename := "test_fixtures/GeodTest-100.dat"

	// Read in the file
	test_cases := read_dot_dat(filename)

	// For each test case, run the direct computation, and compare the results
	for row_num, tc := range test_cases {
		res := geod.DirectCalcAll(tc.lat1, tc.lon1, tc.azi1, tc.s12_dist)

		// Check each result. Note that there should not be any NaN results from this file
		if !almost_equal(res.LatDeg, tc.lat2, float64EqualityThreshold) {
			t.Errorf("Row: %d -- Direct() Lat2 = %v; want %v", row_num, res.LatDeg, tc.lat2)
		}

		if !almost_equal(res.LonDeg, tc.lon2, float64EqualityThreshold) {
			t.Errorf("Row: %d -- Direct() Lon2 = %v; want %v", row_num, res.LonDeg, tc.lon2)
		}

		if !almost_equal(res.AziDeg, tc.azi2, float64EqualityThreshold) {
			t.Errorf("Row: %d -- Direct() AziDeg = %v; want %v", row_num, res.AziDeg, tc.azi2)
		}

		if !almost_equal(res.ReducedLengthM, tc.m12, float64EqualityThreshold) {
			t.Errorf("Row: %d -- Direct() ReducedLengthM = %v; want %v", row_num, res.ReducedLengthM, tc.m12)
		}

		if !almost_equal(res.S12M2, tc.s12_area, 8.8e2) {
			t.Errorf("Row: %d -- Direct() S12M2 = %v; want %v", row_num, res.S12M2, tc.s12_area)
		}

		if !almost_equal(res.A12Deg, tc.a12, float64EqualityThreshold) {
			t.Errorf("Row: %d -- Direct() A12Deg = %v; want %v", row_num, res.A12Deg, tc.a12)
		}

	}
}

func BenchmarkDirect100(b *testing.B) {
	// Create a geod with which to do the computations
	geod := Wgs84()

	// Get the set of 100 test cases from the test_fixtures folder
	filename := "test_fixtures/GeodTest-100.dat"

	// Read in the file
	test_cases := read_dot_dat(filename)

	for i := 0; i < b.N; i++ {
		for _, tc := range test_cases {
			geod.DirectCalcLatLonAzi(tc.lat1, tc.lon1, tc.azi1, tc.s12_dist)
		}
	}
}

func file_exists(filename string) bool {
	_, err := os.Open(filename)
	return !errors.Is(err, os.ErrNotExist)
}

func Test_file_exists(t *testing.T) {
	testCases := []struct {
		desc     string
		filename string
		want     bool
	}{
		{
			desc:     "yes",
			filename: "test_fixtures/GeodTest-100.dat",
			want:     true,
		},
		{
			desc:     "no",
			filename: "somemadeupfile.txt",
			want:     false,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := file_exists(tC.filename)

			if tC.want != got {
				t.Errorf("file_exists() = %v; want %v", got, tC.want)
			}
		})
	}
}

func TestDirectAll(t *testing.T) {
	// Create a geod with which to do the computations
	// geod := Wgs84()

	// Get the set of 100 test cases from the test_fixtures folder
	filename := "test_fixtures/GeodTest-100.dat"

	// Check if the file is there
	if !file_exists(filename) {
		t.SkipNow()
	}

	// If not, skip this test
}

func TestInverseLength(t *testing.T) {
	geod := Wgs84()

	testCases := []struct {
		desc       string
		lat1, lon1 float64
		lat2, lon2 float64
		wants12    float64
	}{
		{
			desc:    "short line bug",
			lat1:    36.493349428792,
			lon1:    0.0,
			lat2:    36.49334942879201,
			lon2:    0.0000008,
			wants12: 0.072,
		},
		{
			desc:    "volatile sbet12a bug test 1",
			lat1:    88.202499451857,
			lon1:    0.0,
			lat2:    -88.202499451857,
			lon2:    179.981022032992859592,
			wants12: 20003898.214,
		},
		{
			desc:    "volatile sbet12a bug test 2",
			lat1:    89.333123580033,
			lon1:    0.0,
			lat2:    -89.333123580032997687,
			lon2:    179.99295812360148422,
			wants12: 20003926.881,
		},
		{
			desc:    "volatile x bug",
			lat1:    56.320923501171,
			lon1:    0.0,
			lat2:    -56.320923501171,
			lon2:    179.664747671772880215,
			wants12: 19993558.287,
		},
		{
			desc:    "tol1_ bug",
			lat1:    52.784459512564,
			lon1:    0.0,
			lat2:    -52.784459512563990912,
			lon2:    179.634407464943777557,
			wants12: 19991596.095,
		},
		{
			desc:    "bet2 = -bet1 bug",
			lat1:    48.522876735459,
			lon1:    0.0,
			lat2:    -48.52287673545898293,
			lon2:    179.599720456223079643,
			wants12: 19989144.774,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			s12 := geod.InverseCalcDistance(tC.lat1, tC.lon1, tC.lat2, tC.lon2)

			if !almost_equal(tC.wants12, s12, 0.5e-3) {
				t.Errorf("InverseCalcDistance() s12 = %v; want %v", s12, tC.wants12)
			}
		})
	}
}

func BenchmarkInverseLength(b *testing.B) {
	geod := Wgs84()
	benchmarks := []struct {
		desc       string
		lat1, lon1 float64
		lat2, lon2 float64
	}{
		{
			desc: "short line bug",
			lat1: 36.493349428792,
			lon1: 0.0,
			lat2: 36.49334942879201,
			lon2: 0.0000008,
		},
		{
			desc: "volatile sbet12a bug test 1",
			lat1: 88.202499451857,
			lon1: 0.0,
			lat2: -88.202499451857,
			lon2: 179.981022032992859592,
		},
		{
			desc: "volatile sbet12a bug test 2",
			lat1: 89.333123580033,
			lon1: 0.0,
			lat2: -89.333123580032997687,
			lon2: 179.99295812360148422,
		},
		{
			desc: "volatile x bug",
			lat1: 56.320923501171,
			lon1: 0.0,
			lat2: -56.320923501171,
			lon2: 179.664747671772880215,
		},
		{
			desc: "tol1_ bug",
			lat1: 52.784459512564,
			lon1: 0.0,
			lat2: -52.784459512563990912,
			lon2: 179.634407464943777557,
		},
		{
			desc: "bet2 = -bet1 bug",
			lat1: 48.522876735459,
			lon1: 0.0,
			lat2: -48.52287673545898293,
			lon2: 179.599720456223079643,
		},
	}

	for _, bm := range benchmarks {
		b.Run(bm.desc, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				geod.DirectCalcLatLonAzi(bm.lat1, bm.lon1, bm.lat2, bm.lon2)
			}
		})
	}
}

func TestInverse100(t *testing.T) {
	// Create a geod with which to do the computations
	geod := Wgs84()

	// Get the set of 100 test cases from the test_fixtures folder
	filename := "test_fixtures/GeodTest-100.dat"

	// Read in the file
	test_cases := read_dot_dat(filename)

	// For each test case, run the direct computation, and compare the results
	for row_num, tc := range test_cases {
		res := geod.InverseCalcAll(tc.lat1, tc.lon1, tc.lat2, tc.lon2)

		// Check each result. Note that there should not be any NaN results from this file
		if !almost_equal(res.DistanceM, tc.s12_dist, float64EqualityThreshold) {
			t.Errorf("Row: %d -- Inverse() DistanceM = %v; want %v", row_num, res.DistanceM, tc.s12_dist)
		}

		if !almost_equal(res.Azimuth1Deg, tc.azi1, 1e-3) {
			t.Errorf("Row: %d -- Inverse() Azimuth1Deg = %v; want %v", row_num, res.Azimuth1Deg, tc.azi1)
		}

		if !almost_equal(res.Azimuth2Deg, tc.azi2, 1e-3) {
			t.Errorf("Row: %d -- Inverse() Azimuth2Deg = %v; want %v", row_num, res.Azimuth2Deg, tc.azi2)
		}

		if !almost_equal(res.ArcLengthDeg, tc.a12, float64EqualityThreshold) {
			t.Errorf("Row: %d -- Inverse() ArcLengthDeg = %v; want %v", row_num, res.ArcLengthDeg, tc.a12)
		}

		if !almost_equal(res.ReducedLengthM, tc.m12, float64EqualityThreshold) {
			t.Errorf("Row: %d -- Inverse() ReducedLengthM = %v; want %v", row_num, res.ReducedLengthM, tc.m12)
		}

		if !almost_equal(res.S12M2, tc.s12_area, 1e3) {
			t.Errorf("Row: %d -- Inverse() S12M2 = %v; want %v", row_num, res.S12M2, tc.s12_area)
		}

	}
}
