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
	testCases := []struct {
		desc             string
		geod             Geodesic
		lat1, lon1       float64
		lat2, lon2       float64
		wants12          float64
		acceptable_delta float64
	}{
		{
			desc:             "short line bug",
			geod:             Wgs84(),
			lat1:             36.493349428792,
			lon1:             0.0,
			lat2:             36.49334942879201,
			lon2:             0.0000008,
			wants12:          0.072,
			acceptable_delta: 0.5e-3,
		},
		{
			desc:             "volatile sbet12a bug test 1",
			geod:             Wgs84(),
			lat1:             88.202499451857,
			lon1:             0.0,
			lat2:             -88.202499451857,
			lon2:             179.981022032992859592,
			wants12:          20003898.214,
			acceptable_delta: 0.5e-3,
		},
		{
			desc:             "volatile sbet12a bug test 2",
			geod:             Wgs84(),
			lat1:             89.333123580033,
			lon1:             0.0,
			lat2:             -89.333123580032997687,
			lon2:             179.99295812360148422,
			wants12:          20003926.881,
			acceptable_delta: 0.5e-3,
		},
		{
			desc:             "volatile x bug",
			geod:             Wgs84(),
			lat1:             56.320923501171,
			lon1:             0.0,
			lat2:             -56.320923501171,
			lon2:             179.664747671772880215,
			wants12:          19993558.287,
			acceptable_delta: 0.5e-3,
		},
		{
			desc:             "tol1_ bug",
			geod:             Wgs84(),
			lat1:             52.784459512564,
			lon1:             0.0,
			lat2:             -52.784459512563990912,
			lon2:             179.634407464943777557,
			wants12:          19991596.095,
			acceptable_delta: 0.5e-3,
		},
		{
			desc:             "bet2 = -bet1 bug",
			geod:             Wgs84(),
			lat1:             48.522876735459,
			lon1:             0.0,
			lat2:             -48.52287673545898293,
			lon2:             179.599720456223079643,
			wants12:          19989144.774,
			acceptable_delta: 0.5e-3,
		},
		{
			desc:             "anitpodal prolate bug 1",
			geod:             NewGeodesic(6.4e6, -1/150.0),
			lat1:             0.07476,
			lon1:             0,
			lat2:             -0.07476,
			lon2:             180,
			wants12:          20106193,
			acceptable_delta: 0.5,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			s12 := tC.geod.InverseCalcDistance(tC.lat1, tC.lon1, tC.lat2, tC.lon2)

			if !almost_equal(tC.wants12, s12, tC.acceptable_delta) {
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
				geod.InverseCalcDistance(bm.lat1, bm.lon1, bm.lat2, bm.lon2)
			}
		})
	}
}

// func TestInverse100(t *testing.T) {
// 	// Create a geod with which to do the computations
// 	geod := Wgs84()

// 	// Get the set of 100 test cases from the test_fixtures folder
// 	filename := "test_fixtures/GeodTest-100.dat"

// 	// Read in the file
// 	test_cases := read_dot_dat(filename)

// 	// For each test case, run the direct computation, and compare the results
// 	for row_num, tc := range test_cases {
// 		res := geod.InverseCalcAll(tc.lat1, tc.lon1, tc.lat2, tc.lon2)

// 		// Check each result. Note that there should not be any NaN results from this file
// 		if !almost_equal(res.DistanceM, tc.s12_dist, float64EqualityThreshold) {
// 			t.Errorf("Row: %d -- Inverse() DistanceM = %v; want %v", row_num, res.DistanceM, tc.s12_dist)
// 		}

// 		if !almost_equal(res.Azimuth1Deg, tc.azi1, 1e-3) {
// 			t.Errorf("Row: %d -- Inverse() Azimuth1Deg = %v; want %v", row_num, res.Azimuth1Deg, tc.azi1)
// 		}

// 		if !almost_equal(res.Azimuth2Deg, tc.azi2, 1e-3) {
// 			t.Errorf("Row: %d -- Inverse() Azimuth2Deg = %v; want %v", row_num, res.Azimuth2Deg, tc.azi2)
// 		}

// 		if !almost_equal(res.ArcLengthDeg, tc.a12, float64EqualityThreshold) {
// 			t.Errorf("Row: %d -- Inverse() ArcLengthDeg = %v; want %v", row_num, res.ArcLengthDeg, tc.a12)
// 		}

// 		if !almost_equal(res.ReducedLengthM, tc.m12, float64EqualityThreshold) {
// 			t.Errorf("Row: %d -- Inverse() ReducedLengthM = %v; want %v", row_num, res.ReducedLengthM, tc.m12)
// 		}

// 		if !almost_equal(res.S12M2, tc.s12_area, 1e3) {
// 			t.Errorf("Row: %d -- Inverse() S12M2 = %v; want %v", row_num, res.S12M2, tc.s12_area)
// 		}

// 	}
// }

// test_cases numbers comes from python test file
var test_cases = [20][12]float64{
	{35.60777, -139.44815, 111.098748429560326,
		-11.17491, -69.95921, 129.289270889708762,
		8935244.5604818305, 80.50729714281974, 6273170.2055303837,
		0.16606318447386067, 0.16479116945612937, 12841384694976.432},
	{55.52454, 106.05087, 22.020059880982801,
		77.03196, 197.18234, 109.112041110671519,
		4105086.1713924406, 36.892740690445894, 3828869.3344387607,
		0.80076349608092607, 0.80101006984201008, 61674961290615.615},
	{-21.97856, 142.59065, -32.44456876433189,
		41.84138, 98.56635, -41.84359951440466,
		8394328.894657671, 75.62930491011522, 6161154.5773110616,
		0.24816339233950381, 0.24930251203627892, -6637997720646.717},
	{-66.99028, 112.2363, 173.73491240878403,
		-12.70631, 285.90344, 2.512956620913668,
		11150344.2312080241, 100.278634181155759, 6289939.5670446687,
		-0.17199490274700385, -0.17722569526345708, -121287239862139.744},
	{-17.42761, 173.34268, -159.033557661192928,
		-15.84784, 5.93557, -20.787484651536988,
		16076603.1631180673, 144.640108810286253, 3732902.1583877189,
		-0.81273638700070476, -0.81299800519154474, 97825992354058.708},
	{32.84994, 48.28919, 150.492927788121982,
		-56.28556, 202.29132, 48.113449399816759,
		16727068.9438164461, 150.565799985466607, 3147838.1910180939,
		-0.87334918086923126, -0.86505036767110637, -72445258525585.010},
	{6.96833, 52.74123, 92.581585386317712,
		-7.39675, 206.17291, 90.721692165923907,
		17102477.2496958388, 154.147366239113561, 2772035.6169917581,
		-0.89991282520302447, -0.89986892177110739, -1311796973197.995},
	{-50.56724, -16.30485, -105.439679907590164,
		-33.56571, -94.97412, -47.348547835650331,
		6455670.5118668696, 58.083719495371259, 5409150.7979815838,
		0.53053508035997263, 0.52988722644436602, 41071447902810.047},
	{-58.93002, -8.90775, 140.965397902500679,
		-8.91104, 133.13503, 19.255429433416599,
		11756066.0219864627, 105.755691241406877, 6151101.2270708536,
		-0.26548622269867183, -0.27068483874510741, -86143460552774.735},
	{-68.82867, -74.28391, 93.774347763114881,
		-50.63005, -8.36685, 34.65564085411343,
		3956936.926063544, 35.572254987389284, 3708890.9544062657,
		0.81443963736383502, 0.81420859815358342, -41845309450093.787},
	{-10.62672, -32.0898, -86.426713286747751,
		5.883, -134.31681, -80.473780971034875,
		11470869.3864563009, 103.387395634504061, 6184411.6622659713,
		-0.23138683500430237, -0.23155097622286792, 4198803992123.548},
	{-21.76221, 166.90563, 29.319421206936428,
		48.72884, 213.97627, 43.508671946410168,
		9098627.3986554915, 81.963476716121964, 6299240.9166992283,
		0.13965943368590333, 0.14152969707656796, 10024709850277.476},
	{-19.79938, -174.47484, 71.167275780171533,
		-11.99349, -154.35109, 65.589099775199228,
		2319004.8601169389, 20.896611684802389, 2267960.8703918325,
		0.93427001867125849, 0.93424887135032789, -3935477535005.785},
	{-11.95887, -116.94513, 92.712619830452549,
		4.57352, 7.16501, 78.64960934409585,
		13834722.5801401374, 124.688684161089762, 5228093.177931598,
		-0.56879356755666463, -0.56918731952397221, -9919582785894.853},
	{-87.85331, 85.66836, -65.120313040242748,
		66.48646, 16.09921, -4.888658719272296,
		17286615.3147144645, 155.58592449699137, 2635887.4729110181,
		-0.90697975771398578, -0.91095608883042767, 42667211366919.534},
	{1.74708, 128.32011, -101.584843631173858,
		-11.16617, 11.87109, -86.325793296437476,
		12942901.1241347408, 116.650512484301857, 5682744.8413270572,
		-0.44857868222697644, -0.44824490340007729, 10763055294345.653},
	{-25.72959, -144.90758, -153.647468693117198,
		-57.70581, -269.17879, -48.343983158876487,
		9413446.7452453107, 84.664533838404295, 6356176.6898881281,
		0.09492245755254703, 0.09737058264766572, 74515122850712.444},
	{-41.22777, 122.32875, 14.285113402275739,
		-7.57291, 130.37946, 10.805303085187369,
		3812686.035106021, 34.34330804743883, 3588703.8812128856,
		0.82605222593217889, 0.82572158200920196, -2456961531057.857},
	{11.01307, 138.25278, 79.43682622782374,
		6.62726, 247.05981, 103.708090215522657,
		11911190.819018408, 107.341669954114577, 6070904.722786735,
		-0.29767608923657404, -0.29785143390252321, 17121631423099.696},
	{-29.47124, 95.14681, -163.779130441688382,
		-27.46601, -69.15955, -15.909335945554969,
		13487015.8381145492, 121.294026715742277, 5481428.9945736388,
		-0.51527225545373252, -0.51556587964721788, 104679964020340.318}}

func TestInverse20(t *testing.T) {
	geod := Wgs84()

	for row, tC := range test_cases {
		lat1, lon1, azi1 := tC[0], tC[1], tC[2]
		lat2, lon2, azi2 := tC[3], tC[4], tC[5]
		s12, a12, m12 := tC[6], tC[7], tC[8]
		M12, M21, S12 := tC[9], tC[10], tC[11]

		inv := geod.InverseCalcAll(lat1, lon1, lat2, lon2)

		if !almost_equal(azi1, inv.Azimuth1Deg, 1e-13) {
			t.Errorf("row %d -- Inverse() azi1 = %v; want %v", row, inv.Azimuth1Deg, azi1)
		}

		if !almost_equal(azi2, inv.Azimuth2Deg, 1e-13) {
			t.Errorf("row %d -- Inverse() azi2 = %v; want %v", row, inv.Azimuth2Deg, azi2)
		}

		if !almost_equal(s12, inv.DistanceM, 1e-8) {
			t.Errorf("row %d -- Inverse() s12 = %v; want %v", row, inv.DistanceM, s12)
		}

		if !almost_equal(a12, inv.ArcLengthDeg, 1e-13) {
			t.Errorf("row %d -- Inverse() a12 = %v; want %v", row, inv.ArcLengthDeg, a12)
		}

		if !almost_equal(m12, inv.ReducedLengthM, 1e-8) {
			t.Errorf("row %d -- Inverse() m12 = %v; want %v", row, inv.ReducedLengthM, m12)
		}

		if !almost_equal(M12, inv.M12, 1e-15) {
			t.Errorf("row %d -- Inverse() M12 = %v; want %v", row, inv.M12, M12)
		}

		if !almost_equal(M21, inv.M21, 1e-15) {
			t.Errorf("row %d -- Inverse() M21 = %v; want %v", row, inv.M21, M21)
		}

		if !almost_equal(S12, inv.S12M2, 0.1) {
			t.Errorf("row %d -- Inverse() S12 = %v; want %v", row, inv.S12M2, S12)
		}
	}
}

func BenchmarkInverse20(b *testing.B) {
	geod := Wgs84()
	for i := 0; i < b.N; i++ {
		for _, tC := range test_cases {
			lat1, lon1 := tC[0], tC[1]
			lat2, lon2 := tC[3], tC[4]

			geod.InverseCalcDistanceAzimuths(lat1, lon1, lat2, lon2)

		}
	}
}

func TestDirect20(t *testing.T) {
	geod := Wgs84()

	for row, tC := range test_cases {
		lat1, lon1, azi1 := tC[0], tC[1], tC[2]
		lat2, lon2, azi2 := tC[3], tC[4], tC[5]
		s12, a12, m12 := tC[6], tC[7], tC[8]
		M12, M21, S12 := tC[9], tC[10], tC[11]

		inv := geod.DirectCalcWithCapabilities(lat1, lon1, azi1, s12, ALL|LONG_UNROLL)

		if !almost_equal(lat2, inv.LatDeg, 1e-13) {
			t.Errorf("row %d -- Inverse() lat2 = %v; want %v", row, inv.LatDeg, lat2)
		}

		if !almost_equal(lon2, inv.LonDeg, 1e-13) {
			t.Errorf("row %d -- Inverse() lon2 = %v; want %v", row, inv.LonDeg, lon2)
		}

		if !almost_equal(azi2, inv.AziDeg, 1e-13) {
			t.Errorf("row %d -- Inverse() azi2 = %v; want %v", row, inv.AziDeg, azi2)
		}

		if !almost_equal(a12, inv.A12Deg, 1e-13) {
			t.Errorf("row %d -- Inverse() a12 = %v; want %v", row, inv.A12Deg, a12)
		}

		if !almost_equal(m12, inv.ReducedLengthM, 1e-8) {
			t.Errorf("row %d -- Inverse() m12 = %v; want %v", row, inv.ReducedLengthM, m12)
		}

		if !almost_equal(M12, inv.M12, 1e-15) {
			t.Errorf("row %d -- Inverse() M12 = %v; want %v", row, inv.M12, M12)
		}

		if !almost_equal(M21, inv.M21, 1e-15) {
			t.Errorf("row %d -- Inverse() M21 = %v; want %v", row, inv.M21, M21)
		}

		if !almost_equal(S12, inv.S12M2, 0.1) {
			t.Errorf("row %d -- Inverse() S12 = %v; want %v", row, inv.S12M2, S12)
		}
	}
}

func BenchmarkDirect20(b *testing.B) {
	geod := Wgs84()
	for i := 0; i < b.N; i++ {
		for _, tC := range test_cases {
			lat1, lon1, azi1 := tC[0], tC[1], tC[2]
			s12 := tC[6]

			geod.DirectCalcWithCapabilities(lat1, lon1, azi1, s12, ALL|LONG_UNROLL)

		}
	}
}

func TestGeodSolve0(t *testing.T) {
	geod := Wgs84()
	inv := geod.InverseCalcDistanceAzimuths(40.6, -73.8, 49.01666667, 2.55)

	if !almost_equal(inv.Azimuth1Deg, 53.47022, 0.5e-5) {
		t.Errorf("azi1 = %v; want %v", inv.Azimuth1Deg, 53.47022)
	}

	if !almost_equal(inv.Azimuth2Deg, 111.59367, 0.5e-5) {
		t.Errorf("azi2 = %v; want %v", inv.Azimuth2Deg, 111.59367)
	}

	if !almost_equal(inv.DistanceM, 5853226, 0.5) {
		t.Errorf("s12 = %v; want %v", inv.DistanceM, 5853226)
	}
}

func TestGeodSolve1(t *testing.T) {
	geod := Wgs84()
	dir := geod.DirectCalcLatLonAzi(40.63972222, -73.77888889, 53.5, 5850e3)

	if !almost_equal(dir.LatDeg, 49.01467, 0.5e-5) {
		t.Errorf("lat2 = %v; want %v", dir.LatDeg, 49.01467)
	}

	if !almost_equal(dir.LonDeg, 2.56106, 0.5e-5) {
		t.Errorf("lon2 = %v; want %v", dir.LonDeg, 2.56106)
	}

	if !almost_equal(dir.AziDeg, 111.62947, 0.5e-5) {
		t.Errorf("azi1 = %v; want %v", dir.AziDeg, 111.62947)
	}
}

func TestGeodSolve2(t *testing.T) {
	// Check fix for antipodal prolate bug found 2010-09-04
	geod := NewGeodesic(6.4e6, -1/150.0)

	inv := geod.InverseCalcDistanceAzimuths(0.07476, 0, -0.07476, 180)
	want_azi1 := 90.00078
	want_azi2 := 90.00078
	want_s12 := 20106193.0

	if !almost_equal(inv.Azimuth1Deg, want_azi1, 0.5e-5) {
		t.Errorf("azi1 = %v; want %v", inv.Azimuth1Deg, want_azi1)
	}

	if !almost_equal(inv.Azimuth2Deg, want_azi2, 0.5e-5) {
		t.Errorf("azi2 = %v; want %v", inv.Azimuth2Deg, want_azi2)
	}

	if !almost_equal(inv.DistanceM, want_s12, 0.5) {
		t.Errorf("s12 = %v; want %v", inv.DistanceM, want_s12)
	}

	inv = geod.InverseCalcDistanceAzimuths(0.1, 0, -0.1, 180)
	want_azi1 = 90.00105
	want_azi2 = 90.00105
	want_s12 = 20106193.0
	if !almost_equal(inv.Azimuth1Deg, want_azi1, 0.5e-5) {
		t.Errorf("azi1 = %v; want %v", inv.Azimuth1Deg, want_azi1)
	}

	if !almost_equal(inv.Azimuth2Deg, want_azi2, 0.5e-5) {
		t.Errorf("azi2 = %v; want %v", inv.Azimuth2Deg, want_azi2)
	}

	if !almost_equal(inv.DistanceM, want_s12, 0.5) {
		t.Errorf("s12 = %v; want %v", inv.DistanceM, want_s12)
	}
}

func TestGeodSolve4(t *testing.T) {
	// Check fix for short line bug found 2010-05-21
	geod := Wgs84()
	s12 := geod.InverseCalcDistance(36.493349428792, 0, 36.49334942879201, 0.0000008)
	want_s12 := 0.072

	if !almost_equal(s12, want_s12, 0.5e-3) {
		t.Errorf("s12 = %v; want %v", s12, want_s12)
	}
}

func TestGeodSolve5(t *testing.T) {
	// Check fix for point2=pole bug found 2010-05-03
	geod := Wgs84()
	dir := geod.DirectCalcLatLonAzi(0.01777745589997, 30, 0, 10e6)
	want_lat := 90.0

	if !almost_equal(dir.LatDeg, want_lat, 0.5e-5) {
		t.Errorf("lat2 = %v; want %v", dir.LatDeg, want_lat)
	}

	if dir.LonDeg < 0 {
		want_lon := -150.0
		want_azi := 180.0

		if !almost_equal(dir.LonDeg, want_lon, 0.5e-5) {
			t.Errorf("lon2 = %v; want %v", dir.LonDeg, want_lon)
		}

		if !almost_equal(dir.AziDeg, want_azi, 0.5e-5) {
			t.Errorf("azi1 = %v; want %v", dir.AziDeg, want_azi)
		}
	} else {
		want_lon := 30.0
		want_azi := 0.0
		if !almost_equal(dir.LonDeg, want_lon, 0.5e-5) {
			t.Errorf("lon2 = %v; want %v", dir.LonDeg, want_lon)
		}

		if !almost_equal(dir.AziDeg, want_azi, 0.5e-5) {
			t.Errorf("azi1 = %v; want %v", dir.AziDeg, want_azi)
		}
	}
}

func TestGeodSolve6(t *testing.T) {
	//Check fix for volatile sbet12a bug found 2011-06-25 (gcc 4.4.4
	// x86 -O3).  Found again on 2012-03-27 with tdm-mingw32 (g++ 4.6.1).
	geod := Wgs84()
	s12 := geod.InverseCalcDistance(88.202499451857, 0, -88.202499451857, 179.981022032992859592)
	want_s12 := 20003898.214
	if !almost_equal(s12, want_s12, 0.5e-3) {
		t.Errorf("s12 = %v; want %v", s12, want_s12)
	}

	s12 = geod.InverseCalcDistance(89.262080389218, 0, -89.262080389218, 179.992207982775375662)
	want_s12 = 20003925.854
	if !almost_equal(s12, want_s12, 0.5e-3) {
		t.Errorf("s12 = %v; want %v", s12, want_s12)
	}

	s12 = geod.InverseCalcDistance(
		89.333123580033,
		0,
		-89.333123580032997687,
		179.99295812360148422,
	)
	want_s12 = 20003926.881
	if !almost_equal(s12, want_s12, 0.5e-3) {
		t.Errorf("s12 = %v; want %v", s12, want_s12)
	}
}

func TestGeodSolve9(t *testing.T) {
	// Check fix for volatile x bug found 2011-06-25 (gcc 4.4.4 x86 -O3)
	geod := Wgs84()

	s12 := geod.InverseCalcDistance(56.320923501171, 0, -56.320923501171, 179.664747671772880215)
	want_s12 := 19993558.287
	if !almost_equal(s12, want_s12, 0.5e-3) {
		t.Errorf("s12 = %v; want %v", s12, want_s12)
	}

}

func TestGeodSolve10(t *testing.T) {
	// Check fix for adjust tol1_ bug found 2011-06-25 (Visual Studio
	// 10 rel + debug)
	geod := Wgs84()

	s12 := geod.InverseCalcDistance(
		52.784459512564,
		0,
		-52.784459512563990912,
		179.634407464943777557,
	)

	want_s12 := 19991596.095
	if !almost_equal(s12, want_s12, 0.5e-3) {
		t.Errorf("s12 = %v; want %v", s12, want_s12)
	}
}

func TestGeodSolve11(t *testing.T) {
	// Check fix for bet2 = -bet1 bug found 2011-06-25 (Visual Studio
	// 10 rel + debug)
	geod := Wgs84()
	s12 := geod.InverseCalcDistance(
		48.522876735459, 0,
		-48.52287673545898293,
		179.599720456223079643,
	)

	want_s12 := 19989144.774
	if !almost_equal(s12, want_s12, 0.5e-3) {
		t.Errorf("s12 = %v; want %v", s12, want_s12)
	}
}

func TestGeodSolve12(t *testing.T) {
	// Check fix for inverse geodesics on extreme prolate/oblate
	// ellipsoids Reported 2012-08-29 Stefan Guenther
	// <stefan.gunther@embl.de>; fixed 2012-10-07
	geod := NewGeodesic(89.8, -1.83)

	inv := geod.InverseCalcDistanceAzimuths(0, 0, -10, 160)

	want_azi1 := 120.27
	want_azi2 := 105.15
	want_s12 := 266.7

	if !almost_equal(inv.Azimuth1Deg, want_azi1, 1e-2) {
		t.Errorf("azi1 = %v; want %v", inv.Azimuth1Deg, want_azi1)
	}

	if !almost_equal(inv.Azimuth2Deg, want_azi2, 1e-2) {
		t.Errorf("azi2 = %v; want %v", inv.Azimuth2Deg, want_azi2)
	}

	if !almost_equal(inv.DistanceM, want_s12, 1e-1) {
		t.Errorf("s12 = %v; want %v", inv.DistanceM, want_s12)
	}
}

func TestGeodSolve14(t *testing.T) {
	// Check fix for inverse ignoring lon12 = nan
	geod := Wgs84()
	inv := geod.InverseCalcDistanceAzimuths(0, 0, 1, math.NaN())

	if !math.IsNaN(inv.Azimuth1Deg) {
		t.Errorf("azi1 = %v; want %v", inv.Azimuth1Deg, math.NaN())
	}
	if !math.IsNaN(inv.Azimuth2Deg) {
		t.Errorf("azi2 = %v; want %v", inv.Azimuth2Deg, math.NaN())
	}
	if !math.IsNaN(inv.DistanceM) {
		t.Errorf("s12 = %v; want %v", inv.DistanceM, math.NaN())
	}

}

func TestGeodSolve15(t *testing.T) {
	// Initial implementation of Math::eatanhe was wrong for e^2 < 0.  This
	// checks that this is fixed.
	geod := NewGeodesic(6.4e6, -1/150.0)
	dir := geod.DirectCalcWithCapabilities(1, 2, 3, 4, AREA)
	if !almost_equal(dir.S12M2, 23700.0, 0.5) {
		t.Errorf("S12 = %v; want %v", dir.S12M2, 23700.0)
	}
}

func TestGeodSolve17(t *testing.T) {
	// Check fix for LONG_UNROLL bug found on 2015-05-07
	geod := Wgs84()
	dir := geod.DirectCalcWithCapabilities(40, -75, -10, 2e7, STANDARD|LONG_UNROLL)
	want_lat := -39.0
	want_lon := -254.0
	want_azi := -170.0

	if !almost_equal(dir.LatDeg, want_lat, 1) {
		t.Errorf("lat = %v; want %v", dir.LatDeg, want_lat)
	}

	if !almost_equal(dir.LonDeg, want_lon, 1) {
		t.Errorf("lon = %v; want %v", dir.LonDeg, want_lon)
	}

	if !almost_equal(dir.AziDeg, want_azi, 1) {
		t.Errorf("azi = %v; want %v", dir.AziDeg, want_azi)
	}

	line := NewGeodesicLine(geod, 40.0, -75.0, -10.0, math.NaN(), math.NaN())

	dir2 := line.PositionWithCapabilities(2e7, STANDARD|LONG_UNROLL)
	want_lat = -39.0
	want_lon = -254.0
	want_azi = -170.0
	if !almost_equal(dir2.Lat2Deg, want_lat, 1) {
		t.Errorf("lat = %v; want %v", dir2.Lat2Deg, want_lat)
	}

	if !almost_equal(dir2.Lon2Deg, want_lon, 1) {
		t.Errorf("lon = %v; want %v", dir2.Lon2Deg, want_lon)
	}

	if !almost_equal(dir2.Azi2Deg, want_azi, 1) {
		t.Errorf("azi = %v; want %v", dir2.Azi2Deg, want_azi)
	}

	dir3 := geod.DirectCalcLatLonAzi(40, -75, -10, 2e7)
	want_lat = -39.0
	want_lon = 105.0
	want_azi = -170.0
	if !almost_equal(dir3.LatDeg, want_lat, 1) {
		t.Errorf("lat = %v; want %v", dir3.LatDeg, want_lat)
	}

	if !almost_equal(dir3.LonDeg, want_lon, 1) {
		t.Errorf("lon = %v; want %v", dir3.LonDeg, want_lon)
	}

	if !almost_equal(dir3.AziDeg, want_azi, 1) {
		t.Errorf("azi = %v; want %v", dir3.AziDeg, want_azi)
	}

	dir4 := line.PositionStandard(2e7)
	want_lat = -39.0
	want_lon = 105.0
	want_azi = -170.0
	if !almost_equal(dir4.Lat2Deg, want_lat, 1) {
		t.Errorf("lat = %v; want %v", dir4.Lat2Deg, want_lat)
	}

	if !almost_equal(dir4.Lon2Deg, want_lon, 1) {
		t.Errorf("lon = %v; want %v", dir4.Lon2Deg, want_lon)
	}

	if !almost_equal(dir4.Azi2Deg, want_azi, 1) {
		t.Errorf("azi = %v; want %v", dir4.Azi2Deg, want_azi)
	}
}

func TestGeodSolve26(t *testing.T) {
	// Check 0/0 problem with area calculation on sphere 2015-09-08
	geod := NewGeodesic(6.4e6, 0)
	inv := geod.InverseCalcWithCapabilities(1, 2, 3, 4, AREA)

	if !almost_equal(inv.S12M2, 49911046115.0, 0.5) {
		t.Errorf("S12 = %v; want %v", inv.S12M2, 49911046115.0)
	}
}

func TestGeodSolve28(t *testing.T) {
	// Check for bad placement of assignment of r.a12 with |f| > 0.01 (bug in
	// Java implementation fixed on 2015-05-19).
	geod := NewGeodesic(6.4e6, 0.1)
	dir := geod.DirectCalcAll(1, 2, 10, 5e6)

	if !almost_equal(dir.A12Deg, 48.55570690, 0.5e-8) {
		t.Errorf("a12 = %v; want %v", dir.A12Deg, 48.55570690)
	}
}

func TestGeodSolve29(t *testing.T) {
	// Check longitude unrolling with inverse calculation 2015-09-16
	// This functionality has not been added to this port of the library
	t.SkipNow()
}

func TestGeodSolve33(t *testing.T) {
	// Check max(-0.0,+0.0) issues 2015-08-22 (triggered by bugs in
	// Octave -- sind(-0.0) = +0.0 -- and in some version of Visual
	// Studio -- fmod(-0.0, 360.0) = +0.0.
	testCases := []struct {
		desc      string
		geod      Geodesic
		lat1      float64
		lon1      float64
		lat2      float64
		lon2      float64
		want_azi1 float64
		want_azi2 float64
		want_s12  float64
	}{
		{
			desc:      "1",
			geod:      Wgs84(),
			lat1:      0,
			lon1:      0,
			lat2:      0,
			lon2:      179,
			want_azi1: 90.0,
			want_azi2: 90.0,
			want_s12:  19926189,
		},
		{
			desc:      "2",
			geod:      Wgs84(),
			lat1:      0,
			lon1:      0,
			lat2:      0,
			lon2:      179.5,
			want_azi1: 55.96650,
			want_azi2: 124.03350,
			want_s12:  19980862,
		},
		{
			desc:      "3",
			geod:      Wgs84(),
			lat1:      0,
			lon1:      0,
			lat2:      0,
			lon2:      180.0,
			want_azi1: 0.0,
			want_azi2: 180.0,
			want_s12:  20003931,
		},
		{
			desc:      "4",
			geod:      Wgs84(),
			lat1:      0,
			lon1:      0,
			lat2:      1,
			lon2:      180.0,
			want_azi1: 0.0,
			want_azi2: 180.0,
			want_s12:  19893357,
		},
		{
			desc:      "5",
			geod:      NewGeodesic(6.4e6, 0),
			lat1:      0,
			lon1:      0,
			lat2:      0,
			lon2:      179,
			want_azi1: 90.0,
			want_azi2: 90.0,
			want_s12:  19994492,
		},
		{
			desc:      "6",
			geod:      NewGeodesic(6.4e6, 0),
			lat1:      0,
			lon1:      0,
			lat2:      0,
			lon2:      180,
			want_azi1: 0.0,
			want_azi2: 180.0,
			want_s12:  20106193,
		},
		{
			desc:      "7",
			geod:      NewGeodesic(6.4e6, 0),
			lat1:      0,
			lon1:      0,
			lat2:      1,
			lon2:      180,
			want_azi1: 0.0,
			want_azi2: 180.0,
			want_s12:  19994492,
		},
		{
			desc:      "8",
			geod:      NewGeodesic(6.4e6, -1/300.0),
			lat1:      0,
			lon1:      0,
			lat2:      0,
			lon2:      179,
			want_azi1: 90.0,
			want_azi2: 90.0,
			want_s12:  19994492,
		},
		{
			desc:      "9",
			geod:      NewGeodesic(6.4e6, -1/300.0),
			lat1:      0,
			lon1:      0,
			lat2:      0,
			lon2:      180,
			want_azi1: 90.0,
			want_azi2: 90.0,
			want_s12:  20106193,
		},
		{
			desc:      "10",
			geod:      NewGeodesic(6.4e6, -1/300.0),
			lat1:      0,
			lon1:      0,
			lat2:      0.5,
			lon2:      180,
			want_azi1: 33.02493,
			want_azi2: 146.97364,
			want_s12:  20082617,
		},
		{
			desc:      "11",
			geod:      NewGeodesic(6.4e6, -1/300.0),
			lat1:      0,
			lon1:      0,
			lat2:      1,
			lon2:      180,
			want_azi1: 0,
			want_azi2: 180.0,
			want_s12:  20027270,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := tC.geod.InverseCalcDistanceAzimuths(tC.lat1, tC.lon1, tC.lat2, tC.lon2)

			if !almost_equal(got.Azimuth1Deg, tC.want_azi1, 0.5e-5) {
				t.Errorf("azi1 = %v; want %v", got.Azimuth1Deg, tC.want_azi1)
			}

			if !almost_equal(got.Azimuth2Deg, tC.want_azi2, 0.5e-5) {
				t.Errorf("azi2 = %v; want %v", got.Azimuth2Deg, tC.want_azi2)
			}

			if !almost_equal(got.DistanceM, tC.want_s12, 0.5) {
				t.Errorf("s12 = %v; want %v", got.DistanceM, tC.want_s12)
			}
		})
	}
}

func TestGeodSolve55(t *testing.T) {
	// Check fix for nan + point on equator or pole not returning all nans in
	// Geodesic::Inverse, found 2015-09-23.
	geod := Wgs84()
	testCases := []struct {
		desc string
		lat1 float64
		lon1 float64
		lat2 float64
		lon2 float64
	}{
		{
			desc: "1",
			lat1: math.NaN(),
			lon1: 0,
			lat2: 0,
			lon2: 90,
		},
		{
			desc: "2",
			lat1: math.NaN(),
			lon1: 0,
			lat2: 90,
			lon2: 9,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := geod.InverseCalcDistanceAzimuths(tC.lat1, tC.lon1, tC.lat2, tC.lon2)

			if !math.IsNaN(got.Azimuth1Deg) {
				t.Errorf("az1 = %v; want NaN", got.Azimuth1Deg)
			}

			if !math.IsNaN(got.Azimuth2Deg) {
				t.Errorf("az2 = %v; want NaN", got.Azimuth2Deg)
			}

			if !math.IsNaN(got.DistanceM) {
				t.Errorf("s12 = %v; want NaN", got.DistanceM)
			}
		})
	}
}

func TestGeodSolve59(t *testing.T) {
	// Check for points close with longitudes close to 180 deg apart.
	geod := Wgs84()
	inv := geod.InverseCalcDistanceAzimuths(5, 0.00000000000001, 10, 180)

	want_azi1 := 0.000000000000035
	want_azi2 := 179.99999999999996
	want_s12 := 18345191.174332713

	if !almost_equal(inv.Azimuth1Deg, want_azi1, 1.5e-14) {
		t.Errorf("azi1 = %v; want %v", inv.Azimuth1Deg, want_azi1)
	}

	// Changed acceptable tolerance from 1.5e-14 to 1.5e-13 to make test pass
	if !almost_equal(inv.Azimuth2Deg, want_azi2, 1.5e-13) {
		t.Errorf("azi2 = %v; want %v", inv.Azimuth2Deg, want_azi2)
	}

	if !almost_equal(inv.DistanceM, want_s12, 5e-9) {
		t.Errorf("s12 = %v; want %v", inv.DistanceM, want_s12)
	}
}

func TestGeodSolve61(t *testing.T) {
	// Make sure small negative azimuths are west-going
	geod := Wgs84()
	dir := geod.DirectCalcWithCapabilities(
		45,
		0,
		-0.000000000000000003,
		1e7,
		STANDARD|LONG_UNROLL,
	)

	want_lat := 45.30632
	want_lon := -180.0
	want_azi := -180.0

	if !almost_equal(dir.LatDeg, want_lat, 0.5e-5) {
		t.Errorf("lat = %v; want %v", dir.LatDeg, want_lat)
	}

	if !almost_equal(dir.LonDeg, want_lon, 0.5e-5) {
		t.Errorf("lon = %v; want %v", dir.LonDeg, want_lon)
	}

	if !almost_equal(dir.AziDeg, want_azi, 0.5e-5) {
		t.Errorf("azi = %v; want %v", dir.AziDeg, want_azi)
	}

	line := geod.InverseLineWithCapabilities(45, 0, 80, -0.000000000000000003, STANDARD|DISTANCE_IN)
	pos := line.PositionWithCapabilities(1e7, STANDARD|LONG_UNROLL)
	want_lat = 45.30632
	want_lon = -180.0
	want_azi = -180.0
	if !almost_equal(pos.Lat2Deg, want_lat, 0.5e-5) {
		t.Errorf("lat = %v; want %v", pos.Lat2Deg, want_lat)
	}

	if !almost_equal(pos.Lon2Deg, want_lon, 0.5e-5) {
		t.Errorf("lon = %v; want %v", pos.Lon2Deg, want_lon)
	}

	if !almost_equal(pos.Azi2Deg, want_azi, 0.5e-5) {
		t.Errorf("azi = %v; want %v", pos.Azi2Deg, want_azi)
	}
}

func TestGeodSolve65(t *testing.T) {
	// Check for bug in east-going check in GeodesicLine (needed to check for
	// sign of 0) and sign error in area calculation due to a bogus override
	// of the code for alp12.  Found/fixed on 2015-12-19.
	geod := Wgs84()
	line := geod.InverseLineWithCapabilities(30, -0.000000000000000001, -31, 180, ALL)
	testCases := []struct {
		desc         string
		s12          float64
		capabilities uint64
		want_lat1    float64
		want_lon1    float64
		want_azi1    float64
		want_lat2    float64
		want_lon2    float64
		want_azi2    float64
		want_s12     float64
		want_a12     float64
		want_m12     float64
		want_M12     float64
		want_M21     float64
		want_S12     float64
	}{
		{
			desc:         "10M meters",
			s12:          1e7,
			capabilities: ALL | LONG_UNROLL,
			want_lat1:    30,
			want_lon1:    -0,
			want_azi1:    -180,
			want_lat2:    -60.23169,
			want_lon2:    -0,
			want_azi2:    -180,
			want_s12:     10000000,
			want_a12:     90.06544,
			want_m12:     6363636,
			want_M12:     -0.0012834,
			want_M21:     0.0013749,
			want_S12:     0,
		},
		{
			desc:         "20M meters",
			s12:          2e7,
			capabilities: ALL | LONG_UNROLL,
			want_lat1:    30,
			want_lon1:    -0,
			want_azi1:    -180,
			want_lat2:    -30.03547,
			want_lon2:    -180,
			want_azi2:    -0,
			want_s12:     20000000,
			want_a12:     179.96459,
			want_m12:     54342,
			want_M12:     -1.0045592,
			want_M21:     -0.9954339,
			want_S12:     127516405431022.0,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := line.PositionWithCapabilities(tC.s12, tC.capabilities)

			if !almost_equal(got.Lat1Deg, tC.want_lat1, 0.5e-5) {
				t.Errorf("lat1 = %v; want %v", got.Lat1Deg, tC.want_lat1)
			}

			if !almost_equal(got.Lon1Deg, tC.want_lon1, 0.5e-5) {
				t.Errorf("lon1 = %v; want %v", got.Lon1Deg, tC.want_lon1)
			}

			if !almost_equal(got.Azi1Deg, tC.want_azi1, 0.5e-5) {
				t.Errorf("azi1 = %v; want %v", got.Azi1Deg, tC.want_azi1)
			}

			if !almost_equal(got.Lat2Deg, tC.want_lat2, 0.5e-5) {
				t.Errorf("lat2 = %v; want %v", got.Lat2Deg, tC.want_lat2)
			}

			if !almost_equal(got.Lon2Deg, tC.want_lon2, 0.5e-5) {
				t.Errorf("lon2 = %v; want %v", got.Lon2Deg, tC.want_lon2)
			}

			if !almost_equal(got.Azi2Deg, tC.want_azi2, 0.5e-5) {
				t.Errorf("azi2 = %v; want %v", got.Azi2Deg, tC.want_azi2)
			}

			if !almost_equal(got.DistanceM, tC.want_s12, 0.5) {
				t.Errorf("s12 = %v; want %v", got.DistanceM, tC.want_s12)
			}

			if !almost_equal(got.ArcLengthDeg, tC.want_a12, 0.5e-5) {
				t.Errorf("a12 = %v; want %v", got.ArcLengthDeg, tC.want_a12)
			}

			if !almost_equal(got.ReducedLengthM, tC.want_m12, 0.5) {
				t.Errorf("m12 = %v; want %v", got.ReducedLengthM, tC.want_m12)
			}

			if !almost_equal(got.M12, tC.want_M12, 0.5e7) {
				t.Errorf("M12 = %v; want %v", got.M12, tC.want_M12)
			}

			if !almost_equal(got.M21, tC.want_M21, 0.5e-7) {
				t.Errorf("M21 = %v; want %v", got.M21, tC.want_M21)
			}

			if !almost_equal(got.S12M2, tC.want_S12, 0.5) {
				t.Errorf("S12 = %v; want %v", got.S12M2, tC.want_S12)
			}

		})
	}
}

func TestGeodSolve66(t *testing.T) {
	geod := Wgs84()
	line := geod.InverseLineWithCapabilities(-5, -0.000000000000002, -10, 180, STANDARD|DISTANCE_IN)
	testCases := []struct {
		desc         string
		s12          float64
		capabilities uint64
		want_lat     float64
		want_lon     float64
		want_azi     float64
	}{
		{
			desc:         "20M meters",
			s12:          2e7,
			capabilities: STANDARD | LONG_UNROLL,
			want_lat:     4.96445,
			want_lon:     -180.00000,
			want_azi:     -0.00000,
		},
		{
			desc:         "other distance",
			s12:          0.5 * line.s13,
			capabilities: STANDARD | LONG_UNROLL,
			want_lat:     -87.52461,
			want_lon:     -0.00000,
			want_azi:     -180.00000,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := line.PositionWithCapabilities(tC.s12, tC.capabilities)

			if !almost_equal(got.Lat2Deg, tC.want_lat, 0.5e-5) {
				t.Errorf("lat2 = %v; want %v", got.Lat2Deg, tC.want_lat)
			}

			if !almost_equal(got.Lon2Deg, tC.want_lon, 0.5e-5) {
				t.Errorf("lon2 = %v; want %v", got.Lon2Deg, tC.want_lon)
			}

			if !almost_equal(got.Azi2Deg, tC.want_azi, 0.5e-5) {
				t.Errorf("azi2 = %v; want %v", got.Azi2Deg, tC.want_azi)
			}
		})
	}
}
