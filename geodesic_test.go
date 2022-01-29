package geographiclibgo

import "testing"

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
