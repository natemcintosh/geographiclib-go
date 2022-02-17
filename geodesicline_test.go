package geographiclibgo

import (
	"math"
	"testing"
)

func TestNewGeodesicLine(t *testing.T) {
	geod := Wgs84()
	gl := NewGeodesicLine(geod, 0.0, 0.0, 0.0, math.NaN(), math.NaN())

	if !f64_equals(gl.a, 6378137.0) {
		t.Errorf("a = %v; want %v", gl.a, 6378137.0)
	}

	if !f64_equals(gl.f, 0.0033528106647474805) {
		t.Errorf("f = %v; want %v", gl.f, 0.0033528106647474805)
	}
	if !f64_equals(gl._b, 6356752.314245179) {
		t.Errorf("_b = %v; want %v", gl._b, 6356752.314245179)
	}
	if !f64_equals(gl._c2, 40589732499314.76) {
		t.Errorf("_c2 = %v; want %v", gl._c2, 40589732499314.76)
	}
	if !f64_equals(gl._f1, 0.9966471893352525) {
		t.Errorf("_f1 = %v; want %v", gl._f1, 0.9966471893352525)
	}
	if gl.caps != uint64(36747) {
		t.Errorf("caps = %v; want %v", gl.caps, 36747)
	}
	if !f64_equals(gl.lat1, 0.0) {
		t.Errorf("lat1 = %v; want %v", gl.lat1, 0.0)
	}
	if !f64_equals(gl.lon1, 0.0) {
		t.Errorf("lon1 = %v; want %v", gl.lon1, 0.0)
	}
	if !f64_equals(gl.azi1, 0.0) {
		t.Errorf("azi1 = %v; want %v", gl.azi1, 0.0)
	}
	if !f64_equals(gl.salp1, 0.0) {
		t.Errorf("salp1 = %v; want %v", gl.salp1, 0.0)
	}
	if !f64_equals(gl.calp1, 1.0) {
		t.Errorf("calp1 = %v; want %v", gl.calp1, 1.0)
	}
	if !f64_equals(gl._dn1, 1.0) {
		t.Errorf("_dn1 = %v; want %v", gl._dn1, 1.0)
	}
	if !f64_equals(gl._salp0, 0.0) {
		t.Errorf("_salp0 = %v; want %v", gl._salp0, 0.0)
	}
	if !f64_equals(gl._calp0, 1.0) {
		t.Errorf("_calp0 = %v; want %v", gl._calp0, 1.0)
	}
	if !f64_equals(gl._ssig1, 0.0) {
		t.Errorf("_ssig1 = %v; want %v", gl._ssig1, 0.0)
	}
	if !f64_equals(gl._somg1, 0.0) {
		t.Errorf("_somg1 = %v; want %v", gl._somg1, 0.0)
	}
	if !f64_equals(gl._csig1, 1.0) {
		t.Errorf("_csig1 = %v; want %v", gl._csig1, 1.0)
	}
	if !f64_equals(gl._comg1, 1.0) {
		t.Errorf("_comg1 = %v; want %v", gl._comg1, 1.0)
	}
	if !f64_equals(gl._k2, geod.ep2) {
		t.Errorf("_k2 = %v; want %v", gl._k2, geod.ep2)
	}
	if !f64_equals(gl.s13, math.NaN()) {
		t.Errorf("s13 = %v; want %v", gl.s13, math.NaN())
	}
	if !f64_equals(gl.a13, math.NaN()) {
		t.Errorf("a13 = %v; want %v", gl.a13, math.NaN())
	}
}

func BenchmarkNewGeodesicLine(b *testing.B) {
	geod := Wgs84()
	for i := 0; i < b.N; i++ {
		NewGeodesicLine(geod, 0.0, 0.0, 0.0, math.NaN(), math.NaN())
	}
}

func Test_gen_position(t *testing.T) {
	geod := Wgs84()
	gl := NewGeodesicLine(geod, 0.0, 0.0, 10.0, math.NaN(), math.NaN())
	a12, lat2, lon2, azi2, s12, m12, M12, M21, S12 := gl._gen_position(false, 150.0, 3979)

	if !f64_equals(a12, 0.0013520059461334633) {
		t.Errorf("a12 = %v; want %v", a12, 0.0013520059461334633)
	}

	if !f64_equals(lat2, 0.0013359451088740494) {
		t.Errorf("lat2 = %v; want %v", lat2, 0.0013359451088740494)
	}
	if !f64_equals(lon2, 0.00023398621812867812) {
		t.Errorf("lon2 = %v; want %v", lon2, 0.00023398621812867812)
	}
	if !f64_equals(azi2, 10.000000002727887) {
		t.Errorf("azi2 = %v; want %v", azi2, 10.000000002727887)
	}
	if !f64_equals(s12, 150.0) {
		t.Errorf("s12 = %v; want %v", s12, 150.0)
	}
	if !f64_equals(math.NaN(), m12) {
		t.Errorf("m12 = %v; want %v", m12, math.NaN())
	}
	if !f64_equals(math.NaN(), M12) {
		t.Errorf("M12 = %v; want %v", M12, math.NaN())
	}
	if !f64_equals(math.NaN(), M21) {
		t.Errorf("M21 = %v; want %v", M21, math.NaN())
	}
	if !f64_equals(math.NaN(), S12) {
		t.Errorf("S12 = %v; want %v", S12, math.NaN())
	}
}

func Benchmark_gen_position(b *testing.B) {
	geod := Wgs84()
	gl := NewGeodesicLine(geod, 0.0, 0.0, 10.0, math.NaN(), math.NaN())
	for i := 0; i < b.N; i++ {
		gl._gen_position(false, 150.0, 3979)
	}
}
