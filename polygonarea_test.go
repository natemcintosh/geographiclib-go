package geographiclibgo

import (
	"math"
	"testing"
)

func planimeter(points [][2]float64) PolygonResult {
	polygon := NewPolygonArea(Wgs84(), false)
	for _, p := range points {
		polygon.AddPoint(p[0], p[1])
	}
	return polygon.Compute(false, true)
}

func polylength(points [][2]float64) PolygonResult {
	polyline := NewPolygonArea(Wgs84(), true)
	for _, p := range points {
		polyline.AddPoint(p[0], p[1])
	}
	return polyline.Compute(false, true)
}

func TestPlanimeter0(t *testing.T) {
	// Check fix for pole-encircling bug found 2011-03-16
	points := [][2]float64{{89, 0}, {89, 90}, {89, 180}, {89, 270}}
	got := planimeter(points)

	want_perimeter := 631819.8745
	want_area := 24952305678.0

	if !almost_equal(got.Perimeter, want_perimeter, 1e-4) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}

	if !almost_equal(got.Area, want_area, 1) {
		t.Errorf("area = %v; want %v", got.Area, want_area)
	}

	points = [][2]float64{{-89, 0}, {-89, 90}, {-89, 180}, {-89, 270}}
	got = planimeter(points)

	want_perimeter = 631819.8745
	want_area = -24952305678.0

	if !almost_equal(got.Perimeter, want_perimeter, 1e-4) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}

	if !almost_equal(got.Area, want_area, 1) {
		t.Errorf("area = %v; want %v", got.Area, want_area)
	}

	points = [][2]float64{{0, -1}, {-1, 0}, {0, 1}, {1, 0}}
	got = planimeter(points)

	want_perimeter = 627598.2731
	want_area = 24619419146.0

	if !almost_equal(got.Perimeter, want_perimeter, 1e-4) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}

	if !almost_equal(got.Area, want_area, 1) {
		t.Errorf("area = %v; want %v", got.Area, want_area)
	}

	points = [][2]float64{{90, 0}, {0, 0}, {0, 90}}
	got = planimeter(points)

	want_perimeter = 30022685.0
	want_area = 63758202715511.0

	if !almost_equal(got.Perimeter, want_perimeter, 1) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}

	if !almost_equal(got.Area, want_area, 1) {
		t.Errorf("area = %v; want %v", got.Area, want_area)
	}

	polyline_got := polylength(points)
	want_perimeter = 20020719.0

	if !almost_equal(polyline_got.Perimeter, want_perimeter, 1) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}

	if !math.IsNaN(polyline_got.Area) {
		t.Errorf("area = %v; want NaN", polyline_got.Area)
	}

}

func TestPlanimeter5(t *testing.T) {
	// Check fix for Planimeter pole crossing bug found 2011-06-24
	points := [][2]float64{{89, 0.1}, {89, 90.1}, {89, -179.9}}
	got := planimeter(points)
	want_perimeter := 539297.0
	want_area := 12476152838.5

	if !almost_equal(got.Perimeter, want_perimeter, 1) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}

	if !almost_equal(got.Area, want_area, 1) {
		t.Errorf("area = %v; want %v", got.Area, want_area)
	}

}

func TestPlanimeter6(t *testing.T) {
	// Check fix for Planimeter lon12 rounding bug found 2012-12-03
	testCases := []struct {
		desc           string
		points         [][2]float64
		want_area      float64
		want_perimeter float64
	}{
		{
			desc:           "1",
			points:         [][2]float64{{9, -0.00000000000001}, {9, 180}, {9, 0}},
			want_area:      0,
			want_perimeter: 36026861,
		},
		{
			desc:           "2",
			points:         [][2]float64{{9, 0.00000000000001}, {9, 0}, {9, 180}},
			want_area:      0,
			want_perimeter: 36026861,
		},
		{
			desc:           "3",
			points:         [][2]float64{{9, 0.00000000000001}, {9, 180}, {9, 0}},
			want_area:      0,
			want_perimeter: 36026861,
		},
		{
			desc:           "4",
			points:         [][2]float64{{9, -0.00000000000001}, {9, 0}, {9, 180}},
			want_area:      0,
			want_perimeter: 36026861,
		},
	}
	for _, tC := range testCases {
		t.Run(tC.desc, func(t *testing.T) {
			got := planimeter(tC.points)

			if !almost_equal(got.Perimeter, tC.want_perimeter, 1) {
				t.Errorf("perimeter = %v; want %v", got.Perimeter, tC.want_perimeter)
			}

			if !almost_equal(got.Area, tC.want_area, 1) {
				t.Errorf("area = %v; want %v", got.Area, tC.want_area)
			}
		})
	}
}

func TestPlanimeter12(t *testing.T) {
	// Area of arctic circle (not really -- adjunct to rhumb-area test)
	points := [][2]float64{{66.562222222, 0}, {66.562222222, 180}}
	got := planimeter(points)

	want_perimeter := 10465729.0
	want_area := 0.0

	if !almost_equal(got.Perimeter, want_perimeter, 1) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}

	if !almost_equal(got.Area, want_area, 1) {
		t.Errorf("area = %v; want %v", got.Area, want_area)
	}

}

func BenchmarkPlanimeter12(b *testing.B) {
	points := [][2]float64{{66.562222222, 0}, {66.562222222, 180}}
	for i := 0; i < b.N; i++ {
		planimeter(points)
	}
}

func TestPlanimeter13(t *testing.T) {
	// Check encircling pole twice
	points := [][2]float64{{89, -360}, {89, -240}, {89, -120}, {89, 0}, {89, 120}, {89, 240}}
	got := planimeter(points)

	want_perimeter := 1160741.0
	want_area := 32415230256.0

	if !almost_equal(got.Perimeter, want_perimeter, 1) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}

	if !almost_equal(got.Area, want_area, 1) {
		t.Errorf("area = %v; want %v", got.Area, want_area)
	}

}

func BenchmarkPlanimeter13(b *testing.B) {
	points := [][2]float64{{89, -360}, {89, -240}, {89, -120}, {89, 0}, {89, 120}, {89, 240}}
	for i := 0; i < b.N; i++ {
		planimeter(points)
	}
}

func TestPlanimeter15(t *testing.T) {
	// Coverage tests, includes Planimeter15 - Planimeter18 (combinations of
	// reverse and sign) + calls to testpoint, testedge.
	lat := []float64{2, 1, 3}
	lon := []float64{1, 2, 3}
	r := 18454562325.45119
	a0 := 510065621724088.5093 // ellipsoid area

	polygon := NewPolygonArea(Wgs84(), false)
	polygon.AddPoint(lat[0], lon[0])
	polygon.AddPoint(lat[1], lon[1])
	got := polygon.TestPoint(lat[2], lon[2], false, true)
	if !almost_equal(got.Area, r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.TestPoint(lat[2], lon[2], false, false)
	if !almost_equal(got.Area, r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.TestPoint(lat[2], lon[2], true, true)
	if !almost_equal(got.Area, -r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.TestPoint(lat[2], lon[2], true, false)
	if !almost_equal(got.Area, a0-r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	inv := Wgs84().InverseCalcAll(lat[1], lon[1], lat[2], lat[2])
	azi1 := inv.Azimuth1Deg
	s12 := inv.DistanceM
	got = polygon.TestEdge(azi1, s12, false, true)
	if !almost_equal(got.Area, r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.TestEdge(azi1, s12, false, false)
	if !almost_equal(got.Area, r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.TestEdge(azi1, s12, true, true)
	if !almost_equal(got.Area, -r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.TestEdge(azi1, s12, true, false)
	if !almost_equal(got.Area, a0-r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	polygon.AddPoint(lat[2], lon[2])
	got = polygon.Compute(false, true)
	if !almost_equal(got.Area, r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.Compute(false, false)
	if !almost_equal(got.Area, r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.Compute(true, true)
	if !almost_equal(got.Area, -r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}

	got = polygon.Compute(true, false)
	if !almost_equal(got.Area, a0-r, 0.5) {
		t.Errorf("area = %v; want %v", got.Area, r)
	}
}

func TestPlanimeter19(t *testing.T) {
	// Coverage tests, includes Planimeter19 - Planimeter20 (degenerate
	// polygons) + extra cases.
	polygon := NewPolygonArea(Wgs84(), false)
	got := polygon.Compute(false, true)
	if got.Area != 0.0 {
		t.Errorf("area = %v; want 0", got.Area)
	}
	if got.Perimeter != 0.0 {
		t.Errorf("perimeter = %v; want 0", got.Perimeter)
	}

	polygon.TestPoint(1, 1, false, true)
	if got.Area != 0.0 {
		t.Errorf("area = %v; want 0", got.Area)
	}
	if got.Perimeter != 0.0 {
		t.Errorf("perimeter = %v; want 0", got.Perimeter)
	}

	got = polygon.TestEdge(90, 1000, false, true)
	if !math.IsNaN(got.Area) {
		t.Errorf("area = %v; want NaN", got.Area)
	}
	if !math.IsNaN(got.Perimeter) {
		t.Errorf("perimeter = %v; want NaN", got.Perimeter)
	}

	polygon.AddPoint(1, 1)
	got = polygon.Compute(false, true)
	if got.Area != 0.0 {
		t.Errorf("area = %v; want 0", got.Area)
	}
	if got.Perimeter != 0.0 {
		t.Errorf("perimeter = %v; want 0", got.Perimeter)
	}

	polyline := NewPolygonArea(Wgs84(), true)
	got = polyline.Compute(false, true)
	if got.Perimeter != 0.0 {
		t.Errorf("perimeter = %v; want 0", got.Perimeter)
	}

	got = polyline.TestPoint(1, 1, false, true)
	if got.Perimeter != 0.0 {
		t.Errorf("perimeter = %v; want 0", got.Perimeter)
	}

	got = polyline.TestEdge(90, 1000, false, true)
	if !math.IsNaN(got.Perimeter) {
		t.Errorf("perimeter = %v; want NaN", got.Perimeter)
	}

	polyline.AddPoint(1, 1)
	got = polyline.Compute(false, true)
	if got.Perimeter != 0.0 {
		t.Errorf("perimeter = %v; want 0", got.Perimeter)
	}

	polygon.AddPoint(1, 1)
	got = polyline.TestEdge(90, 1000, false, true)
	if !almost_equal(got.Perimeter, 1000, 1e-10) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, 1000)
	}

	got = polyline.TestPoint(2, 2, false, true)
	want_perimeter := 156876.149
	if !almost_equal(got.Perimeter, want_perimeter, 0.5e-3) {
		t.Errorf("perimeter = %v; want %v", got.Perimeter, want_perimeter)
	}
}

func TestPlanimeter21(t *testing.T) {
	// Some test to add code coverage: multiple circlings of pole (includes
	// Planimeter21 - Planimeter28) + invocations via testpoint and testedge.
	var a PolygonResult
	lat := 45.0
	azi := 39.2144607176828184218
	s := 8420705.40957178156285
	r := 39433884866571.4277   // Area for one circuit
	a0 := 510065621724088.5093 // Ellipsoid area

	polygon := NewPolygonArea(Wgs84(), false)
	polygon.AddPoint(lat, 60)
	polygon.AddPoint(lat, 180)
	polygon.AddPoint(lat, -60)
	polygon.AddPoint(lat, 60)
	polygon.AddPoint(lat, 180)
	polygon.AddPoint(lat, -60)

	ivals := [2]int{3, 4}
	for _, i := range ivals {
		polygon.AddPoint(lat, 60)
		polygon.AddPoint(lat, 180)
		a = polygon.TestPoint(lat, -60, false, true)
		if !almost_equal(a.Area, float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.TestPoint(lat, -60, false, false)
		if !almost_equal(a.Area, float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.TestPoint(lat, -60, true, true)
		if !almost_equal(a.Area, -float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.TestPoint(lat, -60, true, false)
		if !almost_equal(a.Area, -float64(i)*r+a0, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.TestEdge(azi, s, false, true)

		if !almost_equal(a.Area, float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.TestEdge(azi, s, false, false)
		if !almost_equal(a.Area, float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.TestEdge(azi, s, true, true)
		if !almost_equal(a.Area, -float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.TestEdge(azi, s, true, false)
		if !almost_equal(a.Area, -float64(i)*r+a0, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		polygon.AddPoint(lat, -60)
		a = polygon.Compute(false, true)

		if !almost_equal(a.Area, float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.Compute(false, false)

		if !almost_equal(a.Area, float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.Compute(true, true)
		if !almost_equal(a.Area, -float64(i)*r, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}

		a = polygon.Compute(true, false)
		if !almost_equal(a.Area, -float64(i)*r+a0, 0.5) {
			t.Errorf("i=%v; area = %v; want %v", i, a.Area, float64(i)*r)
		}
	}
}

func TestPlanimeter29(t *testing.T) {
	// Check fix to transitdirect vs transit zero handling inconsistency

	var a PolygonResult

	polygon := NewPolygonArea(Wgs84(), false)
	polygon.Clear()
	polygon.AddPoint(0, 0)
	polygon.AddEdge(90, 1000)
	polygon.AddEdge(0, 1000)
	polygon.AddEdge(-90, 1000)
	a = polygon.Compute(false, true)
	// The area should be 1e6.  Prior to the fix it was 1e6 - A/2, where
	// A = ellipsoid area.
	want_area := 1000000.0
	if !almost_equal(a.Area, want_area, 0.01) {
		t.Errorf("area = %v; want %v", a.Area, want_area)
	}
}
