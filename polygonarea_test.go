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
