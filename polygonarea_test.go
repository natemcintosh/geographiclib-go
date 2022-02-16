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
