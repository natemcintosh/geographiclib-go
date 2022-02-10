package geographiclibgo

import "math"

type PolygonArea struct {
	Earth         Geodesic
	Polyline      bool
	Area0_M2      float64
	_mask         uint64
	_areasum      Accumulator
	_perimetersum Accumulator
	Num           float64
	Lat1_Deg      float64
	_lat0_deg     float64
	Lon1_Deg      float64
	_lon0_deg     float64
	_crossings    int
}

func NewPolygonArea(earth Geodesic, is_polyline bool) PolygonArea {
	// The total area of the ellipsoid in meter^2 (readonly)
	area0 := 4 * math.Pi * earth._c2

	var ternary_opt uint64
	if is_polyline {
		ternary_opt = EMPTY
	} else {
		ternary_opt = AREA | LONG_UNROLL
	}
	mask := LATITUDE | LONGITUDE | DISTANCE | ternary_opt

	areasum := Accumulator{}

	perimetersum := Accumulator{}

	// The current number of points in the polygon
	Num := 0.0

	// The current latitude [degrees]
	Lat1_Deg := math.NaN()

	// The current longitude [degrees]
	Lon1_Deg := math.NaN()

	p := PolygonArea{
		Earth:         earth,
		Polyline:      is_polyline,
		Area0_M2:      area0,
		_mask:         mask,
		_areasum:      areasum,
		_perimetersum: perimetersum,
		Num:           Num,
		Lat1_Deg:      Lat1_Deg,
		_lat0_deg:     Lat1_Deg,
		Lon1_Deg:      Lon1_Deg,
		_lon0_deg:     Lon1_Deg,
		_crossings:    0,
	}

	p.Clear()

	return p
}

// Clear resets to an empty polygon
func (p *PolygonArea) Clear() {
	p.Num = 0
	p._crossings = 0
	if !p.Polyline {
		p._areasum.SetS(0)
	}
	p._perimetersum.SetS(0)
	p._lat0_deg = math.NaN()
	p._lon0_deg = math.NaN()
	p.Lat1_Deg = math.NaN()
	p.Lon1_Deg = math.NaN()
}

// transit counts crossings of prime meridian for AddPoint
func (p PolygonArea) transit(lon1_deg, lon2_deg float64) int {
	// Return 1 or -1 if crossing prime meridian in east or west direction.
	// Otherwise return zero.
	// Compute lon12 the same way as Geodesic::Inverse.
	lon1 := Ang_normalize(lon1_deg)
	lon2 := Ang_normalize(lon2_deg)
	lon12, _ := Ang_diff(lon1, lon2)

	if (lon1 <= 0) && (lon2 > 0) && (lon12 > 0) {
		return 1
	} else {
		if (lon2 <= 0) && (lon1 > 0) && (lon12 < 0) {
			return -1
		} else {
			return 0
		}
	}
}

// AddPoint adds the next vertex to the polygon. This adds an edge from the current
// vertex to the new vertex
func (p *PolygonArea) AddPoint(lat_deg, lon_deg float64) {
	if p.Num == 0 {
		p._lat0_deg = lat_deg
		p.Lat1_Deg = lat_deg

		p._lon0_deg = lon_deg
		p.Lon1_Deg = lon_deg
	} else {
		_, s12, _, _, _, _, _, _, _, S12 := p.Earth._gen_inverse(
			p.Lat1_Deg, p.Lon1_Deg, lat_deg, lon_deg, p._mask,
		)
		p._perimetersum.Add(s12)

		if !p.Polyline {
			p._areasum.Add(S12)
			p._crossings += p.transit(p.Lon1_Deg, lon_deg)
		}
		p.Lat1_Deg = lat_deg
		p.Lon1_Deg = lon_deg
	}
	p.Num += 1
}
