package geographiclibgo

import (
	"math"
)

type PolygonArea struct {
	Earth         Geodesic
	Polyline      bool
	Area0_M2      float64
	_mask         uint64
	_areasum      Accumulator
	_perimetersum Accumulator
	Num           int
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
	Num := 0

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
func (p *PolygonArea) transit(lon1_deg, lon2_deg float64) int {
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

// transit_direct counts the crossings of the prime meridian for AddEdge
func (p *PolygonArea) transit_direct(lon1_deg, lon2_deg float64) int {
	lon1 := math.Mod(lon1_deg, 720.0)
	lon2 := math.Mod(lon2_deg, 720.0)

	var first_num int
	if (lon2 <= 0 && lon2 > -360) || lon2 > 360 {
		first_num = 1
	} else {
		first_num = 0
	}

	var second_num int
	if (lon1 <= 0 && lon1 > -360) || lon1 > 360 {
		second_num = 1
	} else {
		second_num = 0
	}

	return first_num - second_num
}

// AddEdge adds the next edge to the polygon. This specifies the new vertex in terms of
// the edge from the current vertex.
// INPUTS:
// - azi_deg: the azimuth at the current point in degrees
// - s_m: the length of the edge in meters
func (p *PolygonArea) AddEdge(azi_deg, s_m float64) {
	if p.Num != 0 {
		_, lat, lon, _, _, _, _, _, S12, _ := p.Earth._gen_direct(
			p.Lat1_Deg, p.Lon1_Deg, azi_deg, false, s_m, p._mask,
		)
		p._perimetersum.Add(s_m)

		if !p.Polyline {
			p._areasum.Add(S12)
			p._crossings += p.transit_direct(p.Lon1_Deg, lon)
		}

		p.Lat1_Deg = lat
		p.Lon1_Deg = lon
		p.Num += 1
	}
}

// area_reduce_A reduce accumulator area to allowed range
func (p *PolygonArea) area_reduce_A(
	area Accumulator,
	area0 float64,
	crossings int,
	reverse bool,
	sign bool,
) float64 {
	area.Remainder(area0)

	if crossings&1 != 0 {
		var toadd float64
		if area.Sum(0.0) < 0 {
			toadd = 1
		} else {
			toadd = -1
		}
		area.Add(toadd * area0 / 2)
	}

	// area is with the clockwise sense.  If !reverse convert to
	// counter-clockwise convention.
	if !reverse {
		area.Negate()
	}

	// If sign put area in (-area0/2, area0/2], else put area in [0, area0)
	if sign {
		if area.Sum(0.0) > area0/2 {
			area.Add(-area0)
		} else if area.Sum(0.0) <= (-area0 / 2) {
			area.Add(+area0)
		}

	} else {
		if area.Sum(0.0) >= area0 {
			area.Add(-area0)
		} else if area.Sum(0.0) < 0 {
			area.Add(+area0)
		}
	}

	return 0 + area.Sum(0.0)
}

// area_reduce_B reduce double area to allowed range
func (p *PolygonArea) area_reduce_B(
	area float64,
	area0 float64,
	crossings int,
	reverse bool,
	sign bool,
) float64 {
	area = math.Remainder(area, area0)

	if crossings&1 != 0 {
		if area < 0 {
			area += 1 * area0 / 2
		} else {
			area += -1 * area0 / 2
		}
	}

	// area is with the clockwise sense.  If !reverse convert to
	// counter-clockwise convention.
	if !reverse {
		area *= -1
	}

	// If sign put area in (-area0/2, area0/2], else put area in [0, area0)
	if sign {
		if area > area0/2 {
			area -= area0
		} else if area <= -area0/2 {
			area += area0
		}
	} else {
		if area >= area0 {
			area -= area0
		} else if area < 0 {
			area += area0
		}
	}
	return 0 + area
}

type PolygonResult struct {
	Num       int     // the number of vertices in the polygon
	Perimeter float64 // the perimeter fo the polygon or th elength of the polyline [meters]
	Area      float64 // the area of the polygon [meters^2]
}

// Compute the properties of the polygon
//
// INPUTS:
//
// - reverse: if true then clockwise (instead of counter-clockwise) traversal counts as
// a positive area
//
// - sign: if true then return a signed result for the area if the polygon is traversed
// in the "wrong" direction instead of returning the area for the rest of the earth
//
// Arbitrarily complex polygons are allowed.  In the case of
// self-intersecting polygons the area is accumulated "algebraically",
// e.g., the areas of the 2 loops in a figure-8 polygon will partially
// cancel.
//
// If the object is a polygon (and not a polyline), the perimeter
// includes the length of a final edge connecting the current point to
// the initial point.  If the object is a polyline, then area is nan.
//
// More points can be added to the polygon after this call.
func (p *PolygonArea) Compute(reverse, sign bool) PolygonResult {
	if p.Num < 2 {
		var area float64
		if p.Polyline {
			area = math.NaN()
		} else {
			area = 0.0
		}
		return PolygonResult{Num: p.Num, Perimeter: 0.0, Area: area}
	}

	if p.Polyline {
		return PolygonResult{Num: p.Num, Perimeter: p._perimetersum.Sum(0.0), Area: math.NaN()}
	}

	_, s12, _, _, _, _, _, _, _, S12 := p.Earth._gen_inverse(
		p.Lat1_Deg,
		p.Lon1_Deg,
		p._lat0_deg,
		p._lon0_deg,
		p._mask,
	)

	perimeter := p._perimetersum.Sum(s12)

	tempsum := p._areasum
	tempsum.Add(S12)

	crossings := p._crossings + p.transit(p.Lon1_Deg, p._lon0_deg)
	area := p.area_reduce_A(tempsum, p.Area0_M2, crossings, reverse, sign)

	return PolygonResult{Num: p.Num, Perimeter: perimeter, Area: area}
}

// TestPoint return the results assuming a tentative final test point is added;
// however, the data for the test point is not saved. This lets you report
// a running result for the perimeter and area as the user moves the mouse
// cursor. Ordinary floating point arithmetic is used to accumulate the
// data for the test point; thus the area and perimeter returned are less
// accurate than if AddPoint and Compute are used.
//
// INPUTS:
//
// - lat_deg: the latitude of the test point [degrees]. Should be in the range [-90, 90]
//
// - lon_deg: the longitude of the test point [degrees]
//
// - reverse: if true then clockwise (instead of counter-clockwise) traversal counts as
// a positive area
//
// - sign: if true then return a signed result for the area if
// the polygon is traversed in the "wrong" direction instead of returning
// the area for the rest of the earth
func (p *PolygonArea) TestPoint(
	lat_deg float64,
	lon_deg float64,
	reverse bool,
	sign bool,
) PolygonResult {

	if p.Num == 0 {
		var area float64
		if p.Polyline {
			area = math.NaN()
		} else {
			area = 0
		}
		return PolygonResult{Num: 1, Perimeter: 0, Area: area}
	}

	perimeter := p._perimetersum.Sum(0.0)
	var tempsum float64
	if p.Polyline {
		tempsum = 0
	} else {
		tempsum = p._areasum.Sum(0.0)
	}
	crossings := p._crossings
	num := p.Num + 1

	var end_point int
	if p.Polyline {
		end_point = 1
	} else {
		end_point = 2
	}

	for i := 0; i < end_point; i++ {
		var this_lat1 float64
		var this_lon1 float64
		var this_lat0 float64
		var this_lon0 float64
		if i == 0 {
			this_lat1 = p.Lat1_Deg
			this_lon1 = p.Lon1_Deg
			this_lat0 = lat_deg
			this_lon0 = lon_deg
		} else {
			this_lat1 = lat_deg
			this_lon1 = lon_deg
			this_lat0 = p._lat0_deg
			this_lon0 = p._lon0_deg
		}
		_, s12, _, _, _, _, _, _, _, S12 := p.Earth._gen_inverse(
			this_lat1, this_lon1, this_lat0, this_lon0, p._mask,
		)
		perimeter += s12
		if !p.Polyline {
			tempsum += S12
			var arg1 float64
			var arg2 float64
			if i == 0 {
				arg1 = p.Lon1_Deg
				arg2 = lon_deg
			} else {
				arg1 = lon_deg
				arg2 = p._lon0_deg
			}
			crossings += p.transit(arg1, arg2)
		}
	}

	if p.Polyline {
		return PolygonResult{Num: num, Perimeter: perimeter, Area: math.NaN()}
	}

	area := p.area_reduce_B(tempsum, p.Area0_M2, crossings, reverse, sign)
	return PolygonResult{Num: num, Perimeter: perimeter, Area: area}
}

// TestEdge returns the results assuming a tentative final test point is added via an
// azimuth and distance; however, the data for the test point is not saved.
// This lets you report a running result for the perimeter and area as the
// user moves the mouse cursor. Ordinary floating point arithmetic is used
// to accumulate the data for the test point; thus the area and perimeter
// returned are less accurate than if AddPoint and Compute are used.
//
// INPUTS:
//
// - azi_deg: azimuth at the current point [degrees]
//
// - s_m: distance from the current point to the final test point [meters]
//
// - reverse: if true then clockwise (instead of counter-clockwise) traversal counts as
// a positive area
//
// - sign: if true then return a signed result for the area if
// the polygon is traversed in the "wrong" direction instead of returning
// the area for the rest of the earth
func (p *PolygonArea) TestEdge(
	azi_deg float64,
	s_m float64,
	reverse bool,
	sign bool,
) PolygonResult {
	// We don't have a starting point
	if p.Num == 0 {
		return PolygonResult{Num: 0, Perimeter: math.NaN(), Area: math.NaN()}
	}

	num := p.Num + 1
	perimeter := p._perimetersum.Sum(0.0) + s_m

	if p.Polyline {
		return PolygonResult{Num: num, Perimeter: perimeter, Area: math.NaN()}
	}

	tempsum := p._areasum.Sum(0.0)
	crossings := p._crossings

	_, lat, lon, _, _, _, _, _, S12, _ := p.Earth._gen_direct(
		p.Lat1_Deg, p.Lon1_Deg, azi_deg, false, s_m, p._mask,
	)

	tempsum += S12
	crossings += p.transit_direct(p.Lon1_Deg, lon)

	_, s12, _, _, _, _, _, _, _, S12 := p.Earth._gen_inverse(
		lat, lon, p._lat0_deg, p._lon0_deg, p._mask,
	)
	perimeter += s12
	tempsum += S12
	crossings += p.transit(lon, p._lon0_deg)

	area := p.area_reduce_B(tempsum, p.Area0_M2, crossings, reverse, sign)
	return PolygonResult{Num: num, Perimeter: perimeter, Area: area}

}
