package geographiclibgo

// #cgo CXXFLAGS: -std=c++14
// #cgo LDFLAGS: -lm
//
// #include "geodesic.h"
import "C"

// CGeodesic is an object that can perform geodesic operations
// on a given spheroid.
type CGeodesic struct {
	cRepr         C.struct_geod_geodesic
	a             float64
	f             float64
	sphere_radius float64
}

// NewCGeodesic creates a CGeodesic from a radius and flattening.
func NewCGeodesic(radius float64, flattening float64) *CGeodesic {
	minorAxis := radius - radius*flattening
	s := &CGeodesic{
		a:             radius,
		f:             flattening,
		sphere_radius: (radius*2 + minorAxis) / 3,
	}
	C.geod_init(&s.cRepr, C.double(radius), C.double(flattening))
	return s
}

// CWgs84 represents the c wrapper of the default WGS84 ellipsoid
func CWgs84() *CGeodesic {
	return NewCGeodesic(6378137, 1/298.257223563)
}

func (g *CGeodesic) EqualtorialRadius() float64 {
	return g.a
}

func (g *CGeodesic) Flattening() float64 {
	return g.f
}

// DirectCalcLatLon gets the lat and lon of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g *CGeodesic) DirectCalcLatLon(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLon {
	capabilities := LATITUDE | LONGITUDE

	res := g.DirectCalcWithCapabilities(lat1_deg, lon1_deg, azi1_deg, s12_m, capabilities)

	return LatLon{LatDeg: res.LatDeg, LonDeg: res.LonDeg}
}

// DirectCalcLatLonAzi gets the lat, lon, and azimuth of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g *CGeodesic) DirectCalcLatLonAzi(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAzi {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH

	res := g.DirectCalcWithCapabilities(lat1_deg, lon1_deg, azi1_deg, s12_m, capabilities)

	return LatLonAzi{LatDeg: res.LatDeg, LonDeg: res.LonDeg, AziDeg: res.AziDeg}
}

// DirectCalcLatLonAziReducedLength gets the lat, lon, azimuth, and reduced length of geodesic
// of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g *CGeodesic) DirectCalcLatLonAziReducedLength(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAziReducedLength {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH | REDUCEDLENGTH

	res := g.DirectCalcWithCapabilities(lat1_deg, lon1_deg, azi1_deg, s12_m, capabilities)

	return LatLonAziReducedLength{
		LatDeg:         res.LatDeg,
		LonDeg:         res.LonDeg,
		AziDeg:         res.AziDeg,
		ReducedLengthM: res.ReducedLengthM,
	}
}

// DirectCalcLatLonAziGeodesicScales gets the lat, lon, azimuth, and geodesic scales,
// based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g *CGeodesic) DirectCalcLatLonAziGeodesicScales(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAziGeodesicScales {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH | GEODESICSCALE

	res := g.DirectCalcWithCapabilities(lat1_deg, lon1_deg, azi1_deg, s12_m, capabilities)

	return LatLonAziGeodesicScales{
		LatDeg: res.LatDeg,
		LonDeg: res.LonDeg,
		AziDeg: res.AziDeg,
		M12:    res.M12,
		M21:    res.M21,
	}
}

// DirectCalcLatLonAziReducedLengthGeodesicScales gets the lat, lon, azimuth, reduced length,
// and geodesic scales based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g *CGeodesic) DirectCalcLatLonAziReducedLengthGeodesicScales(
	lat1_deg, lon1_deg, azi1_deg, s12_m float64,
) LatLonAziReducedLengthGeodesicScales {
	capabilities := LATITUDE | LONGITUDE | AZIMUTH | REDUCEDLENGTH | GEODESICSCALE

	res := g.DirectCalcWithCapabilities(lat1_deg, lon1_deg, azi1_deg, s12_m, capabilities)

	return LatLonAziReducedLengthGeodesicScales{
		LatDeg:         res.LatDeg,
		LonDeg:         res.LonDeg,
		AziDeg:         res.AziDeg,
		ReducedLengthM: res.ReducedLengthM,
		M12:            res.M12,
		M21:            res.M21,
	}
}

// DirectCalcAll calculates everything possible for the direct method. Takes inputs
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g *CGeodesic) DirectCalcAll(lat1_deg, lon1_deg, azi1_deg, s12_m float64) AllDirectResults {

	capabilities := ALL

	res := g.DirectCalcWithCapabilities(lat1_deg, lon1_deg, azi1_deg, s12_m, capabilities)

	return AllDirectResults{
		LatDeg:         res.LatDeg,
		LonDeg:         res.LonDeg,
		AziDeg:         res.AziDeg,
		ReducedLengthM: res.ReducedLengthM,
		M12:            res.M12,
		M21:            res.M21,
		S12M2:          res.S12M2,
		A12Deg:         res.A12Deg,
	}
}

// DirectCalcWithCapabilities takes inputs
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - capabilities - One or more of the capabilities constant as defined in the file
//     geodesiccapability.go. Usually, they are OR'd together, e.g. LATITUDE | LONGITUDE
func (g *CGeodesic) DirectCalcWithCapabilities(
	lat1_deg float64,
	lon1_deg float64,
	azi1_deg float64,
	s12_m float64,
	capabilities uint64,
) AllDirectResults {
	var retLatDeg, retLonDeg, retAziDeg, ret_s12M, retReducedLength, retM12, retM21, retS12M2, retA12Deg C.double

	C.geod_gendirect(
		&g.cRepr,
		C.double(lat1_deg),
		C.double(lon1_deg),
		C.double(azi1_deg),
		C.unsigned(capabilities),
		C.double(s12_m),
		&retLatDeg,
		&retLonDeg,
		&retAziDeg,
		&ret_s12M,
		&retReducedLength,
		&retM12,
		&retM21,
		&retS12M2,
	)

	return AllDirectResults{
		LatDeg:         float64(retLatDeg),
		LonDeg:         float64(retLonDeg),
		AziDeg:         float64(retAziDeg),
		ReducedLengthM: float64(retReducedLength),
		M12:            float64(retM12),
		M21:            float64(retM21),
		S12M2:          float64(retS12M2),
		A12Deg:         float64(retA12Deg),
	}
}

// InverseCalcDistance returns the distance from point 1 to point 2 in meters. Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g *CGeodesic) InverseCalcDistance(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) float64 {
	capabilities := DISTANCE

	res := g.InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return res.DistanceM
}

// InverseCalcDistanceArcLength returns the distance from one point to the next, and the
// arc length between the points. Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g *CGeodesic) InverseCalcDistanceArcLength(lat1_deg, lon1_deg, lat2_deg, lon2_deg float64) DistanceArcLength {
	capabilities := DISTANCE

	res := g.InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return DistanceArcLength{DistanceM: res.DistanceM, ArcLengthDeg: res.ArcLengthDeg}
}

// InverseCalcAzimuthsArcLength returns the azimuth at point 1, the azimuth at point 2,
// and the arc length between the points. Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g *CGeodesic) InverseCalcAzimuthsArcLength(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) AzimuthsArcLength {
	capabilities := AZIMUTH

	res := g.InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return AzimuthsArcLength{
		Azimuth1Deg:  res.Azimuth1Deg,
		Azimuth2Deg:  res.Azimuth2Deg,
		ArcLengthDeg: res.ArcLengthDeg,
	}
}

// InverseCalcDistanceAzimuths returns the distance from one point to the next, and the
// azimuths. Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g *CGeodesic) InverseCalcDistanceAzimuths(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) DistanceAzimuths {
	capabilities := DISTANCE | AZIMUTH

	res := g.InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return DistanceAzimuths{
		DistanceM:   res.DistanceM,
		Azimuth1Deg: res.Azimuth1Deg,
		Azimuth2Deg: res.Azimuth2Deg,
	}
}

// InverseCalcDistanceAzimuthsArcLength returns the distance from one point to the next,
// the azimuth at point 1, the azimuth at point 2, and the arc length between the points.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g *CGeodesic) InverseCalcDistanceAzimuthsArcLength(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) DistanceAzimuthsArcLength {
	capabilities := DISTANCE | AZIMUTH

	res := g.InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return DistanceAzimuthsArcLength{
		DistanceM:    res.DistanceM,
		Azimuth1Deg:  res.Azimuth1Deg,
		Azimuth2Deg:  res.Azimuth2Deg,
		ArcLengthDeg: res.ArcLengthDeg,
	}
}

// InverseCalcDistanceAzimuthsArcLengthReducedLength returns the distance from one point
// to the next, the azimuth at point 1, the azimuth at point 2, the arc length
// between the points, and the reduceed length of the geodesic.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g *CGeodesic) InverseCalcDistanceAzimuthsArcLengthReducedLength(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) DistanceAzimuthsArcLengthReducedLength {
	capabilities := DISTANCE | AZIMUTH | REDUCEDLENGTH

	res := g.InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return DistanceAzimuthsArcLengthReducedLength{
		DistanceM:      res.DistanceM,
		Azimuth1Deg:    res.Azimuth1Deg,
		Azimuth2Deg:    res.Azimuth2Deg,
		ArcLengthDeg:   res.ArcLengthDeg,
		ReducedLengthM: res.ReducedLengthM,
	}
}

// InverseCalcDistanceAzimuthsArcLengthReducedLengthScales returns everything described
// by the `DistanceAzimuthsArcLengthReducedLengthScales` type.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g *CGeodesic) InverseCalcDistanceAzimuthsArcLengthReducedLengthScales(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) DistanceAzimuthsArcLengthReducedLengthScales {
	capabilities := DISTANCE | AZIMUTH | REDUCEDLENGTH | GEODESICSCALE

	res := g.InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return DistanceAzimuthsArcLengthReducedLengthScales{
		DistanceM:      res.DistanceM,
		Azimuth1Deg:    res.Azimuth1Deg,
		Azimuth2Deg:    res.Azimuth2Deg,
		ArcLengthDeg:   res.ArcLengthDeg,
		ReducedLengthM: res.ReducedLengthM,
		M12:            res.M12,
		M21:            res.M21,
	}
}

// InverseCalcAll returns everything described in the `AllInverseResults` results type.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
func (g *CGeodesic) InverseCalcAll(
	lat1_deg, lon1_deg, lat2_deg, lon2_deg float64,
) AllInverseResults {
	capabilities := DISTANCE | AZIMUTH | REDUCEDLENGTH | GEODESICSCALE | AREA

	res := g.InverseCalcWithCapabilities(lat1_deg, lon1_deg, lat2_deg, lon2_deg, capabilities)

	return AllInverseResults{
		DistanceM:      res.DistanceM,
		Azimuth1Deg:    res.Azimuth1Deg,
		Azimuth2Deg:    res.Azimuth2Deg,
		ArcLengthDeg:   res.ArcLengthDeg,
		ReducedLengthM: res.ReducedLengthM,
		M12:            res.M12,
		M21:            res.M21,
		S12M2:          res.S12M2,
	}
}

// InverseCalcWithCapabilities allows the user to specify which capabilites they wish to use.
// This function is useful if you want some other subset of capabilities than those offered
// by the other InverseCalc...() methods.
// Takes inputs
// - lat1_deg latitude of point 1 [degrees].
// - lon1_deg longitude of point 1 [degrees].
// - lat2_deg latitude of point 2 [degrees].
// - lon2_deg longitude of point 2 [degrees].
// - capabilities - One or more of the capabilities constant as defined in the file
//     geodesiccapability.go. Usually, they are OR'd together, e.g. LATITUDE | LONGITUDE
func (g *CGeodesic) InverseCalcWithCapabilities(
	lat1_deg float64,
	lon1_deg float64,
	lat2_deg float64,
	lon2_deg float64,
	capabilities uint64,
) AllInverseResults {
	var rets12M, retAzi1Deg, retAzi2Deg, retReducedLengthM, retM12, retM21, retS12M2 C.double

	a12 := C.geod_geninverse(
		&g.cRepr,
		C.double(lat1_deg),
		C.double(lon1_deg),
		C.double(lat2_deg),
		C.double(lon2_deg),
		&rets12M,
		&retAzi1Deg,
		&retAzi2Deg,
		&retReducedLengthM,
		&retM12,
		&retM21,
		&retS12M2,
	)

	return AllInverseResults{
		DistanceM:      float64(rets12M),
		Azimuth1Deg:    float64(retAzi1Deg),
		Azimuth2Deg:    float64(retAzi2Deg),
		ArcLengthDeg:   float64(a12),
		ReducedLengthM: float64(retReducedLengthM),
		M12:            float64(retM12),
		M21:            float64(retM21),
		S12M2:          float64(retS12M2),
	}
}

// AreaAndPerimeter computes the area and perimeter of a polygon on a given spheroid.
// The points must never be duplicated (i.e. do not include the "final" point of a Polygon LinearRing).
// Area is in meter^2, Perimeter is in meters.
func (s *CGeodesic) AreaAndPerimeter(points [][2]float64) (area float64, perimeter float64) {
	lats := make([]C.double, len(points))
	lngs := make([]C.double, len(points))
	for i, p := range points {
		lats[i] = C.double(p[0])
		lngs[i] = C.double(p[1])
	}
	var areaDouble, perimeterDouble C.double
	C.geod_polygonarea(
		&s.cRepr,
		&lats[0],
		&lngs[0],
		C.int(len(points)),
		&areaDouble,
		&perimeterDouble,
	)
	return float64(areaDouble), float64(perimeterDouble)
}

// Project returns computes the location of the projected point.
//
// Using the direct geodesic problem from GeographicLib (Karney 2013).
func (s *CGeodesic) Project(
	lat_in_deg float64,
	lon_in_deg float64,
	distance float64,
	azimuth_deg float64,
) (lat_deg, lon_deg float64) {
	var lat, lng C.double

	C.geod_direct(
		&s.cRepr,
		C.double(lat_in_deg),
		C.double(lon_in_deg),
		C.double(azimuth_deg),
		C.double(distance),
		&lat,
		&lng,
		nil,
	)

	return float64(lat), float64(lng)
}
