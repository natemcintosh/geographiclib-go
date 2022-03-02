package geographiclibgo

// #include "geodesic.h"
import "C"

import (
	"math"

	"github.com/golang/geo/s1"
	"github.com/golang/geo/s2"
)

// CGeodesic is an object that can perform geodesic operations
// on a given spheroid.
type CGeodesic struct {
	cRepr        C.struct_geod_geodesic
	Radius       float64
	Flattening   float64
	SphereRadius float64
}

// NewCGeodesic creates a CGeodesic from a radius and flattening.
func NewCGeodesic(radius float64, flattening float64) *CGeodesic {
	minorAxis := radius - radius*flattening
	s := &CGeodesic{
		Radius:       radius,
		Flattening:   flattening,
		SphereRadius: (radius*2 + minorAxis) / 3,
	}
	C.geod_init(&s.cRepr, C.double(radius), C.double(flattening))
	return s
}

// CWgs84 represents the c wrapper of the default WGS84 ellipsoid
func CWgs84() *CGeodesic {
	return NewCGeodesic(6378137, 1/298.257223563)
}

// DirectCalcLatLon gets the lat and lon of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g CGeodesic) DirectCalcLatLon(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLon {
	var retLatDeg, retLonDeg, placeholder C.double
	capabilities := LATITUDE | LONGITUDE

	C.geod_gendirect(
		&g.cRepr,
		C.double(lat1_deg),
		C.double(lon1_deg),
		C.double(azi1_deg),
		C.unsigned(capabilities),
		C.double(s12_m),
		&retLatDeg,
		&retLonDeg,
		&placeholder,
		&placeholder,
		&placeholder,
		&placeholder,
		&placeholder,
		&placeholder,
	)
	lat2 := float64(retLatDeg)
	lon2 := float64(retLonDeg)
	return LatLon{LatDeg: lat2, LonDeg: lon2}
}

// DirectCalcLatLonAzi gets the lat, lon, and azimuth of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g CGeodesic) DirectCalcLatLonAzi(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAzi {
	var retLatDeg, retLonDeg, retAziDeg, placeholder C.double
	capabilities := LATITUDE | LONGITUDE | AZIMUTH

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
		&placeholder,
		&placeholder,
		&placeholder,
		&placeholder,
		&placeholder,
	)
	lat2 := float64(retLatDeg)
	lon2 := float64(retLonDeg)
	azi2 := float64(retAziDeg)
	return LatLonAzi{LatDeg: lat2, LonDeg: lon2, AziDeg: azi2}
}

// DirectCalcAll calculates everything possible for the direct method. Takes inputs
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g CGeodesic) DirectCalcAll(lat1_deg, lon1_deg, azi1_deg, s12_m float64) AllDirectResults {
	var retLatDeg, retLonDeg, retAziDeg, retReducedLength, retM12, retM21, retS12M2, retA12Deg C.double
	capabilities := ALL

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
		&retReducedLength,
		&retM12,
		&retM21,
		&retS12M2,
		&retA12Deg,
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

// Inverse solves the geodetic inverse problem on the given spheroid
// (https://en.wikipedia.org/wiki/Geodesy#Geodetic_problems).
// Returns s12 (distance in meters), az1 (azimuth at point 1) and az2 (azimuth at point 2).
func (s *CGeodesic) Inverse(a, b s2.LatLng) (s12, az1, az2 float64) {
	var retS12, retAZ1, retAZ2 C.double
	C.geod_inverse(
		&s.cRepr,
		C.double(a.Lat.Degrees()),
		C.double(a.Lng.Degrees()),
		C.double(b.Lat.Degrees()),
		C.double(b.Lng.Degrees()),
		&retS12,
		&retAZ1,
		&retAZ2,
	)
	return float64(retS12), float64(retAZ1), float64(retAZ2)
}

// AreaAndPerimeter computes the area and perimeter of a polygon on a given spheroid.
// The points must never be duplicated (i.e. do not include the "final" point of a Polygon LinearRing).
// Area is in meter^2, Perimeter is in meters.
func (s *CGeodesic) AreaAndPerimeter(points []s2.Point) (area float64, perimeter float64) {
	lats := make([]C.double, len(points))
	lngs := make([]C.double, len(points))
	for i, p := range points {
		latlng := s2.LatLngFromPoint(p)
		lats[i] = C.double(latlng.Lat.Degrees())
		lngs[i] = C.double(latlng.Lng.Degrees())
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
func (s *CGeodesic) Project(point s2.LatLng, distance float64, azimuth s1.Angle) s2.LatLng {
	var lat, lng C.double

	C.geod_direct(
		&s.cRepr,
		C.double(point.Lat.Degrees()),
		C.double(point.Lng.Degrees()),
		C.double(azimuth*180.0/math.Pi),
		C.double(distance),
		&lat,
		&lng,
		nil,
	)

	return s2.LatLngFromDegrees(float64(lat), float64(lng))
}
