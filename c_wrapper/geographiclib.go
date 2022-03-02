// Copyright 2020 The Cockroach Authors.
//
// Use of this software is governed by the Business Source License
// included in the file licenses/BSL.txt.
//
// As of the Change Date specified in that file, in accordance with
// the Business Source License, use of this software will be governed
// by the Apache License, Version 2.0, included in the file
// licenses/APL.txt.
//
// This module is included only for testing purposes only.

// Package c_wrapper is a wrapper around the GeographicLib library.
package c_wrapper

// #cgo CXXFLAGS: -std=c++14
// #cgo LDFLAGS: -lm
//
// #include "geodesic.h"
// #include "geographiclib.h"
import "C"

import (
	"math"

	"github.com/golang/geo/s1"
	"github.com/golang/geo/s2"
)

// Spheroid is an object that can perform geodesic operations
// on a given spheroid.
type Geodesic struct {
	cRepr        C.struct_geod_geodesic
	Radius       float64
	Flattening   float64
	SphereRadius float64
}

// NewSpheroid creates a spheroid from a radius and flattening.
func NewGeodesic(radius float64, flattening float64) *Geodesic {
	minorAxis := radius - radius*flattening
	s := &Geodesic{
		Radius:       radius,
		Flattening:   flattening,
		SphereRadius: (radius*2 + minorAxis) / 3,
	}
	C.geod_init(&s.cRepr, C.double(radius), C.double(flattening))
	return s
}

// Wgs84 represents the default WGS84 ellipsoid.
func Wgs84() *Geodesic {
	return NewGeodesic(6378137, 1/298.257223563)
}

// LatLon represents latitude and longitude of a point. All units in degrees
type LatLon struct {
	LatDeg, LonDeg float64
}

// DirectCalcLatLon gets the lat and lon of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcLatLon(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLon {
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

// LatLonAzi represents latitude, longitude, and azimuth of a point. All units in degrees
type LatLonAzi struct {
	LatDeg, LonDeg, AziDeg float64
}

// DirectCalcLatLonAzi gets the lat, lon, and azimuth of the second point, based on input
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcLatLonAzi(lat1_deg, lon1_deg, azi1_deg, s12_m float64) LatLonAzi {
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

// AllDirectResults contains all information that can be computed from the direct method
// latitude, longitude, azimuth, reduced length, geodesic scales, area under the geodesic,
// and arc length between point 1 and point 2
type AllDirectResults struct {
	LatDeg         float64 // Latitude [degrees]
	LonDeg         float64 // Longitude [degrees]
	AziDeg         float64 // Azimuth [degrees]
	ReducedLengthM float64 // Reduced length of the geodesic [meters]
	M12            float64 // Geodesic scale of point 2 relative to point 1 [dimensionless]
	M21            float64 // Geodesic scale of point 1 relative to point 2 [dimensionless]
	S12M2          float64 // Area under the geodesic [meters^2]
	A12Deg         float64 // Arc length between point 1 and point 2 [degrees]
}

// DirectCalcAll calculates everything possible for the direct method. Takes inputs
//   - lat1_deg - Latitude of 1st point [degrees] [-90.,90.]
//   - lon1_deg - Longitude of 1st point [degrees] [-180., 180.]
//   - azi1_deg - Azimuth at 1st point [degrees] [-180., 180.]
//   - s12_m - Distance from 1st to 2nd point [meters] Value may be negative
func (g Geodesic) DirectCalcAll(lat1_deg, lon1_deg, azi1_deg, s12_m float64) AllDirectResults {
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
func (s *Geodesic) Inverse(a, b s2.LatLng) (s12, az1, az2 float64) {
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

// InverseBatch computes the sum of the length of the lines represented
// by the line of points.
// This is intended for use for LineStrings. LinearRings/Polygons should use "AreaAndPerimeter".
// Returns the sum of the s12 (distance in meters) units.
func (s *Geodesic) InverseBatch(points []s2.Point) float64 {
	lats := make([]C.double, len(points))
	lngs := make([]C.double, len(points))
	for i, p := range points {
		latlng := s2.LatLngFromPoint(p)
		lats[i] = C.double(latlng.Lat.Degrees())
		lngs[i] = C.double(latlng.Lng.Degrees())
	}
	var result C.double
	C.CR_GEOGRAPHICLIB_InverseBatch(
		&s.cRepr,
		&lats[0],
		&lngs[0],
		C.int(len(points)),
		&result,
	)
	return float64(result)
}

// AreaAndPerimeter computes the area and perimeter of a polygon on a given spheroid.
// The points must never be duplicated (i.e. do not include the "final" point of a Polygon LinearRing).
// Area is in meter^2, Perimeter is in meters.
func (s *Geodesic) AreaAndPerimeter(points []s2.Point) (area float64, perimeter float64) {
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
func (s *Geodesic) Project(point s2.LatLng, distance float64, azimuth s1.Angle) s2.LatLng {
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
