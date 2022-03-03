package geographiclibgo

// #cgo CXXFLAGS: -std=c++14
// #cgo LDFLAGS: -lm
//
// #include "geodesic.h"
import "C"

type CGeodesicLine struct {
	cRepr C.struct_geod_geodesicline
}

// // NewCGeodesicLine creates a GeodesicLine, with `caps` of STANDARD | DISTANCE_IN
func NewCGeodesicLine(
	geod *CGeodesic,
	lat1, lon1, azi1 float64,
) *CGeodesicLine {

	capabilities := STANDARD | DISTANCE_IN

	l := CGeodesicLine{}
	C.geod_lineinit(
		&l.cRepr,
		&geod.cRepr,
		C.double(lat1),
		C.double(lon1),
		C.double(azi1),
		C.unsigned(capabilities),
	)

	return &l
}
