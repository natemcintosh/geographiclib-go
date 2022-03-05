package geographiclibgo

// #cgo CXXFLAGS: -std=c++14
// #cgo LDFLAGS: -lm
//
// #include "geodesic.h"
import "C"

type CGeodesicLine struct {
	cRepr C.struct_geod_geodesicline
	lat1  float64
	lon1  float64
	azi1  float64
}

// NewCGeodesicLine creates a GeodesicLine, with `caps` of STANDARD | DISTANCE_IN
func NewCGeodesicLine(
	geod *CGeodesic,
	lat1, lon1, azi1 float64,
) *CGeodesicLine {

	capabilities := STANDARD | DISTANCE_IN

	// Note that the C.struct_geod_geodesicline is not initialized here, but rather
	// in geod_lineinit
	l := CGeodesicLine{lat1: lat1, lon1: lon1, azi1: azi1}
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

// NewCGeodesicLineWithCapability is the same as NewCGeodesicLine but the user specifies a
// `capabilities` field.
func NewCGeodesicLineWithCapability(
	geod *CGeodesic,
	lat1, lon1, azi1 float64,
	capabilities uint64,
) *CGeodesicLine {

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

// PositionStandard finds the position on the line given s12_m [meters]. It uses the
// STANDARD capabilities, and returns a PositionResultStandard struct
func (g *CGeodesicLine) PositionStandard(s12_m float64) PositionResultStandard {

	var retLat2Deg, retLon2Deg, retAzi2Deg C.double

	C.geod_position(&g.cRepr, C.double(s12_m), &retLat2Deg, &retLon2Deg, &retAzi2Deg)

	return PositionResultStandard{
		Lat1Deg:   g.lat1,
		Lon1Deg:   g.lon1,
		Azi1Deg:   g.azi1,
		Lat2Deg:   float64(retLat2Deg),
		Lon2Deg:   float64(retLon2Deg),
		Azi2Deg:   float64(retAzi2Deg),
		DistanceM: s12_m,
	}
}

// PositionWithCapabilities finds the position on the line given s12_m [meters]. It uses
// whatever capabilities are handed in. Any results not asked for with the capabilities
// will be math.NaN()
func (g *CGeodesicLine) PositionWithCapabilities(s12_m float64, capabilities uint64) PositionResult {

	var retLat2Deg, retLon2Deg, retAzi2Deg, retDistanceM, retReducedLength, retM12, retM21, retS12 C.double

	retArcLenDeg := C.geod_genposition(
		&g.cRepr,
		C.unsigned(capabilities),
		C.double(s12_m),
		&retLat2Deg,
		&retLon2Deg,
		&retAzi2Deg,
		&retDistanceM,
		&retReducedLength,
		&retM12,
		&retM21,
		&retS12,
	)

	outlon1 := g.lon1
	if capabilities&LONG_UNROLL != 0 {
		outlon1 = ang_normalize(g.lon1)
	}

	return PositionResult{
		Lat1Deg:        g.lat1,
		Lon1Deg:        outlon1,
		Azi1Deg:        g.azi1,
		Lat2Deg:        float64(retLat2Deg),
		Lon2Deg:        float64(retLon2Deg),
		Azi2Deg:        float64(retAzi2Deg),
		DistanceM:      float64(retDistanceM),
		ArcLengthDeg:   float64(retArcLenDeg),
		ReducedLengthM: float64(retReducedLength),
		M12:            float64(retM12),
		M21:            float64(retM21),
		S12M2:          float64(retS12),
	}
}
