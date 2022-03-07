package geographiclibgo

const (
	_CAP_NONE uint64 = 0
	_CAP_C1   uint64 = 1 << 0

	_CAP_C1p  uint64 = 1 << 1
	_CAP_C2   uint64 = 1 << 2
	_CAP_C3   uint64 = 1 << 3
	_CAP_C4   uint64 = 1 << 4
	_CAP_ALL  uint64 = 0x1F
	_CAP_MASK uint64 = _CAP_ALL
	_OUT_ALL  uint64 = 0x7F80

	// Includes LONG_UNROLL
	OUT_MASK      uint64 = 0xFF80
	EMPTY         uint64 = 0
	LATITUDE      uint64 = 1<<7 | _CAP_NONE
	LONGITUDE     uint64 = 1<<8 | _CAP_C3
	AZIMUTH       uint64 = 1<<9 | _CAP_NONE
	DISTANCE      uint64 = 1<<10 | _CAP_C1
	STANDARD      uint64 = LATITUDE | LONGITUDE | AZIMUTH | DISTANCE
	DISTANCE_IN   uint64 = 1<<11 | _CAP_C1 | _CAP_C1p
	REDUCEDLENGTH uint64 = 1<<12 | _CAP_C1 | _CAP_C2
	GEODESICSCALE uint64 = 1<<13 | _CAP_C1 | _CAP_C2
	AREA          uint64 = 1<<14 | _CAP_C4
	LONG_UNROLL   uint64 = 1 << 15

	// Does not include LONG_UNROLL
	ALL uint64 = _OUT_ALL | _CAP_ALL
)
