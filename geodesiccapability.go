package geographiclibgo

const _CAP_NONE uint64 = 0
const _CAP_C1 uint64 = 1 << 0

const _CAP_C1p uint64 = 1 << 1
const _CAP_C2 uint64 = 1 << 2
const _CAP_C3 uint64 = 1 << 3
const _CAP_C4 uint64 = 1 << 4
const _CAP_ALL uint64 = 0x1F
const _CAP_MASK uint64 = _CAP_ALL
const _OUT_ALL uint64 = 0x7F80

// Includes LONG_UNROLL
const OUT_MASK uint64 = 0xFF80
const EMPTY uint64 = 0
const LATITUDE uint64 = 1<<7 | _CAP_NONE
const LONGITUDE uint64 = 1<<8 | _CAP_C3
const AZIMUTH uint64 = 1<<9 | _CAP_NONE
const DISTANCE uint64 = 1<<10 | _CAP_C1
const STANDARD uint64 = LATITUDE | LONGITUDE | AZIMUTH | DISTANCE
const DISTANCE_IN uint64 = 1<<11 | _CAP_C1 | _CAP_C1p
const REDUCEDLENGTH uint64 = 1<<12 | _CAP_C1 | _CAP_C2
const GEODESICSCALE uint64 = 1<<13 | _CAP_C1 | _CAP_C2
const AREA uint64 = 1<<14 | _CAP_C4
const LONG_UNROLL uint64 = 1 << 15

// Does not include LONG_UNROLL
const ALL uint64 = _OUT_ALL | _CAP_ALL
