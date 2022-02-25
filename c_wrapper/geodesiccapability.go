package c_wrapper

const CAP_NONE uint64 = 0
const CAP_C1 uint64 = 1 << 0

const CAP_C1p uint64 = 1 << 1
const CAP_C2 uint64 = 1 << 2
const CAP_C3 uint64 = 1 << 3
const CAP_C4 uint64 = 1 << 4
const CAP_ALL uint64 = 0x1F
const CAP_MASK uint64 = CAP_ALL
const OUT_ALL uint64 = 0x7F80

// Includes LONG_UNROLL
const OUT_MASK uint64 = 0xFF80
const EMPTY uint64 = 0
const LATITUDE uint64 = 1<<7 | CAP_NONE
const LONGITUDE uint64 = 1<<8 | CAP_C3
const AZIMUTH uint64 = 1<<9 | CAP_NONE
const DISTANCE uint64 = 1<<10 | CAP_C1
const STANDARD uint64 = LATITUDE | LONGITUDE | AZIMUTH | DISTANCE
const DISTANCE_IN uint64 = 1<<11 | CAP_C1 | CAP_C1p
const REDUCEDLENGTH uint64 = 1<<12 | CAP_C1 | CAP_C2
const GEODESICSCALE uint64 = 1<<13 | CAP_C1 | CAP_C2
const AREA uint64 = 1<<14 | CAP_C4
const LONG_UNROLL uint64 = 1 << 15

// Does not include LONG_UNROLL
const ALL uint64 = OUT_ALL | CAP_ALL
