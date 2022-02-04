# geographiclib-go
### Author: Nathan McIntosh
A golang port of geographiclib. For an official wrapper of the C functions, see 
[this repo](https://pkg.go.dev/github.com/ruiaylin/pgparser/types/geo/geographiclib).
A version of that wrapper is included in this package, in `geographiclibgo/c_wrapper`.
The more official wrapper does not provide all of the functionality of the original C
library. This repository aims to provide all functionality in both the C wrapper and the
go implementation.

### Aims
 - Mimic the [rust port](https://github.com/georust/geographiclib-rs) as closely as possible
 - Test as extensively as possible
 - Provide useful benchmarks

### Progress
- [X] Geomath
- [X] Geomath tests
- [X] Geodesic Capability constants
- [X] Geodisic Direct
- [X] Geodisic Direct tests
- [X] Geodisic line
- [X] Geodisic line tests
- [X] Geodisic Inverse
- [ ] Geodisic Inverse tests
- [ ] C Wrapper
- [ ] Benchmarks comparing go implementation with C wrapper

#### General Notes
There is a lot of probably unnecessary math in `geomath.go`. First make sure the tests pass, then replace them with standard library functions.
