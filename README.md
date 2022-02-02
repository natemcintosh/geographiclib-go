# geographiclib-go
### Author: Nathan McIntosh
A golang port of geographiclib. For a wrapper of the C functions, see [this repo](https://pkg.go.dev/github.com/ruiaylin/pgparser/types/geo/geographiclib).

### Aims
 - Mimic the [rust port](https://github.com/georust/geographiclib-rs) as closely as possible
 - Test as extensively as possible
 - Provide useful benchmarks

### Progress
- [X] Geomath
- [X] Geomath tests
- [X] Geodesic Capability constants
- [ ] Geodisic
- [ ] Geodisic tests
- [X] Geodisic line
- [X] Geodisic line tests

#### General Notes
There is a lot of probably unnecessary math in `geomath.go`. First make sure the tests pass, then replace them with standard library functions.
