# geographiclib-go
Author: Nathan McIntosh

### About
A golang port of geographiclib, as well as a C wrapper using cgo. This repository aims to provide
all functionality in both the C wrapper and the go implementation.

### Aims
 - Mimic the [rust port](https://github.com/georust/geographiclib-rs) and the [python port](https://pypi.org/project/geographiclib/) as closely as possible
 - Test as extensively as possible
 - Provide useful benchmarks comparing the go implementation to the c-wrapper version

### Progress
- [X] Go translation
    - [X] Geomath
    - [X] Geomath tests
    - [X] Geodesic Capability constants
    - [X] Geodisic Direct
    - [X] Geodisic Direct tests
    - [X] Geodisic line
    - [X] Geodisic line tests
    - [X] Geodisic Inverse
    - [X] Geodisic Inverse tests
    - [X] Polygon Area
    - [X] Polygon Area tests
- [ ] C Wrapper
    - [X] Geodisic Direct
    - [X] Geodisic Direct tests
    - [ ] Geodisic line
    - [ ] Geodisic line tests
    - [X] Geodisic Inverse
    - [X] Geodisic Inverse tests
    - [ ] Polygon Area
    - [ ] Polygon Area tests
- [ ] Benchmarks comparing go implementation with C wrapper
    - [X] Have benchmarks
    - [ ] Create a summary of them for README
- [X] Add DirectAndInverse interface
    - [ ] Explain it in README
- [ ] Add description of the direct and inverse problem, with pictures
- [ ] Add Examples to README
    - [ ] Examples of direct case on Earth
    - [ ] Examples of direct case on Mars
    - [ ] Examples of inverse case
    - [ ] Examples of calculating an area
- [ ] Consider switching from having units in names to using github.com/golang/geo/
which uses types for units. It would perhaps involve more allocations; need to do some
testing.
