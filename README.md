# geographiclib-go
Author: Nathan McIntosh

### About
A golang port of [geographiclib](https://geographiclib.sourceforge.io/).

### Aims
 - Mimic the [rust port](https://github.com/georust/geographiclib-rs), the [python port](https://pypi.org/project/geographiclib/), and the [java port](https://github.com/geographiclib/geographiclib-java) as closely as possible
 - Test as extensively as possible. Match all the rust, python, and java tests

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
