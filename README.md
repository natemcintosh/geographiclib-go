# geographiclib-go
Author: Nathan McIntosh

## About
A golang port of [geographiclib](https://geographiclib.sourceforge.io/).

## Aims
 - Mimic the [rust port](https://github.com/georust/geographiclib-rs), the [python port](https://pypi.org/project/geographiclib/), and the [java port](https://github.com/geographiclib/geographiclib-java) as closely as possible
 - Test as extensively as possible. Match all the rust, python, and java tests

## Outline:
- [Background: Geodesics on an Ellipsoid](#background-geodesics-on-an-ellipsoid)
- [Short Explanation](#short-explanation-of-library)
- [Long Explanation](#long-explanation-of-library)
- [Library Interface](#the-library-interface)
- [Examples](#examples)
- [Progress](#progress)

## Short Explanation of Library
This library can calculate the following:
- Given a latitude/longitude point, an angle from due North, and a distance, calculate the latitude/longitude of the new point. This is calculated with any function starting with `DirectCalc...()`
- Given two latitude/longitude points, calculate the distance between them, and the angles formed from due North to the line connecting the two points. This is calculated with any function starting with `InverseCalc...()`
- Given a set of points or edges that form a polygon, calculate the area of said polygon. This is done by calling `NewPolygonArea()`, adding the points, and finally calling the `Compute()` method to get both the area and the perimeter of the polygon.
- Given a set of points or edges that form a polyline (a set of connected lines), calculate the perimeter of the line. This is done by calling `NewPolygonArea()` with `is_polyline` set to true, adding the points, and finally calling the `Compute()` method to get the length of the lines.

## Long Explanation of Library
This section is copied from the [python documentation](https://geographiclib.sourceforge.io/Python/doc/geodesics.html)
### Introduction
Consider an ellipsoid of revolution with equatorial radius $a$, polar semi-axis $b$, and flattening $f = (a − b)/a$ . Points on the surface of the ellipsoid are characterized by their latitude $\phi$ and longitude $\lambda$. (Note that latitude here means the geographical latitude, the angle between the normal to the ellipsoid and the equatorial plane).

The shortest path between two points on the ellipsoid at $(\phi_1, \lambda_1)$ and $(\phi_2, \lambda_2)$ is called the geodesic. Its length is $s_{12}$ and the geodesic from point 1 to point 2 has forward azimuths $\alpha_1$ and $\alpha_2$ at the two end points. In this figure, we have $\lambda_{12} = \lambda_2 − \lambda_1$.

![Geodesic](assets/geog.svg)
A geodesic can be extended indefinitely by requiring that any sufficiently small segment is a shortest path; geodesics are also the straightest curves on the surface.

### Solution of Geodesic Problems
Traditionally two geodesic problems are considered:

- the direct problem — given $\phi_1$, $\lambda_1$, $\alpha_1$, $s_{12}$, determine $\phi_2$, $\lambda_2$, and $\alpha_2$; this is solved by the various functions starting with `DirectCalc...`

- the inverse problem — given $\phi_1$, $\lambda_1$, $\phi_2$, $\lambda_2$, determine $s_{12}$, $\alpha_1$, and $\alpha_2$; this is solved by the various functions starting with `InverseCalc...`

### Aditional Properties
The routines also calculate several other quantities of interest

- $S_{12}$ is the area between the geodesic from point 1 to point 2 and the equator; i.e., it is the area, measured counter-clockwise, of the quadrilateral with corners $(\phi_1,\lambda_1), (0,\lambda_1), (0,\lambda_2)$, and $(\phi_2,\lambda_2)$. It is given in $\text{meters}^2$.

- $m_{12}$, the reduced length of the geodesic is defined such that if the initial azimuth is perturbed by $\delta\alpha_1$ (radians) then the second point is displaced by $m_{12}\ \delta\alpha_1$ in the direction perpendicular to the geodesic. $m_{12}$ is given in meters. On a curved surface the reduced length obeys a symmetry relation, $m_{12} + m_{21} = 0$. On a flat surface, we have $m_{12} = s_{12}$.

- $M_{12}$ and $M_{21}$ are geodesic scales. If two geodesics are parallel at point 1 and separated by a small distance $\delta t$, then they are separated by a distance $M_{12} \delta t$ at point 2. $M_{21}$ is defined similarly (with the geodesics being parallel to one another at point 2). $M_{12}$ and $M_{21}$ are dimensionless quantities. On a flat surface, we have $M_{12} = M_{21} = 1$.

- $\sigma_{12}$ is the arc length on the auxiliary sphere. This is a construct for converting the problem to one in spherical trigonometry. The spherical arc length from one equator crossing to the next is always $180 \degree$.

If points 1, 2, and 3 lie on a single geodesic, then the following addition rules hold:

- $s_{13} = s_{12} + s_{23}$
- $\sigma_{13} = \sigma_{12} + \sigma_{23}$
- $S_{13} = S_{12} + S_{23}$
- $m_{13} = m_{12}M_{23} + m_{23}M_{21}$
- $M_{13} = M_{12}M_{23} − (1 − M_{12}M_{21}) m_{23}/m_{12}$
- $M_{31} = M_{32}M_{21} − (1 − M_{23}M_{32}) m_{12}/m_{23}$


### Multiple Shortest Geodesics
The shortest distance found by solving the inverse problem is (obviously) uniquely defined. However, in a few special cases there are multiple azimuths which yield the same shortest distance. Here is a catalog of those cases:

- $\phi_1 = −\phi_2$ (with neither point at a pole). If $\alpha_1 = \alpha_2$, the geodesic is unique. Otherwise there are two geodesics and the second one is obtained by setting 
    - $[\alpha1,\alpha2] \leftarrow [\alpha2,\alpha1]$
    - $[M_{12},M_{21}] \leftarrow [M_{21},M_{12}]$
    - $S12 \leftarrow −S12$

    (This occurs when the longitude difference is near $\pm 180 \degree$ for oblate ellipsoids.)
- $\lambda_2 = \lambda_1 \pm 180 \degree$ (with neither point at a pole). If $\alpha_1 = 0 \degree$ or $\pm 180 \degree$, the geodesic is unique. Otherwise there are two geodesics and the second one is obtained by setting
    - $[\alpha1,\alpha2] \leftarrow [−\alpha1,−\alpha2]$
    - $S12 \leftarrow −S12$
    
    (This occurs when $\phi_2$ is near $−\phi_1$ for prolate ellipsoids.)
- Points 1 and 2 at opposite poles. There are infinitely many geodesics which can be generated by setting $[\alpha1,\alpha2] \leftarrow [\alpha1,\alpha2] + [\delta,−\delta]$, for arbitrary $\delta$. (For spheres, this prescription applies when points 1 and 2 are antipodal.)
- $s_{12} = 0$ (coincident points). There are infinitely many geodesics which can be generated by setting $[\alpha1,\alpha2] \leftarrow [\alpha1,\alpha2] + [\delta,\delta]$, for arbitrary $\delta$.


### Area of a Polygon
The area of a geodesic polygon can be determined by summing $S_{12}$ for successive edges of the polygon ($S_{12}$ is negated so that clockwise traversal of a polygon gives a positive area). However, if the polygon encircles a pole, the sum must be adjusted by $\pm A/2$, where $A$ is the area of the full ellipsoid, with the sign chosen to place the result in $(-A/2, A/2]$.

## The library Interface
This section is mostly copied from the [python documentation](https://geographiclib.sourceforge.io/Python/doc/interface.html), and altered where necessary.
### The Units
All angles (latitude, longitude, azimuth, arc length) are measured in degrees with latitudes increasing northwards, longitudes increasing eastwards, and azimuths measured clockwise from north. For a point at a pole, the azimuth is defined by keeping the longitude fixed, writing $\phi = \pm (90 \degree − \varepsilon)$, and taking the limit $\varepsilon → 0+$

### Geodesic Dictionary
The results returned by `DirectCalc...()` and `InverseCalc...()` are structs that may contain some or all of the following:

- **LatDeg** = $φ_1$, latitude of a point (degrees)
- **LonDeg** = $λ_1$, longitude of a point (degrees)
- **AziDeg** = $α_1$, azimuth of line at point (degrees)
- **DistanceM** = $s_{12}$, distance from 1 to 2 (meters)
- **ArcLengthDeg** = $σ_{12}$, arc length on auxiliary sphere from 1 to 2 (degrees)
- **ReducedLengthM** = $m_{12}$, reduced length of geodesic (meters)
- **M12** = $M_{12}$, geodesic scale at 2 relative to 1 (dimensionless)
- **M21** = $M_{21}$, geodesic scale at 1 relative to 2 (dimensionless)
- **S12M2** = $S_{12}$, area between geodesic and equator (meters2)

### Outmaks and caps
There are a number of functions provided by the API. Each returns a struct with the fields calculated. A method that returns a struct with fewer fields than another means fewer calculations were made. *For faster computation, call methods that return only the fields you need.*

You can also specify your own subset of fields by using the `DirectCalcWithCapabilities()` or `InverseCalcWithCapabilities()` and passing in a `capabilities` (sometimes called `caps`) field. `capabilities` are obtained by `or`’ing together the following values

- `EMPTY`, no capabilities, no output
- `LATITUDE`, compute latitude, lat2
- `LONGITUDE`, compute longitude, lon2
- `AZIMUTH`, compute azimuths, azi1 and azi2
- `DISTANCE`, compute distance, s12
- `STANDARD`, all of the above
- `DISTANCE_IN`, allow s12 to be used as input in the direct problem
- `REDUCEDLENGTH`, compute reduced length, m12
- `GEODESICSCALE`, compute geodesic scales, M12 and M21
- `AREA`, compute area, S12
- `ALL`, all of the above;
- `LONG_UNROLL`, unroll longitudes

`DISTANCE_IN` is a capability provided to the GeodesicLine constructor. It allows the position on the line to specified in terms of distance. (Without this, the position can only be specified in terms of the arc length.) This only makes sense in the caps parameter.

`LONG_UNROLL` controls the treatment of longitude. If it is not set then the lon1 and lon2 fields are both reduced to the range $[−180\degree, 180\degree]$. If it is set, then lon1 is as given in the function call and (lon2 − lon1) determines how many times and in what sense the geodesic has encircled the ellipsoid. This only makes sense in the outmask parameter.

### Restrictions on Parameters
- Latitudes must lie in $[−90\degree, 90\degree]$. Latitudes outside this range are replaced by NaNs.
- The distance $s_{12}$ is unrestricted. This allows geodesics to wrap around the ellipsoid. Such geodesics are no longer shortest paths. However they retain the property that they are the straightest curves on the surface.
- Similarly, the spherical arc length $a_{12}$ is unrestricted.
- Longitudes and azimuths are unrestricted; internally these are exactly reduced to the range $[−180\degree, 180\degree]$; but see also the LONG_UNROLL bit.
- The equatorial radius a and the polar semi-axis b must both be positive and finite (this implies that $−\infty < f < 1)$.
- The flattening $f$ should satisfy $f \in [−1/50,1/50]$ in order to retain full accuracy. This condition holds for most applications in geodesy.

Reasonably accurate results can be obtained for $−0.2 \le f \le 0.2$. Here is a table of the approximate maximum error (expressed as a distance) for an ellipsoid with the same equatorial radius as the WGS84 ellipsoid and different values of the flattening.

| $\text{abs}(f)$ | $\text{error}$  |
|--------|--------|
| 0.003  | 15 nm  |
| 0.01   | 25 nm  |
| 0.02   | 30 nm  |
| 0.05   | 10 μm  |
| 0.1    | 1.5 mm |
| 0.2    | 300 mm |


Here 1 nm = 1 nanometer = $10^{−9}$ m (not 1 nautical mile!)

## Examples
```go
package geographiclibgo_test

import (
	"fmt"

	geodesic "github.com/natemcintosh/geographiclib-go"
)

func Example() {
	// Create a struct representing the Earth
	earth := geodesic.Wgs84()

	// If I start in the middle of Times Square in New York City, and head due West for
	// 1000km, where will I end up?
	NY_lat := 40.757954
	NY_lon := -73.985548
	ended_up_at := earth.DirectCalcLatLon(NY_lat, NY_lon, -90, 1000e3)
	fmt.Println("Ended up at", ended_up_at)

	// Looking on a map, we have ended up in Lafayette Township, a little ways outside
	// Indianapolis

	// Now let's do the inverse calculation. What is the distance from New York City to
	// Chicago?
	CHI_lat := 41.882609
	CHI_lon := -87.621978
	NYC_to_CHI_dist := earth.InverseCalcDistance(NY_lat, NY_lon, CHI_lat, CHI_lon)
	fmt.Printf("Distance from NYC to CHI is %0.1f m\n", NYC_to_CHI_dist)

	// Wyoming is a fairly rectangular state. Taking its four corners as its boundaries,
	// what is its area?
	WY_corners_lat_lon := [][2]float64{
		{40.997958, -111.046710},
		{45.001311, -111.055200},
		{44.997380, -104.057699},
		{41.001432, -104.053249},
	}
	polygon := geodesic.NewPolygonArea(earth, false)
	for _, corner := range WY_corners_lat_lon {
		polygon.AddPoint(corner[0], corner[1])
	}
	polygon_result := polygon.Compute(true, false)
	fmt.Printf("Wyoming area is approximatley %0.0f m^2\n", polygon_result.Area)

	// What is the approximate perimeter of Wyoming?
	fmt.Printf("Wyoming perimeter is approximately %0.0f m\n", polygon_result.Perimeter)

	// But we don't only have to do calculations on Earth. Let's try some on Mars!
	mars_equatorial_radius_m := 3396.2e3
	mars_flattening_factor := 5.0304e-3
	mars := geodesic.NewGeodesic(mars_equatorial_radius_m, mars_flattening_factor)

	// What is the distance from Olympus Mons (18.65, -133.8) to the Curiosity Rover's
	// landing site (-4.47, 137.42)?
	mars_distance := mars.InverseCalcDistance(18.65, -133.8, -4.47, 137.42)
	fmt.Printf("Olympus Mons to Curiosity landing site is %0.0f m", mars_distance)

	// Output:
	// Ended up at {40.15431701948773 -85.75720579845405}
	// Distance from NYC to CHI is 1147311.9 m
	// Wyoming area is approximatley 253282066939 m^2
	// Wyoming perimeter is approximately 2028472 m
	// Olympus Mons to Curiosity landing site is 5348380 m
}
```


## Progress
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
- [ ] Consider switching from having units in names to using github.com/golang/geo/ which uses types for units. It would perhaps involve more allocations; need to do some testing.