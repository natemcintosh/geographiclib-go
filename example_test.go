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
