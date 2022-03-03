package geographiclibgo

// func TestAreaAndPerimeter(t *testing.T) {
// 	testCases := []struct {
// 		desc      string
// 		spheroid  *CGeodesic
// 		points    []s2.Point
// 		area      float64
// 		perimeter float64
// 	}{
// 		{
// 			desc:     "WKT from Wikipedia",
// 			spheroid: CWgs84(),
// 			points: []s2.Point{
// 				s2.PointFromLatLng(s2.LatLngFromDegrees(40, 40)),
// 				s2.PointFromLatLng(s2.LatLngFromDegrees(45, 20)),
// 				s2.PointFromLatLng(s2.LatLngFromDegrees(30, 45)),
// 			},
// 			area:      6.91638769184179e+11,
// 			perimeter: 5.6770339842410665e+06,
// 		},
// 	}

// 	for _, tc := range testCases {
// 		t.Run(tc.desc, func(t *testing.T) {
// 			area, perimeter := tc.spheroid.AreaAndPerimeter(tc.points)
// 			require.Equal(t, tc.area, area)
// 			require.Equal(t, tc.perimeter, perimeter)
// 		})
// 	}
// }

// func TestProject(t *testing.T) {
// 	testCases := []struct {
// 		desc     string
// 		spheroid *CGeodesic
// 		point    s2.LatLng
// 		distance float64
// 		azimuth  float64
// 		project  s2.LatLng
// 	}{
// 		{
// 			desc:     "{0,0} project to 100000, radians(45.0) on WGS84Spheroid",
// 			spheroid: CWgs84(),
// 			point:    s2.LatLng{Lat: 0, Lng: 0},
// 			distance: 100000,
// 			azimuth:  45 * math.Pi / 180.0,
// 			project:  s2.LatLng{Lat: 0.011160897716439782, Lng: 0.011086872969072624},
// 		},
// 	}

// 	for _, tc := range testCases {
// 		t.Run(tc.desc, func(t *testing.T) {
// 			project := tc.spheroid.Project(tc.point, tc.distance, s1.Angle(tc.azimuth))
// 			require.Equal(t, tc.project, project)
// 		})
// 	}
// }
