import ee
ee.Authenticate()
ee.Initialize()

test = ee.Image("NASA/NASADEM_HGT/001").get("title").getInfo()
if test == "NASADEM: NASA NASADEM Digital Elevation 30m": print('\nAuthorisation successful.')
else: print("Issues with authentication, seek assistance...")