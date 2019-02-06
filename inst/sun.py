import grass.script as gs
import numpy as np

# day 91 is april 1

## spring campaign apr 27 to may 3; fall campaign oct 2 to oct 14
days = range(117, 124) + range(275,288)  

startTime = 4
endTime = 21
# magic number 6 means every 10 minutes (i.e., 6 times per hour)
times = np.linspace(startTime, endTime, 6*(endTime-startTime), endpoint=False) # every 10 minutes




# compute the raster and output as tiff
# one day (every 10 min from 4-21) is approx 3 hours and 20 GB, so 60 hrs and 400 GB for the entire thing
for d in days:
    for t in times:
        outfile = "vj_irradiance_d" + str(d) + "_t" + str(int(t)) + "_" + str(int(t % 1 * 100)) + ".tif"
        gs.run_command("r.sun", overwrite=True, elevation="FilledDEM@PERMANENT", horizon_basename = horbase, horizon_step = hstep, day=d, time=t, glob_rad="vj_irradiance_tmp", aspect="vjaspect@PERMANENT", slope="vjslope@PERMANENT", long="vjosa_long@PERMANENT",  civil_time=2)
        gs.run_command("r.out.gdal", flags='c', overwrite=True, input="vj_irradiance_tmp@PERMANENT", output=outfile, format="GTiff", type="Float32", createopt="PROFILE=GeoTIFF")
