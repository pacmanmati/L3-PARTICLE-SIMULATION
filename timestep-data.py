import subprocess
import sys
datapoints = int(sys.argv[1])
timestep = 2e-6
for x in range(1, datapoints+1):
    timestep /= 2
    rc = subprocess.call([
        'bin/sol2/sol2', '0', '100.0', str(timestep)] + '0 0 0 0 0 0 4 3 0 0 0 0 0 5 3 4 0 0 0 0 3'.split(' '))
    
