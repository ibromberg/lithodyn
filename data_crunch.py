import os # for navigating directories as needed
import numpy as np

filename = "gtopo1" # without the .xyz extension
# initiate x y z 
x = []
y = []
z = []
R = 6.3781e6 # radius of the earth in meters

# open & read file, make each entry in data_all a separate set of data
with open(filename+".xyz","r") as f:
    data_all = f.read()
data_all = data_all.split('\n\n') 

# loop thru each set of data in the file
for i in range(0,len(data_all)):   
    # process the data
    data2process = data_all[i]
    data2process = data2process.split() # removes all whitespace & puts each 
                                        # text as its own entry
     
# split the numbers into coords based on index
lat = data2process[0:len(data2process):3]
long = data2process[1:len(data2process):3]
r0 = data2process[2:len(data2process):3]
           
# turn strings to floats
lat = [float(i) for i in lat]   # latitude
long = [float(i) for i in long] # longitutde
r0 = [float(i) for i in r0]     # elevation

# Lat from 15 - 45
# Long 70 - 110
lattrim = []
longtrim = []
r0trim = []

for i in range(0,len(lat)):
    if (lat[i]>=15 and lat[i]<=45 and long[i]>=70 and long[i]<=110):
        lattrim.append(lat[i])
        longtrim.append(long[i])
        r0trim.append(r0[i])

r0 = r0trim
lat = lattrim
long = longtrim

# add radius of earth to r0
r = [(i+R) for i in r0] # elevation plus radius of earth for total radius

# convert to cartesian
phi = [(90-i) for i in lat]
theta = long




# convert degrees to radians?
phi = [i*np.pi/180 for i in phi]
theta = [i*np.pi/180 for i in phi]

for i in range (0,len(lat)):
    x.append( r[i] * np.sin(phi[i]) * np.cos(theta[i]) )
    y.append( r[i] * np.sin(phi[i]) * np.sin(theta[i]) )
    z.append( r[i] * np.cos(phi[i]) )
    
# convert x y z floats back to strings for writing
x = [str(i) for i in x] # 
y = [str(i) for i in y] # 
z = [str(i) for i in z] # 

newfilename = filename + "_cartesian.txt"

with open(newfilename,"w") as f:    
    for i in range (0,len(x)):
        f.write(x[i]+"\t\t"+y[i]+"\t\t"+z[i]+"\n")