import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from numpy.linalg import inv


#snowball functions provided:
def temp_warm(lats_in):
   '''
   Create a temperature profile for modern day "warm" earth.


   Parameters
   ----------
   lats_in : Numpy array
       Array of latitudes in degrees where temperature is required.
       0 corresponds to the south pole, 180 to the north.


   Returns
   -------
   temp : Numpy array
       Temperature in Celcius.
   '''


   # Set initial temperature curve
   T_warm = np.array([-47, -19, -11, 1, 9, 14, 19, 23, 25, 25,
                      23, 19, 14, 9, 1, -11, -19, -47])


   # Get base grid:
   npoints = T_warm.size
   dlat = 180 / npoints  # Latitude spacing.
   lats = np.linspace(dlat/2., 180-dlat/2., npoints)  # Lat cell centers.


   # Fit a parabola to the above values
   coeffs = np.polyfit(lats, T_warm, 2)


   # Now, return fitting sampled at "lats".
   temp = coeffs[2] + coeffs[1]*lats_in + coeffs[0] * lats_in**2


   return temp




def insolation(S0, lats):
   '''
   Given a solar constant (`S0`), calculate average annual, longitude-averaged
   insolation values as a function of latitude.
   Insolation is returned at position `lats` in units of W/m^2.


   Parameters
   ----------
   S0 : float
       Solar constant (1370 for typical Earth conditions.)
   lats : Numpy array
       Latitudes to output insolation. Following the grid standards set in
       the diffusion program, polar angle is defined from the south pole.
       In other words, 0 is the south pole, 180 the north.

   Returns
   -------
   insolation : numpy array
       Insolation returned over the input latitudes.
   '''


   # Constants:
   max_tilt = 23.5   # tilt of earth in degrees


   # Create an array to hold insolation:
   insolation = np.zeros(lats.size)


   #  Daily rotation of earth reduces solar constant by distributing the sun
   #  energy all along a zonal band
   dlong = 0.01  # Use 1/100 of a degree in summing over latitudes
   angle = np.cos(np.pi/180. * np.arange(0, 360, dlong))
   angle[angle < 0] = 0
   total_solar = S0 * angle.sum()
   S0_avg = total_solar / (360/dlong)


   # Accumulate normalized insolation through a year.
   # Start with the spin axis tilt for every day in 1 year:
   tilt = [max_tilt * np.cos(2.0*np.pi*day/365) for day in range(365)]


   # Apply to each latitude zone:
   for i, lat in enumerate(lats):
       # Get solar zenith; do not let it go past 180. Convert to latitude.
       zen = lat - 90. + tilt
       zen[zen > 90] = 90
       # Use zenith angle to calculate insolation as function of latitude.
       insolation[i] = S0_avg * np.sum(np.cos(np.pi/180. * zen)) / 365.


   # Average over entire year; multiply by S0 amplitude:
   insolation = S0_avg * insolation / 365


   return insolation


# PART 1


#declaring constants
So = 1370 #solar flux W/m^2
C = 4.2e6 #heat capacity of seawater J/(m^-3 K^-1)
rho = 1020 #density of seawater kg/m^3
r_e = 6378e3 #earth radius in m
lam = 100 #thermal diffusivity m^2/s
dz = 50 #mixed-layer depth of ocean m
sigma = 5.67e-8 #stefan_boltzman constant J/(m^2 s K^-4)
e = 1.0 #emisivity
a_g = 0.3 #albedo of water/ground
a_i = 0.6 #albedo of ice/snow


N = 18 #spacing (needs to be 18 for 10° each latitude band)
#implicit diffuse solver
I = np.eye(N) #identity matrix I


K = np.zeros([N,N])
B = np.zeros([N,N])


#filling in rows for the K matrix
for i in range(N):
   for j in range(N):
       if i == j: #diagonals
           K[i,i] = -2 #all the diagonals should be -2 as per the pattern
       if j == i-1: #left of diagonal
           K[i, j] = 1 #cell to the left of diagonal should be 1
       if j == i+1: #right of diagonal
           K[i, j] = 1 #cell to the right of diagonal should be 1


#boundary conditions for K matrix
K[0,0] = -2
K[0,1] = 2
K[N-1,-1] = -2
K[N-1,-2] = 2
#print(K)


#filling in rows for the B matrix
for i in range(N):
   for j in range(N):
       if j == i-1: #left of diagonal
           B[i, j] = -1 #cell to the left of diagonal should be -1
       if j == i+1: #right of diagonal
           B[i, j] = 1 #cell to the right of diagonal should be 1
#boundary conditions for B matrix
B[0,0] = 0
B[0,1] = 0
B[N-1, -1] = 0
B[N-1, -2] = 0
#print(B) #chat did it work


dlat = 180/ N #lat grid spacing based on npoints (grid spacing should be 10° with N=18 points)
dy = r_e * dlat * (np.pi / 180) #converting dlat to distance dy
lat = np.linspace(dlat/2., 180-dlat/2.,N) #latitude grid
#print(lat)
dt = (3600*24*365.25)*10 # 10 year timestep
total_t = dt * 10000 # 100,000 years
L = I - (lam*dt*K)/dy**2 # dy**2 to scale for latitude
#print(K)

#INITIAL CONDITION
T = temp_warm(lat)

#BASIC DIFFUSION
#timestep T1
t=0
#initial let T1 = T
T1 = T
while t < total_t:
   t = t + dt #update t to advance one timestep
   #animation for understanding
   '''
   plt.clf()
   plt.plot(lat, T, label='Initial Condition', color='blue')#plotting the inital condition temperature
   plt.plot(lat, T1, label='Basic Diffusion', color='red') #plotting basic diffusion
   #plt.plot(lat, Tf,'--', label='Diffusion + Spherical Correction', color='gold') #plotting basic diffusion
   plt.xlabel('Latitude (°)')
   plt.ylabel('Temperature (°C)')
   plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
   plt.grid()
   plt.legend()
   plt.show()
   plt.pause(0.1)
   '''
   T1 = np.linalg.solve(L, T1)


#ADDING SPHERICAL CORRECTION
#need to account for latitudinal spacing for A, converting lat to radians for the calculation
A = np.pi * ((r_e + dz)**2 - r_e**2) * np.sin(lat*(np.pi/180)) #Axz as defined in lecture (convert lat to radians)
#print(A)

def spherical(lat, T, lam):
    '''
    calculates equilibrium temperature profile with speherical correction

    parameters:
    lat: latitude vector
    T: initial condition temperature profile
    lam: thermal diffusivity

    returns:
    Tf: equilibrium temperature profile with spherical correction
    '''
    #time step Tf and average it with spherical correction
    t=0
    #initial let Tf = T
    Tf = T
    while t < total_t:
        t = t + dt #update t to advance one timestep
        ''' animation for understanding
        plt.clf()
        plt.plot(lat, T, label='Initial Condition', color='blue')#plotting the inital condition temperature
        plt.plot(lat, T1, label='Basic Diffusion', color='red') #plotting basic diffusion
        plt.plot(lat, Tf,'--', label='Diffusion + Spherical Correction', color='gold') #plotting basic diffusion
        plt.xlabel('Latitude (°)')
        plt.ylabel('Temperature (°C)')
        plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
        plt.grid()
        plt.legend()
        plt.show()
        plt.pause(0.1)
        '''
        #f changes with each timestep
        f = (lam * dt)/(A*4*(dy**2)) * (np.dot(B,Tf)) * (np.dot(B,A)) #dividng by 4(dy^2) for scaling (as in the lecture slides)
        #print(f)
        #solving for the diffusion with the spherical correction
        Tf = np.linalg.solve(L, (Tf+f))
    return Tf

#SPHERICAL AND RADIATIVE CORRECTION
def sphere_rad(lat, T, lam, rho, C, So, a_g, e):
    '''
    calculates equilibrium temperature profile with speherical and readative correction

    parameters;
    lat: latitude vector
    T: initial condition temperature profile
    lam: thermal diffusivity
    rho: denisty of seawater
    C: heat capacity of seawater
    So: solar flux
    a_g: ground albedo,
    e: emisivity

    returns:
    Tfr: equilibrium temperature profile with the spherical and radiative correction
    '''
    #time step Tf and average it with spherical correction
    t=0
    time_list = []
    #initial let Tfr = T
    Tfr = T
    while t < total_t:
        time_list.append(t)
        t = t + dt #update t to advance one timestep
        #animation for understanding
        '''
        plt.clf()
        plt.plot(lat, T, label='Initial Condition', color='blue')#plotting the inital condition temperature
        plt.plot(lat, T1, label='Basic Diffusion', color='red') #plotting basic diffusion
        plt.plot(lat, Tf, label='Diffusion + Spherical Correction', color='gold') #plotting spherical correction
        plt.plot(lat, Tfr,'--', label='Diffusion + Spherical Correction + Radiative', color='green') #plotting spherical and radiative correction
        plt.xlabel('Latitude (°)')
        plt.ylabel('Temperature (°C)')
        plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
        plt.grid()
        plt.legend()
        plt.show()
        plt.pause(0.1)
        '''
        #f changing with timestep, will then update Tfd to be different
        f = (lam * dt)/(A*4*(dy**2)) * (np.dot(B,Tfr)) * (np.dot(B,A)) #dividng by 4(dy^2) for scaling (as in the lecture slides)
        #print(f)

        #radiative forcing term using ground albedo of 0.3
        d = (dt/(rho * C * dz)) * (insolation(So, lat) * (1-a_g) - (e*sigma * ((Tfr+273.15)**4))) #convert temp to kelvin
        #solving for the diffusion with the spherical correction
        Tfr = np.linalg.solve(L, Tfr + f + d)
    return Tfr




#insol = insolation(So, lat)
#print(t_warm)
def plot_earth(lat, T, T1, Tf, Tfr):

    '''
    plots the equilibrium temperature profile for initial condition, basic diffusion, 
    spherical correction, and spherical+radiative correction

    parameters:
    lat: latitudes
    T: initial conditions temperature profile
    T1: equilibrium temperature profile with basic diffusion
    Tf: equilibrium temperature profile with spherical correction
    Tfr:equilibrium temperature profile with spherical and radiative correction
    '''
    plt.title('Snowball Earth Task 1')
    plt.plot(lat, T, label='Initial Condition', color='blue')#plotting the inital condition temperature
    plt.plot(lat, T1, label='Basic Diffusion', color='red') #plotting basic diffusion
    plt.plot(lat, Tf, label='Diffusion + Spherical Correction', color='gold') #plotting basic diffusion
    plt.plot(lat, Tfr, label='Diffusion + Spherical Correction + Radiative', color='green') #plotting spherical and radiative correction
    plt.xlabel('Latitude (°)')
    plt.ylabel('Temperature (°C)')
    plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
    plt.grid()
    plt.legend()
    plt.show()
    plt.pause(10)
    return


#Task 1 lines to run
Tfr = sphere_rad(lat, T, lam, rho, C, So, a_g, e)
Tf = spherical(lat, T, lam)
#plot_earth(lat, T, T1, Tf, Tfr) #uncomment to reproduce figure 1

#Task 2
#showing effect of ranging lamda values on Tfr
#using values of 0, 25, 50, 75, 100, 125, 150
'''uncomment to reproduce figure 2
plt.title('Snowball earth varying thermal diffusivity from 0-150')
plt.plot(lat, sphere_rad(lat, T, 0, rho, C, So, a_g, e), label='lam=0', color='blue')
plt.plot(lat, sphere_rad(lat, T, 25, rho, C, So, a_g, e), label='lam=25', color='cyan')
plt.plot(lat, sphere_rad(lat, T, 50, rho, C, So, a_g, e), label='lam=50', color='teal')
plt.plot(lat, sphere_rad(lat, T, 75, rho, C, So, a_g, e), label='lam=75', color='orange')
plt.plot(lat, sphere_rad(lat, T, 100, rho, C, So, a_g, e), label='lam=100', color='red')
plt.plot(lat, sphere_rad(lat, T, 125, rho, C, So, a_g, e), label='lam=125', color='magenta')
plt.plot(lat, sphere_rad(lat, T, 150, rho, C, So, a_g, e), label='lam=150', color='grey')
plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
plt.grid()
plt.legend()
plt.show()
plt.pause(10)
'''

#showing effect of ranging emisivity from 0.1-1 on Tfr
#using values of 0, 0.2, 0.4, 0.6, 0.8, 1
'''uncomment to reproduce figure 3
plt.title('Snowball earth varying emisivity from 0.1-1')
plt.plot(lat, sphere_rad(lat, T, lam, rho, C, So, a_g, 0.1), label='e=0.1', color='blue')
plt.plot(lat, sphere_rad(lat, T, lam, rho, C, So, a_g, 0.2), label='e=0.2', color='cyan')
plt.plot(lat, sphere_rad(lat, T, lam, rho, C, So, a_g, 0.4), label='e=0.4', color='teal')
plt.plot(lat, sphere_rad(lat, T, lam, rho, C, So, a_g, 0.6), label='e=0.6', color='orange')
plt.plot(lat, sphere_rad(lat, T, lam, rho, C, So, a_g, 0.8), label='e=0.8', color='red')
plt.plot(lat, sphere_rad(lat, T, lam, rho, C, So, a_g, 1), label='e=1', color='magenta')
plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
plt.grid()
plt.legend()
plt.show()
plt.pause(10)
'''

#trying out some different conditions and comparing to temp_warm (T)
'''uncomment to reproduce figure 4
plt.title('Snowball earth varying emisivity from 0.1-1')
plt.plot(lat, T, label='temp_warm', color='black')
plt.plot(lat, sphere_rad(lat, T, 75, rho, C, So, a_g, 0.4), label='e=0.4, lam = 75', color='blue')
plt.plot(lat, sphere_rad(lat, T, 50, rho, C, So, a_g, 0.3), label='e=0.3, lam=50', color='cyan')
plt.plot(lat, sphere_rad(lat, T, 25, rho, C, So, a_g, 0.6), label='e=0.6, lam=25', color='teal')
plt.plot(lat, sphere_rad(lat, T, 25, rho, C, So, a_g, 0.7), label='e=0.7, lam=25', color='orange')
plt.plot(lat, sphere_rad(lat, T, 5, rho, C, So, a_g, 0.65), label='e=0.65, lam=5', color='red')
plt.plot(lat, sphere_rad(lat, T, 10, rho, C, So, a_g, 0.62), label='e=0.62, lam=10', color='magenta')
plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
plt.grid()
plt.legend()
plt.show()
plt.pause(10)
'''


#Task 3 
#assigning values from the previous task to emisvity and thermal diffusivity
e = 0.7
lam = 25

#setting initial conditions for T_hot and T_cold
# Set initial temperature curve 
T_hot = np.full(18, 60)
T_cold = np.full(18, -60)

# Fit a parabola to the above values
coeffs = np.polyfit(lat, T_hot, 2)

# Now, return fitting sampled at "lats".
temp = coeffs[2] + coeffs[1]*lat + coeffs[0] * lat**2

#redefining our functino to account for dynamic albedo
#SPHERICAL AND RADIATIVE CORRECTION with dynamic albedo
def sphere_rad_dynamic(lat, T, lam, rho, C, So, e):
    '''
    calculates equilibrium temperature  with speherical and readative correction

    parameters;
    lat: latitude vector
    T: initial condition temperature
    lam: thermal diffusivity
    rho: denisty of seawater
    C: heat capacity of seawater
    So: solar flux
    e: emisivity

    returns:
    Tfr: equilibrium temperature with spherical and radiative correction, also using dynamic albedo
    '''
    #time step Tf and average it with spherical correction
    t=0
    time_list = []
    albedo = np.empty_like(T, dtype=float) #declaring variable for albedo which will be changed based on whether ground is frozen or solid
    #initial let Tfr = T
    Tfr = T
    while t < total_t:
        time_list.append(t)
        t = t + dt #update t to advance one timestep

        #f changing with timestep, will then update Tfr to be different
        f = (lam * dt)/(A*4*(dy**2)) * (np.dot(B,Tfr)) * (np.dot(B,A)) #dividng by 4(dy^2) for scaling (as in the lecture slides)
        #print(f)

        #dynamic albedo
        # Update albedo based on conditions:
        loc_ice = Tfr <= -10 #mask for temp values below -10°C
        albedo[loc_ice] = a_i
        albedo[~loc_ice] = a_g

        #radiative forcing term using albedo
        d = (dt/(rho * C * dz)) * (insolation(So, lat) * (1-albedo) - (e*sigma * ((Tfr+273.15)**4))) #convert temp to kelvin
        #solving for the diffusion with the spherical correction
        Tfr = np.linalg.solve(L, Tfr + f + d)
        #animation for understanding
        plt.clf()
        plt.plot(lat, Tfr, color='black')
        plt.xlabel('Latitude (°)')
        plt.ylabel('Temperature (°C)')
        plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
        plt.grid()
        plt.legend()
        plt.show()
        plt.pause(0.1)
        
    return Tfr

#uncomment to reproduce figure 5
'''
plt.title('Dynamic Albedo compared to Warm Earth')
plt.plot(lat, sphere_rad_dynamic(lat, T_hot, lam, rho, C, So, e), label='Dynamic, Hot Earth', color='magenta')
plt.plot(lat, sphere_rad_dynamic(lat, T_cold, lam, rho, C, So, e), label='Dynamic, Cold Earth ', color='blue', ls='--')
plt.plot(lat, T, label='Warm Earth (Initial)', color='black')
plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
plt.grid()
plt.legend()
plt.show()
plt.pause(10)
'''

#starting with initial warm earth conditions and setting albedo to 0.6 to represent a flash freeze

#uncomment to reproduce figure 6
'''
plt.title('"Flash freeze" compared to Warm Earth')
plt.plot(lat, sphere_rad(lat, T, lam, rho, C, So, 0.6, e), label='Flash Freeze', color='cyan')
plt.plot(lat, T, label='Warm Earth (Initial)', color='black')
plt.xlabel('Latitude (°)')
plt.ylabel('Temperature (°C)')
plt.xticks([0,30,60,90,120,150,180],['-90','-60','-30','0','30','60','90']) #making the labels match lattitudes
plt.grid()
plt.legend()
plt.show()
plt.pause(10)
'''

#TASK 4

#initial conditions
gamma =0.4
T_solarforcing = T_cold
albedo=np.empty_like(T, dtype=float)
temps_list = []#initializing list to store the globally averaged temperature values
gamma_list = []

while gamma < 1.4:
    insol = gamma * insolation(So, lat)
    t=0
    while t < total_t:
        t = t + dt #update t to advance one timestep
        #f changing with timestep, will then update Tfr to be different
        f = (lam * dt)/(A*4*(dy**2)) * (np.dot(B,T_solarforcing)) * (np.dot(B,A)) #dividng by 4(dy^2) for scaling (as in the lecture slides)
        #dynamic albedo
        # Update albedo based on conditions:
        loc_ice = T_solarforcing <= -10 #mask for temp values below -10°C
        albedo[loc_ice] = a_i
        albedo[~loc_ice] = a_g
        #radiative forcing term using albedo
        d = (dt/(rho * C * dz)) * (insol * (1-albedo) - (e*sigma * ((T_solarforcing+273.15)**4))) #convert temp to kelvin
        #solving for the diffusion with the spherical correction, using
        T_solarforcing = np.linalg.solve(L, T_solarforcing + f + d)
    gamma_list.append(gamma)
    temps_list.append(np.sum(T_solarforcing*A)/np.sum(A)) #globally averaged weighted temperature 
    gamma = gamma + 0.05

while gamma > 0.4:
    insol = gamma * insolation(So, lat)
    t=0
    while t < total_t:
        t = t + dt #update t to advance one timestep
        #f changing with timestep, will then update T_solarforcing to be different
        f = (lam * dt)/(A*4*(dy**2)) * (np.dot(B,T_solarforcing)) * (np.dot(B,A)) #dividng by 4(dy^2) for scaling (as in the lecture slides)
        #dynamic albedo
        # Update albedo based on conditions:
        loc_ice = T_solarforcing <= -10 #mask for temp values below -10°C
        albedo[loc_ice] = a_i
        albedo[~loc_ice] = a_g
        #radiative forcing term using albedo
        d = (dt/(rho * C * dz)) * (insol * (1-albedo) - (e*sigma * ((T_solarforcing+273.15)**4))) #convert temp to kelvin
        #solving for the diffusion with the spherical correction, using
        T_solarforcing = np.linalg.solve(L, T_solarforcing + f + d)
    gamma_list.append(gamma)
    temps_list.append(np.sum(T_solarforcing*A)/np.sum(A)) #globally averaged weighted temperature 
    gamma = gamma - 0.05

#uncommment to reproduce figure 6
'''
plt.title('Average Global Temperature vs gamma')
plt.plot(gamma_list, temps_list, color='black')
plt.xlabel('Gamma')
plt.ylabel('Temperature (°C)')
plt.grid()
plt.legend()
plt.show()
plt.pause(10)
'''
