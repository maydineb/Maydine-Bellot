import matplotlib.pyplot as plt
from scipy.integrate import ode #for ODE45
import numpy as np

#Part 1 Implementing the Euler and rk45 solver

#starting by defining the equaion used for competition function:
def lotvolt1(t, N, a, b, c, d):
    '''
    Computes the lokta-volterra equations for two species with populations N at time t
    
    parameters:
    t = time
    N = np array of initial populations for the two species [N1, N2]
    a : reproduction rate of population 1
    b : scales the impact of species 2 on species 1
    c : reproduction rate of population 2
    d : scales the impact of species 1 on species 2

    returns: lotvolt competition equations for the given input values
    '''

    dNdT = np.zeros(2) #initializes an array of 2 variables which will be updated corresponding to N1 and N2 at time t
    dNdT[0] = ( ((a * N[0])*(1 - N[0])) - (b * N[0] * N[1]) )
    dNdT[1] = ( ((c * N[1])*(1 - N[1])) - (d * N[1] * N[0]) )
    return dNdT

#defining equation used for the predator prey function:
def predprey(t, N, a, b, c, d):
    '''
    Computes the predator prey rates for two species with populations N at time t
    parameters:
    t = time
    N = np array of initial populations for the two species [N1, N2]
    a : reproduction rate of population 1
    b : scales the impact of species 2 on species 1
    c : reproduction rate of population 2
    d : scales the impact of species 1 on species 2
    returns: predator prey equations for the two species
    '''
    dNdT = np.zeros(2) #initializes an array of 2 variables which will be updated coreesponding to N1 and. N2 at time t
    dNdT[0] = (( a * N[0]) - (b * N[0] * N[1]))
    dNdT[1] = (-(c * N[1]) + (d * N[0] * N[1]))
    return dNdT

#Implementing the Euler method for competition model
def competition_euler(d_t, yr, a, b, c, d, N):
    '''
    uses eulers method to integrate the competition Lotka-Volterra 
    Population Dyanmics equations for 2 species

    parameters:
    d_t : the time step in years
    yr : the number of years simulation is run for
    a : reproduction rate of population 1
    b : scales the impact of species 2 on species 1
    c : reproduction rate of population 2
    d : scales the impact of species 1 on species 2
    N: vector of the 2 species initial population (N1,N2)
    '''

    t = 0 #sets the time variable used in this function to 0 (units are years)
    years = [t] #first element of years array set to inital time (0)
    pop1 = [N[0]] #first element of population arrray for species 1 is set to initial population (N1)
    pop2 = [N[1]] #first element of population array for species 2 is set to initial population (N2)
    #print(pop1)
    #print(pop2)

    while t < yr: #iterate until the time reaches yr
        t = t + d_t #update time to one more timestep ahead
        years.append(t) #add new time to years array

        #updating N to be the new N with dt added to it
        N = N + (d_t * lotvolt1(t, N, a, b, c, d))
        pop1.append(N[0]) #add current population of species 1 to array of species 1 population 
        pop2.append(N[1]) # same thing but for pop2
        #print(pop1)
        #print(pop2)

    N1 = N[0]#reassigning for ease of plotting for my brain :-)
    N2 = N[1]#^
        
    #print the population of species 1 and 2 after the simulation
    #print(f"Species 1 current population: {N1}")
    #print(f"Species 2 current population: {N2}")

    #create plot of the species population over time with these inital conditions
    plt.title(f"Euler Competition Method Species 1 & 2 Population over {yr} years\ntimestep={d_t}, N1={pop1[0]}, N2={pop2[0]}, a={a}, b={b}, c={c}, d={d}", fontsize=12)
    plt.plot(years, pop1, color= 'black', label='Species 1')
    plt.plot(years, pop2, color = 'blue', label = 'Species 2')
    plt.xlabel('Years')
    plt.ylabel('Population')
    plt.legend()
    plt.grid()
    plt.show()
    print(pop1)
    print(pop2)

#Implementing the euler method for the predador prey model:
def predator_euler(d_t, yr, a, b, c, d, N):
    '''
    uses eulers method to integrate the predator-prey version of the Lotka-Volterra 
    Population Dyanmics equations for 2 species

    parameters:
    d_t : the time step in years
    yr : the number of years simulation is run for
    a : reproduction rate of population 1
    b : scales the impact of species 2 on species 1
    c : reproduction rate of population 2
    d : scales the impact of species 1 on species 2
    N: vector of the 2 species initial population (N1,N2)
    '''

    t = 0 #sets the time variable used in this function to 0 (units are years)
    years = [t] #first element of years array set to inital time (0)
    pop1 = [N[0]] #first element of population arrray for species 1 is set to initial population (N1)
    pop2 = [N[1]] #first element of population array for species 2 is set to initial population (N2)
    
    while t < yr: #iterate until the time reaches yr
        t = t + d_t #update time to one more timestep ahead
        years.append(t) #add new time to years array
        
        N = N + (d_t * predprey(t, N, a, b, c, d)) #updates the current population of species 1 to include next timestep
        pop1.append(N[0]) #add current population of species 1 to array of species 1 population 
        pop2.append(N[1]) #add current population of species 2 to array of specieis 2 population
    
    #storing 2 parts of the N vector as N1 and N2:
    N1 = N[0]
    N2 = N[1]
    #print the population of species 1 and 2 after the simulation
    print(f"Species 1 current population: {N1}")
    print(f"Species 2 current population: {N2}")

    #create plot of the species population over time with these inital conditions
    plt.title(f"Euler Predator Prey Model Species 1 & 2 Population over {yr} years\ntimestep={d_t}, N1={pop1[0]}, N2={pop2[0]}, a={a}, b={b}, c={c}, d={d}", fontsize=12)
    plt.plot(years, pop1, color= 'black', label='Species 1')
    plt.plot(years, pop2, color = 'blue', label = 'Species 2')
    plt.xlabel('Years')
    plt.ylabel('Population')
    plt.legend()
    plt.grid()
    plt.show()
    #print(pop1) #expects a list of population for species 1
    #print(pop2) #expects a list of population for species 1
    #print(years) #expects a list of the time elements in order from 0 to yr at time interval d_t


#rk45 sovler for competition model
def rk45_competition(d_t, yr, a, b, c, d, N):
    '''
    uses rk45 solver for the competition Lotka-Volterra 
    Population Dyanmics equations for 2 species

    parameters:
    d_t : the time step in years
    yr : the number of years simulation is run for
    a : reproduction rate of population 1
    b : scales the impact of species 2 on species 1
    c : reproduction rate of population 2
    d : scales the impact of species 1 on species 2
    N: vector of the 2 species initial population (N1,N2)
    '''
    #Set up times to integrate to and initial conditions:
    dt = d_t# 0.1
    t1 = yr
    Ninitial = np.array(N)
    #Create the integrator object
    r = ode(lotvolt1).set_integrator('dopri5')
    #parameters used for the lab report rk45 example
    r.set_initial_value(Ninitial,0).set_f_params(a,b,c,d)
    #Set up our output
    Ns = [Ninitial] #array of populations of species 1 and species 2
    yrs = []
    yrs.append(0)
    while r.successful() and r.t < t1:
        r.integrate(r.t+dt)
        Ns.append(r.y)
        yrs.append(r.t)

    #plotting the rk45 for both species
    Ns = np.array(Ns)
    plt.title(f"RK45 Competition Model Species 1 & 2 Population over {yr} years\ntimestep={d_t}, N1={Ns[0,0]}, N2={Ns[0,1]}, a={a}, b={b}, c={c}, d={d}", fontsize=12)
    plt.plot(yrs,Ns[:,0], color='black', label='Species 1')
    plt.plot(yrs,Ns[:,1], color='blue', label='Species 2')
    plt.xlabel('Years')
    plt.ylabel('Population')
    plt.grid()
    plt.legend()
    plt.show()
    print(Ns.shape)

#rk45 sovler for predator prey model
def rk45_predator(d_t, yr, a, b, c, d, N):
    '''
    uses rk45 solver for the predator prey 
    Population Dyanmics equations for 2 species

    parameters:
    d_t : the time step in years
    yr : the number of years simulation is run for
    a : reproduction rate of population 1
    b : scales the impact of species 2 on species 1
    c : reproduction rate of population 2
    d : scales the impact of species 1 on species 2
    N: vector of the 2 species initial population (N1,N2)
    '''
    #Set up times to integrate to and initial conditions:
    dt = d_t# 0.1
    t1 = yr
    Ninitial = np.array(N)
    #Create the integrator object
    r = ode(predprey).set_integrator('dopri5')
    #parameters used for the lab report rk45 example
    r.set_initial_value(Ninitial,0).set_f_params(a,b,c,d)
    #Set up our output
    Ns = [Ninitial] #array of populations of species 1 and species 2
    yrs = []
    yrs.append(0)
    while r.successful() and r.t < t1:
        r.integrate(r.t+dt)
        Ns.append(r.y)
        yrs.append(r.t)
 
    #plotting the rk45 for both species
    Ns = np.array(Ns)
    plt.title(f"RK45 Predator Prey Model Species 1 & 2 Population over {yr} years\ntimestep={d_t}, N1={Ns[0,0]}, N2={Ns[0,1]}, a={a}, b={b}, c={c}, d={d}", fontsize=12)
    plt.plot(yrs,Ns[:,0], color='black', label='Species 1')
    plt.plot(yrs,Ns[:,1], color='blue', label='Species 2')
    plt.xlabel('Years')
    plt.ylabel('Population')
    plt.grid()
    plt.legend()
    plt.show()

    #plotting the rk45 for species 1 vs species 2 population
    Ns = np.array(Ns)
    plt.title(f"RK45 Predator Prey Model Species 1 vs Species 2 Population over {yr} years\ntimestep={d_t}, N1={Ns[0,0]}, N2={Ns[0,1]}, a={a}, b={b}, c={c}, d={d}", fontsize=12)
    plt.plot(Ns[:,0], Ns[:,1], color='black')
    plt.xlabel('Species 1 Population')
    plt.ylabel('Species 2 Population')
    plt.grid()
    plt.legend()
    plt.show()

#part 4 Lorenz
#defining the equations

def lorenz(t, N, sigma, R, b):
    '''
    lorenz equations for the given input variables
    parameters:
    R: rayleigh number for the lorenz equations
    N: array of values x, y, z
    returns:
    lorenz equations for the input x, y, and z
    '''

    #assigning each element of the array to x y and z
    x = N[0]
    y = N[1]
    z = N[2]

    ddt = np.zeros(3) #initializes array of 3 values, 0 for now
    ddt[0] = -(sigma * x) + (sigma * y)
    ddt[1] = -(x*z) + (R*x) - y
    ddt[2] = (x*y) - (b*z)
    return ddt


def rk45_lorenz(dt, Ttot, N, R):
    '''
    uese rk45 solver for lorenz equations
    parameters
    dt : the time step
    t: time simulation is run for
    N: vector of the initial positions x, y, and z
    R: raleigh number for the lorenz equations
    '''
    #setting sigma and b values because they are going to be the same for each run
    sigma = 10
    b = 8/3

    dims = [N] #x,y,z positions array
    t = []
    t.append(0)
    t1 = Ttot #total time
    r = ode(lorenz).set_integrator('dopri5')
    r.set_initial_value(N,0).set_f_params(sigma, R, b)
    while r.successful() and r.t < t1:
            r.integrate(r.t+dt)
            dims.append(r.y)
            t.append(r.t)

    dims = np.array(dims)
    x = dims[:,0]
    y = dims[:,1]
    z = dims[:,2]
    #print(np.shape(t)) #checking to make sure shape of t and dims are expected and the same lengths
    #print(np.shape(dims))  
    #plotting
    #fig = plt.figure()
    #ax = fig.add_subplot(projection='3d')
    #ax.plot(x, y, z])
    #print(x)
    #print(y)
    #print(z) #printed these to see if the following matched up with the last number of the array for each variable
    print(f"For initial conditions x={x[0]}, y={y[0]}, z={z[0]}")
    print(f"Final X value: {x[-1]}")
    print(f"Final Y value: {y[-1]}")
    print(f"Final Z value: {z[-1]}")
    print(f"Final time: {t[-1]}") #testing to make sure time went up to expected value



#the following were used to produce the graphs used in the report
#uncomment to reproduce the plots

#PART 1
rk45_competition(1, 100, 1,2,1,3, np.array([0.3,0.6])) #figure 1
competition_euler(1, 100, 1,2,1,3, np.array([0.3,0.6])) #figure 2
#predator_euler(0.05, 100, 1,2,1,3, np.array([0.3,0.6])) #figure 3
#predator_euler(0.005, 100, 1,2,1,3, np.array([0.3,0.6])) #figure 4
#rk45_predator(0.05, 100, 1,2,1,3, np.array([0.3,0.6])) #figure 5

#PART 2
#competition_euler(1, 100, 1, 1, 2, 2, np.array([0.6, 0.4])) #figure 6
#competition_euler(1, 100, 1, 1, 2, 2, np.array([0.2, 0.8])) #figure 7

#PART 3
#the following were used to produce the graphs used in the report
#rk45_predator(0.05, 10, 4, 6, 0, 6, np.array([0.3, 0.6])) #figures 8 and 9
#rk45_predator(0.05, 100, 4, 6, 1, 2, np.array([0.3, 0.6])) #figures 10 and 11
#rk45_predator(0.05, 10, 30, 90, 75, 80, np.array([0.5, 0.5])) #figures 12 and 13

#PART 4
#Used for results in part 4
#rk45_lorenz(5, 100, np.array([1,1,1]), 23) 
#rk45_lorenz(5, 100, np.array([10,1,1]), 23) 
#rk45_lorenz(5, 100, np.array([1,10,1]), 23) 
#rk45_lorenz(5, 100, np.array([1,1,10]), 23)

#rk45_lorenz(5, 100, np.array([1,1,1]), 25) 
#rk45_lorenz(5, 100, np.array([10,1,1]), 25) 
#rk45_lorenz(5, 100, np.array([1,10,1]), 25) 
#rk45_lorenz(5, 100, np.array([1,1,10]), 25) 