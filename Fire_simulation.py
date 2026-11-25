#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.colors import ListedColormap

'''
note to self:
1 barren/burnt
2 still vegetation
3 burning

4 newly burning/ about to be burning
'''

#2. declaring constants
nx, ny = 50, 50 #number of cells the x and y directions
prob_spread = 0.8 #chance to spread to adjacent cells
prob_bare = 0.2 #chance of cell to start as bare patch/ vaccinated
prob_start = 0.003 #chance of cell to start on fire/infected

#3. initializing variables
forest = np.zeros([ny,nx]) + 2 #create an initial grid, set all values to 2
forest[2,1] = 3 #set the center cell to burning/infected

#4. use random number generator to set bare cells and burning cells

isbare = np.random.rand(ny, nx) #creates an array of randomly generated numbers between [0,1)
#test for sanity print(isbare)
isbare = isbare < prob_bare #turn into array of true/false values
forest[isbare] = 1 #change the true values of forest grid isbare to barren (1)
#while prob_bare is set to 0, none of the cells will start as burning

isburning = np.random.rand(ny, nx) #creates an array of randomly generated numbers between [0,1)
isburning = isburning < prob_start #turn into array of true/false values
forest[isburning] = 3 #change the true values of forest grid isburning to burning (3)
#while prob_start is set to 0, none of the cells witll start as burning

#5. allow fires to spread (** for now assuming p_spread = 1.0...)
#before:
print("before\n",forest)
#generating plot
forest_cmap = ListedColormap(['tan','darkgreen', 'firebrick']) #colormap
plt.close('all')
fig, ax = plt.subplots(1,1) #figure
ax.pcolor(forest, cmap=forest_cmap, vmin=1, vmax=3)
plt.pause(2)
plt.show()

#defining neighbors function to return the "neighbors" of a cell
def neighbors(i, j):
    '''
    parameters:
    i: column of cell
    j: row of cell

    returns:
    the indicies of the neighboring cells,
    which includes the cells above, bellow,
    to the left, and to the right of the 
    input cell
    '''
    return (i+1,j), (i-1,j), (i,j+1), (i,j-1) #returns the 4 neighboring cells (n1,n2,n3,n4) (top, bottom, right, left)

game = True
#actual start of the simulation:
while game:
    #loop through in the 'x' direction (j):
    for j in range(nx):
        #loop through in the 'y' direction (i):
       for i in range(ny):
            if forest[i,j] == 3: #if cell is actively burning (3)...
                n1,n2,n3,n4 = neighbors(i,j)
                #checking to see if any of the neighboring cells are still vegetation:
              #first check to see if the index is in the right range for each of the neighboring cells
               #N1
              #print(i,j,n1[0] < nx,)
                if ((n1[1] < nx) and (n1[1] >= 0)) and ((n1[0] < ny) and (n1[0]>=0)):
                    #if n1 i index is in the right range (not an edge case), then we can check it
                    if forest[n1]==2:
                        isspread = np.random.rand() # a random number between [0,1)
                        if isspread < prob_spread:
                            forest[n1] = 4 #set n1 cell to 4 so we can keep track of the newly burning cells
                #N2
                if (((n2[1] < nx) and (n2[1] >= 0)) and ((n2[0] < ny) and (n2[0]>=0))):
                    if forest[n2]==2:
                        isspread = np.random.rand() # a random number between [0,1)
                        if isspread < prob_spread:
                            forest[n2] = 4 #set n2 cell...
                #N3
                if ((n3[1] < nx) and (n3[1] >= 0)) and ((n3[0] < ny) and (n3[0]>=0)):
                    if forest[n3]==2:
                        isspread = np.random.rand() # a random number between [0,1)
                        if isspread < prob_spread:
                            forest[n3]=4 #set n3 cell...
                if ((n4[1] < nx) and (n4[1] >= 0)) and ((n4[0] < ny) and (n4[0]>=0)):
                    if forest[n4]==2:
                        isspread = np.random.rand() # a random number between [0,1)
                        if isspread < prob_spread:
                            forest[n4]=4 #set n4 cell...

    #SET CURRENTLY BURNING CELLS TO BURNT:
    #looping through in the x direction (j)
    for j in range(nx):
        #looping through in the y direcion (i)
        for i in range(ny):
            #if the forest is currently burning, set it to burnt
            if forest[i,j] == 3:
                forest[i,j] = 1 #burnt :P

    #SET NEWLY BURNING CELLS TO BURNING
    #looping throgh in the x direction
    for j in range(nx):
        #looping through in the y direction
        for i in range(ny):
            #if newly burning, set it to burning
            if forest[i,j] == 4:
                forest[i,j] = 3 #burning :O
    print("newly burning\n",forest)

    #generating plot
    forest_cmap = ListedColormap(['tan','darkgreen', 'firebrick']) #colormap
    plt.close('all')
    fig, ax = plt.subplots(1,1) #figure
    ax.pcolor(forest, cmap=forest_cmap, vmin=1, vmax=3)
    plt.pause(2)
    plt.show()


    #continue?
    game = False
    #... for now ...
    #loop through in the 'x' direction (j):
    for j in range(nx):
        #loop through in the 'y' direction (i):
       for i in range(ny):
           if forest[i,j] == 3:
            game = True #jk

#print statistics for the run
print('% Barren from burning = ', (np.sum(forest == 1)/(nx * ny)*100)) #number of burning cells divided by total number of cells
print('% Vegetated = ', (np.sum(forest == 2)/(nx*ny))*100) # number of vegetated cells divided by total number of cells
    

#if you want to save plt.savefig
