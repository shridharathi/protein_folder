"""
Shridhar Athinarayanan

This file includes the main functions required for the Metropolis-Monte Carlo algorithm
for finding the ideal configuration of a polymer using a spring potential energy function.
"""

import math
from math import exp
from random import random,randint,gauss, uniform



"""
FUNCTION: move_an_atom
This function chooses a random atom to move (eventually, it will be checked if this random move lowers the spring 
potential energy of the system).

PARAMETERS - x: array of x-coordinates of the protein
             y: array of y-coordinates of the protein
             z: array of z-coordinates of the protein
             radius: distance to move the particle
             rand_atom: index of the atom to move
             
RETURN - x, y, z-coordinates of potentially moved atom
"""
def move_an_atom(x, y, z, radius, rand_atom):
    theta = random()*360.0

    x_step = radius * math.cos(theta)
    y_step = radius * math.sin(theta)
    z_step = uniform(-radius, radius)

    newX = x[rand_atom] + x_step
    newY = y[rand_atom] + y_step
    newZ = z[rand_atom] + z_step

    return newX, newY, newZ


"""
FUNCTION: get_spring_constant
This function will return the spring constant of a polymer from theory.

PARAMETERS - file: path to file of protein coordinates
             kbt: Boltzmann constant times temperature

RETURN - (3*kbt)/dist_squared: spring constant from literature
"""
def get_spring_constant(x, y, z, kbt):
    begin_x, begin_y, begin_z = x[0], y[0], z[0]
    end_x, end_y, end_z =  x[len(x)-1], y[len(x)-1], z[len(x)-1]
    dist_squared = ((end_x - begin_x) ** 2) + ((end_y - begin_y) ** 2) + ((end_z - begin_z) ** 2)
    return (3*kbt)/dist_squared


"""
FUNCTION: get_spring_potential_energy
This function will return the spring potential energy between two atoms (treat the distance between two atoms as a
spring).

PARAMETERS - x1: x-coordinate of the first atom
             x2: x-coordinate of the second atom
             y1: y-coordinate of the first atom
             y2: y-coordinate of the second atom
             z1: z-coordinate of the first atom
             z2: z-coordinate of the second atom
             k: "spring" constant
             
RETURN - (k/2)*(dist**2): spring potential energy between two atoms
"""
def get_spring_potential_energy(x1, x2, y1, y2, z1, z2, k):
    dist = math.sqrt(((x2-x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))
    return (k/2)*(dist**2)


"""
FUNCTION: avg_dist
This function will return the average bond length given a protein specified by its atom coordinates.

PARAMETERS - x: array of x-coordinates of the protein
             y: array of y-coordinates of the protein
             z: array of z-coordinates of the protein

RETURN - num/len(x): average bond length
"""
def avg_dist(x, y, z):
    num = 0
    for i in range(len(x)-1):
        num += math.sqrt(((x[i+1]-x[i])**2) + ((y[i+1]-y[i])**2) + ((z[i+1]-z[i])**2))
    return num/len(x)


"""
FUNCTION: delta_energy
This function will find the change in spring potential energy which results from moving an atom. Either the atom is
the first or last atom in the chain (only one "spring" has been stretched or compressed) or the atom is in the
middle of the chain (two "springs" will be stretched or compressed). 

PARAMETERS - atom: index of the atom to be moved
             x: array of x-coordinates of the protein
             y: array of y-coordinates of the protein
             z: array of z-coordinates of the protein
             newCoords: coordinates of the randomly chosen atom which has been moved
             k: spring constant
             
RETURN - deltaE: change in energy before and after the atom has been moved
"""
def delta_energy(atom, x, y, z, newCoords, k):
    deltaE = 0
    if atom > 0 and atom < len(x)-1:
        newLeftEnergy = get_spring_potential_energy(newCoords[0], x[atom-1], newCoords[1], y[atom-1], newCoords[2], z[atom-1], k)
        oldLeftEnergy = get_spring_potential_energy(x[atom], x[atom-1], y[atom], y[atom-1], z[atom], z[atom-1], k)
        newRightEnergy = get_spring_potential_energy(newCoords[0], x[atom+1], newCoords[1], y[atom+1], newCoords[2], z[atom+1], k)
        oldRightEnergy = get_spring_potential_energy(x[atom], x[atom+1], y[atom], y[atom+1], z[atom], z[atom+1], k)
        deltaE = (newLeftEnergy + newRightEnergy) - (oldLeftEnergy + oldRightEnergy)
    elif atom == 0:
        newEnergy = get_spring_potential_energy(newCoords[0], x[atom+1], newCoords[1], y[atom+1], newCoords[2], z[atom+1], k)
        oldEnergy = get_spring_potential_energy(x[atom], x[atom+1], y[atom], y[atom+1], z[atom], z[atom+1], k)
        deltaE = newEnergy - oldEnergy
    elif atom == len(x)-1:
        newEnergy = get_spring_potential_energy(newCoords[0], x[atom-1], newCoords[1], y[atom-1], newCoords[2], z[atom-1], k)
        oldEnergy = get_spring_potential_energy(x[atom], x[atom-1], y[atom], y[atom-1], z[atom], z[atom-1], k)
        deltaE = newEnergy - oldEnergy
    return deltaE


"""
FUNCTION: accept_move
This function will return True to accept a move and False to reject a move according to the Metropolis Criterion.
This function was pulled from the assn2 starter code and modified.

PARAMETERS - deltaE: change in energy before and after the atom has been moved
             kbT: Boltzmann constant times temperature
RETURN - True/False: move accepted or not
"""
def accept_move(deltaE, kbT):
    if deltaE < 0:
        return True
    prob = exp(-(deltaE/kbT))
    rand = random()
    if rand < prob:
        return True
    return False


"""
FUNCTION: change_atom_coordinates
This function will change the atom coordinates of a particular atom to new coordinates defined by newCoords

PARAMETERS - x: array of x-coordinates of the protein
             y: array of y-coordinates of the protein
             z: array of z-coordinates of the protein
             newCoords: coordinates of the randomly chosen atom which has been moved
             atom: index of the atom to be moved
"""
def change_atom_coordinates(x, y, z, newCoords, atom):
    x[atom] = round(newCoords[0], 3)
    y[atom] = round(newCoords[1], 3)
    z[atom] = round(newCoords[2], 3)


"""
FUNCTION: monte_carlo
This is the main Monte Carlo-Metropolis method

PARAMETERS - iters: how many steps of the Monte Carlo process desired
             k: spring constant of "springs" between atoms
             kBt: Boltzmann constant * temperature
             rad: distance to move an atom
             x: array of x-coordinates of the protein
             y: array of y-coordinates of the protein
             z: array of z-coordinates of the protein
"""
def monte_carlo(iters, k, kBt, rad, x, y, z):
 #main Monte Carlo-Metropolis looping
    for i in range(iters):
        rand_atom = randint(0, len(x)-1) #choose an atom
        new_coords = move_an_atom(x, y, z, rad, rand_atom) #returns new x, y, z coordinates

        deltaE = delta_energy(rand_atom, x, y, z, new_coords, k)

        if accept_move(deltaE, kBt):
            change_atom_coordinates(x, y, z, new_coords, rand_atom)



