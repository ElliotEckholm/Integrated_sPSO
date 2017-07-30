# Integrated_sPSO

Please note this code requires installation of GSL.

Manual on Standard PSO
Author: Elliot Eckholm

Introduction

	Particle Swarm Optimization (PSO) algorithms, have garnered a lot of attention in recent decades for their vast amount of applicability, flexibility and simplicity. PSO is a stochastic algorithm that is based on the behavioral properties of swarms. Since the first developed PSO [add citation] by [add authors] in [add date] there has been a lot of focus on designing the best performing PSO and thus many variants have been produced over the years. In this manual, we will be focusing on the Standard PSO variation. The manual is organized as follows, I will briefly explain how the general PSO Algorithm works and a brief outline of the sPSO Algorithm, as well as how to implement the sPSO code into one’s own code.

You can find the stand alone sPSO code here, please note I am not a contributor to this stand alone code.

It should also be noted that there is a much more comprehensive paper on how sPSO works is located here.


How does the general PSO Algorithm work?

	The general PSO algorithm works by having virtual particles that make up a swarm, fly around a predefined search space and “look” for the most optimal location. If you think of each particle as a bee, then the most optimal location would be the location where there is the highest concentration of flowers for example, and the field of flowers in which the bees fly in would be called the search space. Each bee constantly checks their current location and calculates how “good” that location is. “Good” in this case means the higher the concentration of flowers, the better the location. So each bee keeps track of it’s personal best location it has found so far on it’s own, and is able to constantly compare it’s personal best location with it’s current location to see if it is better. The bee also compares it’s personal best with the swarm’s best location. The swarm’s best location is simply the best of all the personal best locations found by any bee. If a bee is flying around and happens to find a new personal best location, and then compares this new personal best location to the swarm’s best location, then the swarm’s best location is replaced by this bee’s new personal best. So at any given time, one bee’s personal best is also the swarm’s best. The bee’s keep flying around searching for a better location until the swarm’s best has reach a certain criterion, such as having at least thirty flowers in a two inch by two inch area for example. PSO can be generalized to any number of dimensions, not just the four dimensions that bees are limited to, and can be applied to any problem that requires finding an optimal solution. Such as finding the global max of a function in two dimensional space, or even finding the most optimal design parameters for an antenna [add citation].


Standard PSO (sPSO) Pseudo-Code:

Define Search Space (i.e set min and max boundaries in each dimension)
Initialize swarm with random initial positions and velocities inside the search space
FOR each Iteration:
FOR each Particle: Randomize the order of updating the Particles for each iteration
IF Particle is outside search space:  set Particle’s velocity = -0.5 Particle’s velocity or just don’t evaluate fitness
ELSE: Evaluate Fitness
Compare current Particle’s fitness with 2 other randomly selected Particles’ Fitnesses
Find particle with best fitness in this group of 3 particles
IF Fitness of particle with best fitness in the group < Fitness of  lBest:
lBest = This Particle’s position
IF Fitness of Particle’s current position < Fitness of Personal Best:
Personal Best = Particle’s current position
IF Fitness lBest < Fitness of  gBest:
gBest = lBest
Update Particle’s velocity based on the gravitational center 
Update Particle’s position



How to Implement sPSO into Your own code:

	
	




