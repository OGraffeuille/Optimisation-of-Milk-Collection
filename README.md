This repository contains the code for a set partitioning algorithm developed in Julia, as part of my project on "The Optimisation of Milk Collection". This was my honors project as a Engineering Science student at the University of Auckland. 

The collection of milk is the first stage of the milk supply chain. The process of milk collection using a fleet of tankers can be modelled as a Vehicle Routing Problem (VRP), which produces a schedule when solved. The collection of milk has additional challenges including capacity constraints, multiple farm time windows and stochastic milk supply. Stochastic milk supply can result in an event called overflow, where tankers are not being able to collect all the milk from a farm. This incurs recourse costs.

This set partitioning algorithm is an optimisation model developed to generate robust tanker schedules. It is a multi-stage set partitioning algorithm, which considers stochastic milk supply and the risk of overflow when pricing individual tours in order to develop robust schedules.

###set_partitioning 2.0.jl
This file contains the code for the set partitioning algorithm. This is a three stage process detailed in the associated report. This file also contains functions to carry out the experiments detailed in the associated report. 

###data 2.0.jl
This file contains the code that creates the parameters for the set_partitioning 2.0.jl script. This file gathers data from csv files in the data folder, and performs transformations on these to gather all the parameters that allows the set_partitioning script to run.
 
###data
This folder contains data representing a problem instance of the milk collection problem. 

###output
this folder contains the outputted schedules of the set_partitioning script. Whenever a schedule is generated by this script, it creates two files, a beginning file and a trigger file. These contain information about the schedule in a format that is readable by our AnyLogic simulation, and allows us to validate and easily visualise our generated schedule. The julia_trigger file informs the simulation of how tanker routes start, while the julia_farms file informs the simulation of how tanker routes continue. This folder also contains experiments.xlsx, a file which contains the collated experimental results. 