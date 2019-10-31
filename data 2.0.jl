# Data parameters
shift_name = "night2.csv"
load_path = "C:/Users/Oli/Documents/ENGSCI 700/Julia/data/"
save_path = "C:/Users/Oli/Documents/ENGSCI 700/Julia/output/"




# ===== Start of data organisation =====

# Importing data into dataframes
d_df = CSV.read(load_path * "distmatrix_julia.csv")
factory_df = CSV.read(load_path * "factory.csv")
depot_df = CSV.read(load_path * "depot.csv")
farm_df = CSV.read(load_path * "farm_" * shift_name)
tanker_df = CSV.read(load_path * "tanker.csv")

# Finding indexes and number of factories, depots, farms
factories = factory_df.matrixPosition  # List of indices of factories (within d matrix)
depots = depot_df.matrixPosition
farms = farm_df.matrixPosition
tankers = tanker_df.ID

# Making the indexes start at 1 instead of 0
factories = round.(Int, factories + ones(length(factories)))
depots = round.(Int, depots + ones(length(depots)))
farms = round.(Int, farms + ones(length(farms)))
tankers = round.(Int, tankers + ones(length(tankers)))

num_factories = length(factories)
num_depots = length(depots)
num_farms = length(farms)
num_tankers = length(tankers)

#  Distance matrix data (extract releva xant distances)
d = d_df[[factories ; depots ; farms] , [factories ; depots ; farms]]
num_nodes = size(d)[1]

# Update farms list to only include farms in our model
farms = Array(num_factories+num_depots+1 : num_factories+num_depots+num_farms)

# Node data
farm_factory = [zeros(num_factories) ; zeros(num_depots) ; farm_df.factoryPosition]
farm_factory = round.(Int, farm_factory + ones(length(farm_factory)))
milk_volume_mean = [zeros(num_factories) ; zeros(num_depots) ; farm_df.milkPlanned]
milk_volume_mean_original = milk_volume_mean
milk_volume_std = 5.4*sqrt.(milk_volume_mean)
milk_volume_std_original = milk_volume_std
farm_size_limit = [zeros(num_factories) ; zeros(num_depots) ; farm_df.vehicleLimit]
time_windows = [factory_df.startTimeWindow  factory_df.endTimeWindow ;
                depot_df.startTimeWindow    depot_df.endTimeWindow   ;
                farm_df.startTimeWindow     farm_df.endTimeWindow    ]

# Tanker data
tanker_class_vals = tanker_df.class
tanker_sizes = tanker_df.truckSize
tanker_depot = tanker_df.homeDepotPosition
tanker_depot = round.(Int, tanker_depot + ones(length(tanker_depot)))

# Tanker class data
tanker_class_indices = findfirst.(isequal.(unique(tanker_class_vals)), [tanker_class_vals])
tankers = 1:length(tanker_class_vals)
tanker_classes = 1:length(tanker_class_indices)

tanker_class_sizes = tanker_df.truckSize[tanker_class_indices]
tanker_class_capacities = tanker_df.capacity[tanker_class_indices]
tanker_class_count = zeros(length(tanker_class_indices))
for tanker_class in tanker_classes
    tanker_class_count[tanker_class] = count(tanker_class_vals .== tanker_class)
end

loading_time_fixed = tanker_df.fixedLoadTime[tanker_class_indices]
loading_time_variable = tanker_df.serviceTime[tanker_class_indices]

# Time to empty tanker, depends on factory
unloading_time_fixed = factory_df.fixedServiceTime
unloading_time_variable = factory_df.serviceTime

# Tanker travel speed (km/min)
tanker_speed = 80/60 #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
time_matrix = d./tanker_speed

# Cost data
cost_per_km = 3  # Cost per kilometer of travel, for each tanker
cost_per_min = 45/60  # Cost per minute of driver time
cost_late_arrival = 10  # Cost for being late to a farm (per minute)

# The volume of milk collected in a tour is normally distributed. Tours created
# will collect at most this many standard deviations less than the expected milk
# volume, as a way of determining how many farms to include in a tour, and limiting overflow risk
acceptable_overflow_risk = 1.5

# Factor that scales the cost of penalty tankers so they cost more
penalty_coef = 2

latest_tanker_release = 720  # No tanker is allowed to commence a tour after this time

descrete_time_duration = 30  # Duration of each time period we will use to constrain tours to ensure that they don't overlap
earliest_tanker_release = minimum(time_windows[depots, 1])  # No tanker is allowed to commence a tour before this time
latest_tanker_return = maximum(time_windows[depots, 2])  # No tanker is allowed to return to depot after this time
num_time_rows = ceil(Int, (latest_tanker_return - earliest_tanker_release) / descrete_time_duration)  # Number of rows in time constraints model

# We avoid tours longer than this, to avoid factorial computation time of optimal ordering
max_tour_length = 8
perms = []
for i = 1:max_tour_length
    append!(perms, [collect(permutations(1:i;))])
end

# Variable to keep track of solved tight solutions to avoid recomputing them
blacklist = []

# McColl's tour information, outdated
mccoll_tours    = [[[43],
                    [39, 57, 37],
                    [28, 17, 59, 18, 15, 60, 20, 71],
                    [16, 70, 7],
                    [75, 25, 73, 13, 19],
                    [42, 46, 48, 50, 47],
                    [33, 36, 35, 32, 31],
                    [65, 34, 23, 68, 72],
                    [77, 74, 24, 27, 6],
                    [11, 14, 10, 12, 21, 9],
                    [63, 8,  61, 67, 66, 76, 62],
                    [40, 41, 52, 51],
                    [64, 30, 29, 69],
                    [55, 54, 45, 56, 38],
                    [53, 44, 58, 49],
                    [26, 22]],

                   [[17, 22, 18, 24],
                    [31, 38],
                    [23, 20, 19, 21, 12, 14, 7],
                    [32, 26, 27, 47, 33],
                    [36, 35, 37, 50],
                    [6,  16],
                    [48, 15, 43, 46],
                    [34, 30, 49, 29],
                    [42, 39],
                    [41, 40],
                    [28, 25],
                    [45, 44, 10, 11, 13, 9, 8]],

                   [[55, 57],
                    [50, 49, 48, 52, 51],
                    [6,  61, 8,  64, 80, 62],
                    [42, 41, 44, 43],
                    [40, 58],
                    [16, 70, 27, 7],
                    [32, 38, 35, 33],
                    [74, 71, 73, 17, 60, 14, 11, 75],
                    [72, 66, 63, 28, 20],
                    [15, 18, 36, 37, 78, 19],
                    [81, 77, 24, 25],
                    [23, 68, 10, 13, 12, 21, 9],
                    [79, 76, 22, 26, 67],
                    [45, 39, 56],
                    [46, 59, 53, 54, 47],
                    [65, 31, 30, 29, 69]],

                   [[40, 8,  12, 43, 14],
                    [17, 20, 16],
                    [26, 23, 42, 32],
                    [11],
                    [44, 18, 19, 21, 13, 6],
                    [9,  7,  10, 39, 41, 15],
                    [38, 35],
                    [36, 37],
                    [33, 31, 29, 30],
                    [28, 27, 25, 22],
                    [34, 24]]]
mccoll_tour_depots    = [[5, 5, 	4, 	4, 	5, 	5, 	4, 	4, 	5, 	5, 	5, 	5, 	5, 	5, 	5, 	5],
                    [5, 5, 	5, 	5, 	5, 	5, 	5, 	5, 	5, 	5, 	4, 	4],
                    [5, 5, 	5, 	5, 	5, 	5, 	4, 	4, 	4, 	4, 	5, 	5, 	5, 	5, 	5, 	5],
                    [5, 5, 	5, 	5, 	4, 	4, 	5, 	5, 	5, 	5, 	5]]
mccoll_tour_tanker_classes = [[1, 1, 1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	3, 	3, 	3],
                         [1, 1, 1, 	1, 	1, 	1, 	1, 	1, 	3, 	3, 	1, 	1],
                         [1, 1, 1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	1, 	3, 	3],
                         [1, 1, 1, 	1, 	1, 	1, 	1, 	1, 	3, 	3, 	3]]

# McColl's shift information
mccoll_shifts =       [[[43,	1,	39,	57,	37],
                        [28,	17,	59,	18,	15,	60,	20,	71,	2,	16,	70,	7],
                        [75,	25,	73,	13,	19,	2,	42,	46,	48,	50,	47],
                        [33,	36,	35,	32,	31,	2,	65,	34,	23,	68,	72],
                        [77,	74,	24,	27,	6,	2,	11,	14,	10,	12,	21,	9],
                        [63,	8,	61,	67,	66,	76,	62,	2,	40,	41,	52,	51],
                        [64,	30,	29,	69],
                        [55,	54,	45,	56,	38],
                        [53,	44,	58,	49,	1,	26]],
                       [[17,	22,	18,	24,	2,	31,	38],
                        [23,	20,	19,	21,	12,	14,	7,	2,	32,	26,	27,	47,	33],
                        [36,	35,	37,	50,	1,	6,	16],
                        [48,	15,	43,	46,	2,	34,	30,	49,	29],
                        [42,	39,	3,	41,	40],
                        [28,	25,	1,	45,	44,	10,	11,	13,	9,	8]],
                       [[54,	56,	1,	49,	48,	47,	51,	50],
                        [6,	    60,	8,	63,	79,	61,	2,	41,	40,	43,	42],
                        [39,	57,	1,	16,	69,	27,	7],
                        [32,	37,	34,	33,	2,	73,	70,	72,	17,	59,	14,	11,	74],
                        [71,	65,	62,	28,	20,	2,	15,	18,	35,	36,	77,	19],
                        [80,	76,	24,	25,	2,	23,	67,	10,	13,	12,	21,	9],
                        [78,	75,	22,	26,	66,	2,	44,	38,	55],
                        [45,	58,	52,	53,	46],
                        [64,	31,	30,	29,	68]],
                       [[40,	8,	12,	43,	14,	2,	17,	20,	16],
                        [26,	23,	42,	32,	1,	11],
                        [44,	18,	19,	21,	13,	6,	2,	9,	7,	10,	39,	41,	15],
                        [38,	35,	3,	36,	37],
                        [33,	31,	29,	30],
                        [28,	27,	25,	22,	1,	34,	24]]]

mccoll_factories =     [[1, 2, 1, 2, 2, 1, 2, 1, 2],
                        [1, 1, 2, 1, 3, 2],
                        [1, 1, 2, 2, 2, 2, 1, 1, 2],
                        [2, 2, 2, 3, 1, 1]]

mccoll_depots =        [[5,	4,	5,	4,	5,	5,	5,	5,	5],
                        [5,	5,	5,	5,	5,	4],
                        [5,	5,	5,	4,	4,	5,	5,	5,	5],
                        [5,	5,	4,	5,	5,	5]]

mccoll_tanker_classes =[[1, 1,	1,	1,	1,	1,	2,	2,	2],
                        [1,	1,	2,	1,	2,	1],
                        [1,	1,	1,	1,	1,	1,	1,	2,	2],
                        [1,	1,	1,	2,	2,	2]]

mccoll_tanker_classes_orig =[[1, 1,	1,	1,	1,	1,	1,	2,	2],
                        [1,	1,	1,	1,	2,	1],
                        [1,	1,	1,	1,	1,	1,	1,	2,	2],
                        [1,	1,	1,	1,	2,	2]]

## NOTES
# Things to cost:
# cost per distance* = $3/KM
# cost per tanker* = $75/day (assume number of tankers used per day is constant)
# cost per driver* = $45/hour (same w/ overtime)
# cost coefficient for overflow tanker = double for ^
# cost (penalty) per late to farm = $10/minute (choose reasonable number)

# loading time:
# based on speed load profile spreadsheet (per tanker)

# Unloading time:
# based based on depot spreadsheet (per depot)

# Distance matrix format:
# [depots ... factories ... farms]

# Realised milk volume distributions:
# mean = historical realised volume
# std = 5.4 * sqrt(mean [litres])
