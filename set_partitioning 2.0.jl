using Random, Distributions, Combinatorics
using JuMP, GLPK, Gurobi
using CSV, DataFrames

include("data 2.0.jl")

# Prints a 2D matrix into console
function PrintMatrix(M)

    for i = 1:size(M)[1]
        for j = 1:size(M)[2]
            print(M[i,j], " ")
        end
        print("\n")
    end
end

# Finds the norm of a vector
function norm(A)
    return sqrt(sum(A.*A))
end

# Returns an integer between a and b inclusive
function randInt(a, b)
    return Int(floor(a + (b-a+1)*rand()))
end

# Returns a float between a and b
function randFloat(a, b)
    return a + (b-a)*rand()
end

# Returns whether val is in array
function isInArray(val, array)
    for element in array
        if element == val
            return true
        end
    end
    return false
end

# Returns the load time at a farm (depends on tanker type)
function LoadTime(tanker_class, volume)
    time = loading_time_fixed[tanker_class] + volume * loading_time_variable[tanker_class]
    return time
end

# Returns the unloading time at a factory (depends on the factory)
function UnloadTime(factory, volume)
    time = unloading_time_fixed[factory] + volume * unloading_time_variable[factory]
    return time
end

# Analytically calculate the expected cost of a tour
function EvaluateTourCost(tour, factory, depot, tanker_class, start_time, trace=false)

    tour_length = length(tour)

    # Delay start of tour if tanker leaves too early for first farm, or
    # Leave earlier if tanker leaves too late to arrive at first farm on time
    if start_time < time_windows[tour[1], 1] - time_matrix[depot, tour[1]]
        start_time = time_windows[tour[1], 1] - time_matrix[depot, tour[1]]
    elseif start_time > time_windows[tour[1], 2] - time_matrix[depot, tour[1]]
        start_time = time_windows[tour[1], 2] - time_matrix[depot, tour[1]]
    end

    # Variable we will use to keep track of time
    time = start_time

    travel_dist = 0  # Sum of travel distances for main tanker
    driver_time = 0  # Sum of driver time for main tanker
    time_window_lateness = 0  # Sum of late arrivals for breaking time windows

    travel_dist_penalty = 0  # Sum of travel distances for penalty tanker, if overflow occurs
    driver_time_penalty = 0  # Sum of driver time for penalty tanker, if overflow occurs
    time_window_lateness_penalty = 0  # Sum of late arrivals for breaking time windows, if overflow occurs

    tot_wait_time = 0  # Sum of times waiting for time windows

    # Theoretical mean and std of distribution of milk in tanker
    dist_volume_sum = 0
    dist_volume_std = 0

    # Actual milk volume
    volume = 0
    volume_penalty = 0

    p_overflow_cumm = 0
    p_overflow_prev = 0
    prev_node = depot
    prev_farms_set = []
    remaining_nodes_set = copy(tour)

    for node in tour
        append!(prev_farms_set, node)
        dist_volume_sum += milk_volume_mean[node]
        dist_volume_std = norm(milk_volume_std[prev_farms_set])

        #  Probability that tanker will overflow, before this node
        p_overflow_prev = p_overflow_cumm

        #  Probability that tanker will overflow, at this node or before it
        p_overflow_cumm = 1 - cdf(Normal(dist_volume_sum, dist_volume_std), tanker_class_capacities[tanker_class])

        #  Probability that tanker will overflow, at exactly this node
        p_overflow = p_overflow_cumm - p_overflow_prev

        #  Travelling to this node
        travel_dist += (1-p_overflow_prev) * d[prev_node, node]
        driver_time += (1-p_overflow_prev) * time_matrix[prev_node, node]
        travel_dist_penalty += p_overflow_prev * d[prev_node, node]
        driver_time_penalty += p_overflow_prev * time_matrix[prev_node, node]
        time += time_matrix[prev_node, node]

        if isInArray(node, farms)

            #  Checking whether tanker is early and needs to wait for milk collection
            if time < time_windows[node, 1]
                tot_wait_time += (1-p_overflow_prev) * (time_windows[node, 1] - time)
                driver_time += (1-p_overflow_prev) * (time_windows[node, 1] - time)
                time = time_windows[node, 1]
            end

            #  Checking whether tanker is late and needs to be penalised
            if time > time_windows[node, 2]
                time_window_lateness += (1-p_overflow_prev) * (time - time_windows[node, 2])
            end

            #  Simulating ahead to test whether penalty tanker will be late to a farm and needs to be penalised
            temp_time = time
            temp_prev_node = depot
            for temp_node in remaining_nodes_set
                temp_time += time_matrix[temp_prev_node, temp_node]
                if temp_time < time_windows[temp_node, 1]
                    driver_time_penalty += p_overflow * (time_windows[temp_node, 1] - temp_time)
                    temp_time = time_windows[temp_node, 1]
                end
                if temp_time > time_windows[temp_node, 2]
                    time_window_lateness_penalty += p_overflow * (temp_time - time_windows[temp_node, 2])
                end
                temp_time += LoadTime(tanker_class, milk_volume_mean[temp_node])
                if isInArray(temp_node, farms) == false
                    break
                end
            end

            #  Picking up milk
            pickup_volume = min(milk_volume_mean[node], tanker_class_capacities[tanker_class] - volume)
            pickup_volume_penalty = milk_volume_mean[node] - pickup_volume
            volume += pickup_volume
            volume_penalty += pickup_volume_penalty

            load_time = LoadTime(tanker_class, pickup_volume)
            load_time_penalty = LoadTime(tanker_class, pickup_volume_penalty)
            driver_time += (1-p_overflow_prev) * load_time
            driver_time_penalty += p_overflow_cumm * load_time_penalty
            time += load_time

            #  If overflow occurs, send tanker back to factory
            travel_dist += p_overflow * d[node, factory]
            driver_time += p_overflow * time_matrix[node, factory]

            #  If overflow occurs, get new overflow tanker
            travel_dist_penalty += p_overflow * d[depot, node]
            driver_time_penalty += p_overflow * time_matrix[depot, node]

        elseif isInArray(node, factories)

            # Unloading at factory
            unload_time = UnloadTime(node, volume)
            unload_time_penalty = UnloadTime(node, volume_penalty)
            driver_time += unload_time
            driver_time_penalty += p_overflow_cumm * unload_time_penalty
            time += unload_time

            # Resetting milk distribution
            prev_farms_set = []
            dist_volume_sum = 0

            # Resetting actual milk volumes
            volume = 0
            volume_penalty = 0

            p_overflow_cumm = 0
        else
            error_bug
        end

        #  Updating variables
        prev_node = node
        popfirst!(remaining_nodes_set)
    end

    # Travelling back to factory
    travel_dist += (1-p_overflow_cumm) * d[prev_node, factory]
    driver_time += (1-p_overflow_cumm) * time_matrix[prev_node, factory]
    travel_dist_penalty += p_overflow_cumm * d[prev_node, factory]
    driver_time_penalty += p_overflow_cumm * time_matrix[prev_node, factory]

    # Unloading at factory
    unload_time = UnloadTime(factory, volume)
    unload_time_penalty = UnloadTime(factory, volume_penalty)
    driver_time += unload_time
    driver_time_penalty += p_overflow_cumm * unload_time_penalty
    time += unload_time

    # Travelling back to depot
    travel_dist += d[factory, depot]
    driver_time += time_matrix[factory, depot]
    travel_dist_penalty += p_overflow_cumm * d[factory, depot]
    driver_time_penalty += p_overflow_cumm * time_matrix[factory, depot]
    time += time_matrix[factory, depot]

    end_time = time

    #  Check if tanker is too late to return to depot and tour is infeasible
    if end_time > time_windows[depot, 2]
        cost = maxintfloat()
        return cost, start_time, end_time, tot_wait_time
    end

    #  Evaluate total cost
    cost = travel_dist * cost_per_km + driver_time * cost_per_min
    cost += penalty_coef * (travel_dist_penalty * cost_per_km + driver_time_penalty * cost_per_min)
    cost += (time_window_lateness + time_window_lateness_penalty) * cost_late_arrival

    if trace
        println("\nTotal cost for ", tour, " is ", cost)
        println("Travel cost: ",  travel_dist * cost_per_km, " (penalty = ", penalty_coef * travel_dist_penalty * cost_per_km, ")")
        println("Driver cost: ",  driver_time * cost_per_min, " (penalty = ", penalty_coef * driver_time_penalty * cost_per_min, ")")
        println("Lateness: ", time_window_lateness * cost_late_arrival, " (penalty = ", time_window_lateness_penalty * cost_late_arrival, ")")
        println("P(overflow) = ", p_overflow_cumm, " Total waiting time = ", tot_wait_time)
    end

    return cost, start_time, end_time, tot_wait_time
end

# Uses complete enumeration to find the tour (ordering of farms) with the lowest cost
function FindOptimalTour(farm_set, factory, depot, tanker_class, start_time)

    # Initialise variables
    tour_length = length(farm_set)
    tour = farm_set

    optimal_tour = -1
    optimal_cost = maxintfloat()
    optimal_start_time = -1
    optimal_end_time = -1
    optimal_tot_wait_time = -1

    # Iterate through every permutation
    for order in perms[tour_length]

        test_tour = tour[order]
        cost, revised_start_time, end_time, tot_wait_time = EvaluateTourCost(test_tour, factory, depot, tanker_class, start_time)

        # Save tour if it's the best we've found so far
        if cost < optimal_cost
            optimal_tour = test_tour
            optimal_cost = cost
            optimal_start_time = revised_start_time
            optimal_end_time = end_time
            optimal_tot_wait_time = tot_wait_time
        end
    end
    return optimal_tour, optimal_cost, optimal_start_time, optimal_end_time, optimal_tot_wait_time
end

# Uses a GRASP algorithm to return a list of farms to include in the tour
function GenerateFarmSet(tanker_class)

    # Begin at completely random farm
    prev_farm = rand(farms)
    farm_set = [prev_farm]

    milk_sum = milk_volume_mean[prev_farm]
    milk_std = milk_volume_std[prev_farm]

    tour_length = 0

    # Factory used to restrict farm choice if we're partitioning farms by factory
    partition_factory = farm_factory[prev_farm]

    while length(farm_set) < max_tour_length

        # Decide on next farm to visit, randomly, but weighted inversely to distance
        weights = 1 ./ (d[!,prev_farm] + 1e-3*ones(length(d[!,prev_farm])))
        farm = -1

        # Factories and depots can't be added to tour
        for factory in factories
            weights[factory] = 0
        end
        for depot in depots
            weights[depot] = 0
        end

        # Farms already in tour can't be added again
        for farm in farm_set
            weights[farm] = 0
        end

        # Infeasible farms can not be on tour
        for i = 1:length(weights)
            if tanker_class_sizes[tanker_class] > farm_size_limit[i]
                weights[i] = 0
            end
        end

        # Not allowing non-partition tours to be part of solution, if we're partitioning
        if partition
            for i = 1:length(weights)
                if farm_factory[i] != partition_factory
                    weights[i] = 0
                end
            end
        end

        # Randomly select farm
        rand_val = randFloat(0, sum(weights))
        for i = 1:num_nodes
            rand_val -= weights[i]
            if rand_val <= 10^-12
                farm = i
                break
            end
        end

        # Add farm to tour
        append!(farm_set, farm)
        prev_farm = farm

        # If expected milk volume is too large, reject latest farm and return farm set
        milk_sum += milk_volume_mean[farm]
        milk_std = norm(milk_volume_std[farm_set])
        if milk_sum + acceptable_overflow_risk*milk_std > tanker_class_capacities[tanker_class]
            farm_set = farm_set[1:length(farm_set)-1]
            break
        end
    end

    # If we can choose any factory, choose factory based on distance from last farm, otherwise choose partition factory
    factory = -1
    if partition
        factory = partition_factory
    else
        factory_weights = 1 ./ (d[factories,farm_set[end]] + 1e-3*ones(num_factories))
        rand_val = randFloat(0, sum(factory_weights))
        for i = 1:num_factories
            rand_val -= factory_weights[i]
            if rand_val <= 10^-12
                factory = i
                break
            end
        end
    end
    return farm_set, factory
end

# Generates a random tour and the associated cost
function GenerateTour()

    # Begin by randomly choosing a tanker
    tanker = rand(tankers)
    tanker_class = tanker_class_vals[tanker]

    # Depot is dependent on depot
    depot = tanker_depot[tanker]

    # Create random set of farms which will consist our tour
    farm_set, factory = GenerateFarmSet(tanker_class)

    # Generate a round number for departure time of tanker from depot
    start_time = randInt(time_windows[depot, 1], latest_tanker_release)
    start_time -= rem(start_time, descrete_time_duration)

    # Don't solve tour if it already exists
    blacklist_id = vcat(factory, depot, tanker_class, sort(farm_set))
    if isInArray(blacklist_id, blacklist)
        return farm_set, factory, depot, tanker_class, maxintfloat(), -1, -1
    end

    # Compute the optimal tour with this set of farms
    tour, cost, revised_start_time, end_time, tot_wait_time = FindOptimalTour(farm_set, factory, depot, tanker_class, start_time)
    revised_start_time = floor(revised_start_time)

    # Add tour to blacklist if it had no slack (otherwise it might be worth trying different start times)
    if tot_wait_time == 0
        append!(blacklist, [blacklist_id])
    end

    return tour, factory, depot, tanker_class, cost, revised_start_time, end_time
end

# Create set partitioning matrix and corresponding costs
function MakeSetPartitionModel(num_tours)

    # Matrixes used to solve set partitioning IP
    tour_partitions = zeros(Int8, (num_nodes, num_tours))
    time_partitions = zeros(Int8, (num_time_rows, num_tours))
    cost_list = zeros(Float64, num_tours)

    # Lists to keep track of all generated tours
    tour_list = []
    time_list = zeros(Float32, (2, num_tours))
    tanker_class_list = zeros(Int8, num_tours)
    depot_list = zeros(Int8, num_tours)
    factory_list = zeros(Int8, num_tours)

    tour_index = 1
    while tour_index <= num_tours

        # Add column to set partition if feasible
        tour, factory, depot, tanker_class, cost, start_time, end_time = GenerateTour()

        if cost < maxintfloat()

            # Add farms visited in this tour to tour_partitions
            for j = 1:length(tour)
                tour_partitions[tour[j], tour_index] = 1
            end

            # Add time where tanker is busy to time_partitions
            first_time_index = floor(Int, 1 + (start_time - earliest_tanker_release) / descrete_time_duration)# + (tanker_class - 1) * num_time_rows
            last_time_index = floor(Int, 1 + (end_time - earliest_tanker_release) / descrete_time_duration)# + (tanker_class - 1) * num_time_rows
            for j = first_time_index:last_time_index
                time_partitions[j, tour_index] = 1
            end

            # Add cost of tour to relevant matrix
            cost_list[tour_index] = cost

            # Store tour information
            append!(tour_list, [tour])
            time_list[:, tour_index] = [start_time; end_time]
            tanker_class_list[tour_index] = tanker_class
            depot_list[tour_index] = depot
            factory_list[tour_index] = factory

            tour_index += 1
        end

        # Print progress to console at every percentage of completion
        if (mod(tour_index, floor.(num_tours/100)) == 0)
            print("==== Generated ")
            print(tour_index)
            print(" out of ")
            print(num_tours)
            println(" tours.")
        end
    end
    return tour_partitions, time_partitions, tour_list, factory_list, depot_list, tanker_class_list, cost_list, time_list
end

function SolveSetPartitions(tour_partitions, time_partitions, depot_list, tanker_class_list, cost, num_tours)

    # Decide which solver to use
    if solver == "GLPK"
        model = Model(with_optimizer(GLPK.Optimizer))
    else
        model = Model(with_optimizer(Gurobi.Optimizer, TimeLimit=6000, Presolve=0))
    end

    # Inidialise variables, objective and constraints in IP
    @variable(model, 0 <= x[1:num_tours], Int)

    @objective(model, Min, sum(x[j] * cost[j] for j in 1:num_tours))

    RHS = [zeros(Int, num_factories) ; zeros(Int, num_depots) ; ones(Int, num_farms)]
    @constraint(model, [i = 1:num_nodes], sum(tour_partitions[i,j] * x[j] for j in 1:num_tours) == RHS[i])

    @constraint(model, [i = 1:num_time_rows, t in tanker_classes], sum(time_partitions[i,j] * ifelse(tanker_class_list[j] == t,1,0) * x[j] for j in 1:num_tours) <= tanker_class_count[t])

    # Optimise IP
    JuMP.optimize!(model)

    # Collect result
    if (primal_status(model) != MOI.FEASIBLE_POINT)
        return -1
    end
    obj_value = JuMP.objective_value(model)
    println("Objective value: ", obj_value)

    # Determine which tours were visited
    x_value = JuMP.value.(x)
    tours = findall(x_val -> x_val > 1e-3, x_value)
    return tours
end

function SolveTankerAllocation(time_partitions_filtered, depot_list, tanker_class_list, time_list)

    # Decide which solver to use
    if solver == "GLPK"
        model = Model(with_optimizer(GLPK.Optimizer))
    else
        model = Model(with_optimizer(Gurobi.Optimizer, Presolve=0))
    end

    num_accepted_tours = size(time_partitions_filtered)[2]

    @variable(model, 0 <= x[1:num_tankers, 1:num_accepted_tours], Int)  # x(i, j) = 1 if tanker i is selected for tour j
    @variable(model, 0 <= t[1:num_tankers], Int)                        # t(i) = 1 if tanker i has at least one tour

    @objective(model, Min, sum(x) + sum(t))

    # At least one tanker per tour
    @constraint(model, [j = 1:num_accepted_tours], sum(x[i,j] for i in 1:num_tankers) >= 1)

    # At most one time slot per tanker per
    @constraint(model, [i = 1:num_tankers, k = 1:num_time_rows], sum(x[i,j] * time_partitions_filtered[k,j] for j in 1:num_accepted_tours) <= 1)

    # Can only get trucks if that type of truck was assigned to that tour
    @constraint(model, [i = 1:num_tankers, j = 1:num_accepted_tours], x[i,j] <= ifelse(tanker_class_list[j] == tanker_class_vals[i], 1, 0))

    # Enforcing t variable
    @constraint(model, [i = 1:num_tankers], t[i] >= sum(x[i,:]))

    JuMP.optimize!(model)

    obj_value = JuMP.objective_value(model)
    x_value = JuMP.value.(x)

    # Find list of the tours allocated to tanker i, and sort them in order of start time
    tanker_tours_list = []
    for i = 1:num_tankers
        tanker_tours = findall(x_val -> x_val > 1e-3, x_value[i,:])
        order = sort!([1:length(tanker_tours);], by = i -> (time_list[1, tanker_tours[i]]))
        tanker_tours = [tanker_tours[order]]
        append!(tanker_tours_list, tanker_tours)
    end
    return tanker_tours_list
end

# Creates two outputs that contain the information to describe our routes,
# in a format that can be read by the simulation
function CreateOutputFile(tour_list, factory_list, depot_list, tanker_class_list, cost_list, time_list, tanker_tours_list)

    farm_mat = Array{Any,2}(zeros(Int8, (num_farms, 8)))
    trigger_mat = Array{Any,2}(zeros(Int8, (length(tour_list), 4)))

    for tour_index in 1:length(tour_list)

        tour = tour_list[tour_index]
        depot = depot_list[tour_index]
        factory = factory_list[tour_index]

        # Convert from matrix position to simulation ID
        tour -= ones(Int32, length(tour)) * (num_factories + num_depots)
        depot -= num_factories

        last_farm = tour[end]

        for j in 1:length(tour)

            farm = tour[j]
            next_farm = tour[minimum((j+1, length(tour)))]

            # Farm matrix
            farm_mat[farm, 1] = farm - 1
            farm_mat[farm, 2] = farm_df.name[farm]
            farm_mat[farm, 3] = farm_df.x[farm]
            farm_mat[farm, 4] = farm_df.y[farm]
            farm_mat[farm, 5] = farm_df.milkPlanned[farm]
            if farm != last_farm
                farm_mat[farm, 6] = "farm"
                farm_mat[farm, 7] = next_farm - 1
            else
                farm_mat[last_farm, 6] = "factory"
                farm_mat[last_farm, 7] = factory - 1
            end
            farm_mat[farm, 8] = factory - 1
        end
    end

    row = 1
    for tanker in 1:num_tankers

        tour_number = 1
        for tour_index in tanker_tours_list[tanker]

            # Trigger matrix
            trigger_mat[row, 1] = time_list[1, tour_index]
            trigger_mat[row, 2] = tanker - 1
            trigger_mat[row, 3] = tour_list[tour_index][1] - num_factories - num_depots - 1
            trigger_mat[row, 4] = tour_number

            row += 1
            tour_number += 1
        end
    end

    # Save trigger csv file
    trigger_df = DataFrame(trigger_mat)
    names!(trigger_df, Symbol.(["time", "truck", "location", "orderNum"]))
    CSV.write(save_path * "julia_trigger_" * shift_name, trigger_df)

    # Save farms csv file
    farm_output_df = DataFrame(farm_mat)
    names!(farm_output_df, Symbol.(["ID", "name", "x", "y", "milkPlanned", "nextType", "nextID", "factory"]))
    CSV.write(save_path * "julia_farms_" * shift_name, farm_output_df)
end

# Solve a VRP
# - generate 'num tours' tours in set partitioning model
# - solve test_points different times at equal intervals of tours, as tours are being generated
function Optimise(num_tours, test_points=1)

    global blacklist
    blacklist = []

    # Cost of solution at each test point
    tour_costs = zeros(Float32, test_points)
    shift_costs = zeros(Float32, test_points)

    # Time elapsed from start at each test point
    times = zeros(Float32, test_points)
    tic = time()

    # Total number of tours generated at each test point
    tour_nums = Array(round.(Int, num_tours * Array(1:test_points) / test_points))

    # Matrixes used to solve set partitioning IP
    all_tour_partitions = zeros(Int8, (num_nodes, 0))
    all_time_partitions = zeros(Int8, (num_time_rows, 0))
    all_tour_list = []
    all_tanker_class_list = zeros(Int8, 0)
    all_depot_list = zeros(Int8, 0)
    all_factory_list = zeros(Int8, 0)
    all_cost_list = zeros(Float64, 0)
    all_time_list = zeros(Float64, (2, 0) )

    for run_num = 1:test_points

        println("\n\n\n##### Beginning optimisation!")

        ticc = time()

        println("##### Generating matrix information!")
        tours_to_generate = round(Int, num_tours / test_points)
        tour_partitions, time_partitions, tour_list, factory_list, depot_list, tanker_class_list, cost_list, time_list = MakeSetPartitionModel(tours_to_generate)

        tocc = time() - ticc

        # Adding new partitions to problem
        all_tour_partitions = hcat(all_tour_partitions, tour_partitions)
        all_time_partitions = hcat(all_time_partitions, time_partitions)
        append!(all_tour_list, tour_list)
        all_tanker_class_list = vcat(all_tanker_class_list, tanker_class_list)
        all_depot_list = vcat(all_depot_list, depot_list)
        all_factory_list = vcat(all_factory_list, factory_list)
        all_cost_list = vcat(all_cost_list, cost_list)
        all_time_list = hcat(all_time_list, time_list)


        println("##### Solving with ", solver, "!")
        tour_indices = SolveSetPartitions(all_tour_partitions, all_time_partitions, all_depot_list, all_tanker_class_list, all_cost_list, tour_nums[run_num])


        if (tour_indices == -1)
            tour_costs[run_num] = -1
            shift_costs[run_num] = -1
            times[run_num] = tocc
        else


            # Filter chosen tours
            time_partitions_filtered = all_time_partitions[:, tour_indices]

            cost_list = all_cost_list[tour_indices]
            tour_list = all_tour_list[tour_indices]
            time_list = all_time_list[:,tour_indices]
            tanker_class_list = all_tanker_class_list[tour_indices]
            depot_list = all_depot_list[tour_indices]
            factory_list = all_factory_list[tour_indices]

            println("##### Solving tanker allocations!")
            tanker_tours_list = SolveTankerAllocation(time_partitions_filtered, depot_list, tanker_class_list, time_list)

            # Creating solution file
            CreateOutputFile(tour_list, factory_list, depot_list, tanker_class_list, cost_list, time_list, tanker_tours_list)

            # Creating shifts from tour solution
            shift_tours = []
            shift_factories = []
            shift_times = zeros(Float32,(2,0))
            for tanker_tours in tanker_tours_list
                if length(tanker_tours) > 1
                    shift = tour_list[tanker_tours[1]]
                    factory = factory_list[tanker_tours[1]]
                    for tour = 2:length(tanker_tours)
                        shift = vcat(shift, factory)
                        shift = vcat(shift, tour_list[tanker_tours[tour]])
                        factory = factory_list[tanker_tours[tour]]
                    end
                    start_time = time_list[1, tanker_tours[1]]
                    end_time = time_list[2, tanker_tours[length(tanker_tours)]]
                    append!(shift_tours, [shift])
                    append!(shift_factories, factory)
                    shift_times = hcat(shift_times, [start_time ; end_time])
                elseif length(tanker_tours) == 1
                    append!(shift_tours, tour_list[tanker_tours])
                    append!(shift_factories, factory_list[tanker_tours])
                    shift_times = hcat(shift_times, time_list[:,tanker_tours])
                else
                    append!(shift_tours, [[]])
                    append!(shift_factories, -1)
                    shift_times = hcat(shift_times, [-1 ; -1])
                end
            end
            shift_depots = tanker_depot
            shift_tanker_classes = tanker_class_vals

            # Printing and evaluating tour and shift solutions, respectively
            tour_cost = CostTours(tour_list, factory_list, depot_list, tanker_class_list, time_list[1,:], true)
            shift_cost = CostTours(shift_tours, shift_factories, shift_depots, shift_tanker_classes, shift_times[1,:], true)

            times[run_num] = tocc
            tour_costs[run_num] = tour_cost
            shift_costs[run_num] = shift_cost
        end
    end

    for i = 1:test_points
        println(times[i], ", ", tour_costs[i], ", ", shift_costs[i])
    end
    return tour_costs[end], shift_costs[end]
end

# Costs a set of tours, use factories = -1 or start_times = -1 for defaults
function CostTours(tours, factories, depots, tanker_classes, start_times, trace=false)

    # Ensure that we cost with actual mean and standard deviation
    global milk_volume_std
    global milk_volume_std_original
    global milk_volume_mean
    global milk_volume_mean_original

    milk_volume_std_saved = milk_volume_std
    milk_volume_std = milk_volume_std_original
    milk_volume_mean_saved = milk_volume_mean
    milk_volume_mean = milk_volume_mean_original

    costs = zeros(Float32, length(tours))

    for i = 1:length(tours)

        tour = tours[i]

        if length(tour) > 0

            depot = depots[i]
            tanker_class = tanker_classes[i]

            if factories == -1
                factory = farm_factory[tour[end]]
            else
                factory = factories[i]
            end
            if start_times == -1
                start_time = 600
            else
                start_time = start_times[i]
            end

            cost, start_time, end_time = EvaluateTourCost(tour, factory, depot, tanker_class, start_time, trace)

            costs[i] = cost


            if trace
                println("Depot: ", depot, ", Factory: ", factory, ", Tanker class: ", tanker_class, ", Start time: ", start_time)
            end
        end
    end
    if trace
        println("Total cost = ", sum(costs))
    end

    milk_volume_std = milk_volume_std_saved
    milk_volume_mean = milk_volume_mean_saved

    return sum(costs)
end

# Experiment: Tests the performance of solutions found when multiplying
# the variation in milk volume by coefficients
function TestVaryingStd(num_tours, coefs)

    # Storing real milk distribution parameters to restore them later
    global milk_volume_std
    global milk_volume_std_original

    tour_costs = []
    shift_costs = []

    # Solving set partition once for each coefficient
    for coef in coefs
        milk_volume_std = milk_volume_std_original * coef
        tour_cost, shift_cost = Optimise(num_tours)

        append!(tour_costs, tour_cost)
        append!(shift_costs, shift_cost)
    end
    for i in 1:length(coefs)
        println(coefs[i], ", ", tour_costs[i], ", ", shift_costs[i])
    end

    # Restoring milk distribution parameters
    milk_volume_std = milk_volume_std_original
end

# Experiment: Tests the performance of solutions found when setting standard deviation to 0
# and adding coefficients of standard deviations to mean
function TestVaryingMean(num_tours, coefs)

    # Storing real milk distribution parameters to restore them later
    global milk_volume_std
    global milk_volume_std_original
    global milk_volume_mean
    global milk_volume_mean_original

    tour_costs = []
    shift_costs = []

    milk_volume_std = milk_volume_std * 0

    # Solving set partition once for each coefficient
    for coef in coefs
        milk_volume_mean = milk_volume_mean_original + coef * milk_volume_std_original
        tour_cost, shift_cost = Optimise(num_tours)

        append!(tour_costs, tour_cost)
        append!(shift_costs, shift_cost)
    end
    for i in 1:length(coefs)
        println(coefs[i], ", ", tour_costs[i], ", ", shift_costs[i])
    end

    # Restoring milk volume parameters
    milk_volume_mean = milk_volume_mean_original
    milk_volume_std = milk_volume_std_original
end

# Experiment: Tests the effect of stochasticity on tour ordering
function GenerateToursTest(num_tours)
    global acceptable_overflow_risk
    global milk_volume_std
    global milk_volume_std_original

    acceptable_overflow_risk = 0

    num_diff_tours = 0

    percs = [[] for i=1:101]


    for i = 1:num_tours

        # Begin by randomly choosing a tanker
        tanker = rand(tankers)
        tanker_class = tanker_class_vals[tanker]

        # Depot is dependent on depot
        depot = tanker_depot[tanker]

        # Create random set of farms which will consist our tour
        farm_set, factory = GenerateFarmSet(tanker_class)

        # Generate a round number for departure time of tanker from depot
        start_time = randInt(time_windows[depot, 1], latest_tanker_release)
        start_time -= rem(start_time, descrete_time_duration)

        # Compute the optimal tour with this set of farms, with no stochasticity
        milk_volume_std = milk_volume_std * 0
        tour_no, cost, revised_start_time, end_time, tot_wait_time = FindOptimalTour(farm_set, factory, depot, tanker_class, start_time)

        # Compute the optimal tour with this set of farms, with stochasticity
        milk_volume_std = milk_volume_std_original
        cost_no, a, b, c = EvaluateTourCost(tour_no, factory, depot, tanker_class, start_time)
        tour_stoch, cost_stoch, revised_start_time, end_time, tot_wait_time = FindOptimalTour(farm_set, factory, depot, tanker_class, start_time)

        # Compare orderings and print results
        if tour_no != tour_stoch
            num_diff_tours += 1
        end
        milk_perc_unrounded = 100 * sum(milk_volume_mean[farm_set]) / tanker_class_capacities[tanker_class]
        milk_perc = round(Int, milk_perc_unrounded)
        cost_perc = 100 * (cost_no - cost_stoch) / cost_no
        if (milk_perc_unrounded > 80)
            println(milk_perc_unrounded, ", ", cost_perc)
        end
        append!(percs[milk_perc], cost_perc)
    end

    for list in percs
        if length(list) > 0
            println(mean(list))
        else
            println(0)
        end
    end

    println(100 * num_diff_tours / num_tours)

    milk_volume_std = milk_volume_std_original
    acceptable_overflow_risk = 1.5

    return 0
end


# Global optimisation parameters
Random.seed!(0)
solver = "Gurobi"
num_tours = 2500 # How many tours to generate before solving IP
partition = false # Whether to partition farms per fectory (otherwise keep the problem relaxed)

# Run set partitioning model
Optimise(num_tours)
