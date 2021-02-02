
# Create an instance of FAIR used in optimzation and assign it exogenous forcing scenario from NICE.
fair_opt = create_fair_for_opt()
update_param!(fair_opt, :F_exogenous, interpolate_to_annual(bau_model[:radiativeforcing, :forcoth], 10))

# Calculate indices from NICE corresponding to FAIR annual years
nice_decadal_indices = findall((in)(collect(2005:10:2595)), (collect(2005:2595)))

#--------------------------------------------------------------------------#
#------ Carry out optimization for NICE+FAIR with revenue recycling. ------#
#--------------------------------------------------------------------------#

# Create NICE+FAIR with Revenue Recycling objective function with user-specified settings.
fair_recycle_objective = construct_fair_recycle_objective(nice, n_fair_loops, revenue_recycling = true)

# Set upper and lower bounds for global CO₂ tax (upper bound = maximum RICE backstop prices).
upper_bound = maximum(rice_backstop_prices, dims=2)[2:(n_objectives+1)]

# Create an NLopt optimization object for local optimization.
fair_recycle_opt = Opt(local_opt_algorithm, n_objectives)

# Set upper and lower bounds for global CO₂ tax (upper bound = maximum RICE backstop prices).
lower_bounds!(fair_recycle_opt, zeros(n_objectives))
upper_bounds!(fair_recycle_opt, upper_bound)

# Set maximum run time.
maxtime!(fair_recycle_opt, local_stop_time)

# Set convergence tolerance.
ftol_rel!(fair_recycle_opt, local_tolerance)

# Set objective function.
max_objective!(fair_recycle_opt, (x, grad) -> fair_recycle_objective(x))

# Optimize model using local optimization algorithm.
fair_recycle_max_welfare, fair_recycle_opt_tax, fair_recycle_convergence = optimize(fair_recycle_opt, upper_bound .* 0.1)

# ----- Run Model Under Optimal Policy and Save Results ----- #

# Get optimal regional CO₂ mitigation rates from final optimal tax.
fair_recycle_full_opt_tax, fair_recycle_opt_regional_mitigation = mu_from_tax(fair_recycle_opt_tax, rice_backstop_prices)

# Run NICE model under optimal policy.
update_param!(nice, :MIU, fair_recycle_opt_regional_mitigation)
update_param!(nice, :global_carbon_tax, fair_recycle_full_opt_tax)
run(nice)

# Set placeholder

# Now itrate through optimal policy version using FAIR.
for loops = 1:n_fair_loops

    # Run FAIR with annualized NICE CO₂ emissions to calculate temperature.
    update_param!(fair_opt, :E, interpolate_to_annual(nice[:emissions, :E], 10))
    run(fair_opt)

    # Set FAIR temperature projections in NICE.
    set_param!(nice, :TATM, fair_opt[:temperature, :T][nice_decadal_indices])
    run(nice)
end

# Save results (set global optimzation tax to -9999.99 because only doing local optimization).
save_nice_recycle_results(nice, bau_model, fair_recycle_full_opt_tax, ones(length(fair_recycle_full_opt_tax)) .* -9999.99, output_directory, revenue_recycling=true)

#-----------------------------------------------------------------------------#
#------ Carry out optimization for NICE+FAIR without revenue recycling. ------#
#-----------------------------------------------------------------------------#

# Create NICE+FAIR reference objective function with user-specified settings.
fair_reference_objective = construct_fair_recycle_objective(nice, n_fair_loops, revenue_recycling = false)

# Create an NLopt optimization object for local optimization.
fair_reference_opt = Opt(local_opt_algorithm, n_objectives)

# Set upper and lower bounds for global CO₂ tax (upper bound = maximum RICE backstop prices).
lower_bounds!(fair_reference_opt, zeros(n_objectives))
upper_bounds!(fair_reference_opt, upper_bound)

# Set maximum run time.
maxtime!(fair_reference_opt, local_stop_time)

# Set convergence tolerance.
ftol_rel!(fair_reference_opt, local_tolerance)

# Set objective function.
max_objective!(fair_reference_opt, (x, grad) -> fair_reference_objective(x))

# Optimize model using local optimization algorithm.
fair_reference_max_welfare, fair_reference_opt_tax, fair_reference_convergence = optimize(fair_reference_opt, upper_bound .* 0.1)

# ----- Run Model Under Optimal Policy and Save Results ----- #

# Get optimal regional CO₂ mitigation rates from final optimal tax.
fair_reference_full_opt_tax, fair_reference_opt_regional_mitigation = mu_from_tax(fair_reference_opt_tax, rice_backstop_prices)

# Run NICE model under optimal policy.
update_param!(nice, :MIU, fair_reference_opt_regional_mitigation)
update_param!(nice, :global_carbon_tax, fair_reference_full_opt_tax)
run(nice)

# Now itrate through optimal policy version using FAIR.
for loops = 1:n_fair_loops

    # Run FAIR with annualized NICE CO₂ emissions to calculate temperature.
    update_param!(fair_opt, :E, interpolate_to_annual(nice[:emissions, :E], 10))
    run(fair_opt)

    # Set FAIR temperature projections in NICE.
    update_param!(nice, :TATM, fair_opt[:temperature, :T][nice_decadal_indices])
    run(nice)
end

# Save results (set global optimzation tax to -9999.99 because only doing local optimization).
save_nice_recycle_results(nice, bau_model, fair_reference_full_opt_tax, ones(length(fair_recycle_full_opt_tax)) .* -9999.99, output_directory, revenue_recycling=false)
