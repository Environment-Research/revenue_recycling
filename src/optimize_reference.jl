##############################################################################################################
# Optimization Run for NICE Without Revenue Recycling (Reference Case)
##############################################################################################################

#-------------------------------------------------------------------------------------
# Set up first-stage global optimization for reference case without revenue recycling.
#-------------------------------------------------------------------------------------

println("Beginning global optimization for reference case without revenue recycling.")

# Create NICE objective function with user parameter settings and revenue recycling switched off.
reference_objective = construct_nice_recycle_objective(nice, revenue_recycling = false)

# Create an NLopt optimization object.
reference_opt_global = Opt(global_opt_algorithm, n_objectives)

# Set upper and lower bounds for global CO₂ tax (upper bound = maximum RICE backstop prices).
upper_bound = maximum(rice_backstop_prices, dims=2)[2:(n_objectives+1)]

lower_bounds!(reference_opt_global, zeros(n_objectives))
upper_bounds!(reference_opt_global, upper_bound)

# Set max run time.
maxtime!(reference_opt_global, global_stop_time)

# Set convergence tolerance.
ftol_rel!(reference_opt_global, global_tolerance)

# Set objective function.
max_objective!(reference_opt_global, (x, grad) -> reference_objective(x))

# Optimize model using global optimization algorithm, initialized at a fraction of the upper bound.
reference_max_welfare_global, reference_opt_tax_global, reference_convergence_global = optimize(reference_opt_global, upper_bound ./ 100)

# Calculate full tax time series (just to have to cross-check with local solution).
reference_full_global_opt_tax, reference_global_opt_regional_mitigation = mu_from_tax(reference_opt_tax_global, rice_backstop_prices)

println("Convergence result = ", reference_convergence_global)

#--------------------------------------------------------------------------------------
# Set up second-stage local optimization for reference case without revenue recycling.
#--------------------------------------------------------------------------------------

println("Beginning local optimization for reference case without revenue recycling.")

# Create an NLopt optimization object for second-stage local optimization.
reference_opt = Opt(local_opt_algorithm, n_objectives)

# Set upper and lower bounds for global CO₂ tax (upper bound = maximum RICE backstop prices).
lower_bounds!(reference_opt, zeros(n_objectives))
upper_bounds!(reference_opt, upper_bound)

# Set max run time.
maxtime!(reference_opt, local_stop_time)

# Set convergence tolerance.
ftol_rel!(reference_opt, local_tolerance)

# Set objective function.
max_objective!(reference_opt, (x, grad) -> reference_objective(x))

# Optimize model using local optimization algorithm, initialized at global optimization solution.
reference_max_welfare, reference_opt_tax, reference_convergence = optimize(reference_opt, reference_opt_tax_global)

println("Convergence result = ", reference_convergence)

#-----------------------------------------------------------
# Run Reference Model Under Optimal Policy and Save Results
#-----------------------------------------------------------

# Get optimal regional CO₂ mitigation rates from final optimal tax.
reference_full_opt_tax, reference_opt_regional_mitigation = mu_from_tax(reference_opt_tax, rice_backstop_prices)

# Run model under optimal policy.
set_param!(nice, :emissions, :MIU, reference_opt_regional_mitigation)
set_param!(nice, :nice_recycle, :global_carbon_tax, zeros(n_steps))
run(nice)

# Save results.
if save_results == true
    save_nice_recycle_results(nice, bau_model, reference_full_opt_tax, reference_full_global_opt_tax, output_directory, revenue_recycling=false)
end

println("Reference runs complete.")