##############################################################################################################
# Optimization Run for NICE With Revenue Recycling 
##############################################################################################################

#---------------------------------------------------------------------
# Set up first-stage global optimization for NICE + Revenue Recycling.
#---------------------------------------------------------------------

println("Beginning global optimization for NICE + Revenue Recycling.")

# Create NICE + Revenue Recycling objective function with user-specified settings.
recycle_objective = construct_nice_recycle_objective(nice, revenue_recycling = true)

# Create an NLopt optimization object.
recycle_opt_global = Opt(global_opt_algorithm, n_objectives)

# Set upper and lower bounds for global CO₂ tax (upper bound = maximum RICE backstop prices).
upper_bound = maximum(rice_backstop_prices, dims=2)[2:(n_objectives+1)]

lower_bounds!(recycle_opt_global, zeros(n_objectives))
upper_bounds!(recycle_opt_global, upper_bound)

# Set maximum run time.
maxtime!(recycle_opt_global, global_stop_time)

# Set convergence tolerance.
ftol_rel!(recycle_opt_global, global_tolerance)

# Set objective function.
max_objective!(recycle_opt_global, (x, grad) -> recycle_objective(x))

# Optimize model using global optimization algorithm, initialized at a fraction of the upper bound.
recycle_max_welfare_global, recycle_opt_tax_global, recycle_convergence_global = optimize(recycle_opt_global, upper_bound ./ 100)

println("Convergence result = ", recycle_convergence_global)

# Calculate full tax time series (just to have to cross-check with local solution).
recycle_full_global_opt_tax, recycle_global_opt_regional_mitigation = mu_from_tax(recycle_opt_tax_global, rice_backstop_prices)

#----------------------------------------------------------------------
# Set up second-stage local optimization for NICE + Revenue Recycling.
#----------------------------------------------------------------------

println("Beginning local optimization for NICE + Revenue Recycling.")

# Create an NLopt optimization object for second-stage local optimization.
recycle_opt = Opt(local_opt_algorithm, n_objectives)

# Set upper and lower bounds for global CO₂ tax (upper bound = maximum RICE backstop prices).
lower_bounds!(recycle_opt, zeros(n_objectives))
upper_bounds!(recycle_opt, upper_bound)

# Set maximum run time.
maxtime!(recycle_opt, local_stop_time)

# Set convergence tolerance.
ftol_rel!(recycle_opt, local_tolerance)

# Set objective function.
max_objective!(recycle_opt, (x, grad) -> recycle_objective(x))

# Optimize model using local optimization algorithm, initialized at global optimization solution.
recycle_max_welfare, recycle_opt_tax, recycle_convergence = optimize(recycle_opt, recycle_opt_tax_global)

println("Convergence result = ", recycle_convergence)

#---------------------------------------------------
# Run Model Under Optimal Policy and Save Results
#---------------------------------------------------

# Get optimal regional CO₂ mitigation rates from final optimal tax.
recycle_full_opt_tax, recycle_opt_regional_mitigation = mu_from_tax(recycle_opt_tax, rice_backstop_prices)

# Run model under optimal policy.
update_param!(nice, :MIU, recycle_opt_regional_mitigation)
update_param!(nice, :global_carbon_tax, recycle_full_opt_tax)
run(nice)

# Save results.
if save_results == true
    save_nice_recycle_results(nice, bau_model, recycle_full_opt_tax, recycle_full_global_opt_tax, output_directory, revenue_recycling=true)
end

println("Recycle runs complete.")
