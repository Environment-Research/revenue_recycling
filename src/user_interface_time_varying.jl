######################################################################################################################
# This file does a two-stage optimization of NICE (with and without revenue recycling) and saves the results. Users
# can change some of the core model parameters and optimization settings.
######################################################################################################################

# Load required Julia packages.
using NLopt

# Load RICE+AIR source code.
include("MimiNICE_recycle_time_varying.jl")

# ------------------------------------------------------------------------------------------------
# NICE + REVENUE RECYCLING PARAMETERS TO CHANGE
# ------------------------------------------------------------------------------------------------

# Pure rate of time preference.
ρ = 0.015

# Elasticity of marginal utility of consumption.
η = 1.5

# Income elasticity of climate damages (1 = proportional to income, -1 = inversely proportional to income).
damage_elasticity = 1.0

# Share of carbon tax revenue that is lost and cannot be recycled (1 = 100% of revenue lost, 0 = nothing lost)
lost_revenue_share = 0.0

# Shares of regional revenues that are recycled globally as international transfers (1 = 100% of revenue recycled globally).
global_recycle_share = zeros(12)

# Share of recycled carbon tax revenue that each region-quintile pair receives (row = region, column = quintile)
quintile_recycle_share = ones(12,5) .* 0.2

# Should the time-varying elasticity values only change across the range of GDP values from the studies?
# true = limit calculations to study gdp range, false = allow calculations for 0 to +Inf GDP.
bound_gdp_elasticity = false

# Quintile income distribution scenario (options = "constant", "SSP1", "SSP2", "SSP3", "SSP4", or "SSP5")
quintile_income_scenario = "constant"

# Type of slope for regression analysis (options are :central, :steeper, :flatter, :percentile)
regression_slope_type = :central

# Value if setting the regression slope type to a specific elasticity percentile value.
elasticity_percentile = 0.9

# ------------------------------------------------------------------------------------------------
# CHOICES ABOUT YOUR ANALYSIS & OPTIMZATION
# ------------------------------------------------------------------------------------------------

# Name of folder to store your results in (a folder will be created with this name).
results_folder = "My Results"

# Do you also want to perform a reference case optimization run (with no revenue recycling)?
run_reference_case = true

# Number of 10-year timesteps to find optimal carbon tax for (after which model assumes full decarbonization).
n_objectives = 45

# Global optimization algorithm (:symbol) used for initial result. See options at http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
global_opt_algorithm = :GN_DIRECT_L

# Local optimization algorithm (:symbol) used to "polish" global optimum solution.
local_opt_algorithm = :LN_SBPLX

# Maximum time in seconds to run global optimization (in case optimization does not converge).
global_stop_time = 600

# Maximum time in seconds to run local optimization (in case optimization does not converge).
local_stop_time = 300

# Relative tolerance criteria for global optimization convergence (will stop if |Δf| / |f| < tolerance from one iteration to the next.)
global_tolerance = 1e-8

# Relative tolerance criteria for global optimization convergence (will stop if |Δf| / |f| < tolerance from one iteration to the next.)
local_tolerance = 1e-12


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# RUN EVERYTHING & SAVE KEY RESULTS
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Create a folder to save results.
output_directory = joinpath((@__DIR__), "..", "results", results_folder)
mkpath(output_directory)

# Load quintile income distribution scenario data.
income_distribution_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", quintile_income_scenario*"_quintile_distributions_consumption.csv")))

# Clean up and organize time-varying income distribution data into required NICE format (time × regions × quintiles).
income_distributions = get_quintile_income_shares(income_distribution_raw)

#---------------------------------------------------------------------------------------------------
# Create a baseline version of the model without CO₂ mitigation.
#---------------------------------------------------------------------------------------------------

# This includes the user-specifications but has no CO₂ mitigation policy (will be used to calculte global CO₂ policy).
bau_model = create_nice_recycle(slope_type=regression_slope_type, percentile=elasticity_percentile)
n_steps   = length(dim_keys(bau_model, :time))

# Run model and extract some generic values needed to set up optimizations below.
run(bau_model)
rice_backstop_prices  = bau_model[:emissions, :pbacktime] .* 1000

# Now set user-specified parameter settings.
update_param!(bau_model, :MIU, zeros(n_steps, 12))
update_param!(bau_model, :damage_elasticity, damage_elasticity)
update_param!(bau_model, :quintile_income_shares, income_distributions)
update_param!(bau_model, :recycle_share, quintile_recycle_share)
update_param!(bau_model, :lost_revenue_share, lost_revenue_share)
update_param!(bau_model, :global_recycle_share, global_recycle_share)
update_param!(bau_model, :global_carbon_tax, zeros(n_steps))
update_param!(bau_model, :rho, ρ)
update_param!(bau_model, :eta, η)

# If selected, allow elasticities to be calculated for all GDP values (not just those observed in studies).
if bound_gdp_elasticity == false
    update_param!(bau_model, :min_study_gdp, 1e-10)
    update_param!(bau_model, :max_study_gdp, +Inf)
end

run(bau_model)



##############################################################################################################
# Optimization Run for NICE + Revenue Recycling
##############################################################################################################

# Create an instance of NICE with revenue recycling.
nice_rev_recycle = create_nice_recycle(slope_type=regression_slope_type, percentile=elasticity_percentile)

# Set user-specified parameter settings.
update_param!(nice_rev_recycle, :damage_elasticity, damage_elasticity)
update_param!(nice_rev_recycle, :quintile_income_shares, income_distributions)
update_param!(nice_rev_recycle, :recycle_share, quintile_recycle_share)
update_param!(nice_rev_recycle, :lost_revenue_share, lost_revenue_share)
update_param!(nice_rev_recycle, :global_recycle_share, global_recycle_share)
update_param!(nice_rev_recycle, :rho, ρ)
update_param!(nice_rev_recycle, :eta, η)

# If selected, allow elasticities to be calculated for all GDP values (not just those observed in studies).
if bound_gdp_elasticity == false
    update_param!(nice_rev_recycle, :min_study_gdp, 1e-10)
    update_param!(nice_rev_recycle, :max_study_gdp, +Inf)
end

#---------------------------------------------------------------------
# Set up first-stage global optimization for NICE + Revenue Recycling.
#---------------------------------------------------------------------

println("Beginning global optimization for NICE + Revenue Recycling.")

# Create NICE + Revenue Recycling objective function with user-specified settings.
recycle_objective = construct_nice_recycle_objective(nice_rev_recycle, revenue_recycling = true)

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
update_param!(nice_rev_recycle, :MIU, recycle_opt_regional_mitigation)
update_param!(nice_rev_recycle, :global_carbon_tax, recycle_full_opt_tax)
run(nice_rev_recycle)

# Save results.
save_nice_recycle_results(nice_rev_recycle, bau_model, recycle_full_opt_tax, recycle_full_global_opt_tax, output_directory, revenue_recycling=true)



##############################################################################################################
# Optimization Run for NICE Without Revenue Recycling (Reference Case)
##############################################################################################################

if run_reference_case == true

    # Create an instance of NICE to run without revenue recycling.
    nice_reference = create_nice_recycle(slope_type=regression_slope_type, percentile=elasticity_percentile)

    # Set user-specified parameter settings.
    update_param!(nice_reference, :damage_elasticity, damage_elasticity)
    update_param!(nice_reference, :quintile_income_shares, income_distributions)
    update_param!(nice_reference, :recycle_share, quintile_recycle_share)
    update_param!(nice_reference, :lost_revenue_share, lost_revenue_share)
    update_param!(nice_reference, :global_recycle_share, global_recycle_share)
    update_param!(nice_reference, :rho, ρ)
    update_param!(nice_reference, :eta, η)

    # If selected, allow elasticities to be calculated for all GDP values (not just those observed in studies).
    if bound_gdp_elasticity == false
        update_param!(nice_reference, :min_study_gdp, 1e-10)
        update_param!(nice_reference, :max_study_gdp, +Inf)
    end

    #-------------------------------------------------------------------------------------
    # Set up first-stage global optimization for reference case without revenue recycling.
    #-------------------------------------------------------------------------------------

    println("Beginning global optimization for reference case without revenue recycling.")

    # Create NICE objective function with user parameter settings and revenue recycling switched off.
    reference_objective = construct_nice_recycle_objective(nice_reference, revenue_recycling = false)

    # Create an NLopt optimization object.
    reference_opt_global = Opt(global_opt_algorithm, n_objectives)

    # Set upper and lower bounds for global CO₂ tax (upper bound = maximum RICE backstop prices).
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
    update_param!(nice_reference, :MIU, reference_opt_regional_mitigation)
    update_param!(nice_reference, :global_carbon_tax, zeros(n_steps))
    run(nice_reference)

    # Save results.
    save_nice_recycle_results(nice_reference, bau_model, reference_full_opt_tax, reference_full_global_opt_tax, output_directory, revenue_recycling=false)
end

println("Model runs complete.")
