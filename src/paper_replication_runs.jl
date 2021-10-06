######################################################################################################################
# This file replicates the NICE model runs (with and without revenue recycling) presented in Budolfson et al. (2021),
# "Climate Action With Revenue Recycling Has Benefits For Poverty, Inequality, And Wellbeing," Nature Climate Change.
######################################################################################################################

# Activate the project for the paper and make sure all packages we need
# are installed.
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# Load required Julia packages.
using Interpolations, MimiFAIR13, NLopt

# Load NICE+recycling source code.
include("MimiNICE_recycle_time_varying.jl")

# ------------------------------------------------------------------------------------------------
# NICE + REVENUE RECYCLING PARAMETERS TO CHANGE
# ------------------------------------------------------------------------------------------------

# Pure rate of time preference.
ρ =  0.015

# Elasticity of marginal utility of consumption.
η =  1.5

# Income elasticity of climate damages (1 = proportional to income, -1 = inversely proportional to income).
damage_elasticity = 1.0

# Share of recycled carbon tax revenue that each region-quintile pair receives (row = region, column = quintile)
recycle_share = ones(12,5) .* 0.2

# Should the time-varying elasticity values only change across the range of GDP values from the studies?
# true = limit calculations to study gdp range, false = allow calculations for 0 to +Inf GDP.
bound_gdp_elasticity = false

# Quintile income distribution scenario (options = "constant", "lessInequality", "moreInequality", "SSP1", "SSP2", "SSP3", "SSP4", or "SSP5")
quintile_income_scenario = "constant"

# Do you also want to perform a reference case optimization run (with no revenue recycling)?
run_reference_case = true

# Name of folder to store your results in (a folder will be created with this name).
results_folder = "base_case"

# run and SAVE, or just run
save_results = true

# ------------------------------------------------------------------------------------------------
# CHOICES ABOUT YOUR ANALYSIS & OPTIMZATION
# ------------------------------------------------------------------------------------------------

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
# RUN EVERYTHING & SAVE KEY RESULTS
# ------------------------------------------------------------------------------------------------

# ------------------------------
# Baseline Run
# ------------------------------
# Note: model parameters not controlled above are set in the "instantiate_model_in_interface.jl file"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------
# eta=2 and rho=0.1%
# ------------------------------
ρ =  0.001; η =  2
results_folder = "eta_2_rho_01"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------
# Damage Elasticity = 0.5
# ------------------------------
ρ =  0.015; η =  1.5
damage_elasticity = 0.5
results_folder = "damage_elasticity_05"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------
# Damage Elasticity = -1
# ------------------------------
damage_elasticity = -1.0
results_folder = "damage_elasticity_neg1"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# -------------------------------------
# Bounded GDP in Elasticity Regression
# -------------------------------------
damage_elasticity = 1.0
bound_gdp_elasticity = true
results_folder = "bounded_GDP"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# -----------------------------------------------------------
# Bottom Quintile Gets 40% of Revenue and Rest Get 15% Each
# -----------------------------------------------------------
bound_gdp_elasticity = false
recycle_share = [ones(12,1) .* 0.4 ones(12,4) * 0.15]
results_folder = "bottom_recycle_40"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))

# -------------------------------------------------
# No Tax Revenue Recycled to the Bottom Quintile
# -------------------------------------------------
recycle_share = [ones(12,1) .* 0.0 ones(12,4) * 0.25]
results_folder = "no_recycling_for_bottom_quintile"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))

# ------------------------------------------------
# All Revenue to the Bottom Two Quintiles
# ------------------------------------------------
recycle_share = [ones(12,2) .* 0.5 ones(12,3) * 0.0]
results_folder = "all_revenue_to_bottom_2_quintile"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))

# ------------------------------------------------
# 30% of Tax Revenue is the Administrative Cost
# ------------------------------------------------
recycle_share = ones(12,5) .* 0.2 .* 0.7
results_folder = "admin_burden_70_percent_tax"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))

# ---------------------
# No Climate Damages
# ---------------------
recycle_share = ones(12,5) .* 0.2
results_folder = "no_climate_damages"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
update_param!(nice, :slrmultiplier, bau_model[:sealeveldamages,:slrmultiplier] .* 0)
update_param!(nice, :a1, bau_model[:damages,:a1] .* 0)
update_param!(nice, :a2, bau_model[:damages,:a2] .* 0)
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
update_param!(nice, :slrmultiplier, bau_model[:sealeveldamages,:slrmultiplier] .* 0)
update_param!(nice, :a1, bau_model[:damages,:a1] .* 0)
update_param!(nice, :a2, bau_model[:damages,:a2] .* 0)
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Double Conventional Damages (Normal SLR Damages)
# ------------------------------------------------
results_folder = "double_climate_damages"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
update_param!(nice, :a1, bau_model[:damages,:a1] .* 2)
update_param!(nice, :a2, bau_model[:damages,:a2] .* 2)
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
update_param!(nice, :a1, bau_model[:damages,:a1] .* 2)
update_param!(nice, :a2, bau_model[:damages,:a2] .* 2)
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ---------------------------------------------------------------------------
# Double Conventional Damages (Normal SLR Damages) + Damage Elasticity = -1
# ---------------------------------------------------------------------------
damage_elasticity = -1.0
results_folder = "double_climate_damages_elasticity_neg_1"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
update_param!(nice, :a1, bau_model[:damages,:a1] .* 2)
update_param!(nice, :a2, bau_model[:damages,:a2] .* 2)
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
update_param!(nice, :a1, bau_model[:damages,:a1] .* 2)
update_param!(nice, :a2, bau_model[:damages,:a2] .* 2)
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Inequality Projections (less inequality)
# ------------------------------------------------
damage_elasticity = 1
quintile_income_scenario = "lessInequality"
results_folder = "time_varying_consumption_shares_"*quintile_income_scenario
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Inequality Projections (more inequality)
# ------------------------------------------------
quintile_income_scenario = "moreInequality"
results_folder = "time_varying_consumption_shares_"*quintile_income_scenario
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Inequality Projections (SSP1)
# ------------------------------------------------
quintile_income_scenario = "SSP1"
results_folder = "time_varying_consumption_shares_"*quintile_income_scenario
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Inequality Projections (SSP2)
# ------------------------------------------------
quintile_income_scenario = "SSP2"
results_folder = "time_varying_consumption_shares_"*quintile_income_scenario
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Inequality Projections (SSP3)
# ------------------------------------------------
quintile_income_scenario = "SSP3"
results_folder = "time_varying_consumption_shares_"*quintile_income_scenario
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Inequality Projections (SSP4)
# ------------------------------------------------
quintile_income_scenario = "SSP4"
results_folder = "time_varying_consumption_shares_"*quintile_income_scenario
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Inequality Projections (SSP5)
# ------------------------------------------------
quintile_income_scenario = "SSP5"
results_folder = "time_varying_consumption_shares_"*quintile_income_scenario
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------------------
# Use More Basic Consumption Share Studies in Meta-Regression
# ------------------------------------------------------------
quintile_income_scenario = "constant"
results_folder = "expenditure_share_based_studies"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
meta_intercept, meta_slope = meta_regression(elasticity_studies[(elasticity_studies.Type .== "Direct") .| (elasticity_studies.Type .== "EMI"),:])
update_param!(nice, :elasticity_intercept, meta_intercept)
update_param!(nice, :elasticity_slope, meta_slope)
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
meta_intercept, meta_slope = meta_regression(elasticity_studies[(elasticity_studies.Type .== "Direct") .| (elasticity_studies.Type .== "EMI"),:])
update_param!(nice, :elasticity_intercept, meta_intercept)
update_param!(nice, :elasticity_slope, meta_slope)
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ----------------------------------------------------
# Use CGE and Input-Output Studies in Meta-Regression
# ----------------------------------------------------
results_folder = "CGE_IO_studies"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
meta_intercept, meta_slope = meta_regression(elasticity_studies[(elasticity_studies.Type .== "CGE") .| (elasticity_studies.Type .== "IO"),:])
update_param!(nice, :elasticity_intercept, meta_intercept)
update_param!(nice, :elasticity_slope, meta_slope)
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
meta_intercept, meta_slope = meta_regression(elasticity_studies[(elasticity_studies.Type .== "CGE") .| (elasticity_studies.Type .== "IO"),:])
update_param!(nice, :elasticity_intercept, meta_intercept)
update_param!(nice, :elasticity_slope, meta_slope)
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------
# Optimal Global Revenue Recycling
# ------------------------------------------------
results_folder = "optimal_equal_global_recycling"
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
update_param!(nice, :global_recycle_share, ones(12))
include(joinpath("replication_helpers", "optimize_recycle.jl"))
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
update_param!(nice, :global_recycle_share, ones(12))
include(joinpath("replication_helpers", "optimize_reference.jl"))

# ------------------------------------------------------
# ------------------------------------------------------
# Add iterative optimzation with FAIRv1.3 Climate Model
# ------------------------------------------------------
# ------------------------------------------------------

# Reset model parameters to default values.
ρ = 0.015
η = 1.5
damage_elasticity = 1.0
recycle_share = ones(12,5) .* 0.2
bound_gdp_elasticity = false
quintile_income_scenario = "constant"

# Reset optimization parameters to default values
# Note: FAIR optimziation just uses a local optimization due to computational burden of iterative approach
run_reference_case = true
n_objectives = 45
n_fair_loops = 4
local_opt_algorithm = :LN_SBPLX
local_stop_time = 800
local_tolerance = 1e-11
results_folder = "FAIR_climate_model"

# Initialize BAU and base version of NICE with revenue recucling.
include(joinpath("replication_helpers", "instantiate_model_in_interface.jl"))
include(joinpath("replication_helpers", "optimize_with_FAIR.jl"))

# Reload NICE+recycling model.
include("MimiNICE_recycle_time_varying.jl")

# ------------------------------------------------------
# ------------------------------------------------------
# 2°C Recycling Runs (non-optimized)
# ------------------------------------------------------
# ------------------------------------------------------

# Get baseline instance of the model.
nice = create_nice_recycle()
run(nice)
backstopprices = nice[:emissions,:pbacktime] * 1000;

# Load 2°C carbon tax values and resulting mitiation rates.
twodegCP = DataFrame(load(joinpath(@__DIR__, "..", "data", "two_degree_carbon_tax.csv")))[!,1]
twodegMIU = Matrix(DataFrame(load(joinpath(@__DIR__, "..", "data", "two_degree_carbon_mitigation.csv"))))

# Set 2°C rates and prices and run model.
prices, rates = mu_from_tax(twodegCP[2:end,1], backstopprices);
update_param!(nice, :MIU, rates)
update_param!(nice, :global_carbon_tax, prices)
update_param!(nice, :min_study_gdp, 1e-10)
update_param!(nice, :max_study_gdp, +Inf)
run(nice)

# Get an instance of the BAU no-policy model.
bau_model = create_nice_recycle()
update_param!(bau_model, :MIU, rates * 0)
update_param!(bau_model, :min_study_gdp, 1e-10)
update_param!(bau_model, :max_study_gdp, +Inf)
run(bau_model)

# Save the results
output_directory = joinpath(@__DIR__, "..", "results", "two_degree_constraint")
mkpath(output_directory)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=true)

# Save the no reycle runs.
update_param!(nice, :global_carbon_tax, prices .* 0)
run(nice)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=false)

# Save the BAU runs.
output_directory = joinpath(@__DIR__, "..", "results", "bau_no_policy_at_all")
mkpath(output_directory)
save_nice_recycle_results(bau_model, bau_model, prices * 0, prices * 0, output_directory, revenue_recycling=false)

# ---------------------------------------------------------------------------------------
# Carry out runs with only CGE and IO studies, and L-shaped elasticity-income dependence
# ---------------------------------------------------------------------------------------
nice = create_nice_recycle()
update_param!(nice, :MIU, rates)
update_param!(nice, :global_carbon_tax, prices)

# First, perform a meta-regression based on study results to calculate elasticity vs. ln gdp per capita relationship.
meta_intercept, meta_slope = meta_regression(elasticity_studies[elasticity_studies.Type .!= "Direct",:])
update_param!(nice, :elasticity_intercept, 0.7)
update_param!(nice, :elasticity_slope, 0.0)
run(nice)

# Save recyling results.
output_directory = joinpath(@__DIR__, "..", "results", "two_degree_constraint_CGE_IO_studies_only")
mkpath(output_directory)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=true)

# Save no-recycling results.
update_param!(nice, :global_carbon_tax, prices .* 0)
run(nice)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=false)

println("All done!")
