# Load required packages.
using CSVFiles, DataFrames, Mimi, MimiRICE2010, MimiNICE

# Load helper functions and revenue recycling components being added to the NICE model.
include("helper_functions.jl")
include(joinpath("revenue_recycling_components", "nice_revenue_recycle_component_time_varying.jl"))

#-----------------------#
# ----- Load Data ----- #
#-----------------------#

# Load studies used to calculate elasticities.
elasticity_studies = DataFrame(load(joinpath(@__DIR__, "..", "data", "elasticity_study_data.csv")))

# Load updated UN population projections and convert units to millions of people.
un_population_data = convert(Matrix, DataFrame(load(joinpath(@__DIR__, "..", "data", "UN_medium_population_scenario.csv"), skiplines_begin=3))[:, 3:end]) ./ 1000

# Load quintile income distribution data that remains fixed over time (these shares are an update from the original NICE income distributions).
income_distribution_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "constant_quintile_distributions_consumption.csv")))

# Clean up and organize time-varying income distribution data into required NICE format (time × regions × quintiles).
income_distributions = get_quintile_income_shares(income_distribution_raw)

# -------------------------------------------------
# -------------------------------------------------
# Create function to add revenue recyling to NICE.
# -------------------------------------------------
# -------------------------------------------------

function create_nice_recycle()

    # Initialize an instance of NICE to build in revenue recycling.
    nice_rr = MimiNICE.create_nice()

    # Get number of timesteps across model time horizon.
    n_steps = length(dim_keys(nice_rr, :time))

    # Add in NICE revenue recycling component.
    add_comp!(nice_rr, nice_recycle, after = :nice_neteconomy)

    #-------------------------------------#
    # ----- Assign Model Parameters ----- #
    #-------------------------------------#

    # First, perform a meta-regression based on study results to calculate elasticity vs. ln gdp per capita relationship.
    meta_intercept, meta_slope = meta_regression(elasticity_studies)

    # Note, revenue recycling model now uses quintile consumption values calculated in nice_recycle component (not nice_neteconomy).
    set_param!(nice_rr, :grosseconomy, :l, un_population_data)

    set_param!(nice_rr, :nice_neteconomy, :l, un_population_data)

    set_param!(nice_rr, :nice_recycle, :min_study_gdp, minimum(elasticity_studies.GDP))
    set_param!(nice_rr, :nice_recycle, :max_study_gdp, maximum(elasticity_studies.GDP))
    set_param!(nice_rr, :nice_recycle, :elasticity_intercept, meta_intercept)
    set_param!(nice_rr, :nice_recycle, :elasticity_slope, meta_slope)
    set_param!(nice_rr, :nice_recycle, :regional_population, un_population_data)
    set_param!(nice_rr, :nice_recycle, :quintile_income_shares, income_distributions)
    set_param!(nice_rr, :nice_recycle, :damage_elasticity, 1.0)
    set_param!(nice_rr, :nice_recycle, :recycle_share, ones(12,5).*0.2)
    set_param!(nice_rr, :nice_recycle, :global_carbon_tax, zeros(n_steps))

    set_param!(nice_rr, :nice_welfare, :quintile_pop, un_population_data ./ 5)

    # Create parameter connections (:component => :parameter).
    connect_param!(nice_rr, :nice_recycle => :industrial_emissions, :emissions       => :EIND)
    connect_param!(nice_rr, :nice_recycle => :DAMFRAC,              :damages         => :DAMFRAC)
    connect_param!(nice_rr, :nice_recycle => :ABATEFRAC,            :nice_neteconomy => :ABATEFRAC)
    connect_param!(nice_rr, :nice_recycle => :Y,                    :nice_neteconomy => :Y)
    connect_param!(nice_rr, :nice_recycle => :CPC,                  :nice_neteconomy => :CPC)
    connect_param!(nice_rr, :nice_welfare => :quintile_c,           :nice_recycle    => :qc_post_recycle)

    # Return NICE model with revenue recycling.
    return nice_rr
end
