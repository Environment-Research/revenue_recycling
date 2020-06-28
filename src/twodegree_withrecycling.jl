# create a nice model
cd("C:\\Users\\fdennig\\Documents\\code\\revenue_recyling\\src" )
include("MimiNICE_recycle_time_varying.jl")
nice = create_nice_recycle()
# Set user-specified parameter settings.
set_param!(nice, :nice_recycle, :damage_elasticity, 1)
set_param!(nice, :nice_recycle, :recycle_share, ones(12,5) .* 0.2)
# Load quintile income distribution scenario data.
quintile_income_scenario = "constant"
income_distribution_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", quintile_income_scenario*"_quintile_distributions_consumption.csv")))
income_distributions = get_quintile_income_shares(income_distribution_raw)
set_param!(nice, :nice_recycle, :quintile_income_shares, income_distributions)
run(nice)
backstopprices = nice[:emissions,:pbacktime] * 1000

# load two degree run from previous optimisation 
twodegrees = "C:\\Users\\fdennig\\Dropbox\\ARBEIT\\aRESEARCH\\revenue usa\\model_results_20200625\\2_degree_constraint\\"
twodegCP = DataFrame(load(joinpath(twodegrees*"no_revenue_recycling (optimized)\\global_output\\carbon_tax.csv")))
twodegMIU= DataFrame(load(joinpath(twodegrees*"no_revenue_recycling (optimized)\\regional_output\\regional_co2_mitigation.csv")))
# set two degree rates and prices and run model
prices, rates = mu_from_tax(twodegCP[2:end,1], backstopprices)
set_param!(nice, :emissions, :MIU, rates)
set_param!(nice, :nice_recycle, :global_carbon_tax, prices)
set_param!(nice, :nice_recycle, :min_study_gdp, 1e-10)
set_param!(nice, :nice_recycle, :max_study_gdp, +Inf)
run(nice)

#save the results 
output_directory = "C:\\Users\\fdennig\\Dropbox\\ARBEIT\\aRESEARCH\\revenue usa\\model_results_20200625\\2_degree_constraint\\"
bau_model = create_nice_recycle()
set_param!(bau_model, :emissions,    :MIU, rates * 0)
set_param!(nice, :nice_recycle, :min_study_gdp, 1e-10)
set_param!(nice, :nice_recycle, :max_study_gdp, +Inf)
run(bau_model)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=true)

# save bau
set_param!(nice, :emissions,    :MIU, rates * 0 )
set_param!(nice, :nice_recycle, :global_carbon_tax, prices * 0)
run(nice)
output_directory = "C:\\Users\\fdennig\\Dropbox\\ARBEIT\\aRESEARCH\\revenue usa\\model_results_20200625\\bau_no_policy_at_all\\"
save_nice_recycle_results(nice, nice, prices * 0, prices * 0, output_directory)