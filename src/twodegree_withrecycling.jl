# create a nice model
cd("C:\\Users\\fdennig\\Documents\\code\\revenue_recyling\\src" )
include("MimiNICE_recycle_time_varying.jl")

####################################################################################################
## the baseline run: consumption distribution, kevin elasticity data from 20200625 and 2005PPP pcGDP
####################################################################################################
nice = create_nice_recycle()
run(nice)
backstopprices = nice[:emissions,:pbacktime] * 1000;

# load two degree run from previous optimisation 
twodegrees = "C:\\Users\\fdennig\\Dropbox\\ARBEIT\\aRESEARCH\\revenue usa\\model_results\\2_degree_constraint\\"
twodegCP = DataFrame(load(joinpath(twodegrees*"no_revenue_recycling (optimized)\\global_output\\carbon_tax.csv")));
twodegMIU= DataFrame(load(joinpath(twodegrees*"no_revenue_recycling (optimized)\\regional_output\\regional_co2_mitigation.csv")));
# set two degree rates and prices and run model
prices, rates = mu_from_tax(twodegCP[2:end,1], backstopprices);
set_param!(nice, :emissions, :MIU, rates)
set_param!(nice, :nice_recycle, :global_carbon_tax, prices)
set_param!(nice, :nice_recycle, :min_study_gdp, 1e-10)
set_param!(nice, :nice_recycle, :max_study_gdp, +Inf)
run(nice)

#save the results 
output_directory = "C:\\Users\\fdennig\\Documents\\code\\revenue_recyling\\results\\2_degree_constraint\\"
mkpath(output_directory)
bau_model = create_nice_recycle()
set_param!(bau_model, :emissions,    :MIU, rates * 0)
set_param!(bau_model, :nice_recycle, :min_study_gdp, 1e-10)
set_param!(bau_model, :nice_recycle, :max_study_gdp, +Inf)
run(bau_model)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=true)

#save no reycle runs
set_param!(nice,:nice_recycle,:global_carbon_tax, prices .* 0)
run(nice)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=false)

# save bau
output_directory = "C:\\Users\\fdennig\\Documents\\code\\revenue_recyling\\results\\bau_no_policy_at_all\\"
mkpath(output_directory)
save_nice_recycle_results(bau_model, bau_model, prices * 0, prices * 0, output_directory, revenue_recycling=false)

####################################################################################################
## only CGE and IO studies, and L-shaped elasticity-income dependence
####################################################################################################
nice = create_nice_recycle()
set_param!(nice, :emissions, :MIU, rates)
set_param!(nice, :nice_recycle, :global_carbon_tax, prices)
# First, perform a meta-regression based on study results to calculate elasticity vs. ln gdp per capita relationship.
meta_intercept, meta_slope = meta_regression(elasticity_studies[elasticity_studies.Type .!= "Direct",:])
set_param!(nice, :nice_recycle, :elasticity_intercept, 0.7)
set_param!(nice, :nice_recycle, :elasticity_slope, 0.0)
run(nice)
output_directory = "C:\\Users\\fdennig\\Documents\\code\\revenue_recyling\\results\\2_degree_constraint"
mkpath(output_directory)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=true)
set_param!(nice,:nice_recycle,:global_carbon_tax, prices .* 0)
run(nice)
save_nice_recycle_results(nice, bau_model, prices, prices, output_directory, revenue_recycling=false)