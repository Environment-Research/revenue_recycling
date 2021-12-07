# #--------------------------------------------------------------------------------------------------------------------------
# #--------------------------------------------------------------------------------------------------------------------------
# # This file contains functions used for adding revenue recycling to the Nested Inequalities Climate Economy (NICE) model.
# #--------------------------------------------------------------------------------------------------------------------------
# #--------------------------------------------------------------------------------------------------------------------------


#####################################################################################################################
# CLEAN UP AND ORGANIZE INCOME SHARE DATA INTO NICE FORMAT.
#####################################################################################################################
# Description: This function takes a DataFrame of time varying quintile income distributions (for each NICE region
#              from 2015-2105) and converts it into an array with the dimension format required by NICE (time, region,
#              quintile). It sets the 2005 quintile share to the calculated 2015 value and assumes quintile shares after
#              2105 remain fixed.
#
# Function Arguments:
#
#       income_data   = A DataFrame read in from one of the income share files in the "data" folder of NICE.
#--------------------------------------------------------------------------------------------------------------------

function get_quintile_income_shares(income_data::DataFrame)

    # Allocate an array for cleaned up quintile share data in format required for NICE (time × regions × quintile).
    nice_years = collect(2005:10:2595)
    quintile_shares = zeros(length(nice_years), 12, 5);

    # Create vector of years for SSP quintile data (spans 2015-2105 in ten-year increments).
    ssp_years = unique(income_data[:,:decade])    

    # Create vector of column names for each quintile.
    quintile_names = [:s1, :s2, :s3, :s4, :s5]

    # Sort the data by decade periods.
    sort!(income_data, :decade)

    # Loop through each quintile and decade and assign appropriate income shares to each region.
    for (d, decade) in enumerate(ssp_years)
        for (q, quintile) in enumerate(quintile_names)

            # Find indices with data for a particular region in a given decade.
            decade_index = findall(income_data.decade .== decade)

            # Allocate those particular income share values and convert from percentage to income shares.
            quintile_shares[d+1, :, q] = income_data[decade_index, quintile] ./ 100
        end
    end

    # Set shares in 2005 to calculated 2015 value (estiamted values start in 2015, but NICE initializes in 2005).
    quintile_shares[1,:,:] .= quintile_shares[2,:,:]

    # Get index for first decade without calculated shares (i.e. after 2105).
    missing_index = findfirst(x->x==0, quintile_shares[:,1,1])

    # Set income shares constant (to 2105 value) for all periods after 2105
    for q = 1:5
        quintile_shares[missing_index:end,:,q] = repeat(quintile_shares[missing_index-1,:,q]', length(nice_years) - missing_index +1)
    end

    # Return cleaned up income shares.
    return quintile_shares
end


#####################################################################################################################
# CALCULATE DAMAGE, CO₂ MITIGATION COST, OR CO₂ TAX BURDEN DISTRIBUTIONS ACROSS AN INDIVIDUAL REGION'S QUINTILES.
#####################################################################################################################
# Description: This function will calculate quintile distribution shares for an individual RICE region based
#              on a provided income elasticity.
#
# Function Arguments:
#
#       elasticity    = Income elasticity of climate damages, CO₂ mitigation costs, CO₂ tax burdens, etc.
#       income_shares = A vector of quintile income shares for a given RICE region.
#--------------------------------------------------------------------------------------------------------------------

function regional_quintile_distribution(elasticity, income_shares)

    # Apply elasticity to quintile income shares.
    scaled_shares = income_shares .^ elasticity

    # Allocate empty array for distribution across quintiles resulting from the elasticity.
    updated_quintile_distribution = zeros(5)

    # Loop through each quintile to calculate updated distribution.
    for q in 1:5
        updated_quintile_distribution[q] = scaled_shares[q] ./ sum(scaled_shares[:])
    end

    return updated_quintile_distribution
end


#####################################################################################################################
# PERFORM A REGRESSION.
#####################################################################################################################
# Description: This function will perform a regression of 'Y' against a user-supplied 'Data' term.
#
# Function Arguments:
#
#       Y    = Dependent variable in regression.
#       Data = Independent data (x-terms) in regression.
#--------------------------------------------------------------------------------------------------------------------

function regress(Y, Data)

    # Get size of data.
    n = length(Data[:,1])
    k = length(Data[1,:])

    # Add column of ones.
    X = [ones(n) Data]

    # Calculate beta coefficient.
    x_prime_inv = inv(X' * X)
    betahat = x_prime_inv * X' * Y

    # Calculate residuals.
    u = Y - X * betahat

    # Calculate hetero var-cov matrix.
    X_eps = zeros(n,k+1)
    for i = 1:n
        X_eps[i,:] = u[i] * X[i,:]
    end

    varhetero = x_prime_inv * X_eps' * X_eps * x_prime_inv

    return betahat, varhetero
end


#####################################################################################################################
# PERFORM A META-REGRESSION OVER MULTIPLE ELASTICITY STUDIES.
#####################################################################################################################
# Description: This function will perform a regression (elasticity vs. ln gdp per capita) across the results of
#              multiple studies to provide the coefficients for a time-varying income elasticity of CO₂ price
#              exposure relationship.
#
# Function Arguments:
#
#       Data       = Independent data (x-terms) in regression coming from multiple studies.
#       slope_type = The type of analysis to carry out. Options are :central, :steeper, :flatter, :percentile
#                       :central    = returns central regression result
#                       :steeper    = returns results for the steepest slope in the 95% confidence interval range
#                       :flatter    = returns results for the shallowest slope in the 95% confidence interval range
#                       :percentile = returns results for a horizontal line corresponding to a given elasticity percentile
#       percentile = Calculate this percentile for the elasticity values (defaults to 0.9 for the 90th percentile).
#       z_val      = z-value for a given confidence interval (defaults to 1.96 for 95% interval)
#--------------------------------------------------------------------------------------------------------------------

function meta_regression(Data; slope_type::Symbol=:central, percentile::Float64=0.90, z_val::Float64=1.96)

    # Sort study data by region and country names.
    sort!(Data, [:Region, :CountryA3])
    n = size(Data, 1)

    # Allocate array to store elastcities.
    elasticities = zeros(n)

    # Isolate variables needed for each run (quintile fuel spending and income shares).
    cost_shares   = [:costS1, :costS2, :costS3, :costS4, :costS5]
    exp_shares = [:e1, :e2, :e3, :e4, :e5]

    # Loop over each study and calculate elasticity with a log-log regression.
    for i = 1:n

        # Take logs of variables
        log_exp_share = log.(Vector(Data[i, exp_shares]))
        log_cost_share= log.(Vector(Data[i, cost_shares]))

        # Run log-log regression and store elasticity.
        B, V   = regress(log_cost_share, log_exp_share)
        elasticities[i] = B[2]
    end

    # Calculate log of per capita GDP.
    Data.log_pcGDP = log.(Data[!,:pcGDP])

    # Using these elasticities, run meta-regression (elasticity vs. ln gdp per capita).
    meta_B, meta_V  = regress(elasticities, Data[!,:log_pcGDP])
    meta_intercept  = meta_B[1]
    meta_slope      = meta_B[2]

    # ------------------------------------------------------------------------------------------------------
    # Calculate flattest and steepest slopes based on 95% confidence interval range.
    # ------------------------------------------------------------------------------------------------------

    # First, need 95% CI around regression lines (requires squared errors & "x_spread" -> (x - mean(x))²)

    # Predicted elasticity values.
    elasticity_hat = meta_intercept .+ meta_slope .* Data[!,:log_pcGDP]

    # Calculate squared errors of elasticities and their mean.
    sq_err = (elasticity_hat .- elasticities) .^2
    mse = mean(sq_err)

    # Calculate squared differences between per capita GDP values and their mean.
    x_spread = (Data[!,:log_pcGDP] .- mean(Data[!,:log_pcGDP])).^2

    # Calculate lower and upper confidence intervals.
    lower_ci = elasticity_hat .- z_val * sqrt.(mse*(1/n .+ x_spread ./ sum(x_spread)))
    upper_ci = elasticity_hat .+ z_val * sqrt.(mse*(1/n .+ x_spread ./ sum(x_spread)))

    # Second, For the steepest line connect maximum point of upper CI and minimum of lower CI.

    # Calculate steeper slope and intercept.
    steeper_slope     = (maximum(upper_ci) - minimum(lower_ci)) / (minimum(Data[!,:log_pcGDP]) - maximum(Data[!,:log_pcGDP]))
    steeper_intercept = maximum(upper_ci) - steeper_slope * minimum(Data[!,:log_pcGDP])

    # Also calculate a shallower slope and intercept.
    flatter_slope     = (minimum(upper_ci) - maximum(lower_ci)) / (maximum(Data[!,:log_pcGDP]) - minimum(Data[!,:log_pcGDP]))
    flatter_intercept = minimum(upper_ci) - flatter_slope * maximum(Data[!,:log_pcGDP])

    # ------------------------------------------------------------------------------------------------------
    # Calculate horizontal lines for a given elasticity percentile.
    # ------------------------------------------------------------------------------------------------------

    # Horizontal line, so slope = 0.
    percentile_slope = 0.0

    # Calculate intercept for a given percentile.
    percentile_intercept = quantile(elasticities, percentile)

    # ------------------------------------------------------------------------------------------------------
    # Return slope and intercept depending on user-specification for analysis type.
    # ------------------------------------------------------------------------------------------------------

    if slope_type == :central
        return meta_intercept, meta_slope
    elseif slope_type == :steeper
        return steeper_intercept, steeper_slope
    elseif slope_type == :flatter
        return flatter_intercept, flatter_slope
    elseif slope_type == :percentile
        return percentile_intercept, percentile_slope
    else
        println("Incorrect regression analysis type selected. Options are :central, :steeper, :flatter, :percentile")
    end
end


#######################################################################################################################
# CALCULATE REGIONAL CO₂ MITIGATION
########################################################################################################################
# Description: This function calculates regional CO₂ mitigation levels as a function of a global carbon tax.  It
#              uses the RICE2010 backstop price values and assumes a carbon tax of $0 in period 1.  If the number of
#              tax values is less than the total number of model time periods, the function assumes full decarbonization
#              (e.g. the tax = the backstop price) for all future periods without a specified tax.
#
# Function Arguments:
#
#       opt_tax_vals:    A vector of global carbon tax values being optimized.
#       backstop_prices: The regional backstop prices from RICE2010 (units = dollars).
#       theta2:          The exponent on the abatement cost function (defaults to RICE2010 value).
#----------------------------------------------------------------------------------------------------------------------

function mu_from_tax(opt_tax_vals::Array{Float64,1}, backstop_prices::Array{Float64,2}, ; theta2::Float64=2.8)

    # Set first period to $0, then optimized tax values, and then remaining values to maximum backstop price (i.e. full decarbonization).
    full_CO₂_tax = [0.0; opt_tax_vals; maximum(backstop_prices, dims=2)[(length(opt_tax_vals)+2):end]]

    # Calculate regional mitigation rates from the full tax vector.
    regional_CO₂_abatement = min.((max.(((full_CO₂_tax ./ backstop_prices) .^ (1 / (theta2 - 1.0))), 0.0)), 1.0)

    return full_CO₂_tax, regional_CO₂_abatement
end



#######################################################################################################################
# CREATE OBJEVTIVE FUNCTION FOR NICE WITH REVENUE RECYCLING
########################################################################################################################
# Description: This function creates an objective function that uses an instance of NICE with revenue recycling and
#              user-specified settings. The objective function will take in a vector of global carbon tax values and
#              returns the total economic welfare generated by that specifc climate policy.
#
# Function Arguments:
#
#       m:                 An instance of NICE with revenue recycling (type = a Mimi model).
#       revenue_recycling: A check for whether or not to recycle CO₂ tax revenue (true = recycle, false = no recycling).
#----------------------------------------------------------------------------------------------------------------------

function construct_nice_recycle_objective(m::Model; revenue_recycling::Bool=true)

    # Extract backstop prices from model and convert to dollars.
    run(m)
    rice_backstop_prices = m[:emissions, :pbacktime] .* 1000

    # Find number of timesteps across model time horizon.
    n_steps = length(dim_keys(m, :time))

    # Pre-allocate matrix to store optimal tax and mitigation rates.
    optimal_CO₂_tax        = zeros(n_steps)
    optimal_CO₂_mitigation = zeros(n_steps, 12)

    # Create a function to optimize user-specified model for (i) revenue recycling and (ii) reference case.
    nice_recycle_objective =

        if revenue_recycling == true

            # Create revenue recycling version of function.
            function(opt_tax::Array{Float64,1})
                # Calculate mitigation rate resulting from global carbon tax.
                optimal_CO₂_tax[:], optimal_CO₂_mitigation[:,:] = mu_from_tax(opt_tax, rice_backstop_prices)
                # Update CO₂ mitigation rate and global CO₂ tax in NICE.
                update_param!(m, :MIU, optimal_CO₂_mitigation)
                update_param!(m, :global_carbon_tax, optimal_CO₂_tax)
                run(m)
                # Return total welfare.
                return m[:nice_welfare, :welfare]
            end

        else

            # Create reference version of function without revenue recycling.
            function(opt_tax::Array{Float64,1})
                # Calculate mitigation rate resulting from global carbon tax.
                optimal_CO₂_tax[:], optimal_CO₂_mitigation[:,:] = mu_from_tax(opt_tax, rice_backstop_prices)
                # Update CO₂ mitigation rate in NICE.
                # Set tax in recycling component to $0 yielding no tax revenue (i.e. switch off revenue recycling).
                update_param!(m, :MIU, optimal_CO₂_mitigation)
                update_param!(m, :global_carbon_tax, zeros(n_steps))
                run(m)
                # Return total welfare.
                return m[:nice_welfare, :welfare]
            end
        end

    # Return the objective function.
    return nice_recycle_objective
end



#######################################################################################################################
# CALCULATE GLOBAL CO₂ MITIGATION LEVELS
#######################################################################################################################
# Description: This function takes an instance of NICE with regional CO₂ mitigation and returns the resulting global
#              reduction in CO₂ emissions relative to a baseline case without a climate policy.
#
# Function Arguments:
#
#       m_policy: An instance of NICE with CO₂ policy (type = Mimi model).
#       m_bau:    An instance of NICE with 0% mitigation (no CO₂ policy) for all regions and years (type = Mimi model).
#----------------------------------------------------------------------------------------------------------------------

function get_global_mitigation(m_policy::Model, m_bau::Model)

    # Calculate global industrial CO₂ emissions for NICE with a particular mitigation policy.
    global_emissions_policy = sum(m_policy[:emissions, :EIND], dims=2)

    # Calculate global industrial CO₂ emissions in the bau case with no CO₂ policy.
    global_emissions_base = sum(m_bau[:emissions, :EIND], dims=2)

    # Calculate change in global CO₂ mitigation levels.
    global_mitigation_rates = ((global_emissions_base .- global_emissions_policy) ./ global_emissions_base)

    return global_mitigation_rates
end



#######################################################################################################################
# CREATE AN INSTANCE OF FAIR TO USE IN OPTIMIZATION WITH NICE
#######################################################################################################################
# Description: This function creates a version of FAIR v1.3 that can be used in an optimization with NICE. It removes
#              all non-CO₂ radiative forcing components (instead it will take in the endgoenous CO₂ emissions and
#              exogenous non-CO₂ forcings from NICE). It also modifies FAIR to run from 2005-2595 (rather than its
#              default of 1765-2500).
#----------------------------------------------------------------------------------------------------------------------

function create_fair_for_opt()

    #----- Initial Data and Index Prep -----#

    # Set the default start and end years for FAIR.
    fair_start_year = 1765
    fair_end_year   = 2500

    # Set the default start and end years for NICE.
    nice_start_year = 2005
    nice_end_year   = 2595

    # Set annual time steps for NICE and calculate number of annual years.
    nice_annual_years  = collect(nice_start_year:nice_end_year)
    n_annual_nice      = length(nice_annual_years)

    # Calculate FAIR index for year 2005.
    fair_2005_index = findfirst(x -> x == 2005, fair_start_year:fair_end_year)

    # Load RCP4.5 concentrations and get N₂O concentrations (needed for CO₂ forcing calculations in FAIR).
    rcp45_concentrations = DataFrame(load(joinpath(@__DIR__, "..", "data", "RCP45_concentrations.csv"), skiplines_begin=37))
    rcp45_n2o_raw = rcp45_concentrations.N2O

    # Isolate RCP4.5 N₂O data for 2005-2500 (start of NICE to end of RCP data). Then extend to 2595 assuming concentrations remain fixed at 2500 value.
    rcp45_n2o = zeros(n_annual_nice)
    rcp45_n2o[1:length(2005:2500)] = rcp45_n2o_raw[fair_2005_index:end]
    rcp45_n2o[length(2005:2500):end] .= rcp45_n2o_raw[end]

    #----- Create Modified Version of FAIR for NICE -----#

    # Load and run a default version of FAIR.
    fair = MimiFAIR13.create_fair(rcp_scenario="RCP45")
    run(fair)

    # Get initial CO₂ concentration and CO₂ stock perturbation (for NICE years).
    Cacc_2005 = fair[:co2_cycle, :Cacc][fair_2005_index]
    CO2_2004  = fair[:co2_cycle, :C][(fair_2005_index-1)]

    # Remove all FAIR components other than CO₂, radiative forcing, and temperature (use NICE exogenous forcing values).
    delete!(fair, :ch4_cycle)
    delete!(fair, :n2o_cycle)
    delete!(fair, :other_ghg_cycles)
    delete!(fair, :ch4_rf)
    delete!(fair, :n2o_rf)
    delete!(fair, :other_ghg_rf)
    delete!(fair, :trop_o3_rf)
    delete!(fair, :strat_o3_rf)
    delete!(fair, :aerosol_direct_rf)
    delete!(fair, :aerosol_indirect_rf)
    delete!(fair, :bc_snow_rf)
    delete!(fair, :landuse_rf)
    delete!(fair, :contrails_rf)

    # Reset FAIR time dimension to NICE years (annual version).
    set_dimension!(fair, :time, nice_annual_years)

    # Set placeholder value for CO₂ emissions.
    update_param!(fair, :E, ones(n_annual_nice), update_timesteps=true)

    #Set initial Cacc perturbation.
    update_param!(fair, :Cacc_0, Cacc_2005)

    # Set initial (2004) atmospheric CO₂ concentration.
  #  set_param!(fair, :CO₂_0, CO2_2004)
    set_param!(fair, :co2_cycle, :CO₂_0, :co2cycle_2004val, CO2_2004)
#set_param!(m::Model, comp_name::Symbol, param_name::Symbol, ext_param_name::Symbol, val::Any)

    # Set RCP4.5 N₂O concentrations for CO₂ radiative forcing calculations.
    set_param!(fair, :N₂O, rcp45_n2o)

    # Set placeholder for NICE exogenous forcings.
    update_param!(fair, :F_exogenous, ones(n_annual_nice), update_timesteps=true)

    # Set all other radiative forcings to zero.
    set_param!(fair, :F_volcanic, zeros(n_annual_nice))
    set_param!(fair, :F_solar, zeros(n_annual_nice))
    set_param!(fair, :F_CH₄, zeros(n_annual_nice))
    set_param!(fair, :F_CH₄_H₂O, zeros(n_annual_nice))
    set_param!(fair, :F_N₂O, zeros(n_annual_nice))
    set_param!(fair, :F_other_ghg, zeros(n_annual_nice, 28))
    set_param!(fair, :F_trop_O₃, zeros(n_annual_nice))
    set_param!(fair, :F_strat_O₃, zeros(n_annual_nice))
    set_param!(fair, :F_aerosol_direct, zeros(n_annual_nice))
    set_param!(fair, :F_aerosol_indirect, zeros(n_annual_nice))
    set_param!(fair, :F_bcsnow, zeros(n_annual_nice))
    set_param!(fair, :F_landuse, zeros(n_annual_nice))
    set_param!(fair, :F_contrails, zeros(n_annual_nice))

    # Return modified version of FAIR.
    return fair
end



#######################################################################################################################
# LINEARLY INTERPOLATE MODEL RESULTS TO ANNUAL VALUES
#######################################################################################################################
# Description: This function uses linear interpolation to create an annual time series from the model results.
#
# Function Arguments:
#
#       data    = The non-annual model results to be interpolated
#       spacing = Length of time steps between model output.
#----------------------------------------------------------------------------------------------------------------------

function interpolate_to_annual(data, spacing)

    # Create an interpolation object for the data (assume first and last points are end points, e.g. no interpolation beyond support).
    interp_linear = interpolate(data, BSpline(Linear()))

    # Create points to interpolate for (based on spacing term).
    interp_points = collect(1:(1/spacing):length(data))

    # Carry out interpolation.
    return interp_linear[interp_points]
end



function construct_fair_recycle_objective(nice::Model, n_loops::Int; revenue_recycling::Bool=true)

    # Run user-specified instance of NICE to obtain some values.
    run(nice)

    # Extract backstop prices from NICE and convert to dollars.
    nice_backstop_prices = nice[:emissions, :pbacktime] .* 1000

    # Extract exogenous radiative forcing scenario from NICE and convert to annual values.
    annual_nice_exog_rf = interpolate_to_annual(nice[:radiativeforcing, :forcoth], 10)

    # Get an instance of FAIR to use in iterative optimization.
    fair = create_fair_for_opt()

    # Set FAIR exogenous forcing to annualized NICE exogenous forcing scenario.
    update_param!(fair, :F_exogenous, annual_nice_exog_rf)

    # Find number of timesteps across NICE model time horizon.
    n_nice_decadal = length(dim_keys(nice, :time))

    # Get decadal NICE indices for annual FAIR years (2005-2595)
    nice_decadal_indices = findall((in)(collect(2005:10:2595)), (collect(2005:2595)))

    # Pre-allocate matrix to store optimal tax and mitigation rates at 10-year time steps.
    optimal_CO₂_tax        = zeros(n_nice_decadal)
    optimal_CO₂_mitigation = zeros(n_nice_decadal, 12)

    # Create a function to optimize user-specified model for (i) revenue recycling and (ii) reference case.
    fair_recycle_objective =

        if revenue_recycling == true

            # Create revenue recycling version of function.
            function(opt_tax::Array{Float64,1})

                # Calculate mitigation rate resulting from global carbon tax.
                optimal_CO₂_tax[:], optimal_CO₂_mitigation[:,:] = mu_from_tax(opt_tax, nice_backstop_prices)

                # Update CO₂ mitigation rate and global CO₂ tax in NICE.
                update_param!(nice, :MIU, optimal_CO₂_mitigation)
                update_param!(nice, :global_carbon_tax, optimal_CO₂_tax)
                run(nice)

                for loops = 1:n_loops

                    # Run FAIR with annualized NICE CO₂ emissions to calculate temperature.
                    update_param!(fair, :E, interpolate_to_annual(nice[:emissions, :E], 10))
                    run(fair)

                    # Set FAIR temperature projections in NICE.
                    set_param!(nice, :TATM, fair[:temperature, :T][nice_decadal_indices])
                    run(nice)
                end

                # Return total welfare.
                return nice[:nice_welfare, :welfare]
            end

        else

            # Create reference version of function without revenue recycling.
            function(opt_tax::Array{Float64,1})

                # Calculate mitigation rate resulting from global carbon tax.
                optimal_CO₂_tax[:], optimal_CO₂_mitigation[:,:] = mu_from_tax(opt_tax, nice_backstop_prices)

                # Update CO₂ mitigation rate in NICE.
                # Set tax in recycling component to $0 yielding no tax revenue (i.e. switch off revenue recycling).
                update_param!(nice, :MIU, optimal_CO₂_mitigation)
                update_param!(nice, :global_carbon_tax, zeros(n_nice_decadal))
                run(nice)

                for loops = 1:n_loops

                    # Run FAIR with annualized NICE CO₂ emissions to calculate temperature.
                    update_param!(fair, :E, interpolate_to_annual(nice[:emissions, :E], 10))
                    run(fair)

                    # Set FAIR temperature projections in NICE.
                    set_param!(nice, :TATM, fair[:temperature, :T][nice_decadal_indices])
                    run(nice)
                end

                # Return total welfare.
                return nice[:nice_welfare, :welfare]
            end
        end

    # Return the objective function.
    return fair_recycle_objective
end



#######################################################################################################################
# CREATE RESULT DIRECTORIES AND SAVE SPECIFIC MODEL OUTPUT
#######################################################################################################################
# Description: This function will create a folder directory to store results (dividing model output by global,
#              regional, and quintile levels) for either the revenue recyling or reference case without recycling.
#
# Function Arguments:
#
#       m_policy:                 An instance of NICE with CO₂ policy (type = Mimi model).
#       m_bau:                    An instance of NICE with 0% mitigation (no CO₂ policy) for all regions and years (type = Mimi model).
#       opt_tax:                  The optimal tax policy vector resulting from the second-stage local optimization.
#       opt_tax_global_algorithm: The optimal tax policy vector resulting from the first-stage global optimization (just saved for reference).
#       output_directory:         The directory path to the results folder where a particular set of model output will be saved.
#       revenue_recycling:        A check for whether or not the results recycle CO₂ tax revenue (true = recycle, false = no recycling).
#----------------------------------------------------------------------------------------------------------------------

function save_nice_recycle_results(m_policy::Model, m_bau::Model, opt_tax::Array{Float64,1}, opt_tax_global_algorithm::Array{Float64,1}, output_directory::String; revenue_recycling::Bool=true)

    # Make subdirectory folders to store results with and without revenue recycling.
    if revenue_recycling == true

        global_path   = joinpath(output_directory, "revenue_recycling", "global_output")
        regional_path = joinpath(output_directory, "revenue_recycling", "regional_output")
        quintile_path = joinpath(output_directory, "revenue_recycling", "quintile_output")

        mkpath(global_path)
        mkpath(regional_path)
        mkpath(quintile_path)

    else

        global_path   = joinpath(output_directory, "no_revenue_recycling", "global_output")
        regional_path = joinpath(output_directory, "no_revenue_recycling", "regional_output")
        quintile_path = joinpath(output_directory, "no_revenue_recycling", "quintile_output")

        mkpath(global_path)
        mkpath(regional_path)
        mkpath(quintile_path)
    end

    # Save Global Output.
    save(joinpath(global_path, "global_co2_mitigation.csv"), DataFrame(get_global_mitigation(m_policy, m_bau), :auto))
    save(joinpath(global_path, "temperature.csv"), DataFrame(temperature=m_policy[:climatedynamics, :TATM]))
    save(joinpath(global_path, "rough_carbon_tax_global_algorithm.csv"), DataFrame(global_algorithm_tax=opt_tax_global_algorithm))
    save(joinpath(global_path, "carbon_tax.csv"), DataFrame(tax=opt_tax))
    save(joinpath(global_path, "global_emissions_with_landuse.csv"), DataFrame(emissions=m_policy[:emissions, :E]))
    save(joinpath(global_path, "co2_concentration_ppm.csv"), DataFrame(CO2=m_policy[:co2cycle, :MAT] ./ 2.13))

    # Save Regional Output.
    save(joinpath(regional_path, "gross_output.csv"), DataFrame(m_policy[:grosseconomy, :YGROSS], :auto))
    save(joinpath(regional_path, "nice_net_output.csv"), DataFrame(m_policy[:nice_neteconomy, :Y], :auto))
    save(joinpath(regional_path, "consumption.csv"), DataFrame(m_policy[:nice_neteconomy, :C], :auto))
    save(joinpath(regional_path, "population.csv"), DataFrame(m_policy[:nice_neteconomy, :l], :auto))
    save(joinpath(regional_path, "industrial_co2_emissions.csv"), DataFrame(m_policy[:emissions, :EIND], :auto))
    save(joinpath(regional_path, "damage_cost_share.csv"), DataFrame(m_policy[:nice_neteconomy, :DAMFRAC], :auto))
    save(joinpath(regional_path, "abatement_cost_share.csv"), DataFrame(m_policy[:nice_neteconomy, :ABATEFRAC], :auto))
    save(joinpath(regional_path, "regional_co2_mitigation.csv"), DataFrame(m_policy[:emissions, :MIU], :auto))
    save(joinpath(regional_path, "co2_tax_revenue.csv"), DataFrame(m_policy[:nice_recycle, :tax_revenue], :auto))
    save(joinpath(regional_path, "co2_income_elasticity.csv"), DataFrame(m_policy[:nice_recycle, :CO₂_income_elasticity], :auto))

    # Save Quintile Output.
    save(joinpath(quintile_path, "co2_tax_distribution_Q1.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,1], :auto))
    save(joinpath(quintile_path, "co2_tax_distribution_Q2.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,2], :auto))
    save(joinpath(quintile_path, "co2_tax_distribution_Q3.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,3], :auto))
    save(joinpath(quintile_path, "co2_tax_distribution_Q4.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,4], :auto))
    save(joinpath(quintile_path, "co2_tax_distribution_Q5.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,5], :auto))
    save(joinpath(quintile_path, "base_pc_consumption_Q1.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,1], :auto))
    save(joinpath(quintile_path, "base_pc_consumption_Q2.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,2], :auto))
    save(joinpath(quintile_path, "base_pc_consumption_Q3.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,3], :auto))
    save(joinpath(quintile_path, "base_pc_consumption_Q4.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,4], :auto))
    save(joinpath(quintile_path, "base_pc_consumption_Q5.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,5], :auto))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q1.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,1], :auto))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q2.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,2], :auto))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q3.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,3], :auto))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q4.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,4], :auto))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q5.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,5], :auto))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q1.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,1], :auto))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q2.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,2], :auto))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q3.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,3], :auto))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q4.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,4], :auto))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q5.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,5], :auto))

end
