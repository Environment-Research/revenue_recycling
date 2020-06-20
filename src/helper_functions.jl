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
#       Data = Independent data (x-terms) in regression coming from multiple studies.
#--------------------------------------------------------------------------------------------------------------------

function meta_regression(Data)

    # Sort study data by region and country names.
    sort!(Data, [:Region, :Country])
    n = size(Data, 1)

    # Allocate array to store elastcities.
    elasticities = zeros(n)

    # Isolate variables needed for each run (quintile fuel spending and income shares).
    fuel_shares   = [:Q1FS, :Q2FS, :Q3FS, :Q4FS, :Q5FS]
    income_shares = [:Q1IncS, :Q2IncS, :Q3IncS, :Q4IncS, :Q5IncS]

    # Loop over each study and calculate elasticity with a log-log regression.
    for i = 1:n

        # Take logs of variables
        log_income_share = log.(Vector(Data[i, income_shares]))
        log_fuel_total   = log.(Vector(Data[i, fuel_shares]) .* Vector(Data[i, income_shares]))

        # Run log-log regression and store elasticity.
        B, V   = regress(log_fuel_total, log_income_share)
        elasticities[i] = B[2]
    end

    # Using these elasticities, run meta-regression (elasticity vs. ln gdp per capita).
    meta_B, meta_V  = regress(elasticities, log.(Data.GDP))
    meta_intercept  = meta_B[1]
    meta_slope      = meta_B[2]

    return meta_intercept, meta_slope
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
                set_param!(m, :emissions, :MIU, optimal_CO₂_mitigation)
                set_param!(m, :nice_recycle, :global_carbon_tax, optimal_CO₂_tax)
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
                set_param!(m, :emissions, :MIU, optimal_CO₂_mitigation)
                set_param!(m, :nice_recycle, :global_carbon_tax, zeros(n_steps))
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
    save(joinpath(global_path, "global_co2_mitigation.csv"), DataFrame(get_global_mitigation(m_policy, m_bau)))
    save(joinpath(global_path, "temperature.csv"), DataFrame(temperature=m_policy[:climatedynamics, :TATM]))
    save(joinpath(global_path, "rough_carbon_tax_global_algorithm.csv"), DataFrame(global_algorithm_tax=opt_tax_global_algorithm))
    save(joinpath(global_path, "carbon_tax.csv"), DataFrame(tax=opt_tax))
    save(joinpath(global_path, "global_emissions_with_landuse.csv"), DataFrame(emissions=m_policy[:emissions, :E]))
    save(joinpath(global_path, "co2_concentration_ppm.csv"), DataFrame(CO2=m_policy[:co2cycle, :MAT] ./ 2.13))

    # Save Regional Output.
    save(joinpath(regional_path, "gross_output.csv"), DataFrame(m_policy[:grosseconomy, :YGROSS]))
    save(joinpath(regional_path, "nice_net_output.csv"), DataFrame(m_policy[:nice_neteconomy, :Y]))
    save(joinpath(regional_path, "consumption.csv"), DataFrame(m_policy[:nice_neteconomy, :C]))
    save(joinpath(regional_path, "population.csv"), DataFrame(m_policy[:nice_neteconomy, :l]))
    save(joinpath(regional_path, "industrial_co2_emissions.csv"), DataFrame(m_policy[:emissions, :EIND]))
    save(joinpath(regional_path, "damage_cost_share.csv"), DataFrame(m_policy[:nice_neteconomy, :DAMFRAC]))
    save(joinpath(regional_path, "abatement_cost_share.csv"), DataFrame(m_policy[:nice_neteconomy, :ABATEFRAC]))
    save(joinpath(regional_path, "regional_co2_mitigation.csv"), DataFrame(m_policy[:emissions, :MIU]))
    save(joinpath(regional_path, "co2_tax_revenue.csv"), DataFrame(m_policy[:nice_recycle, :tax_revenue]))
    save(joinpath(regional_path, "co2_income_elasticity.csv"), DataFrame(m_policy[:nice_recycle, :CO₂_income_elasticity]))

    # Save Quintile Output.
    save(joinpath(quintile_path, "co2_tax_distribution_Q1.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,1]))
    save(joinpath(quintile_path, "co2_tax_distribution_Q2.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,2]))
    save(joinpath(quintile_path, "co2_tax_distribution_Q3.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,3]))
    save(joinpath(quintile_path, "co2_tax_distribution_Q4.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,4]))
    save(joinpath(quintile_path, "co2_tax_distribution_Q5.csv"), DataFrame(m_policy[:nice_recycle, :carbon_tax_dist][:,:,5]))
    save(joinpath(quintile_path, "base_pc_consumption_Q1.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,1]))
    save(joinpath(quintile_path, "base_pc_consumption_Q2.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,2]))
    save(joinpath(quintile_path, "base_pc_consumption_Q3.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,3]))
    save(joinpath(quintile_path, "base_pc_consumption_Q4.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,4]))
    save(joinpath(quintile_path, "base_pc_consumption_Q5.csv"), DataFrame(m_policy[:nice_recycle, :qc_base][:,:,5]))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q1.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,1]))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q2.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,2]))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q3.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,3]))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q4.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,4]))
    save(joinpath(quintile_path, "post_damage_abatement_pc_consumption_Q5.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_damage_abatement][:,:,5]))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q1.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,1]))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q2.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,2]))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q3.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,3]))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q4.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,4]))
    save(joinpath(quintile_path, "post_recycle_pc_consumption_Q5.csv"), DataFrame(m_policy[:nice_recycle, :qc_post_recycle][:,:,5]))

end
