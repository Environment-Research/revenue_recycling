# #--------------------------------------------------------------------------------------------------------------------------
# #--------------------------------------------------------------------------------------------------------------------------
# # This file contains functions used for adding revenue recycling to the Nested Inequalities Climate Economy (NICE) model.
# #--------------------------------------------------------------------------------------------------------------------------
# #--------------------------------------------------------------------------------------------------------------------------


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
