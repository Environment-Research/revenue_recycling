@defcomp nice_recycle begin

    # --------------------
    # Model Indices
    # --------------------

    regions                  = Index()                                     # Index for RICE regions.
    quintiles                = Index()                                     # Index for regional income quintiles.

    # --------------------
    # Model Parameters
    # --------------------

    min_study_gdp            = Parameter()                                 # Minimum observed per capita GDP value found in elasticity studies ($/person).
    max_study_gdp            = Parameter()                                 # Maximum observed per capita GDP value found in elasticity studies ($/person).
    elasticity_intercept     = Parameter()                                 # Intercept term to estimate time-varying income elasticity.
    elasticity_slope         = Parameter()                                 # Slope term to estimate time-varying income elasticity.
    elasticity_robustness    = Parameter()                                 # Chooses which elasticity robustness run to perform
    meta_n                   = Parameter()                                 # Sample size of meta analysis (needed for robustness on slope parameters)
    meta_MSE                 = Parameter()                                 # Mean Sq Error in meta analysis (needed for robustness on slope parameters)
    meta_xbar                = Parameter()                                 # Avg. Log GDP per cap in meta analysis (needed for robustness on slope parameters)
    meta_xspread             = Parameter()                                 # Sum of sq. dev. from avg Log GDP per cap (needed for robustness on slope parameters)
    elasticity_lower         = Parameter()                                 # Lower percentile of elasticity estimate for robustness
    elasticity_upper         = Parameter()                                 # Upper percentile of elasticity estimate for robustness
    global_carbon_tax        = Parameter(index=[time])                     # Carbon tax ($/ton CO₂).
    regional_population      = Parameter(index=[time, regions])            # Regional population levels (millions of people).
    ABATEFRAC                = Parameter(index=[time, regions])            # Cost of CO₂ emission reductions as share of gross economic output.
    DAMFRAC                  = Parameter(index=[time, regions])            # Climate damages as share of gross output.
    Y                        = Parameter(index=[time, regions])            # Gross world product net of abatement and damages (trillions 2005 USD yr⁻¹).
    CPC                      = Parameter(index=[time, regions])            # Regional per capita consumption (thousands 2005 USD yr⁻¹)
    industrial_emissions     = Parameter(index=[time, regions])            # Industrial carbon emissions (GtC yr⁻¹).
    quintile_income_shares   = Parameter(index=[regions, quintiles])       # Quintile share of regional income.
    damage_dist              = Parameter(index=[regions, quintiles])       # Quintile share of regional climate damages.
    recycle_share            = Parameter(index=[regions, quintiles])       # Share of carbon tax revenue recycled back to each quintile.
    qc_base                  = Parameter(index=[time, regions, quintiles]) # Pre-damage, pre-abatement cost, pre-tax quintile consumption (thousands 2005 USD yr⁻¹).

    # --------------------
    # Model Variables
    # --------------------

    pc_gdp                   = Variable(index=[time, regions])             # Per capita output net of abatement and damages (2005 USD / person).
    CO₂_income_elasticity    = Variable(index=[time, regions])             # Elasticity of CO₂ price exposure with respect to income.
    tax_revenue              = Variable(index=[time, regions])             # Total carbon tax revenue (2005 USD)
    pc_tax_revenue           = Variable(index=[time, regions])             # Total per capita carbon tax revenue ($1000s/person)
    abatement_cost_dist      = Variable(index=[time, regions, quintiles])  # Time-varying regional CO₂ tax distribution share by quintile.
    carbon_tax_dist          = Variable(index=[time, regions, quintiles])  # Time-varying regional CO₂ abatement cost distribution share by quintile.
    qc_post_damage_abatement = Variable(index=[time, regions, quintiles])  # Post-damage, post-abatement cost per capita quintile consumption (thousands 2005 USD/person yr⁻¹).
    qc_post_tax              = Variable(index=[time, regions, quintiles])  # Quintile per capita consumption after subtracting out carbon tax (thousands 2005 USD/person yr⁻¹).
    qc_post_recycle          = Variable(index=[time, regions, quintiles])  # Quintile per capita consumption after recycling tax back to quintiles (thousands 2005 USD/person yr⁻¹).


    function run_timestep(p, v, d, t)

        for r in d.regions

            # Calculate net per capita income ($/person).
            # Note, Y in $trillions and population in millions, so scale by 1e6.
            v.pc_gdp[t,r] = p.Y[t,r] / p.regional_population[t,r] * 1e6

            # Calculate time-varying income elasticity of CO₂ price exposure (requires pc_gdp units in $/person).
            # Note, hold elasticity constant at boundary value if GDP falls outside the study support range.
            if v.pc_gdp[t,r] < p.min_study_gdp
                # GDP below observed study values.
                v.CO₂_income_elasticity[t,r] = p.elasticity_intercept + p.elasticity_slope * log(p.min_study_gdp)
            elseif v.pc_gdp[t,r] > p.max_study_gdp
                # GDP above observed study values.
                v.CO₂_income_elasticity[t,r] = p.elasticity_intercept + p.elasticity_slope * log(p.max_study_gdp)
            else
                # GDP within observed study values.
                v.CO₂_income_elasticity[t,r] = p.elasticity_intercept + p.elasticity_slope * log(v.pc_gdp[t,r])
            end

            ## Add robustness exercises here where we alter CO₂_income_elasticity
            if elasticity_robustness==1
                v.CO₂_income_elasticity[t,r] = v.CO₂_income_elasticity[t,r] - 1.96*sqrt(p.meta_MSE*(1/p.meta_n + (v.pc_gdp[t,r] - p.meta_xbar)^2/p.meta_xspread)) 
            elseif elasticity_robustness==2
                v.CO₂_income_elasticity[t,r] = v.CO₂_income_elasticity[t,r] + 1.96*sqrt(p.meta_MSE*(1/p.meta_n + (v.pc_gdp[t,r] - p.meta_xbar)^2/p.meta_xspread))
            elseif elasticity_robustness==3
                v.CO₂_income_elasticity[t,r] = p.elasticity_lower
            elseif elasticity_robustness==4
                v.CO₂_income_elasticity[t,r] = p.elasticity_upper
            end

            # Calculate total carbon tax revenue for each region (dollars).
            # Note, emissions in GtC and tax in dollars, so scale by 1e9.
            v.tax_revenue[t,r] = p.industrial_emissions[t,r] * p.global_carbon_tax[t] * 1e9

            # Calculate per capita tax revenue for each region (convert to $1000s/person to match pc consumption units).
            # Note, tax in dollars and population in millions, so scale by 1e9.
            v.pc_tax_revenue[t,r] = v.tax_revenue[t,r] / p.regional_population[t,r] / 1e9

            # Calculate quintile distribution share of CO₂ tax burden and mitigation costs (assume both distributions are equal).
            v.abatement_cost_dist[t,r,:] = regional_quintile_distribution(v.CO₂_income_elasticity[t,r], p.quintile_income_shares[r,:])
            v.carbon_tax_dist[t,r,:]     = regional_quintile_distribution(v.CO₂_income_elasticity[t,r], p.quintile_income_shares[r,:])

            # Create a temporary variable used to calculate NICE baseline quintile consumption (just for convenience).
            temp_C = 5.0 * p.CPC[t,r] * (1.0 + p.DAMFRAC[t,r]) / (1.0 - p.ABATEFRAC[t,r])

            for q in d.quintiles

                # Calculate post-damage, post-abatement cost per capita quintile consumption (bounded below to ensure consumptions don't collapse to zero or go negative).
                # Note, this differs from standard NICE equation because quintile CO₂ abatement cost shares now vary over time.
                v.qc_post_damage_abatement[t,r,q] = max(p.qc_base[t,r,q] - (5.0 * p.CPC[t,r] * p.DAMFRAC[t,r] * p.damage_dist[r,q]) - (temp_C * p.ABATEFRAC[t,r] * v.abatement_cost_dist[t,r,q]), 1e-8)

                # Subtract tax revenue from each quintile based on quintile CO₂ tax burden distributions.
                # Note, per capita tax revenue and consumption should both be in 1000s dollars/person.
                v.qc_post_tax[t,r,q] = v.qc_post_damage_abatement[t,r,q] - (5 * v.pc_tax_revenue[t,r] * v.carbon_tax_dist[t,r,q])

                # Recycle tax revenue by adding shares back to all quintiles (assume recycling shares constant over time).
                v.qc_post_recycle[t,r,q] = v.qc_post_tax[t,r,q] + (5 * v.pc_tax_revenue[t,r] * p.recycle_share[r,q])
            end
        end
    end
end
