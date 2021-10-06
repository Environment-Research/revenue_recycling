## Replication Code for Budolfson et al. (2021), "*Climate Action with Revenue Recycling Has Benefits for Poverty, Inequality, and Wellbeing*."

### How To Install Required Packages

This code runs on [Julia v1.6](https://julialang.org/downloads/) and uses the [Mimi](https://www.mimiframework.org/) framework for building integrated assessment models. To install the required packages, first enter the Julia package manager using the `]` key. Once in the package manager, run the following code:

```julia
add CSVFiles  
add DataFrames  
add Interpolations
add Mimi  
add NLopt
```

Run the following lines of code to install the MimiRICE model. This code will create a connection to the online Mimi model library and then install Mimi-RICE2010. An in-depth description of these steps can be found at the [MimiRICE2010 GitHub Repository](https://github.com/anthofflab/MimiRICE2010.jl).

Still in the Julia pacakge manager, run the following lines (note, these only need to be run once per computer):

```julia
registry add https://github.com/mimiframework/MimiRegistry.git
add MimiRICE2010
```

While still in the package manager, run the following line to install MimiNICE:

```julia
add https://github.com/Environment-Research/MimiNICE.jl.git
```

One of the sensitivity analyses also uses the FAIR climate model. Run the following line to install MimiFAIR v1.3:

```julia
add https://github.com/FrankErrickson/MimiFAIR13.jl.git
```

After this, hit the `backspace` key to exit back to Julia.


### Run NICE with Revenue Recycling and View Results

Next, clone or download the `revenue_recycling` Git repository. Once this is on your computer, set this folder as your working directory in the Julia console and run the following code to get an instance of NICE with revenue recycling:

```julia
# Load the file to create an instance of Mimi-NICE with revenue recycling.
include("src/MimiNICE_recycle_time_varying.jl")

# Create an instance of the model.
m = create_nice_recycle()

# Run the model.
run(m)

# View the projected global temperature anomalies.
m[:climatedynamics, :TATM]

# If you want to see interactive plots of all model output, run the following code. Note that it only works for 2-d data (i.e. time x region), it won't show quintile level plots.
explore(m)
```

### Run Paper Replication Code

To reproduce the results from "*Climate Action with Revenue Recycling Has Benefits for Poverty, Inequality, and Wellbeing*," run:
```julia
include("src/paper_replication_runs.jl")
```
The replication output will be saved in the `results` folder.
