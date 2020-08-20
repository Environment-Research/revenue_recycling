## Code for NICE + Revenue Recycling

### How To Install Required Packages

This code runs on Julia v1.2. Install these packages if you don't have them. To enter the package manager, hit the `]` key. To exit it back to Julia, hit the `backspace` key. Once in the package manager, run the following code:

```julia
add CSVFiles  
add DataFrames  
add ExcelReaders
add Interpolations
add Mimi  
```

Run the following lines of code to install the MimiRICE model. This code will create a connection to the online Mimi model library and then install Mimi-RICE2010. An in-depth description of these steps can be found at https://github.com/anthofflab/MimiRICE2010.jl

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

```julia
# To reproduce the results from "Progressive Revenue Recycling Can Alleviate Poverty, Reduce Inequality, and Improve Wellbeing While Avoiding Dangerous Climate Change," run the following file:

include("src/user_interface_time_varying_for_paper_submission.jl")
```
