push!(LOAD_PATH, "../src")  #make the package available to the script
using Documenter, ClenshawCurtisQuadrature #load the package and Documenter

#generate the documentation
makedocs(
    sitename="ClenshawCurtisQuadrature.jl",
)