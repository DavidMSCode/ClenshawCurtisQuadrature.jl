push!(LOAD_PATH, "../src")  #make the package available to the script
using Documenter, ClenshawCurtisQuadrature #load the package and Documenter

#generate the documentation
makedocs(
    sitename="ClenshawCurtisQuadrature.jl",
    doctest = false,    #FIXME: enable once examples are added to the documenation,

)

deploydocs(
    repo = "github.com/DavidMSCode/ClenshawCurtisQuadrature.jl.git",
)