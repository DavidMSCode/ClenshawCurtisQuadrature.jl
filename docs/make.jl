push!(LOAD_PATH, "../src")  #make the package available to the script
using Documenter, ClenshawCurtisQuadrature #load the package and Documenter

doctest(ClenshawCurtisQuadrature, fix=true)

#generate the documentation
makedocs(
    sitename="ClenshawCurtisQuadrature.jl",
    doctest = true,
    format = Documenter.HTML(assets = ["assets/custom.css"])
)



deploydocs(
    repo = "github.com/DavidMSCode/ClenshawCurtisQuadrature.jl.git",
)