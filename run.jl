using HDF5
using ODE

@everywhere include("Ionization.jl")

@everywhere generate(10_000)
