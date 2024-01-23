module MicrotubuleSpringModel

using StaticArrays
using LinearAlgebra
using LinearAlgebra: norm, normalize
using Parameters
using Base.Threads
using Distributions
using Random
using Quaternions
using Configurations: @option, to_toml
using Setfield: @set
using DelimitedFiles
using CSV
using DataFrames
using ProgressMeter: @showprogress

using GLMakie
using ColorSchemes

export 
    BeadPos,

    Lattice,
    BeadPars,
    SpringConst,

    create_lattice,
    initialise,
    initialise_dimer,
    initialise_PF,
    set_bond_angles,

    bond_indices,
    intra_dimer_index,

    total_energy,

    iterate!,

    MicrotubuleConfig,
    PatchConfig,

    surface_area,
    youngs_modulus,
    microtubule_length,
    deflection_end,
    stiffness,

    save_to_csv,
    save_params,
    load_from_csv,

    plot,
    plot_flat!,
    plot_individual!,
    plot_energy!,
    colorschemes,
    COLORSCHEME,

    burnin


BeadPos = MVector{3, Float64}

include("utils/quaternions.jl")

include("initialise/lattice.jl")
include("initialise/neighbours.jl")
include("initialise/angles.jl")
include("initialise/conf.jl")

include("integrate/internal_forces.jl")
include("integrate/external_forces.jl")
include("integrate/energy.jl")
include("integrate/step.jl")

include("simulate.jl")

include("utils/visuals.jl")
include("utils/measurements.jl")
include("utils/io.jl")

end
