using Logging
using FilePaths; using FilePathsBase: /
using Parameters: @unpack
using ProgressMeter
using CairoMakie
using LaTeXStrings

using StaticArrays
using LinearAlgebra
using LinearAlgebra: normalize, norm
using Distributions
using Colors
using GLMakie
using Configurations
using Setfield
using Revise

using DelimitedFiles, CSV

#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using ColorSchemes
theme = include("scripts/theme.jl")
theme isa Attributes && set_theme!(theme)

