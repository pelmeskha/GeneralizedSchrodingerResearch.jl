coupled_gl_xi(; Δβ, β2, γ, c0, k, ω) = (
    xi1 = (-Δβ - 2*ω/c0 + β2*ω^2 + 2*k) / (β2*c0^2),
    xi2 = -2*γ / (β2*c0^2),
)

common_dependent_args(indep) = (
    δ = (1 - indep.β2*indep.c0*indep.ω) / indep.c0,
    g = 0,
)

"""
    ModelingArgs(independent, dependent)
    ModelingArgs(independent, construct_dependent_args)

Parameters of one analytical solution. `independent` holds the values chosen by
the user, `dependent` holds the values derived from them. Fields are read with
`args.name` from either part.
"""
struct ModelingArgs{I<:NamedTuple, D<:NamedTuple}
    independent::I
    dependent::D

    function ModelingArgs(independent::I, dependent::D) where {I<:NamedTuple, D<:NamedTuple}
        clashes = intersect(keys(independent), keys(dependent))
        isempty(clashes) || throw(ArgumentError("independent and dependent args share keys: $clashes"))
        return new{I,D}(independent, dependent)
    end
end

ModelingArgs(independent::NamedTuple, construct_dependent_args) =
    ModelingArgs(independent, construct_dependent_args(independent))

function Base.getproperty(args::ModelingArgs, name::Symbol)
    independent = getfield(args, :independent)
    dependent = getfield(args, :dependent)
    haskey(independent, name) && return getproperty(independent, name)
    haskey(dependent, name) && return getproperty(dependent, name)
    return getfield(args, name)
end
