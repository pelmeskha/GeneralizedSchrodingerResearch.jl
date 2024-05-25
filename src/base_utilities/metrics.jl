function relative_error_to_amplitude(reference_vector::Vector{<:Number}, vector::Vector{<:Number})
    return maximum( 
        abs.((abs.(reference_vector) - abs.(vector)))
    ) / maximum(abs.(reference_vector)) * 100
end