function parse_to_vector(filename::String)
    file_content = open(filename) do file
        read(file, String)
    end
    str_vector = split(file_content, ", ")
    vector = [eval(Meta.parse(s)) for s in str_vector]
    return vector
end