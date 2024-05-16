"""
-2   1   0   0   1
 1  -2   1   0   0
 0   1  -2   1   0
 0   0   1  -2   1
 1   0   0   1  -2
"""
function create_FD_matrix(n)
    mat = zeros(Int, n, n)
    
    for i = 1:n
        mat[i, i] = -2
        if i > 1
            mat[i, i-1] = 1
        end
        if i < n
            mat[i, i+1] = 1
        end
    end
    
    mat[1, n] = 1
    mat[n, 1] = 1
    
    return mat
end