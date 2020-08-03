function load_data(dataset::String)
    if dataset == "facebook_data"
        return load_facebook_data()
    end
end

function load_facebook_data()
    facebook_adj_matrix = fill(0, (4039, 4039))
    open("data/facebook_combined.txt") do file
        for ln in eachline(file)
            pieces = split(ln,' ', keepempty=false)
            facebook_adj_matrix[(parse(Int64,pieces[1]) + 1), (parse(Int64,pieces[2]) + 1)] = 1
        end
    end

    facebook_adj_matrix = sparse(facebook_adj_matrix)
    facebook_adj_matrix = facebook_adj_matrix + facebook_adj_matrix'
    facebook_adj_matrix = min.(facebook_adj_matrix, 1)
    facebook_adj_matrix = MatrixNetwork(facebook_adj_matrix)
end




