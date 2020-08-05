function load_data(dataset::AbstractString)
    if dataset == "facebook_data"
        return load_facebook_data()
    end
end

function read_txt(filename::AbstractString)
    adj_matrix = fill(0, (4039, 4039))
    open(filename) do file
        for ln in eachline(file)
            pieces = split(ln,' ', keepempty=false)
            adj_matrix[(parse(Int64,pieces[1]) + 1), (parse(Int64,pieces[2]) + 1)] = 1
        end
    end

    adj_matrix = sparse(adj_matrix)
    adj_matrix = adj_matrix + adj_matrix'
    adj_matrix = min.(adj_matrix, 1)
    adj_matrix = MatrixNetwork(adj_matrix)

    return adj_matrix
    
end

# function load_data(dataset::AbstractString)
#     pathname = joinpath(dirname(dirname(@__FILE__)),"data")
#     file = joinpath(pathname, "$(dataset).txt")

#     if isfile(file)
#         return load_facebook_data()
#     else
#         error("Dataset does not exist")
#     end
# end

function add_data(filename::AbstractString)
    filename_begins = findlast(isequal('/'), filename)
    fname_without_directory = filename[filename_begins:length(filename)]
    cp(filename, joinpath("../data/", fname_without_directory), force=true)
end

function load_facebook_data()
    pathname = joinpath(dirname(dirname(@__FILE__)),"data")
    fname = joinpath(pathname, "facebook_data.txt")
    facebook_adj_matrix = fill(0, (4039, 4039))
    open(fname) do file
        for ln in eachline(file)
            pieces = split(ln,' ', keepempty=false)
            facebook_adj_matrix[(parse(Int64,pieces[1]) + 1), (parse(Int64,pieces[2]) + 1)] = 1
        end
    end

    facebook_adj_matrix = sparse(facebook_adj_matrix)
    facebook_adj_matrix = facebook_adj_matrix + facebook_adj_matrix'
    facebook_adj_matrix = min.(facebook_adj_matrix, 1)
    facebook_adj_matrix = MatrixNetwork(facebook_adj_matrix)

    return facebook_adj_matrix
end




