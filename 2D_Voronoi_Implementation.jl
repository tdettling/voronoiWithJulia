using Pkg
# for scatter plot to visualize voronoi
Pkg.add("Plots")
Pkg.add("GR")

# defining a coordinate structure for our points and seeds
# :: is datatype declaration
struct Point
    x::Int
    y::Int
end

# Defining quadtree node
struct Node
    # defining the boundary of the voronoi diagram with two points,
    # one in the upper right and one in the lower left corner.
    boundary::Tuple{Point, Point}

    # recursive data structure
    children::Vector{Node}

    # storing the data inside our nodes.
    # could be empty if node is non-leaf
    dataInNode::Vector{Point}
end

# function to determine if a point/seed falls in our proper boundaries
function in_boundary(point::Point, boundary::Tuple{Point, Point})
    isIn = True

    # making sure x is within the boundaries 
    if point.x < boundary[1].x || point.x > boundary[2].x
        isIn = False
    end

    # making sure y is within the boundaries
    if point.y < boundary[1].y || point.y > boundary[2].y
        isIn = False
    end

    return isIn
end

function findClosestSeed(point::Point)
    winner = seeds[1]
    for seed in seeds
        winDistance = sqrt((winner.x - point.x)^2 + (winner.y - point.y)^2)
        seedDistance = sqrt((seed.x - point.x)^2 + (seed.y - point.y)^2)

        if seedDistance <= winDistance
            winner = seed
        end
    end
    return winner
end

function split(node::Node)
    # calculating new corners
    low_right_X = div(node.boundary[1].x + node.boundary[2].x, 2)
    low_right_Y = div(node.boundary[1].y + node.boundary[2].y, 2)

    low_right = (Point(low_right_X, low_right_Y), node.boundary[2])
    up_left = (node.boundary[1], Point(low_right_X, low_right_Y))
    
    up_right = (Point(low_right_X, node.boundary[1].y), Point(node.boundary[2].x, low_right_Y))
    low_left = (Point(node.boundary[1].x, low_right_Y), Point(low_right_X, node.boundary[2].y))

    # adding the children
    node.children = [(low_right, [], []) , (up_left, [], []), (up_right, [], []), (low_left, [], []) ]

    # adding the data back into the children nodes
    # https://discourse.julialang.org/t/help-design-a-node-for-a-tree/67444/11
    # https://stackoverflow.com/questions/41946007/efficient-and-well-explained-implementation-of-a-quadtree-for-2d-collision-det

    for point in node.dataInNode
        for child in node.children
            closeSeed = findClosestSeed(child.dataInNode)
            insertNode(child, point, closeSeed)
        end
    end
    
    # removing data from parent
    node.dataInNode = []

end

# inserting node
function insertNode(node::Node, point::Point, closestSeed::Point)
    # is empty means that there are only refrences to children
    # not a full region to split
    if isempty(node.children)
        push!(node.dataInNode, point)
        # could switch this to two rather than 4.. depends on the number of seeds.
        if length(node.dataInNode) > 4
            split(node)
        end

    # if its not empty, then it has points 
    # and cannot/should not be subdivded further
    else
        for child in node.children
            closeSeed = findClosestSeed(child.dataInNode)
            insertNode(child, point, closeSeed)
        end
    end
end


# creating a quadtree
function create_tree(boundary::Tuple{Point, Point})
    return Node(boundary, [], [])
end


#########################################################################################################################
## Start of Main
#########################################################################################################################
println("START")
global seeds = [Point(1,1), Point(3,9), Point(7,2)]
global size = 10

#getting boundary from size
up_left_X = div((-1 * size), 2)
up_left_Y = div(size, 2)
low_right_X = div(size, 2)
low_right_Y =  div((-1*size), 2)

boundary = (Point(up_left_X, up_left_Y), Point(low_right_X, low_right_Y))
quadtree = create_tree(boundary)

for seed in seeds
    #seed will be closest to itself
    insertNode(quadtree, seed, seed)
end
println("END")

