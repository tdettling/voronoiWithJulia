using Pkg
# for scatter plot to visualize voronoi
Pkg.add("Plots")
Pkg.add("GR")

# defining a coordinate structure for our points and seeds
# :: is datatype declaration
struct Coordinate
    x::Int
    y::Int
end

# Defining quadtree node
struct QuadtreeNode
    # defining the boundary of the voronoi diagram with two points,
    # one in the upper right and one in the lower left corner.
    boundary::Tuple{Coordinate, Coordinate}

    # recursive data structure
    children::Vector{QuadtreeNode}

    # storing the data inside our nodes.
    # could be empty if node is non-leaf
    dataInNode::Vector{Coordinate}
end

# function to determine if a point/seed falls in our proper boundaries
function in_boundary(coordinate::Coordinate, bound::Tuple{Coordinate, Coordinate})
    # making sure x is within the boundaries 
    xRange = bound[1].x <= coordinate.x <= bound[2].x

    # making sure y is within the boundaries
    yRange = bound[1].y <= coordinate.y <= bound[2].y

    # returning bool
    return xRange && yRange
end

function split(node::QuadtreeNode)
    # splitting grid
    xCenter = div(node.boundary[1].x + node.boundary[2].x, 2)
    yCenter = div(node.boundary[1].y + node.boundary[2].y, 2)

    # calculating new corners
    up_left = (node.boundary[1], Coordinate(xCenter, yCenter))
    up_right = (Coordinate(xCenter, node.boundary[1].y), Coordinate(node.boundary[2].x, yCenter))
    low_left = (Coordinate(node.boundary[1].x, yCenter), Coordinate(xCenter, node.boundary[2].y))
    low_right = (Coordinate(xCenter, yCenter), node.boundary[2])

    # adding the children
    node.children = [
        QuadtreeNode(up_left, [], []),
        QuadtreeNode(up_right, [], []),
        QuadtreeNode(low_left, [], []),
        QuadtreeNode(low_right, [], [])
    ]

    # adding the data back into the children nodes
    for coordinate in node.dataInNode
        for child in node.children
            insertNode(child, coordinate)
        end
    end

    # removing data from parent
    node.dataInNode = []
end

# inserting node
function insertNode(node::QuadtreeNode, coordinate::Coordinate)
    if in_boundary(coordinate, node.boundary)
        # is empty means that we could have to split
        if isempty(node.children)
            push!(node.dataInNode, coordinate)
            # could switch this to two rather than 4.. depends on the number of seeds.
            if length(node.dataInNode) > 4
                split(node)
            end
        # storing data points in leaf
        else
            for child in node.children
                insertNode(child, coordinate)
            end
        end
    end 
end

# creating a quadtree
function create_tree(boundary::Tuple{Coordinate, Coordinate})
    return QuadtreeNode(boundary, [], [])
end


#########################################################################################################################
## Start of Voronoi Definitions 
#########################################################################################################################


function voronoi_seed_gen(tree::QuadtreeNode, boundary::Tuple{Coordinate, Coordinate}, seeds::Vector{Coordinate})
    # if in a given area, there is a single seed
    if isempty(tree.children)
        # calculates center of region?
        push!(seeds, Coordinate(div(tree.boundary[1].x, tree.boundary[2].x),
                                div(tree.boundary[1].y, tree.boundary[2].y)))
    # more than 1 seed in an area
    else
        for child in tree.children
            voronoi_seed_gen(child, boundary, seeds)
        end
    end
end








#testing
# Example usage

println("START")
bounds = (Coordinate(0, 0), Coordinate(100, 100))
quadtree = create_tree(bounds)

points_to_insert = [Coordinate(10, 20), Coordinate(30, 40), Coordinate(80, 90)]
for point in points_to_insert
    insertNode(quadtree, point)
end

voronoi_points = Vector{Coordinate}(undef,0)
voronoi_seed_gen(quadtree, bounds, voronoi_points)

using Plots
gr()

scatter([point.x for point in voronoi_points], [point.y for point in voronoi_points], ratio = 1, legend = false)
scatter!([point.x for point in points_to_insert], [point.y for point in points_to_insert], markersize = 10, markerstrokewidth = 0)

plot!()
println("END")
