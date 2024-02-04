# imports
using Pkg
Pkg.add("BenchmarkTools")
using BenchmarkTools

# defining a coordinate structure for our points and seeds
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
    isIn = true

    # making sure x is within the boundaries 
    if point.x < boundary[1].x || point.x > boundary[2].x
        isIn = false
    end

    # making sure y is within the boundaries
    if point.y < boundary[1].y || point.y > boundary[2].y
        isIn = false
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


function remove_chunk(bound::Tuple{Point, Point}, chunk_size::Tuple{Int, Int})
    # Extract the corners of the chunk
    chunk_up_left = bound[1]
    chunk_low_right = Point(bound[1].x + chunk_size[1] - 1, bound[1].y - chunk_size[2] + 1)

    # Extract the remaining larger area
    remaining_up_left = Point(bound[1].x, bound[1].y - chunk_size[2])
    remaining_low_right = Point(bound[2].x, bound[2].y)

    return ((remaining_up_left, remaining_low_right), (chunk_up_left, chunk_low_right))
end


function split(node::Node)
    bounds = node.boundary

    # Checking to see if we have an odd area to split
    checkBound = oddDivision(bounds)
    boundPrimary = checkBound[1]
#=
    # If area is odd, calculate the closest seed by brute force for the smaller area
    if length(checkBound) > 1
        # Extracting the removed smaller chunk
        (_, removed_chunk) = remove_chunk(boundPrimary, (1, 1))

        # Brute-force computation for the removed smaller chunk
        brute_force_result = brute_force_partition(removed_chunk, seeds)

        # Add the brute-force result to the quadtree
        for point in brute_force_result
            insertNode(node, point, findClosestSeed(point))
        end
    end
=#

    # calculating new corners
    low_right_X = div(boundPrimary[1].x + boundPrimary[2].x, 2)
    low_right_Y = div(boundPrimary[1].y + boundPrimary[2].y, 2)

    low_right = (Point(low_right_X, low_right_Y), boundPrimary[2])
    up_left = (boundPrimary[1], Point(low_right_X, low_right_Y))
    
    up_right = (Point(low_right_X, boundPrimary[1].y), Point(boundPrimary[2].x, low_right_Y))
    low_left = (Point(boundPrimary[1].x, low_right_Y), Point(low_right_X, boundPrimary[2].y))

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

function oddDivision(boundary::Tuple{Point, Point})
    boundarySet = Tuple{Tuple{Point, Point}, Tuple{Point, Point}}

    #placeholders for values
    #new larger area
    new_up_right = boundary[1]
    new_low_left = boundary[2]

    # small area that we cut off of the larger area
    extra_up_right = boundary[1]
    extra_low_left = boundary[2]

    if ((boundary[2].x - boundary[1].x) % 2 == 0) && ((boundary[2].y - boundary[1].y) % 2 == 0)
        boundarySet = (boundary, boundary)
    else
        # horozontal length odd
        if (boundary[2].x - boundary[1].x) % 2 != 0
            new_up_right = (boundary[1].x + 1, boundary[1].y)
            extra_low_left = (boundary[1].x + 1, boundary[2].y)
        end
        # vertical length odd
        if (boundary[2].y - boundary[1].y) % 2 != 0
            new_low_left = (boundary[2].x, boundary[2].y + 1)
            extra_up_right = (boundary[1].x, boundary[2].y - 1)
        end
        #returning two areas to calculate the diagram with
        boundarySet = ((new_up_right, new_low_left), (extra_up_right, extra_low_left))
    end

    return boundarySet
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

# given two bounary sets, make the larges box and return it. 
function maxBoundary(bound1::Tuple{Point, Point}, bound2::Tuple{Point, Point})
    maxBoundary = Tuple{Point, Point}


    # new upper right
    max_X1 = min(bound1[1].x, bound2[1].x)
    max_Y1 = max(bound1[1].y, bound2[1].y)

    # new lower left
    max_X2 = max(bound1[2].x, bound2[2].x)
    max_Y2 = min(bound1[2].y, bound2[2].y)

    maxBoundary = (Point(max_X1, max_Y1), Point(max_X2, max_Y2))
    return maxBoundary
end

# Brute-force function to compute closest seed for each pixel in a given boundary
function brute_force_partition(boundary::Tuple{Point, Point}, seeds::Vector{Point})
    closest_seeds = Vector{Point}()

    # for the entrie x-axis of the boundary
    for x in boundary[1].x:boundary[2].x
        #for the entire y axis of the boundary
        for y in boundary[1].y:boundary[2].y
            curPoint = Point(x, y)
            closest_seed = findClosestSeed(curPoint)
            # add closest seed for that point in a list
            push!(closest_seeds, closest_seed)
        end
    end

    return closest_seeds
end

# creating a quadtree
function create_tree(boundary::Tuple{Point, Point})
    return Node(boundary, [], [])
end

# moving through tree to calculate grid
function populateDiagram(node::Node, boundary::Tuple{Point, Point}, diagram::Matrix{Point})
    # leaf check - then no need to traverse tree. 
    # manually add
    if isempty(node.children)
         # for the entrie x-axis of the boundary
        for x in boundary[1].x:boundary[2].x
            #for the entire y axis of the boundary
            for y in boundary[1].y:boundary[2].y
                curPoint = Point(x, y)
                closest_seed = findClosestSeed(curPoint)
                # populating diagram
                diagram[x,y] = closest_seed
            end
        end
        return diagram
    else
        # non-leaf, stores refrences t children. 
        # need recursivly break down to manually add
        for child in node.children
            newBound = (child.boundary[1], child.boundary[2])
            populateDiagram(child, newBound, diagram)
        end
    end
end

function readFile(file)
    # needs to return x, y, number of seeds, and a list of seeds
end

function printGrid()
    for row_of_grid in eachrow(grid)
        println(row_of_grid)
    end
end

function assignGrid()
    # assigning rows and cols for readability
    rows = size(grid, 1)
    columns = size(grid, 2)


    # looping through grid
    for row in 1:rows
        for column in 1:columns
            # assigning actual coordinate to space
            # -1 for row/column
            newX = voronoi_bound[1].x + (column-1)
            newY = voronoi_bound[1].y - (row-1)
            grid[row, column] = Point(newX, newY)
        end
    end
end

function getBoundary()
    new_bound = (Point(0,0), Point(0,0))

    if size_of_grid % 2 != 0
        # we need to shift grid up to get accurate oddDivision
        up_left_X = div((-1 * size_of_grid), 2) 
        up_left_Y = div(size_of_grid , 2) + 1
        low_right_X = div(size_of_grid, 2) + 1
        low_right_Y =  div((-1*size_of_grid), 2)

        new_bound = (Point(up_left_X, up_left_Y), Point(low_right_X, low_right_Y))
    else
        up_left_X = div((-1 * size_of_grid), 2)
        up_left_Y = div(size_of_grid , 2)
        low_right_X = div(size_of_grid, 2)
        low_right_Y =  div((-1*size_of_grid), 2)

        new_bound = (Point(up_left_X, up_left_Y), Point(low_right_X, low_right_Y))
    end
    return new_bound
end

# main runner for timing
function generateVoronoi(tree::Node, boundary::Tuple{Point, Point}, diagram::Matrix{Point})
    # Traverse the tree and populate voronoi
    grid = populateDiagram(tree, boundary, diagram)
    return grid
end

#########################################################################################################################
## Start of Main
#########################################################################################################################

println("START")

global seeds = [Point(1,1), Point(3,9), Point(7,2)]
global size_of_grid = 6
global grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
global voronoi_bound = getBoundary()

#getting boundary from size_of_grid
#=up_left_X = div((-1 * size_of_grid), 2)
up_left_Y = div(size_of_grid + 1, 2)
low_right_X = div(size_of_grid + 1, 2)
low_right_Y =  div((-1*size_of_grid), 2)

global voronoi_bound = (Point(up_left_X, up_left_Y), Point(low_right_X, low_right_Y))
=#

assignGrid()
quadtree = create_tree(voronoi_bound)

for seed in seeds
    #seed will be closest to itself
    insertNode(quadtree, seed, seed)
end

# Create a 2D array to represent the voronoi diagram
printGrid()

# start of runtime        
sec = @benchmark begin
    diagram = generateVoronoi(quadtree, voronoi_bound, grid)
end

time_seconds = time(sec)/1e9
println("Elapsed Time: ", time_seconds, " seconds")

println("END")
