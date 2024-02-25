# imports
#https://docs.julialang.org/en/v1/base/arrays/#Core.Array
# defining a coordinate structure for our points and seeds
using Pkg
Pkg.add("BenchmarkTools")
Pkg.add("Random")
using BenchmarkTools
using Random

# https://lhendricks.org/econ890/julia/user_defined_types.html
# for struct of point idea
# Struct defines a coordniate Point
# copntainms a Point.x and a Point.y call
struct Point
    x::Int
    y::Int
end

# Defining quadtree node
#Node structure to store the quadtree
mutable struct Node
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
# takes in a Point, and a boundary (Tuple of Points)
# returns a bool of whetehr or not a Point is within the provided boundary 
function in_boundary(point, boundary)
    isIn = true

    # making sure x is within the boundaries 
    if point.x < boundary[1].x || point.x > boundary[2].x ||
       point.y < boundary[1].y || point.y > boundary[2].y
        isIn = false
    end
    return isIn
end

# This function takes in a Point and a List of Seeds (Points), 
# and determines the closest need to that point. 
# returns the closest seed (Point)
function findClosestSeed(point, seeds)
    winner = seeds[1]
    winDistance = sqrt((winner.x - point.x)^2 + (winner.y - point.y)^2)

    for seed in seeds
        seedDistance = sqrt((seed.x - point.x)^2 + (seed.y - point.y)^2)
        if seedDistance <= winDistance
            winner = seed
            winDistance = seedDistance
        end
    end
    return winner
end

# this functiuon determines if all four corners of a grid
# (stored within a node)
# have the same closest seed. 
# returns a bool determining if all 4 corners
# have the same cloestest seed. 
function sameSeedCorners(node, seeds)
    up_left_seed = findClosestSeed(Point(node.boundary[1].x, node.boundary[1].y), seeds)
    up_right_seed = findClosestSeed(Point(node.boundary[2].x, node.boundary[1].y), seeds)
    low_left_seed = findClosestSeed(Point(node.boundary[1].x, node.boundary[2].y), seeds)
    low_right_seed = findClosestSeed(Point(node.boundary[2].x, node.boundary[2].y), seeds)
    isSame = true

    if !(isequal(up_left_seed, up_right_seed) && isequal(up_left_seed, low_left_seed) && isequal(up_left_seed, low_right_seed))
        isSame = false
    end

    return isSame
end


# This function takes in a particular node, a grid (matrix of points), amnd a tuple of Points. 
# when the children of a node are empty, that means we can spliut our area. 
# We split our node into 4 children, and repeat our process of splitting
# We incorperate sameSeedCorners() for optimization
function splitNodes(node, grid, seeds)
    bounds = node.boundary
    # Checking to see if we have an odd area to split
    # calculating new corners
    low_right_X = div(bounds[1].x + bounds[2].x, 2)
    low_right_Y = div(bounds[1].y + bounds[2].y, 2)

    low_right = (Point(low_right_X, low_right_Y), bounds[2])
    up_left = (bounds[1], Point(low_right_X, low_right_Y))
    
    up_right = (Point(low_right_X, bounds[1].y), Point(bounds[2].x, low_right_Y))
    low_left = (Point(bounds[1].x, low_right_Y), Point(low_right_X, bounds[2].y))

    # adding the children
    node.children = [
        Node(low_right, [], []),
        Node(up_left, [], []),
        Node(up_right, [], []),
        Node(low_left, [], [])
    ]

    # assigning new boundaries to the children
    node.children[1].boundary = (Point(low_right_X, low_right_Y), Point(bounds[2].x, bounds[2].y))
    node.children[2].boundary = (Point(bounds[1].x, bounds[1].y), Point(low_right_X, low_right_Y))
    node.children[3].boundary = (Point(low_right_X, bounds[1].y), Point(bounds[2].x, low_right_Y))
    node.children[4].boundary = (Point(bounds[1].x, low_right_Y), Point(low_right_X, bounds[2].y))


    # Check if all four corners have the same closest seed
    if sameSeedCorners(node, seeds)
        # If all four corners have the same closest seed, stop subdividing
        # Set each point in the diagram to the closest seed
        node.children = []
        print("--populaitng diagram--")
        populateDiagram(node, grid, seeds)
        return
    end

    # adding the data back into the children nodes
    for point in node.dataInNode
        #closestSeed = findClosestSeed(point, seeds)
        insertNode(node, point, seeds, grid)
    end

    # removing data from parent
    node.dataInNode = []
end

# fucntion that generates a truple of Points 
# used in testing
# retunrs a Tupole of Points
function generateRandomSeeds(numSeeds, boundary)
    #grabbing max values for generation
    minX = boundary[1].x
    minY = boundary[2].y
    maxX = boundary[2].x
    maxY = boundary[1].y

    seedArray = [Point(0, 0)]

    count = 1
    while count <= numSeeds
        newX = rand(minX:maxX)
        newY = rand(minY:maxY)

        newPoint = Point(newX, newY)
        push!(seedArray, newPoint)
        count +=1
    end

    # delete the placeholder fake seed
    popfirst!(seedArray)
    return seedArray
end


# this function inserts a new node into our tree
function insertNode(node, point, seeds, grid)
    # is empty means that there are only references to children
    # not a full region to split
    if isempty(node.children)
        if length(node.dataInNode) >= 4
            print("--splitting nodes--")
            splitNodes(node, grid, seeds)
        end
        push!(node.dataInNode, point)
        return
    # if it's not empty, then it has points 
    # and cannot/should not be subdivided further
    else        
        # Use the same closestSeed for all children
        for child in node.children
            insertNode(child, point, seeds, grid)
            #insertNode(child, point, seeds)
        end
    end
end


# Brute-force function to compute closest seed for each pixel in a given boundary
# will return a filled out grid 
function brute_force_partition(seedList, gridPoints)
    for row in size(gridPoints, 1)
        for col in size(gridPoints, 2)
            shortestDistance = typemax(Int32)
            closestSeed = nothing 
            for seed in seedList
                seedDistance = distance(seed, gridPoints[row, col])
                if seedDistance < shortestDistance
                    shortestDistance = seedDistance
                    closestSeed = seed
                end
            end
            gridPoints[row, col] = closestSeed
        end
    end
    return gridPoints
end

# simple distance function for code clarity
# input is two Points 
function distance(p1, p2)
    return sqrt((p2.x - p1.x)^2 + (p2.y-p1.y)^2)
end

# creating a quadtree
function create_tree(boundary)
    return Node(boundary, [], [])
end

# fucntion that manually populates our new grid based on the closest seed
# recursive nature
function populateDiagram(node, diagram, seeds)
    if isempty(node.children)
        # Leaf node
        for row in eachindex(view(diagram, 1, :))
            for col in eachindex(view(diagram, :, 1))
                newPoint = diagram[row, col]
                closest_seed = findClosestSeed(newPoint, seeds)
                diagram[row, col] = closest_seed
            end
        end
    else
        # Non-leaf node, recursively process children
        for child in node.children
            populateDiagram(child, diagram, seeds)
        end
    end

    return diagram
end




# file reads in data

# https://www.geeksforgeeks.org/opening-and-reading-a-file-in-julia/

function readFile(filePath)
    # needs to return x, y, number of seeds, and a list of seeds
    seeds = Point[]
    println("opening path")
    #f = open(filePath, "r", encoding="UTF-16 LE")
    f = open(filePath, "r")

    println("setting num seeds")
    numSeeds = parse(Int64, readline(f))

    println("reading actual seeds")
    for i in 1:numSeeds
        line = readline(f)

        # Split the line into parts
        point = split(line)

        # Extract x and y coordinates and convert them to integers
        seedX = parse(Int, point[1])
        seedY = parse(Int, point[2])
        newPoint = Point(seedX, seedY)
        # Create a Point object and push it to the array
        push!(seeds, newPoint)
    end
    close(f)
    return seeds
end

# Provides a clean print of any matrix of points
# used with samll test cases in Dev
function printGrid(grid_to_print)
    for row_of_grid in eachrow(grid_to_print)
        println(row_of_grid)
    end
end


#=
function readFile(filePath)
    # needs to return x, y, number of seeds, and a list of seeds
    seeds = Point[]
    println("opening path")

    # Read file with UTF-16LE encoding
    content = read(filePath, String, encoding="UTF-16LE")

    # Split the content into lines
    lines = split(content, '\n')

    # Extract the number of seeds from the first line
    numSeeds = parse(Int, lines[1])

    println("reading actual seeds")
    for i in 2:eachindex(lines)
        line = lines[i]

        # Split the line into parts
        point = split(line)

        # Extract x and y coordinates and convert them to integers
        seedX = parse(Int, point[1])
        seedY = parse(Int, point[2])
        newPoint = Point(seedX, seedY)
        # Create a Point object and push it to the array
        push!(seeds, newPoint)
    end

    return seeds
end

=#


# used to automatically populate the grid with starting points
function assignGrid(grid, voronoi_bound)
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

# used to get an accruate voronoi "box"
# this function will determine the side lengths of a given area. 
# Based on the area (odd or even), 
# we can "shift" the box foirmed on an x,y corrdniate plane such that
# odd areas shoift up one unit while keeping a box shape on the plane.
# for example, a 4x3 box with an origin of (0,0) will have an upper left
# coordniare of (-2,2) and the lower right bound coordniate will be (2, -1). 
# still a square area, just "shifts" the box on the plane. 
function getBoundary(size_of_grid)
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
function generateVoronoi(tree, diagram, seeds)
    # retun in front originally 
    populateDiagram(tree, diagram, seeds)
    return
end


# actual main runner

function mainMain(filepath)
    println("initalizing boudary")
    size_of_grid = 1024
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)
    println("begin reading file")
    seeds = readFile(filepath)

    #populating the voronoi grid with actual points
    println("populate grid, creating quadtree")
    #assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    println("adding in seeds")
    # Adding in seeds
    for seed in seeds
        # Seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    println("start timing")
    # Start timing
    start_time = time()
    println("generateing digram")
    generateVoronoi(quadtree, grid, seeds)

    # End timing
    end_time = time()

    elapsed_time = end_time - start_time

    println(elapsed_time)
end
#########################################################################################################################
## Start of Main
#########################################################################################################################

if length(ARGS) == 0
    println("No command-line arguments provided.")
else
    filepath = ARGS[1]
    mainMain(filepath)
end



#=

println("START")

#initalization of stuff
seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]
size_of_grid = 20
grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
voronoi_bound = getBoundary(size_of_grid)

#seeds = generateRandomSeeds(10, voronoi_bound)

#populating the voronoi grid with actual points
assignGrid(grid, voronoi_bound)
quadtree = create_tree(voronoi_bound)

#adding in seeds
for seed in seeds
    #seed will be closest to itself
    insertNode(quadtree, seed, seed)
end

#starting
println("*****************************************************************************************")
println(" Starting Grid")
printGrid(grid)
println("*****************************************************************************************")

# start of runtime        
sec = @benchmark begin
    global final_diagram
    final_diagram = generateVoronoi(quadtree, grid, seeds)
end


println("*****************************************************************************************")
println("Final Diagram - Nodes")
println("*****************************************************************************************")
printGrid(final_diagram)
println("*****************************************************************************************")

bruteDiagram = brute_force_partition(seeds, grid)
println("*****************************************************************************************")
println("Final Diagram - Brute")
println("*****************************************************************************************")
printGrid(bruteDiagram)
println("*****************************************************************************************")

time_seconds = time(sec)/1e9
println("*****************************************************************************************")
println("Elapsed Time: ", time_seconds, " seconds")
println("*****************************************************************************************")

println("END")
=#
#listT = readFile("C:/Users/tdett/2024_research/forkedJuliaVoronoi/voronoiWithJulia/TestCase1")
#println(listT)


