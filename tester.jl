using Test
include("2D_Voronoi_Implementation.jl")
include("VoronoiWithThreads.jl")

# Define a test function
function test_voronoi_small1()
    #initalization of stuff
    seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]
    size_of_grid = 7
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seed)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    # You can customize the actual and expected diagrams based on your expectations
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

function test_voronoi_small2()
    #initalization of stuff
    seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]
    size_of_grid = 11
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seed)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    # You can customize the actual and expected diagrams based on your expectations
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

function test_voronoi_small3()
    #initalization of stuff
    seeds = [Point(4, 7), Point(9,1), Point(2,3), Point(0,0)]
    size_of_grid = 25
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seed)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    # You can customize the actual and expected diagrams based on your expectations
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

# Run the test
@testset "Voronoi Diagram Test" begin
    #small batches for comple errors (not real tests)
    test_voronoi_small1()
    test_voronoi_small2()
    test_voronoi_small3()
end
