using Test
include("2D_Voronoi_Implementation.jl")
include("VoronoiWithThreads.jl")

# Define a test function

function test_voronoi_oneSeed()
    #initalization of stuff
    seeds = [Point(-3, 3)]
    size_of_grid = 8
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

function test_voronoi_small1()
    #initalization of stuff
    seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]
    size_of_grid = 8
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
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
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
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
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

function test_voronoi_med1()
    #initalization of stuff
    seeds = [Point(4, 7), Point(9,1), Point(2,3), Point(0,0), Point(12,2), Point(8,4)]

    #seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]

    size_of_grid = 30
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)
    #seeds = generateRandomSeeds(10, voronoi_bound)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

function test_voronoi_med2()
    #initalization of stuff
    #seeds = [Point(4, 7), Point(9,1), Point(2,3), Point(0,0), Point(12,2), Point(8,4)]

    #seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]

    size_of_grid = 30
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)
    seeds = generateRandomSeeds(10, voronoi_bound)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end
function test_voronoi_med3()
    #initalization of stuff
    #seeds = [Point(4, 7), Point(9,1), Point(2,3), Point(0,0), Point(12,2), Point(8,4)]

    #seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]

    size_of_grid = 200
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)
    seeds = generateRandomSeeds(13, voronoi_bound)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

function test_voronoi_lag1()
    #initalization of stuff
    #seeds = [Point(4, 7), Point(9,1), Point(2,3), Point(0,0), Point(12,2), Point(8,4)]

    #seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]

    size_of_grid = 1024
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)
    seeds = generateRandomSeeds(5, voronoi_bound)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

function test_voronoi_lag2()
    #initalization of stuff
    #seeds = [Point(4, 7), Point(9,1), Point(2,3), Point(0,0), Point(12,2), Point(8,4)]

    #seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]

    size_of_grid = 1024
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)
    seeds = generateRandomSeeds(5, voronoi_bound)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end

function test_voronoi_lag3()
    #initalization of stuff
    #seeds = [Point(4, 7), Point(9,1), Point(2,3), Point(0,0), Point(12,2), Point(8,4)]

    #seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]

    size_of_grid = 1024
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    voronoi_bound = getBoundary(size_of_grid)
    seeds = generateRandomSeeds(5, voronoi_bound)

    #populating the voronoi grid with actual points
    assignGrid(grid, voronoi_bound)
    quadtree = create_tree(voronoi_bound)

    #adding in seeds
    for seed in seeds
        #seed will be closest to itself
        insertNode(quadtree, seed, seeds, grid)
    end

    final_diagram = generateVoronoi(quadtree, grid, seeds)
    expected_diagram = brute_force_partition(seeds, grid)
    
    @test final_diagram == expected_diagram
end


@testset "Voronoi Diagram Test" begin
    # single seed case
    test_voronoi_oneSeed()

    #small batches for comple errors (not real tests)
    test_voronoi_small1()
    test_voronoi_small2()
    test_voronoi_small3()

    # testing voronoi and random seed generateRandomSeeds
    test_voronoi_med1()
    test_voronoi_med2()
    test_voronoi_med3()

    # actual feesable area tests
    test_voronoi_lag1()
    test_voronoi_lag2()
    test_voronoi_lag3()
end
