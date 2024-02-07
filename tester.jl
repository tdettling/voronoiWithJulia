using Test
include("2D_Voronoi_Implementation.jl")
include("VoronoiWithThreads.jl")

# Define a test function
function test_voronoi_small1()
    seeds = [Point(-3, 3), Point(3, -3), Point(-3, -3), Point(3,3)]
    size_of_grid = 6
    voronoi_bound = getBoundary()
    grid = Matrix{Point}(undef, size_of_grid + 1, size_of_grid + 1)
    quadtree = create_tree(voronoi_bound)

    for seed in seeds
        # seed will be closest to itself
        insertNode(quadtree, seed, seed)
    end

    # Create a 2D array to represent the voronoi diagram
    final_diagram = generateVoronoi(quadtree, voronoi_bound, grid)

    # You can customize the actual and expected diagrams based on your expectations
    expected_diagram = fill(Point(-3, 3), size(final_diagram))
    
    @test final_diagram == expected_diagram
end

# Run the test
@testset "Voronoi Diagram Test" begin
    test_voronoi_small1()
end
