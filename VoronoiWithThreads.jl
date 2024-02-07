# An approximation to the Voronoi Diagram
# This code creates a discrete approximation 
# to the Voronoi Diagram
# This version has a larger area of 4096 x 4096
# This version uses threads to parallelize the execution of the code

using Base.Threads

const SIZE_AREA_NT = 10
const NSEEDS = 4

global voronoi = zeros(Int32, SIZE_AREA_NT, SIZE_AREA_NT)
global seeds = [(1, 1);
                (1, SIZE_AREA_NT);
                (SIZE_AREA_NT, 1);
                (SIZE_AREA_NT, SIZE_AREA_NT)]

function calculateDiagram(voronoi, seeds, size)
    #@threads
    for j in 1:size
        for i in 1:size
            shortestDistance = typemax(Int32)
            closestSeed = 0
            for k in 1:NSEEDS
                distance = (i - seeds[k][1])^2 + (j - seeds[k][2])^2
                distance = sqrt(distance)
                if distance < shortestDistance
                    shortestDistance = distance
                    closestSeed = k
                end
            end
            voronoi[i, j] = closestSeed
        end
    end
    return voronoi
end

function printGrid(grid_to_print)
    for row_of_grid in eachrow(grid_to_print)
        println(row_of_grid)
    end
end

function toPointGrid(voronoi, seeds)
    newGrid = Tuple{Int, Int}[][]

    for row in 1:size(voronoi, 1)
        row_values = Tuple{Int, Int}[] 
        for col in 1:size(voronoi, 2)
            seedIndex = voronoi[row, col]
            push!(row_values, seeds[seedIndex])
        end
        push!(newGrid, row_values)
    end

    return newGrid
end


newDiagram = calculateDiagram(voronoi, seeds, SIZE_AREA_NT)
grid = toPointGrid(newDiagram, seeds)
printGrid(grid)
print(size(grid))
print(typeof(grid))
