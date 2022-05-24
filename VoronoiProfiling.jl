# An approximation to the Voronoi Diagram
# This code creates a discrete approximation 
# to the Voronoi Diagram
# This version has a larger area of 4096 x 4096
const SIZE_AREA_P = 4096
const NSEEDS = 4

function calculateDiagram(voronoi, seeds)
    for i in 1:SIZE_AREA_P
        for j in 1:SIZE_AREA_P
            shortestDistance = typemax(Int32)
            closestSeed = 0
            for k in 1:NSEEDS
                distance = (i - seeds[k,1])^2 + (j - seeds[k,2])^2
                distance = sqrt(distance)
                if distance < shortestDistance
                    shortestDistance = distance
                    closestSeed = k
                end
            end
            voronoi[i,j] = closestSeed
        end
    end
end

global voronoi = zeros(Int32,SIZE_AREA_P,SIZE_AREA_P)
global seeds = [1 1;1 SIZE_AREA_P;SIZE_AREA_P 1;SIZE_AREA_P SIZE_AREA_P]

#using Profile

#@profile calculateDiagram(voronoi, seeds)
@time calculateDiagram(voronoi, seeds)

@time calculateDiagram(voronoi, seeds)

#print(voronoi)

#Profile.print()
