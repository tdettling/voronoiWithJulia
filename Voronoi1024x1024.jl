# An approximation to the Voronoi Diagram
# This code creates a discrete approximation 
# to the Voronoi Diagram
# This version has a larger area of 1024 x 1024

const SIZE_AREA = 1024
global voronoi = zeros(Int32,SIZE_AREA,SIZE_AREA)

const NSEEDS = 4
global seeds = [1 1;1 SIZE_AREA;SIZE_AREA 1;SIZE_AREA SIZE_AREA]
#= seeds = zeros(Int32,4,2)
seeds[1,1]= 1
seeds[1,2]= 1
seeds[2,1]= 1
seeds[2,2]= SIZE_AREA
seeds[3,1]= SIZE_AREA
seeds[3,2]= 1
seeds[4,1]=SIZE_AREA
seeds[4,2]=SIZE_AREA =#
# print(typeof(seeds))
for i in 1:SIZE_AREA
    for j in 1:SIZE_AREA
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
#print(voronoi)
