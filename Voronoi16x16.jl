# An approximation to the Voronoi Diagram
# This code creates a discrete approximation 
# to the Voronoi Diagram

const SMALL_SIZE_AREA = 16
voronoi = zeros(Int32,SMALL_SIZE_AREA,SMALL_SIZE_AREA)

const NSEEDS = 4
seeds = [1 1;1 SMALL_SIZE_AREA;SMALL_SIZE_AREA 1;SMALL_SIZE_AREA SMALL_SIZE_AREA]
#= seeds = zeros(Int32,4,2)
seeds[1,1]= 1
seeds[1,2]= 1
seeds[2,1]= 1
seeds[2,2]= SMALL_SIZE_AREA
seeds[3,1]= SMALL_SIZE_AREA
seeds[3,2]= 1
seeds[4,1]=SMALL_SIZE_AREA
seeds[4,2]=SMALL_SIZE_AREA =#
# print(typeof(seeds))
for i in 1:SMALL_SIZE_AREA
    for j in 1:SMALL_SIZE_AREA
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
print(voronoi)
