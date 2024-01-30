# An approximation to the Voronoi Diagram
# This code creates a discrete approximation 
# to the Voronoi Diagram
# This version has a larger area of 4096 x 4096
# This version uses threads to parallelize the execution of the code
# These are outputs produced by the execution of this program
# on a Windows machine with an Intel microprocessor with 6 cores
# and 12 Gigabytes of main memory:
# Intel(R) Core(TM) i5-8400 CPU @ 2.80GHz   2.81 GHz
# 12.0 GB (11.8 GB usable)
# In sequential mode:
# C:\Users\Owner\AppData\Local\Programs\Julia-1.7.2\bin\julia.exe --threads 1 VoronoiWithThreads.jl
# 0.293455 seconds (141.52 k allocations: 8.098 MiB, 21.86% compilation time)
# 0.232730 seconds (6 allocations: 688 bytes)
# 1234
# In Parallel - With 6 threads:
# C:\Users\Owner\AppData\Local\Programs\Julia-1.7.2\bin\julia.exe --threads 6 VoronoiWithThreads.jl
# 0.109866 seconds (141.62 k allocations: 8.105 MiB, 62.42% compilation time)
# 0.044498 seconds (32 allocations: 3.703 KiB)
# 1234
# The speedup is almost 6

using Base.Threads
const SIZE_AREA_NT = 4096
const NSEEDS = 4

function calculateDiagram(voronoi, seeds)
    @threads for j in 1:SIZE_AREA_NT
        for i in 1:SIZE_AREA_NT
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

global voronoi = zeros(Int32,SIZE_AREA_NT,SIZE_AREA_NT)
global seeds = [1 1;1 SIZE_AREA_NT;SIZE_AREA_NT 1;SIZE_AREA_NT SIZE_AREA_NT]

#@profile calculateDiagram(voronoi, seeds)
@time calculateDiagram(voronoi, seeds)

@time calculateDiagram(voronoi, seeds)

#print(voronoi)

#Profile.print()
print(voronoi[1,1])
print(voronoi[1,SIZE_AREA_NT])
print(voronoi[SIZE_AREA_NT,1])
print(voronoi[SIZE_AREA_NT,SIZE_AREA_NT])
