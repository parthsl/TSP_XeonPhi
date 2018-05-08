
# PIHC
#### Parallel Iterative Hill Climibing Algorithm for TSP heuristics problem.

#### Project Goal:
Optimisation of PIHC algorithm on various system architectures mainly on IBM POWER, Intel Xeon Phi,NVIDIA GPUs and also on Heterogeneous HPC Cluster. PIHC algorithm should produce output in optimal time with highest accuracy and efficiently utilise HPC systems.

#### Project Outline:
This project is used to solve any Combinatorial Problem or any NP-Hard problem and focus here is made on to the Travelling Salesman Problem. Circuit Desiging, Air-crew scheduling, Vehicle Routing, 3D rendering, and all complex combinatorial task can be mapped to TSP problem. Hence this problem is solved in a novel approach of PIHC which gives better accuracy and time optimal results.

#### TSP: 
Given a set of cities and distance between them, aim is to travel all the given cities in minimal time exactly once and return back to the starting city. This problem is Graph Problem stated as 

> For given Graph G consisting of V vertices and E edges, the problem of finding closed path on the graph starting from a vertex and traveling through all the vertices exactly once and then returning back to the starting vertex is termed as Traveling Salesman Problem.

 - Here problem is solved for Symmetric Graph, such that traveling from vertex A to vertex B is same as traveling from vertex B to vertex A.
 - The input files are in [TSPLIB](https://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/) input format which gives co-ordinates of the cities. Distance between these cities can be find out using [Euclidean Distance Formula](https://en.wikipedia.org/wiki/Euclidean_distance).

Algorithm is combination of two different approaches:

 1.  **Nearest Neighbour(NN):** This is greedy approach to select the nearest city in each hop by iteratively scanning all unvisited city. This approach creates an initial route with 2-approximation Eulerian circuit value. See [NN Aproach](https://en.wikipedia.org/wiki/Nearest_neighbour_algorithm).
 2. **2opt:** For every pair of the edges, find if the swapping these egdes gives the benefit and how much. Pick the best benefit giving swap and swap the edges resulting in optimising tour by the picked swap benefit. See [2-opt method](https://en.wikipedia.org/wiki/2-opt).
One such iteration makes result converge to optimal value, and hence doing so for several iterations unless not getting benefit results in most optimal value in Eulerian Circuit generated.

## BUILD:

```sh
$ make clean 
$ make CC=<compiler-binary-name> OPT="<compiler-optimisation-options>"
```
__Example:__
```sh
$ make CC=gcc OPT="-O3 -maltivec"
```


## Execute:

**Syntax**: 
```sh
$ 2opt <path-to-tsplib-input-file> [number-of-threads]
```
_number-of-threads : Threads for execution (by default to number-of-cpu)._

**Example:** 
```sh
$ ./2opt input/pcb3038.tsp
```

