Contribution can be done by implementing TSP solving meta Heuristic algorithms highly optimized for Intel Xeon processors.

Requirements include:
- Genetic algorithm : Sequential and Parallel implementation of Genetic Algorithm to solve TSP problems to approximate values with minimum time.
- Max-Min Ant System  Optimization : Parallel implementation of MAAS Algorithm which uses MPICH to parallel execute in cluster of any Architecture.
- CUDA+KNL implementation of 2-opt Method : Algorithm to solve TSP using 2-opt algorithm in CUDA and Xeon PHI 2 Knights Landing Architecutre.
- Ant-Colony Optimization : Parallel implementation of ACO on Knights Landing Intel Processors and it's clusters.


Each code should achieve following four optimization levels:
1) Scaler Tuning : Code should effectively use cache to reduce False sharing. Peeled loop should be ignored as much as possible.
2) Vectorization : Code should effectively use AVX-256 and AVX-512 technology.
3) OpenMP : Code should take advantage of Multi-Core systems.
4) Memory : Code should optimize for least usage of disk memory and high use of High Bandwidth Memory.

Any other Strategy, Algorithms and Issues are welcomed.
