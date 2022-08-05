# Stony-Brook-Research

This repository contains the C++ programs I wrote for research at Stony Brook University in summer 2022. 

Every file uses the C++ linear algebra library Eigen. If you want to use my code, make sure you download this library, otherwise you are responsible for finding an alternative way to compute eigenvalues of a matrix. 

Explanations of what each file does (somewhat in chronological order):
1. connected_builder: randomly generates rooted and strongly connected directed graphs according to L=BB^T, then computes their eigenvalues to verify nonnegative and nondecreasing properties
2. all_rooted_iteration: same as above, except it builds all possible rooted and strongly connected directed graphs iteratively
3. paper_definition: randomly generates rooted and strongly connected directed graphs according to L=D-A, where A is out degree, then computes their eigenvalues to show that nonnegative and nondecreasing properties do not hold
4. add_reverse_edge: randomly generates rooted and strongly connected directed graphs according to L=BB^T, then adds edges of opposite directions of existing edges to show that Fiedler value may remain the same or increase
5. vertex_connectivity: simulates Theorem 1.10 (regarding a relationship between vertex connectivity and Fiedler value) in Z.-M. Hong et al. / Linear Algebra and its Applications 579 (2019) 72–88 
6. vc_ford_fulkerson: randomly generates connected undirected graph, then uses the max flow min cut theorem and Ford-Fulkerson process to find vertex connectivity of graph. The correctness of the result is checked by brute force. 
