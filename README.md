# Stony-Brook-Research

This repository contains the C++ programs (and two short python programs) that I wrote for research at Stony Brook University as a member of the CSIRE program in summer 2022. In addition, I wrote an [overleaf document](https://www.overleaf.com/read/vdzrmsjtgrrh) detailing what I did for the program. Specifically, all below files are explained in the overleaf. 

Every C++ file uses the C++ linear algebra library Eigen. If you want to use my code, make sure you download this library, otherwise you are responsible for having an alternative way to compute eigenvalues of a matrix. 

Explanations of the function/purpose of each file (somewhat in chronological order):
1. undirected_laplacian: verifies 7 properties (stated in overleaf) of undirected Laplacian matrices under the definition of L=D-A 
2. incidence_matrix_laplacian: verifies the same 7 properties of undirected Laplacian matrices under the definition of L=BB^T
3. connected_builder: randomly generates rooted and strongly connected directed graphs according to L=BB^T, then computes their eigenvalues to verify nonnegative and nondecreasing properties
4. all_rooted_iteration: same as above, except the program builds all possible rooted and strongly connected directed graphs iteratively
5. paper_definition: randomly generates rooted and strongly connected directed graphs according to L=D-A, where A is out degree, then computes their eigenvalues to show that nonnegative and nondecreasing properties do not hold
6. add_reverse_edge: randomly generates rooted and strongly connected directed graphs according to L=BB^T, then adds edges of opposite directions of existing edges to show that Fiedler value may remain the same or increase
7. vertex_connectivity: simulates Theorem 1.10 (regarding a relationship between vertex connectivity and Fiedler value) in Z.-M. Hong et al. / Linear Algebra and its Applications 579 (2019) 72–88 
8. vc_ford_fulkerson: randomly generates connected undirected graph, then uses the max flow min cut theorem and Ford-Fulkerson process to find vertex connectivity of graph. The correctness of the result is then checked by brute force. 
