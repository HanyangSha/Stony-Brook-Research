import numpy as np
import random as rand

#matrix size
mat_sz = int(input("Enter the matrix size: "))

adj = np.zeros((mat_sz,mat_sz))

#assign random values for adjacency matrix
flag = 1
while (flag):
    for i in range(0,mat_sz):
        for j in range(i+1,mat_sz):
            if (i != j):
                adj[i][j] = (int)(rand.random()*2)
                adj[j][i] = adj[i][j]
                if (adj[i][j] == 1): 
                    flag = 0

# compute degree matrix 
def compute_deg(adj, n):
    d = np.zeros((mat_sz,mat_sz))
    for i in range(0,n): 
        for j in range(0,n):
            d[i][i] += adj[i][j]
    return d 
deg = compute_deg(adj, mat_sz)

#by definition, L = D-A
laplacian = deg-adj

def print_matrices(): 
    print("A: ", adj)
    print("D: ", deg)
    print("L: ", laplacian)

print_matrices()

#compute eigenvalues
eigenvalues, eigenvectors = np.linalg.eig(laplacian)
eigenvalues.sort()

eps = 1e-14 #floating point error threshold 

# checking symmetry
def check_symmetry(l,n): 
    for i in range(0,n):
        for j in range(0,n): 
            if (l[i][j] != l[j][i]): 
                return 0 
    return 1 

#check positive semidefinite
def check_positive_semidefinite(e):
    for i in e:
        if (i < -eps): 
            return 0
    return 1 

#check off diagonal entries are nonpositive
def check_offdiagonal(l,n):
    for i in range(0,n):
        for j in range(0,n):
            if (i != j and l[i][j] > 0): 
                return 0
    return 1 

#check row and colum sum are zero
def check_sums(l,n):
    for i in range(0,n):
        sumr = 0
        sumc = 0
        for j in range(0,n):
            sumr += l[i][j]
            sumc += l[j][i]
        if (sumr != 0 or sumc != 0): 
            return 0
    return 1

#lambda_0 = 0
def check_smallest_eigen(e):
    mn = 1e9
    for i in e:
        mn = min(mn, i)
    if (abs(mn) < eps):
        return 1
    return 0


#check that for a connected graph, the second eigenvalue is positive
visited = np.zeros((mat_sz))
components = []

def dfs(v):
    global visited

    visited[v] = 1
    for i in range(0,mat_sz):
        if (adj[v][i] == 1 and visited[i] == 0):
            dfs(i)

# returns if graph is connected 
def find_components(): 
    global components

    for i in range(0, mat_sz):
        if (visited[i] == 0):
            dfs(i)
            components.insert(len(components), i)
    if (len(components) > 1): # more than 1 component -> graph is not connected 
        return 0
    return 1 

# make the graph connected if it is not currently connected, by connecting disjoint components 
def make_connected(): 
    global adj
    global deg 
    global laplacian
    global eigenvalues 
    global eigenvectors

    for i in range(0,len(components)-1):
        a = components[i]
        b = components[i+1]
        adj[a][b] = 1
        adj[b][a] = 1

    deg = compute_deg(adj,mat_sz)
    laplacian = deg-adj
    eigenvalues, eigenvectors = np.linalg.eig(laplacian)
    eigenvalues.sort()

    print("New: ")
    print_matrices()


def check_2ndsmallest_eigen(e): 
    if (find_components() == 0): # graph is not connected 
        make_connected()
    if (e[1] > -eps): 
        return 1
    return 0

# add one link to create new graph 
def add_link(): 
    adj2 = adj

    cnt = 0
    for i in range(0,mat_sz):
        for j in range(i+1,mat_sz):
            if (i != j and adj[i][j] == 0):
                cnt += 1
    rng = (int)(rand.random()*cnt) 

    cnt = 0
    for i in range(0,mat_sz):
        for j in range(i+1,mat_sz):
            if (i != j and adj[i][j] == 0):
                if (cnt == rng): 
                    adj2[i][j] = 1
                    adj2[j][i] = 1
                else: 
                    cnt += 1

    deg2 = compute_deg(adj2, mat_sz)
    laplacian2 = deg2 - adj2

    e2, _ = np.linalg.eig(laplacian2)
    e2.sort()

    #print(eigenvalues[1], " ", e2[1])
    if (eigenvalues[1] <= e2[1] + eps): 
        return 1 
    return 0 

def check_all_properties():
    print("symmetry: ", check_symmetry(laplacian, mat_sz))
    print("positive semidefinite: ", check_positive_semidefinite(eigenvalues)) 
    print("off diagonal entries are nonpositive: ", check_offdiagonal(laplacian,mat_sz))
    print("every row sum and column sum is zero: ", check_sums(laplacian,mat_sz))
    print("lambda_0 (smallest eigenvalue) = 0: ", check_smallest_eigen(eigenvalues))
    print("For a connected graph, the second eigenvalue (lambda_1) is positive: ", check_2ndsmallest_eigen(eigenvalues))
    print("Adding an edge on an undirected graph produces a graph with equal or greater algebraic connectivity: ", add_link())

check_all_properties()
