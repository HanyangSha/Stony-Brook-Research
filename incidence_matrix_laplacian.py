import numpy as np
import random as rand

eps = 1e-14 #floating point error threshold 

mat_sz = int(input("Enter the matrix size: "))
max_edges = (int)(mat_sz * (mat_sz-1) / 2)
num_edges = (int)(rand.random()*max_edges+1)

print("num edges:", num_edges)

b = np.zeros((mat_sz,num_edges))

for i in range(0,num_edges):
    x = (int)(rand.random()*mat_sz)
    y = 0 
    flag = 1
    while (flag):
        y = (int)(rand.random()*mat_sz)
        if (y != x):
            flag = 0 # get out 
    b[min(x,y)][i] = 1
    b[max(x,y)][i] = -1

'''
f = open("edges.in", "r")

s = f.readline()

mat_sz, num_edges = [int(x) for x in s.split()]
b = np.zeros((mat_sz,num_edges))

for i in range(0,num_edges):
    s = f.readline()
    x, y = [int(x) for x in s.split()] # edge from a to b
    #x -= 1 # let input be 1 indexed 
    #y -= 1
    b[x][i] = 1
    b[y][i] = -1

f.close()
'''

b_T = b.transpose()
laplacian = np.matmul(b, b_T)

print("incidence matrix B: \n", b)
print("laplacian: \n", laplacian)

eigenvalues, _ = np.linalg.eig(laplacian)
eigenvalues.sort()

print(eigenvalues[1])

adj = np.zeros((mat_sz,mat_sz))
for i in range(0,num_edges):
    x = 0
    y = 0
    for j in range(0,mat_sz):
        if (b[j][i] == 1): 
            x = j
        elif (b[j][i] == -1):
            y = j
    adj[x][y] = 1

print("adj: ", adj)


####################

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
    global b 
    global b_T
    global laplacian
    global eigenvalues 

    for i in range(0,len(components)-1):
        new_col = np.zeros((mat_sz,1))
        x = components[i]
        y = components[i+1]
        new_col[x] = 1
        new_col[y] = -1
        b = np.hstack((b,new_col))

    b_T = b.transpose()
    laplacian = np.matmul(b, b_T)
    eigenvalues, _ = np.linalg.eig(laplacian)
    eigenvalues.sort()

    print("New B: ", b)
    print("New L: ", laplacian)

def check_2ndsmallest_eigen(e): 
    if (find_components() == 0): # graph is not connected 
        make_connected()
    if (e[1] > -eps): 
        return 1
    return 0

# add one link to create new graph 
def add_link(b_, e): 
    b2 = b_

    x = (int)(rand.random()*mat_sz)
    row_sum = np.sum(b[x])
    while (row_sum == num_edges or row_sum == -num_edges):
        x = (int)(rand.random()*mat_sz)
        row_sum = np.sum(b[x])

    y = 0
    tag = 1
    while (tag):
        y = (int)(rand.random()*mat_sz)
        for i in range(0,num_edges):
            if (not(b2[x][i] == 1 and b2[y][i] == -1)):
                tag = 0 # get out 
                new_col = np.zeros((mat_sz,1))
                new_col[x] = 1
                new_col[y] = -1
                b2 = np.hstack((b2,new_col))
    
    b_T2 = b2.transpose()
    laplacian2 = np.matmul(b2, b_T2) 
    e2, _ = np.linalg.eig(laplacian2)
    e2.sort()

    #print(eigenvalues[1], " ", e2[1])
    if (e[1] <= e2[1] + eps): 
        return 1 
    return 0 

def check_all_properties():
    print("symmetry: ", check_symmetry(laplacian, mat_sz))
    print("positive semidefinite: ", check_positive_semidefinite(eigenvalues)) 
    print("off diagonal entries are nonpositive: ", check_offdiagonal(laplacian,mat_sz))
    print("every row sum and column sum is zero: ", check_sums(laplacian,mat_sz))
    print("lambda_0 (smallest eigenvalue) = 0: ", check_smallest_eigen(eigenvalues))
    print("For a connected graph, the second eigenvalue (lambda_1) is positive: ", check_2ndsmallest_eigen(eigenvalues))
    print("Adding an edge on an undirected graph produces a graph with equal or greater algebraic connectivity: ", add_link(b, eigenvalues))

check_all_properties()
