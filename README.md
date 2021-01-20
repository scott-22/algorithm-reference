# Algorithms and Data Structures

Personal reference for data structures and algorithms for competitive programming. This would be awful to learn from, and it’s meant to be a refresher only. It’s also perpetually work-in-progress.

The names of algorithms or important data structures will be bolded. The names of problems (not actual names in any scientific literature, just for my personal ease of access) will be in unbolded caps, start and end with *, and be provided at the beginning of an algorithm description.

### Table of Contents
[Fundamental Ideas](#fundamentalideas)  
[Math](#math)  
[Arrays](#arrays)  
[Sequences](#sequences)  
[Recursion](#recursion)  
[Strings](#strings)  
[Graphs](#graphs)  
[Dynamic Programming](#dynamicprogramming)  



<a name="fundamentalideas"/>

# Fundamental Ideas
### Overview:
Note: For a problem, don’t just use the “correct” approach. Think from many perspectives, use multiple algorithm designs, and master general ideas and analysis over memorizing strategies. Understand the rationale for each strategy.

Algorithm := finite sequence of steps that will stop and yield a result (alternatively some def using turing machines)  
Invariants := a property which never changes throughout a certain transformation, such as the steps of an algorithm  
Algorithms can be proved through maintaining an invariant, then proving the invariant holds after each passing through the algorithm. At the end of the algorithm, the desired result should be deduced from the invariant.

Big O notation := _f_ is in class O(_g_) if f(_x_) < c * g(_x_) as _x_ → infinity for some _c_ ∈ R+, assuming _f_ and _g_ are positive functions  
Algorithm speed is evaluated using complexity analysis, seeing how running time is affected as input size grows. The running time is expressed as a function of the input size, and that function is placed in an O complexity class. When f is the sum of multiple functions, it has complexity equal to the fastest growing of them.

Program flow can be changed through goto statements and loops, or functions (being recursively called). Most problems solved with recursion can be solved with a loop, and the math behind what can be solved using what is pretty interesting.

Data structure := a hypothetical structure that stores data in a certain way for reading/modification  
Data structures are usually an implementation of a mathematical concept, or at least can be described mathematically. Each data structure has different properties, which can be used to solve problems, or combined to form new data structures.

<a name="math"/>

# Math
### Overview:
All problems involve some form of math, but some problems are exclusively math. Usually the categories of math covered is number theory, combinatorics, algebra, discrete mathematics, and some other stuff. Even geometry is sometimes used. Moral of the story is get good at math.

Note: some nice primes for hashing: 31, 1e9+7, 1<<31 - 1, 1e11+3, 1<<61 - 1

### Coordinate/Array index compression:
#### \*COORDCOMPRESS*
Instead of storing empty coordinates, only store the ones that contain info. One way to do this is to find all the coordinates, then sort them and put them into an array/dynamic array. Then, **binary search** them to find the index of a coordinate.

### Matrices:
Can be used in fibonacci-like sequences. Will fill in later

#### \*MATRIXMULTIPLY*

### Sets
This heavily borrows from set theory. A set data structure is extremely useful, and can be implemented as a simple bool “taken” array, a binary tree, or a hash table. Some other kinds of sets are outlined below.

#### \*DISJOINT-SET*
A **disjoint set** is just a set with no common elements (or basically a normal set). What makes it unique is that it’s able to keep track of multiple elements and what set each of those elements belong to. Each set of elements S is represented by some “representative” x ∈ S. To merge sets S1 and S2, we simply add the representative of S1 to S2, or vice versa. That’s why it’s best to represent it as a forest, where each parent node is a leader of a set of all nodes that are descended from it. In practice, we implement it as a parent array. A disjoint set must support two operations: finding the representative of an element/node (what set it belongs to) and conducting a union of two sets given an element of each.

In order to find the representative of a node, we need to go up the tree, which could become very deep. We can employ **path compression**, a method where we set the leader of all descendents of a node to the root, reducing the depth of that specific path to 1. This can be done easily with a **DFS**.

Another technique to reduce the running time of a disjoint set is **union by rank**. We keep track of the maximum possible depth of any set, and then when merging two sets, we use the representative from the set with the least maximum depth. With these two modifications, the running time becomes almost an amortized constant (I have no idea why it’s the inverse of the Ackermann function, probably some advanced maths).

### Binary exponentiation:
a^b = (a^(b/2))^2

Just apply that recursively for logarithmic time exponentiation.



<a name="arrays"/>

# Arrays
### Overview:
In the most general sense, there are three categories of things we can do with arrays. The first are queries, or some property we want to find. The second are transformations, or some action we perform on the array. The third is creating an array to model a problem.

There are multiple useful properties to be found in arrays. Common queries include
- a specific value
- a pair that sums or subtracts to a value
- individual elements that sum to a value (k-sum problem) or satisfy some other property
- the sum of a certain subarray
- a subarray that sums to a value
- the minimum/nearest minimum in a subarray
- the longest/shortest suffix/prefix/subarray with a certain property
and more. For queries, it’s important to analyze what each property of the query is, and how that property can be efficiently computed, before considering how to solve it.

Common transformations include
- sorting by some property
- multiplying/dividing/adding/subtracting certain subarrays by a value
- separating the array into categories with O(1) update and O(1) retrieval
and more. The same logic applies to transformations, where a certain property is used to make it efficient.

## Properties/Queries
Use some sort of pattern to prevent searching through every possibility. If there is no pattern, sort it. For a single value, use binary search. For a pair use a two pointers approach. For range queries, use a prefix array or a tree approach. Be able to use other data structures to store relevant information that can solve the problem at hand.

### Searching:
#### \*SEARCH-VALUE*
When searching for a single value, the worst case is to look through each value, giving O(n) time. For an unsorted array, there isn’t any structure, so we can’t do it efficiently. For a sorted array, we can do **binary search** in O(log n) time.

#### \*SEARCH-TUPLE*
When searching for a pair that either adds or subtracts to a value, the worst case is O(n^2) time, pretty intuitively. However, in a sorted array, we can use **two pointers** to do it in O(n) time. This works since we always move the pointers towards each other. In fact, for subsets of length k (i.e, k random values of the array that add or subtracts to a value), we can do it in O(n^(k-1)) time, since we can go through each element, then take the remaining elements k-1 at a time.

There is another approach to this problem, and it can be a bit more generalized. For searching a pair that sums up to s, you can iterate through every element i, then do a **binary search** for s-a[i]. Thus you can find a pair in O(n log n) time. In fact, for subsets of length k, we can do it in O(n^(k-1) * log n) time, which is always log n times slower than the multiple pointers approach, but the idea can be applied to relatively more problems.

#### \*MAX-SUBARRAY-VALUE*
Let’s say we had an array to be split into k contiguous subarrays, and we’d like to maximize the sum of any subarray. The naive solution would be to iterate through and find each of the k subarrays, which even if we used a **PSA** to find the sum, it would still take O(n^(k-1)) time. Rather, we could do a binary search of the maximum value. Let’s say the total sum of the array was z. We would binary search all the values between 1 and z, and set that as the maximum. Then, we can iterate through the **PSA** and mark off whenever we surpass that value, and see if we have more or less groups than the target value k. Thus, we would only need O(n log z) time, which can be much less depending on the value of z.

### Range Queries:
Obviously the worst case for doing a range query is just to go through each of the values.

#### \*RSQ*
For range sum queries, and some others, we can do a **prefix sum array (PSA)**, giving us O(1) query complexity and O(n) update complexity. We can apply this to a 2D array too, making a **2D prefix sum array** that allows calculation of a rectangular sum query using the principle of inclusion exclusion. It gets the sum of the rectangle defined by [0][0] as the top left corner, and [x][y] as the bottom right, which is the 2D equivalent of a prefix sum. Then we can just manipulate the four corners of the rectangle we want to sum. It calculates the query in O(1) complexity, and initializes/updates in O(n^2) complexity.

Another way for range sum queries, and some others, is to use a tree, such as a **binary indexed tree/Fenwick tree**. The hypothetical tree stores certain precomputed range queries that can be added together to fulfill the complete query. It’s stored as an array, so it takes up the same space as a prefix sum array, and allows O(log n) query and update complexity, and O(n log n) initialization, making it really good. The reason it works out is because it uses binary representations of numbers in a smart way. We can apply this to a **2D BIT** as well.

#### \*RMQ*
For range minimum queries, and some others, one way is to use a **sparse table**. This table precomputes a property for a subarray starting at every possible position, with length of all powers of 2. This way, it can either break every range query down into powers of 2 (storing it as binary) for an O(log n) query, or if the subarrays are able to overlap as in the case of a \*RMQ*, an O(1) query. It takes up O(n log n) space. The table can be computed in O(n log n) time by using a kind of dynamic programming, using two lengths of 2^n to make a length of 2^(n+1).

The sparse table can also compute the **lowest common ancestor** of trees.

#### \*RQ-MISC1*
Let’s say we wanted to find a property of all elements in an array except in a certain subarray, let’s say a[i, j]. We can create a prefix and suffix array for that property, and find the property of the subarray a[0, i] with the prefix array, and the property of the subarray a[i, end] with the suffix array, and combine the two ranges.

### Inversion count:
#### \*INVERSION*
An inversion of array A is two indices i and j such that i < j and A[i] > A[j]. Inversions can be counted through brute force in O(n^2) complexity, by iterating through each element, and checking for the number of elements to the left of it (index less than it) and a value greater than it. We can do this faster with a **BIT**.

### Subsequences and Subarrays:
#### \*LIS*
For longest increasing/decreasing subsequences, the naive solution is recursive, and in O(n!) time (I think). However, we can use **DP** for a more efficient solution.

For many problems on subarrays, especially when involving their cumulative sum, it’s helpful to turn them into a **PSA** and then use one of the searching algorithms. For example, take the painter’s partition problem.

#### \*PAINTERS-PARTITION*
The most generalized form of this problem is to split an array into k contiguous subarrays such that
1. The difference between the largest and smallest sums of all subarrays is minimized,
or
2. The largest sum of all subarrays is minimized.

Both of these should be equivalent (have not proved yet), but lead into different methods. 

For 1., after converting the array into a **PSA**, you can simply do the \*SEARCH-TUPLE* problem, and find a tuple of indices using binary search (I don’t think two pointers works here).
For 2., after converting the array into a **PSA**, you can do the \*DISTRIBUTE-VALUE* problem, also using binary search.

## Transformations
Different sorts of transformations can be applied to arrays and subarrays. The minimum complexity is O(n), so we generally try to achieve this with some smart tricks.

### Sorting:
#### \*SORT*
Involves two major components: comparisons and re-arrangement of the array. Minimizing the running time of the entire sort involves minimizing running time of both of them.

Lower bound: n log n comparisons (do proof later)  
In practice, just use the built-in sort function.  
Lots of arrays can be sorted to avoid going through all subsets, and instead using structure to be more efficient.  

### Operations on subarrays:
#### \*RUQ*
When performing a range update query, it’s only possible to do it in O(n) time, since you have to iterate through and modify each item. However, if we need to do this multiple times on multiple subarrays, let’s say m times in total, we can do it in total in O(n) time instead of O(n*m), for certain operations like addition and multiplication. We do this by writing out a sort of **difference array**, and putting the operation and inverse of the operation at the start index and one past the end index. Just like the prefix sum array, we can apply this in two dimensions, for a O(n^2) update at the end.

If the operations are something like adding to subarrays, and there are significantly less operations than array elements, we can pair together the operations and store them in a separate array, as a kind of **coordinate compression**. For example, if we were to add m from index i to j, we would put pairs (i, m) and (j+1, -m) into the other array. At the end, we sort the other array, and iterate through the operations instead of each element.




<a name="sequences"/>

# Sequences
### Overview:
I’ve put any sequence type problem in here that I couldn’t fit anywhere else. I define a sequence type problem to be any problem that involves multiple objects or things that have some sort of ordering (similar to an array, but either with more abstract types or with inputs far too large to store as an actual array).

### Overlapping Intervals:
Many problems involve multiple intervals that either represent a transformation on a number line (e.g. a sort of difference array, but with far too many elements to create an actual array), or a distance, period of time, etc.

#### \*MERGE-INTERVAL*
 Let’s say we had multiple intervals that needed to be merged, represented as an ordered pair (_x_, _y_) that represents the beginning and end of the interval, respectively. The naive way would be to select every pair of intervals, and if they were overlapping, we merge them into a new one. This would work in O(n^2) time.

A better solution (and a common theme) is to sort them by x. Then, we can simply do a transversal through each of the intervals, see if it overlaps the last one, and then merge it with the last one if it does. This works in O(n log n) time.

#### \*TRANSFORM-INTERVAL*
When we need to apply many transformations to an array, we usually use a difference array. We can use a similar idea to apply many transformations along intervals to an sequential or array-like object (such as a number line). Assume we had multiple intervals in the form (_x_, _y_) that add some number s to each of the elements in the interval. The naive solution is to simply iterate through each of the intervals and add them, resulting in O(n*k) time, where n is the number of intervals and k is the length of the sequence.

We can again speed this up by sorting, and by drawing inspiration from the difference array. We express each interval as two ordered pairs: (_x_, _s_) and (_y_+1, -_s_). Then, we sort each pair by the first value, and iterate through each of them while keeping a total sum.

### Scheduling
Many task scheduling problems can be solved by expressing each task as an ordered pair and then sorting them. The pair could be (starting time, duration), or (deadline, duration). This makes them essentially like intervals, but in this case it may be to put the most amount of tasks in a set amount of time, or some other thing.



<a name="recursion"/>

# Recursion
### Overview:
Recursion in CS borrows from the idea from mathematics, being a sort of self referencing, and is implemented with putting function calls on “The Stack”. Given that an algorithm should be finite, recursive functions must have base cases in programming.

There are some problems that are inherently recursive, but most have an iterative solution. The power of recursion lies in it usually being more straightforward, and is especially useful if it can reduce the time complexity of a problem (although the overhead is greater). 

There are some interesting uses of recursion, which are highlighted below.

### Combinatorics:
Many useful combinatorial properties can be obtained through recursion. For example, we can find all permutations or combinations of an array.

#### \*PERMUTATION*
All permutations can be achieved by a function that iterates through all elements of an array, stores that element in a vector and marking it as taken, then calls itself again. Once the vector is the same size as the array, you can process the elements, pop the last element, and return.

#### \*COMBINATION*
All combinations can be achieved by a function that takes the _i_ th element, either takes it and adds it to a vector or skips over it, then recursively calls itself for the _i_+1 th element. Once _i_ reaches the end of the array, then you can process the elements, pop the last element, and return.


### Indirect Recursion
Also called mutual recursion. It occurs where two or more functions recursively call one another instead of themselves directly.


### Backtracking
**Backtracking** is the process of solving a problem with a kind of brute force, recursive approach. It resembles **DFS**, where we start at a beginning state, and want to reach a success state. At each call of the solving function, we take every single possible path to the next state. If we ever reach a failed state, we go back to the last non-failed state and continue with the other paths. It’s slow, but it works.


<a name="strings"/>

# Strings
### Overview:
String hashing :(
<a name="graphs"/>

# Graphs
### Overview:
Borrows from graph theory in math. A special type of connected graph with no cycles is called a “tree”. Multiple trees makes a “forest”.

## Trees
A few definitions to get started.
Tree := a connected graph with any n vertices and n-1 edges  
Root := the node considered to be the starting point of the tree, with depth 0  
Depth := the depth of a node is the distance from that node to the root, while the depth
of a tree is the maximum depth of any node in the tree  
Height := the height of a node is the maximum distance from that node to a leaf node,
while the height of a tree is the height of the root  
Diameter := the longest possible path(s) through a tree, or alternatively, the length of
that path  
Center := either one or two nodes in the tree such that the depth is minimized (note the
	center is thus the centre node(s) of the diameter)  
Radius := the longest possible path(s) starting at a center of the tree

### Center, Diameter, Radius:
These ideas are probably the most important in terms of solving tree-specific problems.

Finding the diameter is pretty simple. Start at a random node _s_ and find the farthest node from it, called _u_. This can be done with **BFS** or **DFS**. Then root the tree at _u_, and find the farthest node _v_, repeating the same process. _u_ and _v_ are then endpoints of a diameter. To find the diameter itself, during the second transversal starting at _u_, keep track of each node’s parent. Then, follow the path backwards from _v_ to _u_.

Once the diameter path is obtained, the center and radius are trivial. Some problems that require knowing the centre, diameter, and radius are below:

#### \*DEPTH-QUERY*
To get the depth of a tree by arbitrarily rooting it any node, let’s say _u_, we need to find the maximum distance from _u_ to any other part of the tree. This distance will always be part of the diameter. To prove this, let p(_a_, _b_) and d(_a_, _b_) be the path and distance from _a_ to _b_, respectively. If the diameter of a tree was p(_m_, _n_), and there was another node _x_ such that d(_u_, _x_) > _max_(_d_(_u_, _m_), _d_(_u_, _n_)), then that means either _d_(_m_, _x_) or _d_(_n_, _x_) is larger than _d_(_m_, _n_), which is a contradiction. Thus the depth of a node _u_ is _max_(_d_(_u_, _m_), _d_(_u_, _n_)). How can we find that value? We could very easily simply do **BFS** three times: first to find one endpoint of the diameter, then to find the degree from that node as well as the other endpoint of the diameter, finally to find the degree from the other endpoint. With multiple queries, we just need to store these degrees for an O(1) query.

Another (and worse) way is to find the distance from the center of the tree to each point, to keep track of parents, and of the endpoints of the diameter.

### Lowest Common Ancestor:
#### \*LCA*
Defined as the common ancestor of two nodes with the furthest degree of separation from the root node. There are two ways to do it.

The **Euler tour** of the tree can find the \*LCA* in O(1) time. It compresses the entire tree into one sequence, and keeps track of the first time and last time you encounter each node during a **DFS**. 

**Binary Lifting** can also find the \*LCA*, but in O() time.

## General Graphs
### Transversal:
There is **BFS** and **DFS**, implemented with a queue and stack respectively. For **DFS**, it’s easier to use recursion and just keep track of the states of all the nodes.

### Queries:
#### \*SHORTEST-PATH*
To get the shortest path in an unweighted graph, we can simply use **BFS**.
Otherwise, we can use the **Dijkstra** algorithm. It is a greedy algorithm similar to **BFS**, but we keep track of the “temporary distance” of each node from the source. At the beginning, we set all of these distances to be infinity. Then, for each iteration, we select the next unvisited node with the least temporary distance, mark it as visited, finalize its distance, and then go through all of its neighbours and see if we can get a shorter temporary distance. That’s why rather than using a normal queue we use a priority queue and sort by weights (or we just iterate if we have a complete graph). This ensures we make the optimal decision at any point in time. We end up potentially adding every edge to the priority queue, getting O(E log E) time for that part. For each next node, we also must remove from the priority queue, getting O(V log E). In total, we have O((E + V) log E) time. Since there are at least most V^2 edges, we can simplify it to O(E log V) time. Note that when E becomes very large (e.g, V^2) it essentially becomes quadratic, and it would be better not to use a priority queue and instead just to iterate for the naive O(V^2) implementation. 

**Dijkstra** fails if we have negative weights, so in that case we use the **Bellman-Ford** algorithm.

If we need to find multiple shortest path queries (shortest path between multiple points), it’s more efficient to use **Floyd-Warshall**. 

#### \*CONNECTED-COMPONENTS*
A connected component of a graph is a set of vertices such that a path exists from each vertex to each other one. A very intuitive way of determining whether two vertices lie on the same connected component is to use **BFS**. If we have multiple connected component queries and we also have a completed graph, we can simply iterate through every node and do a **BFS** if it hasn’t been visited, for O(n) time. An online way of determining connected components is to use a **disjoint set**.

#### \*MINIMUM-SPANNING-TREE*
A minimum spanning tree of a connected undirected graph is a subset of its edges such that  
1. The set of edges and nodes form a tree
2. The total sum of the weights of this subset of edges is minimized
There also exists a minimum spanning forest, for a non-connected undirected graph. There are two algorithms to find a MST, both of which are greedy.

The first is **Prim’s algorithm**, where we start at one node and eventually add other nodes to our tree until we have a MST. As it turns out, it works very similarly to **Dijkstra**. We start with any arbitrary node, and then put all outgoing edges into a priority queue, sorted by weights. Then, at each iteration, we pick the next unvisited node with the least edge weight, mark it as visited (and keep track of the edge too), then add all of it’s neighbours that aren’t visited into the priority queue. That gives us exactly the same complexity as **Dijkstra**: O(E log V) time complexity and O(E) space complexity.  
The second algorithm is **Kruskal’s algorithm**, where we add edges to our minimum spanning forest until it becomes a tree (or not, if the graph is not connected). We must also keep track of \*CONNECTED-COMPONENTS* so as to avoid cycles. We first sort all the edges by their weight. At each iteration, we take the edge with the least weight, given it doesn’t form a cycle. We then merge the connected components of the two endpoints, and repeat. Using a **disjoint set** to keep track of connected components, we have basically constant time queries and merges. In total we get O(E log E + V) = O(E log V) time, and O(E) space, the same as Prim’s.

### Cycles:
We use **DFS** to find cycles as well as their lengths in both directed and undirected graphs.

#### \*CYCLE-DIR*
For a directed graph, we keep a stack of all the vertices we have been to. When we backtrack from a vertex, we pop the stack. Then, there is a cycle if we have access to a vertex that is already in the stack. To find the length of this cycle, we can count the degree of separation between the first and the second occurrences of the vertex in the stack. To differentiate between the vertices we’ve already finished traversing and the ones in the stack, we keep track of each node’s state. We can save ourselves the implementation of the stack by using a recursive **DFS**. We could also use a topological sort (\*TOPOLOGICAL-SORT*) instead, although this doesn’t find cycle length.

#### \*CYCLE-UNDIR*
For an undirected graph, we can easily detect a cycle by running **DFS** and seeing if we ever visit an already visited node. If we wanted to count the length of this cycle, we’d keep a stack just as for a directed graph.

What if we wanted to find all cycles in a graph? If we had an edge list, we could use a **disjoint set** to keep track of all connected components as we iterate through all edges. Whenever we bump into an edge that forms a cycle, we increment a counter. Disclaimer: I do not know if this works for sure, I haven’t proven it.

#### \*CYCLES-DIR*
For a directed graph, we perform a similar operation to finding one cycle. For nodes numbered 1 to _n_, we start a **DFS** starting at every node 1 ≤ _u_ ≤ _n_, and find all cycles that end again on _u_. The reason the cycle must start at _u_ and end at _u_ is because if there was a graph _a_ → _b_ → _c_ → _d_ → _b_, then the cycle would be counted both when we start at node _a_ and when we start on any of the nodes in the cycle. Also, we must ensure we do not count the cycle more than once when we start at _b_, _c_, or _d_. One way to do this is to keep track of the start node _u_, and ensure during the **DFS** that any _v_ we travel to satisfies _v_ > _u_.

### Sorts:
#### \*TOPOLOGICAL-SORT*
A topological sort of a directed graph is a sorting of its vertices such that if there is a path from _u_ to _v_, then _u_ precedes _v_ in the sorting. Thus it’s obvious that such a sorting only exists if the graph is a DAG (directed acyclic graph). There’s two ways to find a topological sort, both in O(n) time.

The first is **Kahn’s algorithm** where we keep track of the indegree of each vertex (we can do this with an array), as well as which vertices have an indegree of 0 (we do this with a queue). We continue to take from the front of the queue and decrement all the nodes it leads into. If at any point the queue is empty but we have not included all of the vertices in the graph, there is a cycle.

The second is through a **DFS** and a stack, and keeping track of all vertices’ states, just like with a normal **DFS**. We pick a random vertex, and begin a **DFS** from that vertex. Whenever from a certain vertex we have exhausted all options, we push that vertex to the stack (e.g., if we reach a vertex with outdegree 0, or if we finish the search through all of a vertex’s child vertices). This ensures that all vertices are pushed to the stack before its parent. If we ever bump into a node that is currently in our path, we know there is a cycle. If we bump into a node that has already been finished, then we can skip over that node.




<a name="dynamicprogramming"/>

# Dynamic Programming
### Overview:
Can be used to solve problems with optimal substructure. If the optimal solution to a problem can be determined by the optimal solution to a subproblem, then a greedy algorithm should work. However, if it can only be determined by the optimal solution to multiple subproblems, then dynamic programming should be used.

In theory, it recursively breaks a problem _P_(_n_) into multiple subproblems
_P1_(_n_-1)… _Pk_(_n_-1), and stores the answers of subproblems _Pi_(_j_) such that _j_ < _n_ in a data structure, to ensure that the same problems are not calculated multiple times.


