# Algorithms and Data Structures

Personal reference for data structures and algorithms for competitive programming. This would be awful to learn from, and it’s meant to be a refresher only. It’s also perpetually work-in-progress.

The names of algorithms or important data structures will be bolded. The names of problems or techniques (not actual names in any scientific literature, just for my personal ease of access) will be in unbolded caps, start and end with *, and be provided at the beginning of an algorithm description.

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

Note: some nice primes for hashing: 31, 1e9+7, 1e9+9, 1<<31 - 1, 1e11+3, 13e15+19

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

#### \*STD::SET/BST*
There’s a lot of really useful uses of the std::set, many of which I’m not very adept at doing. Ultimately, a std::set can do exactly what it’s meant to do: keep track of which items belong in it or not. This can extend to a few uses, from keeping track of the availability of things to the distance of nodes in a graph. There’s four things to keep in mind:
1. std::set::insert() and std::set::erase() can add, remove, and also sort of edit (be removing, updating, then inserting again) items in O(log n) time
2. std::set::lower_bound() can do a search within the set in O(log n) time, so it’s extremely simple to find elements that are greater/less than/equal compared to a certain element
3. std::set::begin() returns an iterator to the least element, std::set::rbegin() returns an iterator to the greatest element, so getting the minimum/maximum of a set is simple
4. std::set&lt;object&gt; can store any object that has the &lt; operator defined (e.g. std::pairs), which can come in handy for storing and updating distances or any other value associated with a node

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

### Searches:
#### \*SEARCH-VALUE*
When searching for a single value, the worst case is to look through each value, giving O(n) time. For an unsorted array, there isn’t any structure, so we can’t do it efficiently. For a sorted array, we can do **binary search** in O(log n) time.

#### \*SEARCH-TUPLE*
When searching for a pair that either adds or subtracts to a value, the worst case is O(n^2) time, pretty intuitively. However, in a sorted array, we can use the **two pointers** approach to do it in O(n) time. This works since we always move the pointers towards each other. In fact, for subsets of length k (i.e, k random values of the array that add or subtracts to a value), we can do it in O(n^(k-1)) time, since we can go through each element, then take the remaining elements k-1 at a time.

There is another approach to this problem, and it can be a bit more generalized. For searching a pair that sums up to s, you can iterate through every element i, then do a **binary search** for s-a[i]. Thus you can find a pair in O(n log n) time. In fact, for subsets of length k, we can do it in O(n^(k-1) * log n) time, which is always log n times slower than the multiple pointers approach, but the idea can be applied to relatively more problems.

#### \*MAXIMUM-SUBARRAY-SUM*
To find the max subarray sum, we can use **DP** with **Kadane’s algorithm**. Note that the answer is trivial if all elements are positive. For each element at index _i_ of the array, we store the maximum subarray that ends at that element, denoted as _m(i)_. We use the substructure _m(i) = max(m(i-1)+a[i], a[i])_.

#### \*INVERSION*
An inversion of array A is two indices i and j such that i < j and A[i] > A[j]. Inversions can be counted through brute force in O(n^2) complexity, by iterating through each element, and checking for the number of elements to the left of it (index less than it) and a value greater than it. We can do this faster with a **BIT**.

### Range Queries:
Obviously the worst case for doing a range query is just to go through each of the values (linear in both array length and query number). There’s often better ways of doing so, either with lazy evaluation or preprocessing.

#### \*RSQ*
For range sum queries, and some others, we can do a **prefix sum array (PSA)**, giving us O(1) query complexity and O(n) update complexity. We can apply this to a 2D array too, making a **2D prefix sum array** that allows calculation of a rectangular sum query using the principle of inclusion exclusion. It gets the sum of the rectangle defined by [0][0] as the top left corner, and [x][y] as the bottom right, which is the 2D equivalent of a prefix sum. Then we can just manipulate the four corners of the rectangle we want to sum. It calculates the query in O(1) complexity, and initializes/updates in O(n^2) complexity.

Another way for range sum queries, and some others, is to use a tree, such as a **binary indexed tree/Fenwick tree (BIT)**. The hypothetical tree stores certain precomputed range queries that can be added together to fulfill the complete query. It’s stored as an array, so it takes up the same space as a prefix sum array, and allows O(log n) query and update complexity, and O(n log n) initialization, making it really good. The reason it works out is because it uses binary representations of numbers in a smart way. We can apply this to a **2D BIT** as well.

#### \*RMQ*
For range minimum queries, and some others, one way is to use a **sparse table**. This table precomputes a property for a subarray starting at every possible position, with length of all powers of 2. This way, it can either break every range query down into powers of 2 (storing it as binary) for an O(log n) query, or if the subarrays are able to overlap as in the case of a \*RMQ*, an O(1) query. It takes up O(n log n) space. The table can be computed in O(n log n) time by using a kind of dynamic programming, using two lengths of 2^n to make a length of 2^(n+1).

The sparse table can also compute the **lowest common ancestor** of trees.

#### \*EXCEPT-SUBARRAY-QUERY*
Let’s say we wanted to find a property of all elements in an array except in a certain subarray, let’s say a[i, j]. It’s kind of like the opposite of a range query. We can create both a prefix and suffix array for that property, and find the property of the subarray a[0, i] with the prefix array, and the property of the subarray a[i, end] with the suffix array, and combine the two ranges.

### Subsequences and Subarrays:
Many subsequence problems can be solved through DP. For that reason, I put them under the DP section. A subsequence _A_ of _B_ is defined as a sequence achieved by deleting 0 or more elements from _B_. Note they don’t have to be contiguous. A subarray is defined as a contiguous subsequence. For many problems on subarrays, especially when involving their cumulative sum, it’s helpful to turn them into a **PSA** and then use one of the searching algorithms.

#### \*MAX-SUBARRAY-VALUE* or \*DISTRIBUTE-VALUE*
Let’s say we had an array to be split into k contiguous subarrays, and we’d like to maximize the sum of any subarray. The naive solution would be to iterate through and find each of the k subarrays, which even if we used a **PSA** to find the sum, it would still take O(n^(k-1)) time. Rather, we could do a binary search of the maximum value. Let’s say the total sum of the array was z. We would binary search all the values between 1 and z, and set it as the maximum. Then, we can iterate through the **PSA** and mark off whenever we surpass that value, and see if we have more or less groups than the target value k. Thus, we would only need O(n log z) time, which can be much less depending on the value of z.

#### \*PAINTERS-PARTITION*
The most generalized form of the painter’s partition problem is to split an array into k contiguous subarrays such that
1. The difference between the largest and smallest sums of all subarrays is minimized,
or
2. The largest sum of all subarrays is minimized.

Both of these should be equivalent, but lead into different methods. TODO: PROVE THIS

For 1., after converting the array into a **PSA**, you can simply do a repeated \*SEARCH-TUPLE* problem, and find a tuple of indices using **binary search** (two pointers doesn’t work here, since we need to find multiple tuples in the array at specific positions). TODO: EXPAND ON THIS ALGORITHM

For 2., after converting the array into a **PSA**, you can do the \*DISTRIBUTE-VALUE* problem, also using binary search. This works in a pseudo-polynomial O(n log s) time, where s is the sum of all elements in the array.

## Transformations
Different sorts of transformations can be applied to arrays and subarrays. The minimum complexity for any transformation is O(n), and O(n log n) for sorting, so we generally try to achieve this with some smart tricks. Lots of arrays can be sorted to avoid going through all subsets, and instead using structure to be more efficient. 

### Sorting:
#### \*SORT*
Involves two major components: comparisons and re-arrangement of the array. Minimizing the running time of the entire sort involves minimizing running time of both of them. In practice, just use the built-in sort function. Some useful types of sorting are below.

**Merge sort** is a recursive sorting algorithm using a divide and conquer approach: to sort each interval, we split it in half, then apply **Merge sort** to both intervals, then combine the two together. The base case is one or two items, which is either trivial or at most a single swap. We can merge two sorted intervals in O(n) time, and we go over the entire interval in an amortized way for each division of all intervals, which occurs O(log n) times, giving us that O(n log n) time.

**Quicksort** is also a divide and conquer algorithm. It’s intuitively recursive, but can be made to work iteratively. It picks an element to be the “pivot”, and moves all elements that are less than it to its left, and all elements greater than it to its right. This step is referred to as partitioning, and will take O(n) time. We now know that the pivot is in its correct position, and we thus apply **Quicksort** to the left and right side. Again, all of the partitioning will take an amortized O(n) time for each iteration. So, how many iterations are there? Well, if we were able to split the array in half each time, we’d get roughly O(log n) iterations, for a total of O(n log n) time. If, however, we picked sub-optimal values (like always picking the least or greatest element in the interval) we’d end up with an O(n^2) solution. That’s why picking the pivot is quite important in quicksorts.

**Bubble sort** is the most well known bad sorting algorithm. It sweeps from left to right, and if two elements are out of order, it swaps them. Each O(n) iteration is only guaranteed to move the greatest element to its place, meaning it comes out to O(n^2).

**Patience sort** is an interesting one, similar to a card game called Patience. It separates into two parts: firstly, it keeps a record of multiple decreasing piles, and then it merges all of them. To sort into piles, every element can only be put on a pile where the last element is greater than it, and if there are multiple possibilities of piles, it gets put on the left-most one. If it cannot be put in any pile, it should form a new pile at the very right. This rule ensures that we create the minimum number of piles, since the top-most value decreases as we go to the left. To find the specific pile a card goes on, we do a binary search of all piles. At the end of this procedure, we have _k_ piles (_1 ≤ k ≤ n_) that are sorted in decreasing order, in O(n log n) time. We now do a \*K-WAY-MERGE* in O(n log n) time to complete the algorithm. TODO: TALK ABOUT K-WAY MERGES.

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

### Sets:
Many problems also end up involving sequences that can be efficiently represented as sets. For a rundown of sets, see the [Math](#math) section.

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
Most string problems can be reduced to other types of problems (often array problems). Some of them though, such as pattern matching, can’t.

### Pattern Matching
#### \*MATCH-STRING*
How would we search for all occurrences of string, let’s say of length _p_, in some other string of length _q_? The naive approach would be to simply loop through for a O(pq) search. However, we can do better with two algorithms, both in O(q) time.

The first quicker way is to do a **Rolling hash**. For this algorithm, we usually do a polynomial hash of a substring, with a certain power _p_ and modulus _m_, both prime numbers. Then, to “roll” the hash we simply subtract the highest degree term, multiply the result by _p_, and finally add the constant term. Then, if this hash matches the hash of our initial string, we have a match!

The second way is to use the **Knuth-Morris-Pratt** or **KMP** algorithm. TODO

#### \*MATCH-PERMUTED-STRINGS*
What if we wanted to match all permutations of a certain string within a larger string, and keep track of the number of distinct permutations we encounter? TODO

<a name="graphs"/>

# Graphs
### Overview:
Borrows from graph theory in math. A special type of connected graph with no cycles is called a “tree”. Multiple trees makes a “forest”.

### Transversal:
There is **BFS** and **DFS**, implemented with a queue and stack respectively. They work for both directed and undirected graphs. **DFS** for transversal is a bit complicated, as it involves keeping track of the state of each vertex to avoid cycles. It’s easier to use recursion (usually).

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

### Center, Diameter, Radius, Depth:
These ideas are probably the most important in terms of solving tree-specific problems.

Finding the diameter is pretty simple. Start at a random node _s_ and find the farthest node from it, called _u_. This can be done with **BFS** or **DFS**. Then root the tree at _u_, and find the farthest node _v_, repeating the same process. _u_ and _v_ are then endpoints of a diameter. To find the diameter itself, during the second transversal starting at _u_, keep track of each node’s parent. Then, follow the path backwards from _v_ to _u_.

Once the diameter path is obtained, the center, radius, and depth are trivial. Some problems that require knowing the centre, diameter, and radius are below:

#### \*DEPTH-QUERY*
To get the depth of a tree by arbitrarily rooting it any node, let’s say _u_, we need to find the maximum distance from _u_ to any other part of the tree. This distance will always be part of the diameter. To prove this, let p(_a_, _b_) and d(_a_, _b_) be the path and distance from _a_ to _b_, respectively. If the diameter of a tree was p(_m_, _n_), and there was another node _x_ such that d(_u_, _x_) > _max_(_d_(_u_, _m_), _d_(_u_, _n_)), then that means either _d_(_m_, _x_) or _d_(_n_, _x_) is larger than _d_(_m_, _n_), which is a contradiction. Thus the depth of a node _u_ is _max_(_d_(_u_, _m_), _d_(_u_, _n_)). How can we find that value? We could very easily simply do **BFS** three times: first to find one endpoint of the diameter, then to find the degree from that node as well as the other endpoint of the diameter, finally to find the degree from the other endpoint. With multiple queries, we just need to store these degrees for an O(1) query.

Another (and worse) way is to find the distance from the center of the tree to each point, to keep track of parents, and of the endpoints of the diameter.

### Lowest Common Ancestor:
#### \*LCA*
Defined as the common ancestor of two nodes with the furthest degree of separation from the root node. There are two ways to do it.

The **Euler tour** of the tree can find the \*LCA* in O(1) time. It compresses the entire tree into one sequence, and keeps track of the first time and last time you encounter each node during a **DFS**. 

**Binary Lifting** can also find the \*LCA*, but in O() time.

## General Graphs

### Paths:
#### \*SHORTEST-PATH*
The shortest path problem is probably the most famous graph problem. Each edge is given a weight, and the weight of a path is defined as the sum of the weights of each composing edge. In an unweighted graph, all edges can be seen as having a weight of 1 unit. The aim is to find the path with the minimum weight. All the algorithms below work for both directed and undirected graphs.

##### SINGLE SOURCE
To get the shortest path in an unweighted graph, we can simply use **BFS**.

Otherwise, we can use the **Dijkstra** algorithm, which works if all weights are positive. It is a greedy algorithm similar to **BFS**, but we keep track of the temporary distance of each node from the source. At the beginning, we set all of these distances to be infinity. Then, for each iteration, we take the node with the minimum temporary distance and finalize its distance. We know it’s the optimal solution since if the path contained any other node with a greater or equal temporary distance, then we’d either have the same or greater distance. After that, we perform “relaxation”, which is just updating the minimum path length to all the other nodes. The actual implementation of that is by visiting that node (and marking it as so), and checking all of its neighbours to see if it can reduce the distance. Then, if it can, we need to add it to some sort of queue that’s able to get the minimum temporary distance. That’s why rather than using a normal queue we use a priority queue and sort by weights (or we just iterate if we have a complete graph). This ensures we make the optimal decision at any point in time. Note that this property is no longer true if there are negative weights. We end up potentially adding every edge to the priority queue, getting O(E log E) time for that part. For each next node, we also must remove from the priority queue, getting O(V log E). In total, we have O((E + V) log E) time. Since there are at least most V^2 edges, we can simplify it to O(E log V) time. Note that when E becomes very large (e.g, V^2) it essentially becomes quadratic, and it would be better not to use a priority queue and instead just to iterate for the naive O(V^2) implementation. 

**Dijkstra** fails if we have negative weights, so in that case we use the **Bellman-Ford** algorithm. It starts by setting all distances to infinity, and the source distance to 0. Then, for relaxation, instead of finding the node with the least temporary distance, it iterates through all edges and updates the distance of a node if there’s a better path. This is guaranteed to relax at least one node. TODO: PROVE THIS. That’s why we perform this step V-1 times, for a total of O(EV) runtime.

There’s another better algorithm that can handle negative weights, called the **Shortest Path Faster Algorithm**, or the **SPFA**. Similarly to the other two, it sets all node distances to infinity, and the source to 0. It also uses a queue (normal queue, not priority queue) to keep track of next nodes. At each iteration, it pops a node from the queue, and then relaxes all of its neighbours. Unlike **Dijkstra** however, whenever it updates a neighbouring node, it adds it to the queue (possibly more than once), since there’s a chance it can update other paths from that node. At each node, the algorithm can potentially push all the other edges into the queue, which means it runs in O(EV) time. However, experimentally, it runs in O(E) time, a better complexity than all of the other algorithms. This allows it to always outperform **Bellman-Ford**, and usually **Dijkstra**.

##### ALL PAIRS
If we need to find multiple shortest path queries (shortest path between multiple points), it’s usually more efficient to use **Floyd-Warshall**. The Floyd-Warshall algorithm is a dynamic programming algorithm, with O(V^2) space complexity and O(V^3) time complexity. It works for negative weights (no negative cycles though) as well. For denser graphs, this beats out running **Dijkstra** V times resulting in O(VE log V) → O(V^3 log V). The DP transition used is that  
_distance(i, j, k) = min(distance(i, j, k-1), distance(i, k-1, k-1)+distance(k-1, j, k-1))_  
where _distance(i, j, k)_ represents the minimum distance from _i_ to _j_ where the path only consists of those two nodes and optionally the nodes from 1 to _k_. It’s quite easy to see why it runs in cubic time, and by iterating through _k_ (similar to the space optimized knapsack) we only need quadratic space.

#### \*BOTTLENECK-PATH*
What’s referred to as the widest path problem defines the bottleneck of a path to be the edge with the least weight that makes up this path. The problem seeks to find the shortest path that maximizes the value of the bottleneck (in other words, the shortest path that maximizes the minimum value of any of the edges it’s made up of), akin to the maximin problem from other fields of math. It’s also often useful to do the opposite: to find the path that minimizes the maximum edge weight in a graph, akin to the minimax problem.

##### SINGLE SOURCE
For the widest path problem, we can slightly modify **Dijkstra**. Instead of keeping track of distance, we’re keeping track of the temporary maximum bottleneck from the source to another node. At each iteration, we take the node _n_ with the highest temporary bottleneck, finalize its bottleneck, and perform relaxation of its neighbours. The relaxation step updates the temporary bottleneck of a _v_ only if _min(w(n, v), b(n)) > b(v)_ where _b(v)_ is the temporary bottleneck of node _v_ (or finalized, in the case of _n_). The time complexity ends up being the same as that of **Dijkstra**, and it similarly works on undirected and directed graphs.

The process is almost exactly the same for the minimax equivalent. However, for each iteration we instead take the node _n_ with the minimum max-edge, and for the relaxation of a neighbour _v_, we update _m(v)_ only if _max(w(n, v), m(n)) < m(v)_ where _m(v)_ is the minimum max-edge of node _v_.

##### ALL PAIRS
One way to find the path with the maximum bottleneck is to create a _maximum_ spanning tree, which is almost exactly the same as the \*MINIMUM-SPANNING-TREE*  problem. This only works on undirected graphs. Once we have our MST, it’s quite simple to find the single path between each two nodes using a transversal in O(n) time.

The minimax equivalent of this is quite literally the \*MINIMUM-SPANNING-TREE* problem.

#### \*KTH-SHORTEST-PATH*
This problem is similar to the \*SHORTEST-PATH* problem, but instead of finding the shortest, we find the k-th shortest (either strictly, i.e. the k-th minimum out of all possible path _lengths_, or non-strictly, i.e. the k-th minimum out all possible _paths_). It can be solved using multiple algorithms, but a simple one (with I think not too great time complexity) is simply a modified **Dijkstra’s**. Rather than storing the shortest distance to get to any node, we store the _k_ shortest distances instead (this can be done as a set). Then, the _i_th time we reach node _p_, we store it as the _i_th greatest distance to get to _p_ if it’s greater than the _i-1_th shortest distance. (note that, since we order nodes in a priority queue, this is guaranteed to be the optimal distance). If we reach a node more than _k_ times we simply discard that value. We thus continue with this algorithm until we either run out of paths, or we reach our destination node more than _k_ times. This takes what I believe to be O(k E log V) time.

##### SECOND-SHORTEST-PATH
There’s a special case to consider for the second shortest path if we know the start node _s_ and the end node _k_. This technique works for the strictly 2nd shortest path. We can compute it by solving \*SHORTEST-PATH* twice: once from the beginning and once from the end (so we have a prefix shortest distance and a suffix shortest distance). This is a very useful technique that’s used often for all sorts of graph problems, not just the second shortest path. Note that if the graph is directed then we need to reverse it before calculating the suffix shortest path. Then, for each edge _<u, v>_, we find _dst(s, u) + <u, v> + dst(v, k)_. If this value is greater than the shortest path, we consider it as a potential 2nd shortest path. TODO: PROVE THIS. If we’re allowed to visit each edge more than once, we can simply perform the same thing with the edge _<v, u>_.

#### \*SHORTEST-PATH-WITH-SECONDARY-RESTRICTION*
We define each edge to have a weight _w_ and a secondary weight _k_.  The shortest path with secondary restriction should find the path with the minimum weight (sum of _w_) that satisfies the fact that the sum of all _k_ should be less than a certain _K_. To solve this, we modify **Dijkstra**. We keep track of the shortest path to each node that has a _k_ sum of exactly _ki_, for 1 ≤ _i_ ≤ _K_. This gives us O(VK) storage complexity, with the same time complexity (but worse runtime than normal **Dijkstra** in practice, given that every edge will likely be used).

Sometimes the problem is defined in another way: we define each edge to have weight _w_ and a secondary weight _k_. The shortest path with the secondary restriction should find the path with the minimum sum of _w_, and then break ties with the second sum of _k_. Here, we modify **Dijkstra** again, where we compare paths by first comparing their _w_ sum, and then their _k_ sum. This has the same time and space complexity as normal.

### Components:
#### \*CONNECTED-COMPONENTS*
A connected component of a graph is a set of vertices such that a path exists from each vertex to each other one. A very intuitive way of determining whether two vertices lie on the same connected component is to use **BFS** or **DFS**. If we have multiple connected component queries and we also have a completed graph, we can simply iterate through every node and do a **BFS** if it hasn’t been visited, for O(n) time. An online way of determining connected components is to use a **disjoint set**.

#### \*ARTICULATION-POINTS* or \*BRIDGE-FINDING*


#### \*MINIMUM-SPANNING-TREE*
A minimum spanning tree (MST) of a connected undirected graph is a subset of its edges such that  
1. The set of edges and nodes form a tree (forms a spanning tree, or ST)
2. The total sum of the weights of this subset of edges is minimized
There also exists a minimum spanning forest, for a non-connected undirected graph. There are two algorithms to find a MST, both of which are greedy.

The first is **Prim’s algorithm**, where we start at one node and eventually add other nodes to our tree until we have a MST. As it turns out, it works very similarly to **Dijkstra**. What’s different, though, is that it works even with negative edge weights (even negative cycles). That’s because it doesn’t accumulate distance, it rather just looks at individual edges. We start with any arbitrary node, and then put all outgoing edges into a priority queue, sorted by weights. Then, at each iteration, we pick the next unvisited node with the least edge weight, mark it as visited (and keep track of the edge too), then add all of it’s neighbours that aren’t visited into the priority queue. That gives us exactly the same complexity as **Dijkstra**: O(E log V) time complexity and O(E) space complexity.

The second algorithm is **Kruskal’s algorithm**, where we add edges to our minimum spanning forest until it becomes a tree (or not, if the graph is not connected). Like **Prim’s**, it’s able to work even if there are negative edges. We must keep track of \*CONNECTED-COMPONENTS* so as to avoid cycles in the ST. We first sort all the edges by their weight. At each iteration, we take the edge with the least weight, given it doesn’t form a cycle. We then merge the connected components of the two endpoints, and repeat. Using a **disjoint set** to keep track of connected components, we have basically constant time queries and merges. In total we get O(E log E + V) = O(E log V) time, and O(E) space, the same as Prim’s.

### Cycles:
We use **DFS** to find cycles as well as their lengths in both directed and undirected graphs.

#### \*CYCLES-DIR*
How do we detect a cycle in a directed graph? We keep a stack of all the vertices we have been to. When we backtrack from a vertex, we pop the stack. Then, there is a cycle if we have access to a vertex that is already in the stack. To find the length of this cycle, we can count the degree of separation between the first and the second occurrences of the vertex in the stack. To differentiate between the vertices we’ve already finished traversing and the ones in the stack, we keep track of each node’s state. We can save ourselves the implementation of the stack by using a recursive **DFS**. We could also use a topological sort (\*TOPOLOGICAL-SORT*) instead, although this doesn’t find cycle length.

To find all cycles in an undirected graph, we perform a similar operation to finding one cycle. For nodes numbered 1 to _n_, we start a **DFS** starting at every node 1 ≤ _u_ ≤ _n_, and find all cycles that end again on _u_. The reason the cycle must start at _u_ and end at _u_ is because if there was a graph _a_ → _b_ → _c_ → _d_ → _b_, then the cycle would be counted both when we start at node _a_ and when we start on any of the nodes in the cycle. Also, we must ensure we do not count the cycle more than once when we start at _b_, _c_, or _d_. One way to do this is to keep track of the start node _u_, and ensure during the **DFS** that any _v_ we travel to satisfies _v_ > _u_.

#### \*CYCLES-UNDIR*
How do we detect a cycle in an undirected graph? We can easily detect a cycle by running **DFS** and seeing if we ever visit an already visited node. If we wanted to count the length of this cycle, we’d keep a stack just as for a directed graph.

What if we wanted to find all cycles in an undirected graph? If we had an edge list, we could use a **disjoint set** to keep track of all connected components as we iterate through all edges. Whenever we bump into an edge that forms a cycle, we increment a counter. TODO: PROVE THIS. I DO NOT KNOW IF THIS WORKS.

### Sorts:
#### \*TOPOLOGICAL-SORT*
A topological sort of a directed graph is a sorting of its vertices such that if there is a path from _u_ to _v_, then _u_ precedes _v_ in the sorting. Thus it’s obvious that such a sorting only exists if the graph is a DAG (directed acyclic graph). There’s two ways to find a topological sort, both in O(n) time.

The first is **Kahn’s algorithm** where we keep track of the indegree of each vertex (we can do this with an array), as well as which vertices have an indegree of 0 (we do this with a queue). We continue to take from the front of the queue and decrement all the nodes it leads into. If at any point the queue is empty but we have not included all of the vertices in the graph, there is a cycle.

The second is through a **DFS** and a stack, and keeping track of all vertices’ states, just like with a normal **DFS**. We pick a random vertex, and begin a **DFS** from that vertex. Whenever from a certain vertex we have exhausted all options, we push that vertex to the stack (e.g., if we reach a vertex with outdegree 0, or if we finish the search through all of a vertex’s child vertices). This ensures that all vertices are pushed to the stack before its parent. If we ever bump into a node that is currently in our path, we know there is a cycle. If we bump into a node that has already been finished, then we can skip over that node.




<a name="dynamicprogramming"/>

# Dynamic Programming
### Overview:
Can be used to solve optimization problems or counting problems that have optimal substructure. If the optimal solution to a problem can be determined by the optimal solution to a subproblem, then a greedy algorithm should work. However, if it can only be determined by the optimal solution to multiple subproblems, then dynamic programming should be used.

In theory, it recursively breaks a problem, and stores the answers of subproblems in a data structure, to ensure that the same problems are not calculated multiple times. Each subproblem is referred to as a DP state, and is uniquely identified by some set of variables _(n1, n2, … nk)_. Usually the first step of solving a DP problem would be to identify the DP state. The equation to get from one state to another (i.e. solving a problem based on its subproblems) is referred to as a DP transition, and is the second step. The formula to solve each problem (aka our DP transition) usually takes the form of a function _P(n1, n2, … nk)_. Each subproblem would have its DP state simplified, usually taking on the form  _P(n1-d1, n2-d2, … nk-dk)_.

There are usually two approaches to achieving DP: top-down and bottom-up. Top-down approaches are recursive in nature, calculating a subproblem when it is needed and storing its value. It usually uses a process called **memoization**, basically remembering the answer to a subproblem when we need to calculate it. Bottom-up approaches are iterative in nature, calculating the answer to every query and storing it in a table. For this reason, this process is called **tabulation**.

Often, there’s a large memory overhead with DP. That’s because for a problem _P(n1, n2, … nk)_ we often need to store every possible value of each variable _ni_ associated with the DP state, reaching polynomial space to the power of _k_. There’s a really common memory optimization though. The way we tabulate or memoize often requires us to iterate through each value of each _ni_. However, if we won’t ever be looking at subproblems more than some constant _c_ away (that is, the subproblem _P(n1-d1, n2-d2, … nk-dk)_ has _di ≤ c_), we can get rid of one factor. Instead of storing all _n1 * n2 * … nk_ values, we store _c * n2 * … nk_ values. That’s because for the calculation of each subproblem, we only ever need _c_ values of _n1_ anyway. Notice we can’t get rid of more than 1 factor, since as we iterate using for loops, once each inner for loop finishes its iteration, it jumps back to the start value. The outer for loop is the only one that iterates “smoothly”, without any resets.

### Knapsacks:
Many problems can be reduced into some form of the knapsack problem. It involves a theoretical knapsack with a weight capacity of _W_. There are _n_ given objects/weights, each with weight _wi_ and value _vi_, for 1 ≤ _i_ ≤ _n_. The objective is to maximize the total sum of values of the objects we can put in the knapsack such that the total weight is within the knapsack’s capacity. It can be easily shown (by counterexample) that greedy solutions for multiple knapsack problems won’t work. The different variants of problems are below.

#### \*0-1-KNAPSACK*
The 0-1 Knapsack problem is pretty straightforward: for each weight, we can either take it or not. As with most knapsack problems, we consider each weight individually. Our DP transition is  
_maxV(w, i) = max(maxV(w, i-1), vi+maxV(w-wi, i-1))_  
where _maxV(w, i)_ is the maximum value attainable for a knapsack of weight _w_ only considering the first _i_ objects. What this does is that for each weight, we consider the case of whether we take it or not, and find the maximum of the two cases. How do we calculate this then? One way is to create a 2D array _maxV[n][W]_, and then for each weight we then iterate through and find the answer from the entries in the previous row. This stores every single optimal answer up to our actual weight capacity _W_, and works in O(_n*W_) complexity for both time and space.

Since we only use entries from the previous row, we can optimize it for space a little. In fact, we can only store a single array _maxV[W]_. Then, we simply update all the values starting from the end, allowing us to have O(_W_) space complexity but with the same time complexity.

If our weights are exceedingly large and our different values are by comparison smaller, we can modify the table a bit. Notice that we can store the minimum weight required to achieve each value, also by iterating through each element. Our DP transition becomes  
_minW(v, i) = min(minW(v, i-1), wi+minW(v-vi, i-1))  
where _minW(v, i)_ is the minimum weight required to achieve a value of exactly _v_ using the first _i_ objects. At the beginning, we would need to initialize all minimum weights to infinity, and set the minimum weight of 0 to always be 0. Then, we can find the answer in O(_n*V_) time (assuming _V_ is the sum of all the values of all objects), and either O(_n*V_) or O(_V_) space, depending on how we initialize the array.

#### \*EXACT-0-1-KNAPSACK*
What if the problem instead asked for the knapsack to be filled exactly to capacity, even if it meant decreasing the value? The modification is actually exceedingly simple: we simply initialize the array to store negative infinity as the max value for each weight, except for 0. Then, we simply perform the knapsack algorithm as usual. It’s not difficult to see why this works.

#### \*UNLIMITED-KNAPSACK*
A common variant is where we have potentially infinite of each object, with the rest of the conditions being identical. How do we solve the problem this time? Let’s write out a DP transition: again, we can either take each object or not, and we look at all the objects one at a time. However, when we take an object, we don’t necessarily decrease the number of objects we look at, since we could take the same object again. This leads to the following transition (with the same variables as above):  
_maxV(w, i) = max(maxV(w, i-1), vi+maxV(w-wi, i))_  
Note the only difference is we have _i_ rather than _i-1_ at the end. And, to solve the problem we use largely the same approach: we can store the max value for all weights up to _W_, but move from the beginning instead of from the end (since we need to access values of the current row, not of the previous row). This runs in the same time and space complexity: O(n*W) and O(W) respectively.

#### \*LIMITED-KNAPSACK*
Yet another common variation is where we have a finite amount _k_ of every type of object, with the same constraints. The naive approach would be to run a 0-1 knapsack, but instead of taking one object, we loop through the amount of objects there are. That is, the DP transition is  
_maxV(w, i) = max(maxV(w, i-1), vi*maxV(w-wi, i-1), 2*vi*maxV(w-2*wi, i-1)…ki*vi*maxV(w-k*wi, i-1))_  
This would run in O(n*k*W) time, and O(W) space. For larger values of _ni_, this could be quite bad.

One optimization is to use a technique called **binary packaging**. Like with many other data structures, we take advantage of the properties of binary numbers and easy bit operations. We first express _k_ as the sum of consecutive powers of 2, and some remainder. In other words, we find integers _p_ and _r_ such that 1+2+...+_2^p_+_r_ = _k_, and where _p_ is as large as possible. Then, for each integer _q_ in this set, we treat that amount as a single weight. In other words, we “add” an object with weight _q*wi_ and value _q*vi_ for _1 ≤ i ≤ n_. So, in total we will have _n*ceil(log2 k)_ objects. Now, we can just treat it like a 0-1 knapsack. This runs in O(n*log k*W) time and O(W) space.

### Sequences
DP can be used to solve multiple problems related to sequences. Many of these sequential DP problems are related to subsequences. A subsequence _A_ of _B_ is defined as a sequence achieved by deleting 0 or more elements from _B_. Note they don’t have to be contiguous. In a lot of the sections I use “sequence” and “array” synonymously.

#### \*LCS*
The longest common subsequence of two arrays _A_ and _B_ is defined as the sequence of longest possible length that is a subsequence of both _A_ and _B_. The naive solution would be to check every single possible subsequence for exponential time. And once again, we can use DP. Let _lcs(X, Y)_ be the LCS of X and Y. Let _x, y_ be the last elements of _X, Y_ respectively. If _x_ = _y_, then we can include them in the LCS. This is pretty trivial. If _x_ != _y_, then either _x_ cannot be in the LCS or _y_ cannot be in the LCS. This is also quite trivial. Thus our recursive formula looks like this:  
_lcs(X, Y) = {  
Ø if X = Ø or Y = Ø,  
lcs(X-x, Y-y) + x if x = y,  
max(lcs(X-x, Y), lcs(X, Y-y)) if x != y  
}_.  
How do we turn this into a DP transition? First, let the length of _A_ be  _a_ and the length of _B_ be _b_. Let _c(i, j)_ be the length of LCS of the first x elements of A and the first y elements of B. We can then modify our recursive formula above to get the following DP transition that works in O(n^2) time and space:  
_c(i, j) = {  
	0 if i = 0 or j = 0,  
	c(i-1, j-1) + 1 if A[i] = B[j],  
	max(c(i-1, j), c(i, j-1)) if A[i] != B[j]  
}_.  
Note by keeping only two rows we can optimize to O(n) space. We can get the LCS itself in O(n) time by stepping back through our table.

#### \*LCS-DISTINCT-ELEMENTS*
If at least one of the arrays has only unique elements, we can actually improve our time complexity to O(n log n) by reducing the problem. We take the array with unique elements and do a reverse mapping of each element to its index. We then iterate through every element in the other array, and if it is present in the map, add the index of itself in the first array to a dynamic array. The indices of the elements in the LCS is simply the \*LIS* of our dynamic array.

#### \*LIS* or \*LDS*
The longest increasing/decreasing subsequence of one array _A_ is defined as the subsequence of _A_ of longest possible length such that the elements are in strictly increasing/decreasing order. Like the \*LCS*, the naive solution is exponential. However, we can use DP for a more efficient solution. Let’s only look at the longest increasing subsequence, since the logic maps trivially for the LDS. Let _lis(x)_ be the LIS that ends at the _x_th index. Thus, we can consider the following DP transition:  
_lis(x) = max(lis(i)) + A[x]_, if _A[i] < A[x]_ where _1 ≤ i < x_.  
This gives us our result in O(n^2) time and we can get the list itself in O(n) time.

However, we can actually do better, using something close to **Patience sort**. We only do the first part without the merging, and end up with multiple piles. How do we use this to find the LIS? It turns out that the index of the pile each element is added to is also the length of the LIS that ends at that element. TODO: PROVE THIS. This allows us to get the LIS in O(n log n) time. If we only wanted a single LIS (rather than every single LIS that ends with each element), note we only need to keep one card from each pile: the minimum one. That means we can represent our piles as a single vector. Beyond that, we follow largely the same procedure: put the card in the leftmost pile if possible, otherwise create a new array.


#### \*EDIT-DISTANCE*
The edit distance of two sequences is defined as the minimum number of edits applied to one to get the other. Different edit distances are possible depending on the legal edits that can be applied to a sequence. The most common edit distance, the Levenshtein distance, allows for the following:
1. Deleting any element of a sequence
2. Inserting any element to any position in a sequence
3. Modifying any element of a sequence
Let _A_ be the first sequence, and _B_ be the second, such that _|A| ≤ |B|_. Since each operation can be reversed, we can look at the edit distance to transform _A_ to _B_ WLOG. Let’s look at the last element of both. Note that at some point during our transformation, we must do something to make the _A_’s last element equal to _B_’s. If they are equal, then there’s nothing needed. Apart from that, we can consider the case there’s an insertion, deletion, or a modification. This gives rise to the DP transition.

Let _lv(i, j)_ be the edit distance of the first _i_ elements of array _A_ and the first _j_ elements of array _B_. When we consider _lv(i, j)_, there’s again the same cases. If _A[i] = B[j]_, then _lv(i, j) = lv(i-1, j-1)_. Otherwise, we can delete the last character of _A_ and compare the edit distance of the two results, insert _B_’s last character to the end of _A_ (equivalent to deleting the last character of _B_) and compare the results, or modify the last character of _A_ to be the same as _B_’s. Finally, if any sequence is empty, then we have no choice but to add/delete all characters from the non-empty one. Our transition is thus  
_lv(i, j) = {  
j if i = 0,  
i if j = 0,  
min(  
	lv(i-1, j) + 1,  
	lv(i, j-1) + 1,  
	lv(i-1, j-1) + cost(i, j)  
) otherwise  
}_  
where _cost(i, j) = 0_ if _A[i] = B[j]_, and _cost(i, j) = 1_ otherwise. The three values correspond to deleting the last character of _A_, appending the last character of _B_, and modifying the last character of _A_ to _B_’s (if needed), then ignoring them. This works in O(n^2) time and space, but an O(n) space optimization should be possible.
