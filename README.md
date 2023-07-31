# Solution to TSP
Many distinct solution to TSP.
## Method
### brute force
Enumerate all the possible path.

### greedy
Find the closet city. 

### dynamic programming
Memorize the state. (I use bits to compress the state)

### branch and bound (TODO)

### simulated annealing
ref: http://www2.stat.duke.edu/~scs/Courses/Stat376/Papers/TemperAnneal/KirkpatrickAnnealScience1983.pdf

#### simulated annealing with 2-opt
Based on SA, a neighbour of a state is produced by reverse a random path within.

### ant colony optimization
ref: M. Dorigo, M. Birattari and T. Stutzle, "Ant colony optimization," in IEEE Computational Intelligence Magazine, vol. 1, no. 4, pp. 28-39, Nov. 2006, doi: 10.1109/MCI.2006.329691.

#### ant colony optimization with 2-opt
After an ant find a path, try to reverse random path, recording the shorest path.

### ant colony system
ref: M. Dorigo and L. M. Gambardella, "Ant colony system: a cooperative learning approach to the traveling salesman problem," in IEEE Transactions on Evolutionary Computation, vol. 1, no. 1, pp. 53-66, April 1997, doi: 10.1109/4235.585892.

## HOW TO USE
normal mode:  
```
$ cd XXX_TSP
$ make dep all clean
$ cat testdata/XXX | ./main
$ make plot
```

debug mode:  
```
$ cd XXX_TSP
$ make debug all clean
$ cat testdata/XXX | ./main
```

record mode (record how the method find the answer, only apply on SA/ACO/ACS):
```
$ cd XXX_TSP
$ make record all clean
$ cat testdata/XXX | ./main
$ python3 visualize.py
```
