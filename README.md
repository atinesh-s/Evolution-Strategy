# Evolution Strategy
Evolution Strategy was developed by Rechenberg and Schwefel at the Technical University of Berlin in the mid-1960s. For the brief history of Evolution Strategy refer [1].

## Benchmark Functions
1. Ackley
2. Rastrigin
3. Rosenbrock
4. Schwefel
5. Sphere

## Requirements
* g++

## Commands
```
$ make
$ ./main function dimension runs
```

## Results
#### /plot 
* `<function>_<dimension>_<run>.txt`: Contains plotting related data
#### /results
* `best.txt`: The Global best value of all the runs
* `time.txt`: Total time required for all the runs
* `measures.txt`: Mean and Standard deviation of the global best value of all the runs

## References
[[1] Hans-Georg Beyer and Hans-Paul Schwefel, Evolution strategies - a comprehensive introduction, Natural Computing 1 (2002), no. 1, 3–52](https://link.springer.com/article/10.1023/A:1015059928466)

[[2] A. E. Eiben and James E Smith, Introduction to Evolutionary Computing, Natural Computing Series (2003), no. 1, Springer-Verlag Berlin Heidelberg, 300](https://www.springer.com/in/book/9783642072857)