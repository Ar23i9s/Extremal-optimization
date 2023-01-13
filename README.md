
# C++ Implementation of Extremal Optimization for Modularity and In-block Nestedness

This code is a C++ program built on **Visual Studio 2022 17.4.4**.

## Requirements
Operating System: Windows

Compiler: MSVC 14.34, igraph C library 0.9.9 (https://igraph.org/c/html/0.9.9/)

## Usage
### Command line
```
<exe_path> <csv_path> <object_flag> <alpha> <recursive> <loop>
```

### Input
1) `exe_path`: The path of the executable file built from the source code. **The file `eo.exe` can be directly used.**
2) `csv_path`: The path of the CSV file that records the adjacency matrix of a network.
Given a 3-node network, the content in the CSV file should be like:
```         
0,1,1
1,0,1
1,1,0 
```
3) `object_flag`: A binary variable decides the optimization target, 0 for modularity-based optimization, 1 for in-block nestedness-based optimization. 
4) `alpha`: The parameter α in the process of extremal optimization. If the optimized value is not improved in `α*N` steps(N-the number of nodes), the optimization process will stop. 
5) `recursive`: The number of times the whole network will be recursively partitioned. In the first recursion, the whole network will be split into two partitions. In the second recursion, each partition will be split into two smaller partitions.
6) `loop`: The number of times the code will run, corresponding to the number of optimized values/partitions.

### output
The optimized results will be saved in a CSV file named by the above input parameters, stored in `csv_path`.
1) For modularity-based optimization(`object_flag = 0`), the results should be like: 
```
<Q1>
<R1>
<C1>
<T1>      
<Q2>
<R2> 
<C2>
<T2>      
...
```
`Qx`: The optimized modularity in the xth loop.

`Rx`: The block index of each row node(country, company, etc.), in the xth loop.

`Cx`: The block index of each col node(product, technology, etc.), in the xth loop.

`Tx`: The running time of the xth loop.

2) For in-block nestedness-based optimization(`object_flag = 1`), the results should be like: 
``` 
<N>    
<I1>    
<R1> 
<C1>
<T1>      
<I2>      
<R2>      
<C2>
<T2>      
...  
``` 
`N`: The value of global nestedness.

`Ix`: The optimized in-block nestedness in the xth loop. 

`Rx`: The block index of each row node(country, company, etc.), in the xth loop.

`Cx`: The block index of each col node(product, technology, etc.), in the xth loop.

`Tx`: The running time of the xth loop. 

 By computing the ratio `N/Ix`, we can determine whether nesting is a local or global attribute.

## Example
#### File structure
```
C:\Users\Administrator\eo_ibn
│  eo.exe
│  matrix.csv
```

#### Command shell
```
C:\Users\Administrator>cd eo_ibn
C:\Users\Administrator\eo_ibn>eo.exe matrix.csv 1 0.5 3 3
```

#### Output
```
C:\Users\Administrator\eo_ibn
│  eo.exe
│  matrix.csv
│  matrix_0.5_3_ibn.txt
```

In matrix_0.5_3_ibn.txt
```
0.07878
0.129972
1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0
running time 0.11s
0.167711
1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0
1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0
running time 0.178s
0.147363
0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1
0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1
running time 0.135s
```
