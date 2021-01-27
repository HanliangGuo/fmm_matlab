# fmm_matlab
This is an implementation of a (toy) 2D Laplace problem using fast multipole method (FMM) in Matlab. 
The algorithm follows L. Greengard & V. Rokhlin, Journal of Computational Physics, 1987
Some variable names are motivated by L. Ying et al., Journal of Computational Physics, 2004

## What does the code do
N source points are randomly distributed in a unit box centered at (0.5, 0.5)
The strengths of the charges are also random.
The fast multipole method is used to evaluate the potential at the target points (collocated with the source points).
The results are compared to direct evaluation as a validation.

## Intent of the project
This code was a self-exercise when I was in grad school at the University of Southern California.
I used this code to learn the algorithm of the fast multipole method (FMM). 
The code is reasonably optimized and one could test the complexity and convergence by changing the inputs. 
I shall stress that this code is only for “proof-of-concept” purpose, as there are plenty of existing (high-performance) FMM libraries that are far superior than this implementation (notably at Courant mathematics and Computing Laboratory).
Nevertheless, this code might be a reasonable starting point for people who are familiar with Matlab and wants to get a hand on FMM.
