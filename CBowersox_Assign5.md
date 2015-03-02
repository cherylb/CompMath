
# Week 5 Assignment
Cheryl Bowersox 
IS 605 - 1
3/1/15

========================================================
## Problem Set 1

```r
# function for converting matrix format from R
convmat <- function(X){
  M <- paste("\\begin{bmatrix}","\n")
  
  i = 1
  n = nrow(X)
  if(is.null(n)){
    n <- length(X)
    while(i <= n){
       M <- paste(M, paste(round(X[i],4), collapse = " & "),"\\\\ \n")
       i = i + 1
    }
  }else{
    while(i <= n){
       M <- paste(M, paste(round(X[i,],4), collapse = " & "),"\\\\ \n")
       i = i + 1
    }
  }
  M <- paste(M, "\\end{bmatrix}","\n", collapse = "")
  return (M)
}

A <- matrix(c(1,1,1,1,0,1,3,4), nrow = 4)

b <- matrix(c(0,8,8,19), nrow = 4)
```
$$

  A =
  \begin{bmatrix} 
 1 & 0 \\ 
 1 & 1 \\ 
 1 & 3 \\ 
 1 & 4 \\ 
 \end{bmatrix} 
\\ 
 
  A^TA = 
  \begin{bmatrix} 
 4 & 8 \\ 
 8 & 26 \\ 
 \end{bmatrix} 
\\
  
  A^T\vec{b} = 
  \begin{bmatrix} 
 35 \\ 
 108 \\ 
 \end{bmatrix} 
\\
  
$$

###  least squares solution and squared error
  

```r
x <- solve(t(A)%*%A) %*% t(A) %*% b

eb <- b - (A %*% x)
E <- sum(eb^2)
```

The least squares solution  for
$$
  \vec{b} = \begin{bmatrix} 
 0 \\ 
 8 \\ 
 8 \\ 
 19 \\ 
 \end{bmatrix} 
\\
  
  
  \hat{x} = 
  \begin{bmatrix} 
 1.15 \\ 
 3.8 \\ 
 \end{bmatrix} 
\\ 
$$  


The resulting squared error $E_b$ is:
$$
  E_b = 38.35
  
$$

### system is consistent and solvable for p

```r
# using p to solve system 
p <- matrix(c(1,5,13,17), nrow = 4)
x <- solve(t(A)%*%A) %*% t(A) %*% p

e <- round(p - (A %*% x),6)

E <- round(sum(e^2), 4)
```

The least squares solution  for
$$
 \vec{p} = \begin{bmatrix} 
 1 \\ 
 5 \\ 
 13 \\ 
 17 \\ 
 \end{bmatrix} 
 is\\
 
 
 \hat{x} = 
  \begin{bmatrix} 
 1 \\ 
 4 \\ 
 \end{bmatrix} 
\\
  
$$
With error vector $\vec{e}$ of zero:
$$
  \vec{e} = 
  \begin{bmatrix} 
 0 \\ 
 0 \\ 
 0 \\ 
 0 \\ 
 \end{bmatrix} 
\\
$$

the zero error indicates that $\hat{x} = \vec{x}$, and $\vec{p} $ is consistent


### e is for orthoganol to p  for C(A)
The error vector $\vec{e}$ is the difference between $\vec{b}$ and $\vec{p}$:
$$
  \vec{e} = \begin{bmatrix} 
 -1 \\ 
 3 \\ 
 -5 \\ 
 2 \\ 
 \end{bmatrix} 
\\
$$

The dot product of $\vec{e} and \vec{p}$ is 0 
e is orthogonal to p


```r
i <- 1
B <- rep(NA,ncol(A))
while (i <= ncol(A)){
    B[i] <- (round(sum(e * A[,i])))
    i <- i +1
}
```
The dot products of $\vec{e}$ across the columns of A are

0, 0

$\vec{e}$ is orthogonal to columns of A



========================================================
## Problem Set 2

