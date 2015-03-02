
# Week 4 Assignment
Cheryl Bowersox 
IS 605 - 1
2/22/15

========================================================
## Problem Set 1

$$

  A = \begin{bmatrix}
  1 & 2 & 3\\
  -1 & 0 & 4\end{bmatrix}\\
  
  
  X = AA^t
  Y = A^tA

$$


```r
library("MASS")
A <- matrix(c(1,-1,2,0,3,4),nrow= 2)

X <- A %*% t(A)
Y <- t(A) %*% A
```
$$

  X = \begin{bmatrix}
  14 & 11\\
  11 & 17\end{bmatrix}\\


  Y = \begin{bmatrix}
  2 & 2 & -1\\
  2 & 4 & 6\\
  -1 & 6 & 25\end{bmatrix}\\

$$

### Eigenvalues and vectors for X, Y

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
```

Eigenvalues and Vectors for X:
$$
  W_{x1} = 
  \begin{bmatrix} 
 0.6576 \\ 
 0.7534 \\ 
 \end{bmatrix} 
     
 
  \lambda_{x1} = 26.6018\\
  
   W_{x2} =
  \begin{bmatrix} 
 \end{bmatrix} 
    
 
  \lambda_{x2} = 4.3982\\
$$  

Eigenvalues and Vectors for Y:
$$
  W_{Y1} = 
  \begin{bmatrix} 
 -0.0186 \\ 
 0.255 \\ 
 0.9668 \\ 
 \end{bmatrix} 
     
 
  \lambda_{Y1} = 26.6018\\
  
   W_{Y2} =
  \begin{bmatrix} 
 -0.6728 \\ 
 -0.7185 \\ 
 0.1766 \\ 
 \end{bmatrix} 
    
 
  \lambda_{21} = 4.3982\\
  
   W_{Y3} =
  \begin{bmatrix} 
 0.7396 \\ 
 -0.6472 \\ 
 0.1849 \\ 
 \end{bmatrix} 
    
 
  \lambda_{31} = 0\\
$$  


### SVD of A

```r
# X = A 
v_x <- eigen(X)$vectors
v_y <- eigen(Y)$vectors
e_x <- eigen(X)$values
e_y <- eigen(Y)$values

U <- svd(A)$u
D <- diag(svd(A)$d)
V <- svd(A)$v
```

### singular vectors are the eigenvectors of X and Y
(1) U and V are eigenvectors(W) of $X$ iff $XW = \lambda_W$

(2) The columns of U are equal to scaler multiples of the eigenvectors of X given by the eigen() function
(3) If the columns of U are a scaled version of the eigenvectors of X, they must have equivalent eigenvalues

$$
  U = 
  \begin{bmatrix} 
 -0.6576 & -0.7534 \\ 
 -0.7534 & 0.6576 \\ 
 \end{bmatrix} 
     
$$

(4) eigenvalues of X as $\lambda$ and the corresponding vectors of U satisfy the conditions of (1) for each column vector of U

**U1**
$$
  W_{u1} = 
  \begin{bmatrix} 
 -0.6576 \\ 
 -0.7534 \\ 
 \end{bmatrix} 
     
 
  \lambda_{x1} = 26.6018
$$

```r
C <- X %*% U[,1]
E <- e_x[1] * U[,1]
```
$$ 
  W_{u1} = 
  \begin{bmatrix} 
 -17.4935 \\ 
 -20.0408 \\ 
 \end{bmatrix} 
    = 
  \begin{bmatrix} 
 -17.4935 \\ 
 -20.0408 \\ 
 \end{bmatrix} 
    =
  W_{x1}\\
$$

**U2**
 
$$
  W_{u2} = 
  \begin{bmatrix} 
 -0.7534 \\ 
 0.6576 \\ 
 \end{bmatrix} 
     
  \lambda_{x2} = 4.3982
$$

```r
C <- X %*% U[,2]
E <- e_x[2] * U[,2]
```
$$ 
  W_{u2} = 
  \begin{bmatrix} 
 -3.3134 \\ 
 2.8923 \\ 
 \end{bmatrix} 
    = 
  \begin{bmatrix} 
 -3.3134 \\ 
 2.8923 \\ 
 \end{bmatrix} 
    =
  W_{x2}\\
$$


(4) Similarly to (2), V contains eigenvectors of Y
(5) eigenvalues of Y as $\lambda$ and the corresponding vectors of V satisfy the conditions of (1) for each column vector of V

**V1**
$$
  W_{v1} = 
  \begin{bmatrix} 
 0.0186 \\ 
 -0.255 \\ 
 -0.9668 \\ 
 \end{bmatrix} 
     

  \lambda_{y1} = 26.6018
$$

```r
C <- Y %*% V[,1]
E <- e_y[1] * V[,1]
```
$$ 
  W_{v1} = 
  \begin{bmatrix} 
 0.4939 \\ 
 -6.7834 \\ 
 -25.7176 \\ 
 \end{bmatrix} 
    = 
  \begin{bmatrix} 
 0.4939 \\ 
 -6.7834 \\ 
 -25.7176 \\ 
 \end{bmatrix} 
    =
  W_{y1}\\
$$

**V2**
 
$$
  W_{v2} = 
  \begin{bmatrix} 
 -0.6728 \\ 
 -0.7185 \\ 
 0.1766 \\ 
 \end{bmatrix} 
     
  \lambda_{v2} = 4.3982
$$

```r
C <- Y %*% V[,2]
E <- e_y[2] * V[,2]
```
$$ 
  W_{v2} = 
  \begin{bmatrix} 
 -2.9591 \\ 
 -3.1599 \\ 
 0.7766 \\ 
 \end{bmatrix} 
    = 
  \begin{bmatrix} 
 -2.9591 \\ 
 -3.1599 \\ 
 0.7766 \\ 
 \end{bmatrix} 
    =
  W_{y2}\\
$$

**V3**

$\lambda_{v3} = 0$

The 3rd eigenvalue of V is essentially zero, resulting in the trivial solution and is not shown.


(6) The non-zero eigenvalues for X and Y are equal to eachother
$$
 \lambda_{Y1} = 26.6018  =  \lambda_{x1} = 26.6018\\
 \lambda_{Y2} = 4.3982  =  \lambda_{x1} = 4.3982\\
$$

(7) the eigenvalues for X and Y are equal to square of non-zero singular values of A

$$ 
  D = 
  \begin{bmatrix} 
 5.1577 & 0 \\ 
 0 & 2.0972 \\ 
 \end{bmatrix} 
\\
  
  
5.1577^2 = 26.6018 = \lambda_1\\
2.0972^2 = 4.3982 = \lambda_2\\
$$


========================================================
## Problem Set 2

### Inverse of a Matrix A

```r
myinverse <- function(A){
  #  function to compute the inverse of well-conditioned, full rank square matrix using cofactors
  n <- nrow(A)
  i = 1
  j = 1
  w = 0
  C = diag(0,n)
  while(i <= n){ #for each row in A
    j = 1
    while(j <= n){  #for each column A
      M <- A[-i,-j]
      s = (-1)^w
      w = w + 1
      C[i,j] <- s * det(M)   # assign to cofactor matrix
      j = j + 1
    }
    i = i + 1
  }
  B <- t(C) * 1/det(A)
return(B)
}   

# example
A <- matrix(c(1,2,3,1,2,4,3,5,3),nrow = 3)
B <- myinverse(A)
(C <- A %*% B)
```

```
##            [,1]       [,2] [,3]
## [1,]  1.000e+00 -4.441e-16    0
## [2,] -3.553e-15  1.000e+00    0
## [3,]  0.000e+00 -3.109e-15    1
```

```r
#rounding to 8 places to compare to I
round(C,8) == diag(1,nrow(A))  #rounded C = I
```

```
##      [,1] [,2] [,3]
## [1,] TRUE TRUE TRUE
## [2,] TRUE TRUE TRUE
## [3,] TRUE TRUE TRUE
```


$$ 
  A = 
  \begin{bmatrix} 
 1 & 1 & 3 \\ 
 2 & 2 & 5 \\ 
 3 & 4 & 3 \\ 
 \end{bmatrix} 
\\
  
  A^-1 = 
  \begin{bmatrix} 
 -14 & 9 & -1 \\ 
 9 & -6 & 1 \\ 
 2 & -1 & 0 \\ 
 \end{bmatrix} 
\\
  
  A %*% A^-1 =
  \begin{bmatrix} 
 1 & 0 & 0 \\ 
 0 & 1 & 0 \\ 
 0 & 0 & 1 \\ 
 \end{bmatrix} 
    =  I
  
  
$$



