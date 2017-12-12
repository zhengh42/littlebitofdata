---
title: "Variance, Covariance, and Correlation"
date: '2017-12-08'
slug: statsbasic
categories: ["data at fingertips"]
tags: ["statistics"]
---

__Variance__

Variance is the difference between when we square the inputs to expectation and when we square the expectation itself.

`$Var(x) =  \frac {1}{n} \sum_{i=1}^{n} (x_{i} - \bar{x})^2 $` =   
`$ E[(x-\bar{x})^2]$` =   
`$ E[x^2 - 2 \times x \times \bar{x} + (\bar{x})^2] $` =  
`$ E[x^2] -2 \times E[X] \times E[\bar{x}] + E[(\bar{x})^2] $` =   
`$ E[x^2] - (E[x])^2$`

Note:   

* `$x$` is a vector (`$n * 1$` matrix).
* `$E[\bar{x}] ==  E[x] == \bar{x}$`

__Covariance__

It measures the variance between two variables.

We can rewrite the variance equation as:

`$Var(x) =E[xx]−E[x]E[x]$`

What if one of the `$x$` is another random variable?", so that we would have:

`$E[xy]−E[x]E[y]$`

which is the definition of covariance between `$x$` and `$y$`: `$Cov(x,y)$`

It can also be written as `$ \frac {1}{n} \sum_{i=1}^{n} (x_{i} - \bar{x})(y_{i} - \bar{y}) $`

Note: `$x$` and `$y$` are both vectors (`$n * 1$` matrix).

__Correlation__

`$Cor(x,y) = \frac {Cov(x,y)} {\sqrt{(Var(x)Var(y))}}$` =
`$\frac {\sum_{i=1}^{n} (x_{i} - \bar{x})(y_{i} - \bar{y})} {\sqrt{ \sum_{i=1}^{n} (x_{i} - \bar{x})^2} \sqrt{\sum_{i=1}^{n}(y_{i} - \bar{y})^2}}  $`

It is the Pearson correlation coefficient between variables `$x$` and `$y$`.

Covariance is just an unstandardized version of correlation. To compute any correlation, we divide the covariance by the standard deviation of both variables to remove units of measurement. So a covariance is just a correlation measured in the units of the original variables. 

Note: `$x$` and `$y$` are both vectors ( `$n * 1$` matrix).

__Covariance matrix__

`$$
\left(\begin{array}{cc} 
s_{1}^2 & s_{12} & ... & s_{1p} \\
s_{21} & s_{2}^2 & ... & s_{2p} \\
... & ... & ... & ... \\
s_{p1} & s_{p2} & ... & s_{p}^2
\end{array}\right)
$$`

* `$X$` is a `$n * p$` matrix.
* `$p$`: number of features 
* `$n$`: number of observations
* `$s_{j}^2$` is the variance of the j-th variable. `$ \frac {1}{n} \sum_{i=1}^{n} (x_{ij} - \bar{x_{j}})^2 $`
* `$s_{jk}$` is the covariance between the j-th and k-th variables. `$ \frac {1}{n} \sum_{i=1}^{n} (x_{ij} - \bar{x_{j}})(x_{ik} - \bar{x_{k}}) $`


In matrix form:  
`$$ S  = \frac {1} {n} Xc^TXc$$`

or

`$$ S  = \frac {1} {n-1} Xc^TXc$$`
 
* `$Xc$`, centered matrix of `$X$`. `$Xc =  X - 1_{n} \bar{X}'  $`.  `$\bar{X}'$` is column means of `$X$`, in the form of a `$1*p$` matrix. `$1_{n}$` is a `$n*1$` matrix.

`$$
\left(\begin{array}{cc} 
x_{11}-\bar{x_{1}} & x_{12}-\bar{x_{2}}  & ... & x_{1p}-\bar{x_{p}}  \\
x_{21}-\bar{x_{1}} & x_{22}-\bar{x_{2}}  & ... & x_{2p}-\bar{x_{p}} \\
... & ... & ... & ... \\
x_{n1}-\bar{x_{1}} & x_{n2}-\bar{x_{2}}  & ... & x_{np}-\bar{x_{p}}
\end{array}\right)
$$`

*  Sometimes (and in R) it is also divided by `$n−1$`, which is a typical way to correct for the bias introduced by using the sample mean instead of the true population mean. 


Calculate covariance matrix in R:

```{r}
S <- cov(X)
```

__Correlation matrix__

`$$
\left(\begin{array}{cc} 
1 & r_{12} & ... & r_{1p} \\
r_{21} & 1 & ... & r_{2p} \\
... & ... & ... & ... \\
r_{p1} & r_{p2} & ... & 1
\end{array}\right)
$$`

where

`$ r_{jk} = \frac {s_{jk}}{s_{j}s_{k}}  $` =
`$\frac {\sum_{i=1}^{n} (x_{ij} - \bar{x_{j}})(x_{ik} - \bar{x_{k}})} {\sqrt{ \sum_{i=1}^{n} (x_{ij} - \bar{x_{j}})^2} \sqrt{\sum_{i=1}^{n}(x_{ik} - \bar{x_{k}})^2}}  $` is the Pearson correlation coefficient between variables `$x_{j}$` and `$x_{k}$`.

In matrix form:  
`$$ R  = \frac {1} {n} Xs^TXs$$`

or

`$$ R  = \frac {1} {n-1} Xs^TXs$$`

* `$Xs = XcD^{-1}$`, where `$D = diag(s_{1}, . . . , s_{p})$` is the diagonal scaling matrix.

`$$
\left(\begin{array}{cc} 
\frac {x_{11}-\bar{x_{1}}}{s_{1}} & \frac {x_{12}-\bar{x_{2}}}{s_{2}} & ... & \frac {x_{1p}-\bar{x_{p}}}{s_{p}} \\
\frac {x_{21}-\bar{x_{1}}}{s_{1}} & \frac {x_{22}-\bar{x_{2}}}{s_{2}} & ... & \frac {x_{2p}-\bar{x_{p}}}{s_{p}} \\
... & ... & ... & ... \\
\frac {x_{21}-\bar{x_{1}}}{s_{1}} & \frac {x_{n2}-\bar{x_{2}}}{s_{2}} & ... & \frac {x_{np}-\bar{x_{p}}}{s_{p}}
\end{array}\right)
$$`
 
Calculate correlation matrix in R:

```{r}
S <- cor(X)
```

Further readings:  
https://www.countbayesie.com/blog/2015/2/21/variance-co-variance-and-correlation
http://www.theanalysisfactor.com/covariance-matrices/  
http://users.stat.umn.edu/~helwig/notes/datamat-Notes.pdf  
【__Linear Algebra__】 http://www4.ncsu.edu/~slrace/LinearAlgebra.pdf 


