---
title: "Volesti tutorial"
author:
- name: Vissarion Fisikopoulos
date: "17 October 2019"
output:
  html_document:
    df_print: paged
    toc: true
    toc_depth: 2
    keep_md: true
  pdf_document: default
params:
  printcode: false  # or set it to true
---

<!--
This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 
-->

***

`volesti` is a `C++` package (with an `R` interface) for computing estimations of volume of polytopes given by a set of points or linear inequalities or Minkowski sum of segments (zonotopes). There are two algorithms for volume estimation and algorithms for sampling, rounding and rotating polytopes.

We can download the `R` package from https://CRAN.R-project.org/package=volesti


```r
# first load the volesti library
library(volesti)
```

```
## Loading required package: Rcpp
```

```r
packageVersion("volesti")
```

```
## [1] '1.0.3'
```

You have access to the documentation of volesti functions like volume computation and sampling.


```r
help("volume")
help("sample_points")
```

#Sampling

Sampling in the square.


```r
library(ggplot2)
P = GenCube(2, 'H')
points1 = sample_points(P, WalkType = "RDHR", walk_step = 1, N=1000)
g<-ggplot(data.frame( x=points1[1,], y=points1[2,] )) + geom_point( aes(x=x, y=y))
g<-g+annotate("path",
   x=cos(seq(0,2*pi,length.out=100)),
   y=sin(seq(0,2*pi,length.out=100)),color="red")+coord_fixed()
plot(g)
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

Naive Monte Carlo fails in (not so) high dimensions.


```r
for (d in 2:20) {
  P = GenCube(d, 'H')

  num_of_points <- 100000
  count_inside <- 0
  tim1 = system.time({ points1 = sample_points(P, WalkType = "RDHR", walk_step = 1, N=num_of_points) })

  for (i in 1:num_of_points) {
    if (norm(points1[,i], type="2") < 1) {
      count_inside <- count_inside + 1
    }
  }
  vol_estimation <- count_inside*2^d/num_of_points
  vol_exact <- pi^(d/2)/gamma(d/2+1) 
  
  cat(d, vol_estimation, vol_exact, abs(vol_estimation- vol_exact)/ vol_exact, "\n")
}
```

```
## 2 3.12408 3.141593 0.005574451 
## 3 4.224 4.18879 0.008405719 
## 4 4.96608 4.934802 0.006338207 
## 5 5.3088 5.263789 0.008551062 
## 6 5.11936 5.167713 0.009356708 
## 7 4.76288 4.724766 0.008066861 
## 8 4.18304 4.058712 0.03063235 
## 9 3.24608 3.298509 0.01589473 
## 10 2.22208 2.550164 0.1286521 
## 11 2.10944 1.884104 0.1195986 
## 12 1.88416 1.335263 0.4110781 
## 13 1.14688 0.9106288 0.2594375 
## 14 0.32768 0.5992645 0.4531964 
## 15 0.65536 0.3814433 0.718106 
## 16 0.65536 0.2353306 1.784848 
## 17 0 0.1409811 1 
## 18 0 0.08214589 1 
## 19 0 0.0466216 1 
## 20 0 0.02580689 1
```

##Sampling via random walks

`volesti` supports 3 types of random walks 

1. Ball walk 

<img src="figs/ball1.png" alt="drawing" width="150"/> 
<img src="figs/ball2.png" alt="drawing" width="150"/> 
<img src="figs/ball3.png" alt="drawing" width="150"/> 
<img src="figs/ball4.png" alt="drawing" width="150"/> 
<img src="figs/ball5.png" alt="drawing" width="150"/> 

2. Random directions hit-and-run

<img src="figs/hnr1.png" alt="drawing" width="150"/> 
<img src="figs/hnr2.png" alt="drawing" width="150"/> 
<img src="figs/hnr3.png" alt="drawing" width="150"/> 

3. Coordinate directions hit-and-run

<img src="figs/cdhr1.png" alt="drawing" width="150"/> 
<img src="figs/cdhr2.png" alt="drawing" width="150"/> 
<img src="figs/cdhr3.png" alt="drawing" width="150"/> 

There are two important parameters `cost per step` and `mixing time` that affects the accuracy and performance of the walks. Below we illustrate this by choosing different walk steps for each walk while sampling on the 100-dimensional cube.



```r
library(ggplot2)
library(volesti)
for (step in c(1,20,100,150)){
  for (walk in c("CDHR", "RDHR", "BW")){
    P <- GenCube(100, 'H')
    points1 <- sample_points(P, WalkType = walk, walk_step = step, N=1000)
    g<-plot(ggplot(data.frame( x=points1[1,], y=points1[2,] )) + geom_point( aes(x=x, y=y, color=walk)) + coord_fixed(xlim = c(-1,1), ylim = c(-1,1)) + ggtitle(sprintf("walk length=%s", step, walk)))
  }
}
```

<img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-1.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-2.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-3.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-4.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-5.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-6.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-7.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-8.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-9.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-10.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-11.png" width="33%" /><img src="volesti_tutorial_files/figure-html/unnamed-chunk-5-12.png" width="33%" />

***

#Volume estimation

Now let's compute our first example. The volume of the 3-dimensional cube. 


```r
library(geometry)

PV <- GenCube(3,'V')
str(PV)
```

```
## Reference class 'Rcpp_Vpolytope' [package "volesti"] with 3 fields
##  $ V        : num [1:8, 1:3] -1 1 -1 1 -1 1 -1 1 -1 -1 ...
##  $ dimension: num 3
##  $ type     : int 2
##  and 16 methods, of which 2 are  possibly relevant:
##    finalize, initialize
```

```r
#P = GenRandVpoly(3, 6, body = 'cube')
tim1 <- system.time({ geom_values = convhulln(PV$V, options = 'FA') })
tim2 <- system.time({ vol2 = volume(PV) })

cat(sprintf("exact vol = %f\napprx vol = %f\nrelative error = %f\n", 
            geom_values$vol, vol2, abs(geom_values$vol-vol2)/geom_values$vol))
```

```
## exact vol = 8.000000
## apprx vol = 7.776406
## relative error = 0.027949
```

Now try a higher dimensional example. By setting the `error` parameter we can control the apporximation of the algorithm.


```r
PH = GenCube(10,'H')
volumes <- list()
for (i in 1:10) {
  volumes[[i]] <- volume(PH, error=1) # default parameters
}
options(digits=10)
summary(as.numeric(volumes))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  987.4002 1007.0985 1041.0874 1038.4261 1054.6985 1113.2900
```


```r
volumes <- list()
for (i in 1:10) {
  volumes[[i]] <- volume(PH, error=0.5)
}
summary(as.numeric(volumes))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 1019.978 1031.679 1039.877 1038.224 1044.439 1059.876
```

Deterministic algorithms for volume are limited to low dimensions (e.g. less than $15$)


```r
library(geometry)

P = GenRandVpoly(15, 30)
# this will return an error about memory allocation, i.e. the dimension is too high for qhull
tim1 <- system.time({ geom_values = convhulln(P$V, options = 'FA') })

#warning: this also takes a lot of time in v1.0.3
print(volume(P))
```

```
## [1] 3.605252155e-10
```

##Volume of Birkhoff polytopes

We now continue with a more interesting example, the 10-th Birkhoff polytope. It is known from https://arxiv.org/pdf/math/0305332.pdf that its volume equals 

$\text{vol}(\mathcal{B}_{10}) = \frac{727291284016786420977508457990121862548823260052557333386607889}{828160860106766855125676318796872729344622463533089422677980721388055739956270293750883504892820848640000000}$

obtained via massive parallel computation.


```r
library(volesti)
P <- fileToMatrix('data/birk10.ine')
exact <- 727291284016786420977508457990121862548823260052557333386607889/828160860106766855125676318796872729344622463533089422677980721388055739956270293750883504892820848640000000 

# warning the following will take around half an hour
#print(volume(P, Algo = 'CG'))
```

Compare our computed estimation with the "normalized" floating point version of  $\text{vol}(\mathcal{B}_{10})$


```r
n <- 10
vol_B10 <- 727291284016786420977508457990121862548823260052557333386607889/828160860106766855125676318796872729344622463533089422677980721388055739956270293750883504892820848640000000
print(vol_B10/(n^(n-1)))
```

```
## [1] 8.782005031e-55
```

***

#Rounding

We generate skinny polytopes, in particular skinny cubes of the form $\{x=(x_1,\dots,x_d)\ |\ x_1\leq 100, x_1\geq-100,x_i\leq 1,x_i\geq-1,x_i\in \mathbb{R}, \text{ for } i=2,\dots,d\}$. 

Our random walks perform poorly on those polytopes espacially as the dimension increases. Note that if we use the `CDHR` walk here is cheating since we take advanatage of the instance structure.


```r
library(ggplot2)

P = GenSkinnyCube(2)
points1 = sample_points(P, WalkType = "CDHR", N=1000)
ggplot(data.frame(x = c(points1[1,]), y = c(points1[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed(ylim = c(-10,10))
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
points1 = sample_points(P, WalkType = "RDHR", N=1000)
ggplot(data.frame(x = c(points1[1,]), y = c(points1[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed(ylim = c(-10,10))
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-12-2.png)<!-- -->


```r
P = GenSkinnyCube(10)
points1 = sample_points(P, WalkType = "CDHR", N=1000)
ggplot(data.frame(x = c(points1[1,]), y = c(points1[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed(xlim = c(-100,100), ylim = c(-10,10))
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
points1 = sample_points(P, WalkType = "RDHR", N=1000)
ggplot(data.frame(x = c(points1[1,]), y = c(points1[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed(xlim = c(-100,100), ylim = c(-10,10))
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

```r
P = GenSkinnyCube(100)
points1 = sample_points(P, WalkType = "RDHR", N=1000)
ggplot(data.frame(x = c(points1[1,]), y = c(points1[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed(xlim = c(-100,100), ylim = c(-10,10))
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-13-3.png)<!-- -->


Then we examine the problem of rounding by sampling in the original and then in the rounded polytope and look at the effect in volume computation.


```r
library(ggplot2)

d <- 10

P = GenSkinnyCube(d)
points1 = sample_points(P, WalkType = "CDHR", N=1000)
ggplot(data.frame(x = c(points1[1,]), y = c(points1[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed(ylim = c(-10,10))
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
P <- rand_rotate(P)$P

points1 = sample_points(P, WalkType = "RDHR", N=1000)
ggplot(data.frame(x = c(points1[1,]), y = c(points1[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed()
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-14-2.png)<!-- -->

```r
exact <- 2^d*100
cat("exact volume                 = ", exact , "\n")
```

```
## exact volume                 =  102400
```

```r
cat("volume estimation (no round) = ", volume(P, WalkType = "RDHR", rounding=FALSE), "\n")
```

```
## volume estimation (no round) =  30625.55588
```

```r
cat("volume estimation (rounding) = ", volume(P, WalkType = "RDHR", rounding=TRUE), "\n")
```

```
## volume estimation (rounding) =  132716.6547
```

```r
# 1st step of rounding
res1 = round_polytope(P)
points2 = sample_points(res1$P, WalkType = "RDHR", N=1000)
ggplot(data.frame(x = c(points2[1,]), y = c(points2[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed()
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-14-3.png)<!-- -->

```r
volesti <- volume(res1$P) * res1$round_value
cat("volume estimation (1st step) = ", volesti, " rel. error=", abs(exact-volesti)/exact,"\n")
```

```
## volume estimation (1st step) =  103961.1625  rel. error= 0.01524572708
```

```r
# 2nd step of rounding
res2 = round_polytope(res1$P)
points2 = sample_points(res2$P, WalkType = "RDHR", N=1000)
ggplot(data.frame(x = c(points2[1,]), y = c(points2[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed()
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-14-4.png)<!-- -->

```r
volesti <- volume(res2$P) * res1$round_value * res2$round_value
cat("volume estimation (2nd step) = ", volesti, " rel. error=", abs(exact-volesti)/exact,"\n")
```

```
## volume estimation (2nd step) =  111074.6912  rel. error= 0.0847137808
```

```r
# 3rd step of rounding
res3 = round_polytope(res2$P)
points2 = sample_points(res3$P, WalkType = "RDHR", N=1000)
ggplot(data.frame(x = c(points2[1,]), y = c(points2[2,])), aes(x=x, y=y)) + geom_point() +labs(x =" ", y = " ")+coord_fixed()
```

![](volesti_tutorial_files/figure-html/unnamed-chunk-14-5.png)<!-- -->

```r
volesti <- volume(res3$P) * res1$round_value * res2$round_value * res3$round_value
cat("volume estimation (3rd step) = ", volesti, " rel. error=", abs(exact-volesti)/exact,"\n")
```

```
## volume estimation (3rd step) =  105218.0717  rel. error= 0.02752023144
```

***

#Integration

We can use sampling and volume estimation to estimate integrals over polyhedral domains. Below there is an example with a degree 2 polynomial over a 3-dimensional cube. 


```r
library(cubature) # load the package "cubature"
f <- function(x) { 2/3 * (2 * x[1]^2 + x[2] + x[3]) + 10 }  # "x" is vector 
adaptIntegrate(f, lowerLimit = c(-1, -1, -1), upperLimit = c(1, 1, 1))$integral
```

```
## [1] 83.55555556
```

```r
# Simple Monte Carlo integration
# https://en.wikipedia.org/wiki/Monte_Carlo_integration
P = GenCube(3, 'H')
num_of_points <- 10000
points1 <- sample_points(P, WalkType = "RDHR", walk_step = 100, N=num_of_points)
int<-0
for (i in 1:num_of_points){
  int <- int + f(points1[,i])
}
V <- volume(P)
print(int*V/num_of_points)
```

```
## [1] 85.93261174
```

# Counting linear extensions

Let $G= (V, E)$ be an acyclic digraph with $V= [n] :=\{1,2, . . . , n\}$. One might want to consider $G$ as a representation of the partially ordered set (poset) $V:i > j$ if and only if there is a directed path from node $i$ to node $j$.A permutation $\pi$ of $[n]$ is called a linear extension of $G$ (or the associated poset $V$) if $\pi^{−1}(i)> \pi^{−1}(j)$ for every edge $i\rightarrow j \in E$.

Let $P_{LE}(G)$ be the polytope in $R^n$ defined by $P_{LE}(G) ={x\in R^n|1\geq x_i \geq 0 for all i=1,2,\dots ,n, and $x_i\geq x_j$ for all directed edges $i\rightarrow j \in E$.

A well known result is $#_{LE}G=vol(P)*n!$.

The following example from https://inf.ethz.ch/personal/fukudak/lect/pclect/notes2016/expoly_order.pdf has $9$ linear extensions. See also https://www.aaai.org/ocs/index.php/AAAI/AAAI18/paper/viewFile/16957/15838 for counting linear extensions in practice.

![graph](figs/graph.png)

We can approximate this number by the following code.


```r
A = matrix(c(-1,0,1,0,0,0,-1,1,0,0,0,-1,0,1,0,0,0,0,-1,1,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,-1,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,-1), ncol=5, nrow=14, byrow=TRUE)
b = c(0,0,0,0,0,0,0,0,0,1,1,1,1,1)
P = Hpolytope$new(A, b)
volume(P,error=0.2)*factorial(5)
```

```
## [1] 8.914854383
```

#Questions

##Volume of intersection of polytopes

How can we compute the volume of the intersection of two V-polytopes? What is the complexity? Can you provide any implementations (deterministic or randomized).

#Volume of Minkowski sums

The Minkowski sum of two sets (e.g. polytopes) A and B in Euclidean space is formed by adding each vector in A to each vector in B, i.e., the set
$${\displaystyle A+B=\{a +b \ |\ a \in A,\ b \in B\}.}$$

1. Describe an algorithm for the volume of the Minkowski sum.  
2. What is the complexity of your algorithm?
3. Implement you algorithm.
