---
title: MTH3045 May Exam 2024 
author: "** 154268 **"
output: pdf_document
fontsize: 11pt
classoption: a4paper
---

# Q1

## Q1(a)

```{r}
Q <- cbind(
  c( 0.856, -0.480, -0.192),
  c(-0.480, -0.600, -0.640), 
  c(-0.192, -0.640,  0.744)
)
R <- cbind(
  c(0.53, 0.00, 0.00),
  c(0.59, 0.85, 0.00),
  c(0.36, 0.57, 0.40)
)
b <- c(3.080, 4.376, 2.945)
```

### Q1(a)i
```{r}
all.equal(crossprod(Q),diag(ncol(Q)))
all(R[lower.tri(R)]==0)
sum(abs(diag(R)))
```

### Q1(a)ii
```{r}
x1 <- backsolve(R, crossprod(Q, b))
x1

```

### Q1(a)iii
```{r}
AtA <- crossprod(R)
AtA
```
### Q1(a)iv
```{r}
x2 <- solve(AtA,b)
x2
```
## Q1(b)

### Q1(b)i
```{r}
f <- function(x) sin(3 *pi * x) + sqrt(1-x^2)
N1 <- 5
a <- -0.6
b <- 0.6
I_gq <- function(N, a, b){
  xw <- pracma::gaussLegendre(N, a, b)
  sum(xw$w * f(xw$x))
}
I1 <- I_gq(N1, a, b)
I1
```

### Q1(b)ii
```{r}
I_approx <- function(N, a, b){
  h <- (b - a) / N
  x <- seq(a, b, length.out = N + 1)
  x_mid <- (x[-1] + x[-length(x)]) / 2
  integral <- h * (0.5 * f(a) + sum(f(x_mid)) + 0.5 * f(b))
  return(integral)
}
N2 <- 10
I2 <- I_approx(N2, a, b)
I2
```

### Q1(b)iii
```{r}
I_exact <- 0.5 * (asin(0.6) + 0.6 * sqrt(1 - 0.6^2) - asin(-0.6) - (-0.6) * sqrt(1 - (-0.6)^2))

relative_error_gaussian <- abs(I1 - I_exact) / abs(I_exact)
relative_error_approx <- abs(I2 - I_exact) / abs(I_exact)

relative_error_gaussian
relative_error_approx

```
Gaussian is better as error is smaller

### Q1(b)iv
```{r}
find_min_nodes <- function(f, a, b, target_error, I_exact) {
  N <- 5
  relative_error <- 1
  while (relative_error > target_error) {
    integral <- I_gq(N, a, b)
    relative_error <- abs(integral - I_exact) / abs(I_exact)
    N <- N + 1
  }
  return(N - 1)
}

min_nodes <- find_min_nodes(f, -0.6, 0.6, 1e-4, I_exact)
min_nodes

```

## Q1(c)

```{r}
y <- c(0.5, 0.5, 0.5, 1.5, 1.8, 2.5, 3.2, 3.9, 4.1, 9.6)
```

### Q1(c)i
```{r}
n <- length(y)
w <- sum((y + 2)^3)
f0 <- function(gamma, y){
  n * log(3) + n * log(gamma) + 2 * sum(log(y + 2)) -gamma * w
}
gamma <- 0.004
f0(gamma, y)
```

### Q1(c)ii
```{r}
f1 <- function(gamma, y){
  n / gamma + w
}
f1(gamma, y)
```
### Q1(c)iii
```{r}
f2 <- function(gamma, y){
  matrix(-n / gamma^2)
}
f2(gamma, y)
```
### Q1(c)iv
```{r}
gamma - f1(gamma, y) / f2(gamma, y)
```

### Q1(c)v
```{r}
gamma_t <- n / w
abs((gamma - gamma_t) / gamma_t)
```
# Q2

## Q2(a)
```{r}
c_fun <- function(z, x) {
  exp(-(outer(z, x, FUN = "-")^2))
}

# Function to compute m(x) and v(x)
mvt <- function(x, z, H, a, n) {
  A <- c_fun(z, z)
  L <- chol(A)
  
  HtAinv <- forwardsolve(t(L), H)
  HtAinvH <- crossprod(HtAinv)
  
  beta_hat <- backsolve(chol(HtAinvH), crossprod(HtAinv, forwardsolve(t(L), a)))
  
  residual <- a - H %*% beta_hat
  
  sigma_hat_sq <- as.numeric(crossprod(forwardsolve(t(L), residual))) / (n - 3)
  
  t_x <- c_fun(z, x)
  
  A_inv_residual <- backsolve(L, forwardsolve(L, residual))
  m_x <- beta_hat + crossprod(t_x, A_inv_residual)
  
  A_inv_tx <- backsolve(L, forwardsolve(L, t_x))
  temp <- 1 - crossprod(t_x, A_inv_tx)
  
  H_inv_HtAinv <- backsolve(chol(HtAinvH), diag(ncol(H)))
  HtAinvH_inv <- H_inv_HtAinv %*% t(H_inv_HtAinv)
  
  v_x <- sigma_hat_sq * (temp + (1 - crossprod(t_x, forwardsolve(t(L), H))) %*% HtAinvH_inv %*% (1 - crossprod(t_x, forwardsolve(t(L), H))))
  
  return(c(m_x, v_x))
}

```
## Q2(b)

```{r}
# Given data
n <- 5
a <- c(0.76, 1.57, 1.14, 0.53, 1.03)
z <- c(0.2, 0.6, 1.0, 1.4, 1.8)
H <- matrix(1, nrow = 5, ncol = 1)  
x <- 0.5
result <- mvt(x, z, H, a, n)
result
```

# Q3

## Q3(a)
```{r}
IN <- function(N){
  x <- cos(pi * (1:N) / (N + 1))
  w <- (pi / (N + 1)) * sin(pi * (1:N) / (N + 1))^2
  I_N <- sum(w * cos(x))
  return(I_N)
}
```

## Q3(b)
```{r}
IN(N2)
```
## Q3(c)
```{r}
I_MC <- function(f, a ,b ,N){
 x <- runif(N, a, b)
 (b - a) * mean(f(x))
}
```

## Q3(d)
```{r}
f <- function(x) cos(x) * sqrt(1 - x^2)
N3 <- 1e2
I_MC1 <- I_MC(f, -1, 1, N3) 
I_MC1
```

## Q3(e)
```{r}
fd <- function(f, x, delta) {
  (f(x - 2 * delta) - 8 * f(x - delta) + 8 * f(x + delta) - f(x + 2 * delta)) / (12 * delta)
}

f_prime <- function(x) {
  -sin(x) * sqrt(1 - x^2) - x * cos(x) / sqrt(1 - x^2)
}

delta <- 1e-3
f_approx_prime <- fd(f, 0, delta)
f_exact_prime <- f_prime(0)

relative_error <- abs(f_approx_prime - f_exact_prime) / abs(f_exact_prime)
relative_error
```


# Q4

```{r}
# y
c(2.2, 2.6, 3.3, 3.3, 4.5, 5.6, 6.4, 7.7, 9.9, 10.4)
# (tau_0, omega_0)
c(1.1, 1.7)
```

