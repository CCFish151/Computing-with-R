---
title： 统计计算第三次作业
author: 2019302100125 赵晨宇
---

# 统计计算第三次作业

2019302100125 赵晨宇



## Exercise 6.3

### 思路分析

检验的功效可以定义为：
$$
\pi(\theta)=P_{\theta}(Reject\; H_0)
$$
一型错误率通过显著性水平$\alpha$的选择控制，较高的功效则意味着较低的二型错误率。因此，在对两个具有相同显著性水平的检验进行比较时，我们感兴趣的是比较他们的功效。

检验的经验功效可以通过Monte Carlo算法来计算：

1. 利用Monte Carlo算法计算t检验的功效
2. 绘制功效曲线，并分析结果



### 代码实现

首先，编写计算t-test功效的函数` Power.t`，其实现如下（或直接利用`Power.t.test`）

```R
Power.t = function(n, mu0 = 500, sigma = 100, mu) { 
    ##计算t检验功效,n代表样本数量
    m <- 1000
    M <- length(mu)
    power <- numeric(M)
    
    for (i in 1:M) {
        mu1 <- mu[i]
        ##计算检验的p值并与显著性水平比较
        pvalues <- replicate(m, expr = {
            x <- rnorm(n, mean = mu1, sd = sigma)
            ttest <- t.test(x, alternative = "greater", mu = mu0)
            ttest$p.value})
        power[i] <- mean(pvalues <= .05)
    }
    return(power)
}
```

之后，调用函数`Power.t`计算不同样本量大小下的功效：

```R
n <- c(seq(10, 50, 10))
mu <- c(seq(450,650,10))
pw <- matrix(0, length(n), length(mu))
##调用函数计算功效
for (i in 1:5) {
    pw[i,] <- Power.t(n[i], mu = mu)}
```

最后，利用R中的画图功能绘制功效曲线，其代码如下：

```R
##绘制功效曲线
limit <- max(pw)
color <- c("black","red","skyblue","green","blue")
type <- c(1,1,2,3,4)
plot(mu, numeric(length(mu)), ylab = "Power" ,ylim = c(0, limit), type = "n")

for (i in 1:5) {
    lines(mu, pw[i,], col = color[i], lty = type[i])
}
##添加图例
legend("topleft",legend = c("n=10", "n=20", "n=30", "n=40", "n=50"),
       lty = type, col = color)
```

### 结果解释

代码运行绘制图象如下显示，其中`n`代表样本数量。

![Power Curves](D:\统计计算\RPractice\作业3\image\PowerCurve .jpeg)

由图像可以看出，随着样本数量的增加，检验的功效随之增加；同时随着`mu`的增加，功效也增加。

可以得出结论，样本数量越多，偏离$ H_0$程度越大，t检验的经验功效越大，即把握度也就越大。



## Exercise 6.4

### 思路分析

由于置信水平即为置信区间能够覆盖 $\theta$ 真值的概率，也即接受原假设的概率，因此，经验的置信水平可通过如下 $Monte\; Carlo$ 算法来实现：

1. 产生第 $j$ 个随机样本 $x_1^{(j)},\;x_2^{(j)},\;...,\;x_n^{(j)}$
2. 计算基于第 $j$ 个样本的置信区间 $ C_j$
3. 计算 $ y_j = I(\theta \in C_j) $
4. 计算经验的置信水平 $\bar y = \frac {1} {m} \Sigma _ {j=1} ^ {m} y_j$



### 代码实现

代码如下

```R
n <- 20
m <- 10000
alpha <- .05
##利用Monte Carlo方法估计置信水平
ucls <- replicate(m, expr = {
    x <- rlnorm(n, mean=1, sd=2) ##生成对数正态分布
    ucls[i] <- (mean(log(x))-1)/sqrt(var(log(x))/n)#构建枢轴量
})
##计算经验置信水平
cov.rate<-cumsum(ucls>qt(alpha, df = n-1))/1:m
#绘制图象
plot(2:m, cov.rate[-1], type = "n", 
     xlab = "样本数", ylab = "置信水平估计", ylim = c(0.8,1.0))
abline(h=0.95, col = "gray")
lines(2:m, cov.rate[-1], lwd = 1.5)
```



### 结果解释

运行结果如下：

![Conf Level](D:\统计计算\RPractice\作业3\image\ConfLevel.jpeg)

图像显示，随着样本数量逐渐增大，经验置信水平趋于稳定，并与名义置信水平保持一致。



## Exercise 6.5

### 思路分析

针对服从 $ \chi ^2(2)$ 分布的随机变量 $X$ ，进行方差和 $ t$ 的置信区间估计，并对经验置信水平进行比较。

方差估计的置信上界为
$$
\frac {(n-1)S^2} {\chi ^2 _{\alpha} (n-1)}
$$
$ t-interval$ 的置信上界为
$$
\frac {\bar x -2} {t_{\alpha}(n-1)}
$$
之后利用 $Monte \; Carlo$ 方法对经验置信水平作出估计。 

### 代码实现

```R
##calculate experience conf.level
expCL = function(n, m, alpha, var = TRUE) {
    #两种不同检验方法的封装
    CI.var <- function(n,alpha) {
        x <- rchisq(n, df = 2)
        return((n-1) * var(x)/qchisq(alpha, df = n-1) > 4)}
    CI.t <- function(n,alpha) {
        x <- rchisq(n, df = 2)
        return((mean(x)-2)/qt(alpha, df = n-1) > 2)
    }
    
    if (!var) UCL<-replicate(m,expr=CI.t(n=20,alpha=.05))
    else  UCL<-replicate(m,expr=CI.var(n=20,alpha=.05))
    
    return(mean(UCL))
}
##调用函数对结果进行计算
ucls.var <- expCL(100,100000,.05)
ucls.t <- expCL(100,100000,.05,var = FALSE)
```



### 结果解释

结果显示，

> `ucls.var = 0.7827`
>
> `ucls.t = 0.89212`

比较可知，此时利用 $t-interval$ 估计的置信区间更加稳健。

## Exercise 6.8

### 思路分析

针对简单两样本有“Count Five”和F-检验的方法，本题即为针对简单正态样本在样本数量不同的情况下，比较两种检验方法的功效。



**Count Five**方法：

针对两组均值和样本量均相同的样本，通过计算一组样本相对于另一个样本极值点的个数，即不在样本范围内的观测值的个数，若每组样本相对另一组都至少有5个极值点，则拒绝等方差的假设。



**F-test**方法：*仅适用于正态总体*

对于正态总体的样本，在方差未知的情况下，我们可以构建检验统计量
$$
F = \frac {S_1^2} {S_2^2}
$$
此时由于检验统计量 $F \sim F(n-1,n-1)$ ，因此可求出其拒绝域，并判断是否接受原假设。



最终，利用 $Monte Carlo$ 方法估计两种方法在不同样本量下的经验功效，并进行比较。

### 代码实现

```R
## var2test提供了两样本方差检验的两种方法
## m:样本数；sigma:样本标准差；
## count5表示默认使用Count Five方法；alpha代表显著性水平
var2test = function(m, sigma, count5 = TRUE, alpha = 0.055) {
    ##count five method
    count5test <- function(x, y) {
        X <- x - mean(x)
        Y <- y - mean(y)
        outx <- sum(X > max(Y)) + sum(X < min(Y))
        outy <- sum(Y > max(X)) + sum(Y < min(X))
        # return 1 (reject) or 0 (do not reject H0)
        return(as.integer(max(c(outx, outy)) > 5))
    }
    
    power <- mean(replicate(m, expr={
        x <- rnorm(20, 0, sigma[1])
        y <- rnorm(20, 0, sigma[2])
        # use F-test
        if (!count5) as.integer(var.test(x, y)$p.value<alpha)
        # use Count Five
        else count5test(x, y)
    }))
    
    return(power)}

m <- c(20, 1000, 100000)
sigma <- c(1, 1.5) 
powerm <- matrix(0, 2, length(m))
##调用函数计算并对功效进行比较
for (i in 1:3) {
    powerm[1,i] <- var2test(m[i], sigma)
    powerm[2,i] <- var2test(m[i], sigma, count5 = FALSE)
}

power <- data.frame(powerm, row.names = c("count5","F-test"))
names(power) <- c("small", "medium", "large")
power
```



### 结果解释

取小样本 `n=20`，中样本 `n=1000`，大样本 `n=100000`进行方差齐性检验，得到结果如下表所示

|  ***Power***   | Small | Medium |  Large  |
| :------------: | :---: | :----: | :-----: |
| **Count Five** | 0.30  | 0.326  | 0.31326 |
|   **F-test**   | 0.15  | 0.412  | 0.41823 |

由表格可以看出，**F-test**方法相对**Count Five**方法在样本数量较多的情况下具有较高的功效，而小样本下**Count Five**方法功效更高可能是由于样本数量过少，对**F-test**产生了较大的影响。



## Exercise 6.9

### 思路分析

根据题目中给出条件，可知在次序统计量条件下，基尼系数的估计方法为：
$$
\hat G =\frac{ 1}{n^2\mu} \Sigma _{i = 1}^{n}(2i-n-1)x_{(i)}
$$
因此，生成服从 $logNormal,\; Uniform,\;Bernoulli $ 分布的随机数，排序后计算。

### 代码实现

```R
library(purrr) #get rbernoulli()
## distribution = lnorm|unif|bernoulli
Gini.cal = function(num.sample, num.repeat, distribution = "unif") {
    gini.sam <- numeric(num.repeat)
    
    gini.sam <- replicate(num.repeat, expr = {
        sample <- switch(distribution,
                         lnorm = rlnorm(num.sample, 0, 1),
                         unif = runif(num.sample),
                         bernoulli = as.integer(
                             rbernoulli(num.sample,0.1)),
                         warning("can't find such distribution")) 
                         #增加分布类型或调整参数请修改这里
        order.sam <- sort(sample)#获得次序统计量
        order <- rep(1:num.sample)
        
        tem <- (2*order-num.sample-1)*order.sam
        sum(tem)/(num.sample^2 * mean(sample))
    })
    
    gini.order <- sort(gini.sam)
    #绘制直方图
    hist(gini.order, main = paste("Histogram", distribution))
    
    result1 <- mean(c(gini.order[round(num.repeat/10)],
                      gini.order[round(num.repeat/10)+1]))
    result2 <- median(gini.order) ##median
    result3 <- mean(gini.sam, na.rm = TRUE)  ##mean
    
    return(list(mean = result3, median = result2, deciles = result1))
}

dis <- c("lnorm","unif","bernoulli")
##按照题目要求计算各分布基尼系数
Gini.hat <- as.data.frame(mapply(Gini.cal, 20, 10000, distribution = dis),
                          row.name = c("mean","median","deciles"))
names(Gini.hat) <- dis
```



### 结果解释

得到各数据如下：

|         | lNorm | Unif  | Bernoulli |
| ------- | ----- | ----- | --------- |
| mean    | 0.481 | 0.322 | 0.886     |
| median  | 0.476 | 0.321 | 0.900     |
| deciles | 0.383 | 0.254 | 0.800     |

$ \hat G$ 在不同样本分布下的密度曲线如下：

<img src="D:\统计计算\RPractice\作业3\image\lnorm.jpeg" alt="lnorm"  />

![unif](D:\统计计算\RPractice\作业3\image\unif.jpeg)

![](D:\统计计算\RPractice\作业3\image\bernoulli.jpeg)



## Project 6.A

### 思路分析

利用 $Monte Carlo$ 算法计算经验一型错误率即计算显著的检验比例 $ \frac {1}{m}\Sigma _{j=1}^{m} I_j$ ，其中 $ I_j = 1$ ，若 $H_0$ 在显著性水平 $\alpha$ 下被拒绝，否则 $I_j = 0$ 。

之后，将一型错误率与名义置信水平 $ \alpha$ 作比较。

### 代码实现

```R
emIerror = function(n,m,alpha,population) {
    ##计算经验I型错误率
    p <- numeric(m) #storage for p-values
    p <- replicate(m, expr = {
        sample <- switch(population,
                         chisq = rchisq(n, df = 1),
                         unif = runif(n,0,2),
                         exp = rexp(n),
                         warning("can't find such population"))
        ttest <- t.test(sample, alternative = "two.sided", mu = 1)
        p <- ttest$p.value
    })

    p.hat <- mean(p < alpha)
    return(p.hat)
}

#example
emIerror(100,1000000,.05,"unif")
```



### 结果解释

当样本服从 $ \chi ^2(1)$ 时，其经验一型错误率为`Ichisq = 0.0654` ;

当样本服从 $\epsilon (1)$ 时，其经验一型错误率为`Iexp = 0.0578` ;

当样本服从 $U(0,1)$ 时，其经验一型错误率为`Iunif = 0.0502` ;

可以看到，均匀分布总体下，t检验的经验一型错误率与名义置信水平较为接近，其他均要明显大于名义置信水平。

## Project 6.B

### 思路分析

#### 问题1

针对二元正态向量，$ \left( \begin{array} { l } { X } \\ { Y } \end{array} \right) \sim N_2 (\left( \begin{array} { l } { 2 } \\ { 1 } \end{array} \right),\:\left( \begin{array} { l l} { 1 }&{0.5} \\ {0.5}&{ 1 } \end{array} \right) )$ ，显然此时 $X,\;Y$ 具有相关性。

此时利用`cor.test`进行检验并计算功效。

#### 问题2

若取服从标准柯西分布随机变量 $X$ ， $Z_1 \sim N(-10,2)$ ， $Z_2 \sim N(4,2)$ ， $Z_3 \sim N(20,1)$ 

$Y = 0.001\;X+Z_1 + Z_2+Z_3$ ，此时显然 $X, \;Y$ 具有相关性，因此通过`cor.test `模拟进行检验并计算功效。

### 代码实现

首先调用`MASS`软件包

```R
library(MASS)
```



针对问题一，可以使用如下代码进行模拟

```R
##Question One
mean <- c(2,1)
sigma <- matrix(c(1,0.5,0.5,1),nrow = 2,ncol = 2)
n <- 20
m <- 10000
alpha <- .05
##生成符合条件随机数
munorm = function(n, mu, Sigma)
{
    d <- length(mu)
    Q <- chol(Sigma) 
    Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
    X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
    return(X)
}
#二元正态随机变量相关性检验功效
corPower1 = function(m, n, mean, sigma, method) {
    power <- mean(replicate(m, expr = {
        sample <- munorm(n, mean, sigma) #生成二元正态分布
        x <- sample[,1]
        y <- sample[,2]
        result <- cor.test(x,y, alternative = "two.sided",
                           method = method, conf.level = 0.95)
        result$p.value <= alpha
    }))
    return(power)
}

pearson.result <- corPower(m, n, mean, sigma, "pearson")
spearman.result <- corPower(m, n, mean, sigma, "spearman")
kendall.result <- corPower(m, n, mean, sigma, "kendall")
```



针对问题2，可用如下代码实现：

```R
##Question two
corPower2 = function(m, n, mean, sigma, method) {
    power <- mean(replicate(m, expr = {
        x <- rcauchy(n)
        z1 <- rnorm(n, -10, sqrt(2))
        z2 <- rnorm(n, 4, sqrt(2))
        z3 <- rnorm(n, 20, 1)
        y <- 0.001*x+z1+z2+z3
        result <- cor.test(x,y, alternative = "two.sided",
                           method = method, conf.level = 0.95)
        result$p.value <= alpha
    }))
    return(power)
}

pearson.result2 <- corPower2(m, n, mean, sigma, "pearson")
spearman.result2 <- corPower2(m, n, mean, sigma, "spearman")
kendall.result2 <- corPower2(m, n, mean, sigma, "kendall")
```



### 结果解释

#### 问题1

代码运行结果为 

```R
> pearson.result1 = 0.6438
> spearman.result1 = 0.5656
> kendall.result1 = 0.5535
```

#### 问题2

代码运行结果为

```R
> pearson.result2 = 0.0569
> spearman.result2 =0.0533
> kendall.result2 = 0.0470
```



综合上述结果，可知非参数检验的经验功效要高于相关性检验的功效。





