#EX6.3
Power.t = function(n, mu0 = 500, sigma = 1000, mu) { 
    ##计算t检验功效,n代表样本数量
    m <- 1000
    M <- length(mu)
    power <- numeric(M)
    
    for (i in 1:M) {
        mu1 <- mu[i]
        ##计算检验的p值并与显著性水平比较
        pvalues <- replicate(m,expr = {
            x <- rnorm(n, mean = mu1, sd = sigma)
            ttest <- t.test(x, alternative = "greater", mu = mu0)
            ttest$p.value})
        power[i] <- mean(pvalues <= .05)
    }
    return(power)
}
##调用函数计算功效
n <- c(seq(10, 50, 10))
mu <- c(seq(450,650,10))
pw <- matrix(0, length(n), length(mu))

for (i in 1:5) {
    pw[i,] <- Power.t(n[i], mu = mu)
}
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


#6.4
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

#6.5
##calculate experience conf.level
expCL = function(n, m, alpha, var = TRUE) {
    
    CI.var <- function(n,alpha) {
        x <- rchisq(n, df = 2)
        return((n-1) * var(x)/4 > qchisq(alpha, df = n-1))}
    CI.t <- function(n,alpha) {
        x <- rchisq(n, df = 2)
        return((mean(x)-2)/sqrt(var(x)/n) > qt(alpha, df = n-1))
    }
    
    if (!var) UCL<-replicate(m,expr=CI.t(n=20,alpha=.05))
    else  UCL<-replicate(m,expr=CI.var(n=20,alpha=.05))
    
    return(mean(UCL))
}

ucls.var <- expCL(100,100000,.05)
ucls.t <- expCL(100,100000,.05,var = FALSE)

#6.8
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

#6.9 Gini Ratio calculate
library(purrr) #get rbernoulli
## distribution = lnorm|unif|bernoulli
Gini.cal = function(num.sample, num.repeat, distribution = "unif") {
    gini.sam <- numeric(num.repeat)
    
    gini.sam <- replicate(num.repeat, expr = {
        sample <- switch(distribution,
                         lnorm = rlnorm(num.sample, 0, 1),
                         unif = runif(num.sample),
                         bernoulli = as.integer(rbernoulli(num.sample, 0.1)),
                         warning("can't find such distribution")) #增加分布类型或调整参数请修改这里
        order.sam <- sort(sample)#获得次序统计量
        order <- rep(1:num.sample)
        
        tem <- (2*order-num.sample-1)*order.sam
        sum(tem)/(num.sample^2 * mean(sample))
    })
    
    gini.order <- sort(gini.sam)
    
    result1 <- mean(c(gini.order[round(num.repeat/10)],gini.order[round(num.repeat/10)+1]))
    result2 <- median(gini.order) ##median
    result3 <- mean(gini.sam, na.rm = TRUE)  ##mean
    
    return(list(mean = result3, median = result2, deciles = result1))
}

dis <- c("lnorm","unif","bernoulli")
##按照题目要求计算各分布基尼系数
Gini.hat <- as.data.frame(mapply(Gini.cal, 20, 10000, distribution = dis, USE.NAMES = TRUE),
                          row.name = c("mean","median","deciles"))
names(Gini.hat) <- dis

#6.A
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
Iunif <- emIerror(100,1000000,.05,"unif")
Ichisq <- emIerror(100,1000000,.05,"chisq")
Iexp <- emIerror(100,1000000,.05,"exp")


#6.B
library(MASS)
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

pearson.result1 <- corPower1(m, n, mean, sigma, "pearson")
spearman.result1 <- corPower1(m, n, mean, sigma, "spearman")
kendall.result1 <- corPower1(m, n, mean, sigma, "kendall")

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




