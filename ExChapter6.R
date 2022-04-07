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

#6.5##?``
##Example6.4中置信上界为7.224
n <- 20
m <- 1000000
alpha <- .05
x <- rchisq(n, df = 2)
UCL = mean(replicate(m,expr = (n-1)*var(x)/qchisq(alpha, df = n-1)))


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
## care = mean(default)|median|deciles
## 添加可选分布参数，将随机嵌入函数中
Gini.cal = function(num.sample, num.repeat, distribution, care = "mean") {
    gini.sam <- numeric(num.repeat)
    m <- length(sample)
    
    order.sam <- sort(sample)#获得次序统计量
    order <- rep(1:m)
    gini.sam <- replicate(num.repeat, expr = {
        tem <- (2*order-m-1)*order.sam
        sum(tem)/(m^2 * mean(sample))
    })
    
    gini.order <- sort(gini.sam)
    
    if (care == "deciles") result <- mean(c(gini.order[round(num.repeat/10)],gini.order[round(num.repeat/10)+1]))
    else if (care == "median") result <- median(gini.order)
    else result <- mean(gini.sam)
    
    return(result)
}









