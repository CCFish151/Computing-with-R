# 3.3
n <- 1200
u <- runif(n) #生成均匀分布随机数
a <- 2
b <- 2
x <- b*(1-u)^(-1/a) #利用反函数生成特定分布随机数
hist(x , prob=TRUE , main=expression(F(x)==1-(b/x)^a))
y<-seq(2,80,0.01)
z <- 8/(y^3)
lines(y, z) #直接绘制密度曲线

#3.4
#生成Rayleigh分布随机数
rRayleigh = function(sigma, n=1000)
{
  u <- runif(n)
  x <- (-2*(sigma)^2*log(u))^0.5
  hist(x, prob=TRUE)
}
#选取一些绘制图形
rRayleigh(1)
rRayleigh(0.5)
rRayleigh(2)

#3.9
# 生成算法函数
selectf = function(n) 
{
  u <- rep(0,n)
  for(i in 1:n)
  {
    u1 <- runif(1,-1,1)
    u2 <- runif(1,-1,1)
    u3 <- runif(1,-1,1)
    if(abs(u3) >= abs(u2) && abs(u3) >= abs(u1))
    {
      u[i] <- u2
    }
    else
    {
      u[i] <- u3
    }
  }
  return(u)
}
x <- selectf(10000) #抽样
hist(x , prob = TRUE) #绘制直方图
#直接绘制密度函数
y <- seq(-1,1,0.001) 
z <- 0.75*(1-y^2)
lines(y,z) #输出对比


#3.11
#生成混合分布随机数
GenerateMixNorm = function(n,p) 
{
  u1 <- rnorm(n,0,1) # N(0,1)
  u2 <- rnorm(n,3,1) # N(3,1)
  x <- rep(0,n)
  #利用sample抽样
  k <- sample(1:2,size = n, replace = TRUE, prob = c(p, 1-p)) 
  #通过循环选取随机数
  for(i in 1:n) 
  {
    if(k[i] == 1){x[i] = u1[i]}
    else{x[i] = u2[i]}
  }
  return(x)
}
#生成混合高斯分布 绘制直方图 直观检验
x <- GenerateMixNorm(100000,0.75) 
hist(x , prob = TRUE)  
# 以步长0.1寻找双峰图像
for(i in seq(0,1,0.1))
{
  x <- GenerateMixNorm(100000,i) 
  hist(x , prob = TRUE)
}


#3.13
n <- 1000
u <- runif(n)
x <- (2/(1-u)^(1/4))-2 
hist(x, prob = TRUE, main = expression(f(x)==64/(x+2)^5)) 
y <- seq(0, 20, 0.01)
lines(y, 64/(y+2)^5)  

#3.14
mv <- c(0,1,2)
cov <- matrix(c(1,-.5,.5,-.5,1,-.5,.5,-.5,1),3,3)
##Choleski分解方法生成多元随机变量
rmvCholeski = function(n,mv,Sigma)
{
  d <- length(mv)
  Q <- chol(Sigma) 
  Z <- matrix(rnorm(n*d), nrow=n, ncol=d)
  X <- Z %*% Q + matrix(mv, n, d, byrow=TRUE)
  X
}

X <- rmvCholeski(1000,mv,cov)
##pairs and plot
pairs(X)
plot(X[,1:2], xlab = "x1", ylab = "x2", pch = 16)
plot(X[,1:3], xlab = "x1", ylab = "x3", pch = 16)
plot(X[,2:3], xlab = "x2", ylab = "x3", pch = 16)

#3.16
##获取数据集
library(bootstrap) 
##将数据标准化并计算协方差矩阵
cov(scale(scor[,1:2]))
cov(scale(scor[,3:5]))
cov(scale(scor))

#3.20
##创建生成函数
simu.cpoigamma = function(lambda,shape,scale,size,t)
  #lambda:齐次Poisson过程速率
  #shape:Gamma分布形状参数;scale:Gamma分布尺度参数
{
  ##生成齐次poisson过程
  rpoisp = function(t0,upper=1000)
  {
    Tn <- rexp(upper,lambda) 
    ## upper表示模拟步数，过少可能无法得到结果
    Sn <- cumsum(Tn)
    n <- min(which(Sn > t0))-1
    return(n)
  }
  pprocess <- replicate(size,rpoisp(t))
  ##生成题目要求的X(t)
  X <- sapply(pprocess,function(n)
  {
    vY <- rgamma(n = 100, shape = shape, scale = scale) 
    #生成独立同分布Gamma随机变量
    sum(vY[1:n])
  })
  ##比较模拟与理论差异
  ##模拟值
  msim <- mean(X)
  vsim <- var(X)
  ##理论值
  mthe <- lambda*t*shape*scale
  vthe <- lambda*t*(shape+1)*shape*scale^2
  show = data.frame(msim,vsim,mthe,vthe)
  return(show)
}
##调整参数进行测试
simu.cpoigamma(2,1,2,10000,10)
simu.cpoigamma(3,5,4,10000,10)
