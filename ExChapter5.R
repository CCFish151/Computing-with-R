#Exercise5.3
m <- 10000
set.seed(0315)
## Simple Monte Carlo
x1 <- runif(m,0,0.5)
g1 <- exp(-x1)
th.hat1 <- 0.5*mean(g1)   #estimate
vg1 <- mean((g1-mean(g1))^2)/m   #variance
## hit-or-miss
x2 <- rexp(m)
g2 <- (x2 <= 0.5)
th.hat2 <- mean(g2)   #estimate
vg2 <- mean((g2-mean(g2))^2)/m   #variance
##比较两种方法方差
(vg2-vg1)/vg2
###g1由均匀分布衍生,g2实际服从二项分布

#Exercise5.6
m <- 10000
set.seed(0315)
##对偶方法计算
x <- runif(m/2)
u <- exp(x)
v <- exp(1-x)
result1 <- cov(u,v)
result2 <- var((u+v))
##简单蒙特卡洛方法
y <- runif(m)
result3 <- var(exp(y))
##计算减少的比例
cat(100*(result3-result2/2)/result3,"%") 

#Exercise5.7 
##MC法估计函数，包含简单和对偶两种方法
MC.I = function(m = 10000 ,antithetic = TRUE)
{
  u <- runif(m/2)   ###先生成一半随机数
  if(!antithetic) v <- runif(m/2)  ###不采用对偶方法独立生成另一半
  else v <- 1-u  ###采用对偶方法时对偶产生另一半
  u <- c(u,v)
  g <- exp(u)
  cdf <- mean(g)
  return(cdf)
}
##调用函数并进行方差比较
MC1 <- MC2 <-  numeric(1000)
for(i in 1:1000)
{
  MC1[i] <- MC.I(antithetic = FALSE)
  MC2[i] <- MC.I()
}
cat("方差减少了",100*(var(MC1)-var(MC2))/var(MC1),"%")

#Exercise5.11
## 从方程出发推导公式

#Exercise5.13
library(bayesmeta)
g = function(x){ (x^2/sqrt(2*pi))*exp(-x^2/2) }
## 通过下面的步骤确定分布具体参数
xr <- seq(0,10,0.1)
### 实际可求出当x=sqrt(2)为最大值点，猜测取mean = scale = 1.4
for(i in seq(1,2,0.1) )
{
  gr <- g(xr)
  f1 <- dnorm(xr,mean = i)
  f2 <- drayleigh(xr, scale= i )
  
  lim = max(c(gr , f1 , f2))
  plot(xr, gr , type= "l" ,ylim = c(0, lim))
  lines(xr, f1, col="red")
  lines(xr, f2, col="blue")
}
###根据图像可得 mean = scale = 1.4时效果最好
### f1(x) = dnorm(x,mean = 1.4) 
### f2(x) = drayleigh(x , scale = 1.4)
image_Ratio = function(i)
{
  gr <- g(xr)
  f1 <- dnorm(xr,mean = i)
  f2 <- drayleigh(xr, scale= i )

  plot(xr, gr/f1,type = "l" ,ylab = "Ratio", col="red")
  lines(xr, gr/f2, col="blue")
}
image_Ratio(1.4)
legend("topright",legend = c("g/f1" ," g/f2"),lty = c(1,1),col = c("red","blue"))
### f2更为"靠近"g，因此猜测其方差相对更小

#5.14
##进行重要抽样计算
theta.hat <- se <- numeric(3)
m <- 10000
### 由于是无穷积分，做对称处理
x <- runif(m) #using g
fg <- g(x)
theta.hat[1] <- 0.5-mean(fg)
se[1] <- sd(fg)

x <- rnorm(m,mean = 1.4) #using f1
fg <- (g(x)/(dnorm(x,mean = 1.4)))*(x>1)
theta.hat[2] <- mean(fg , na.rm = TRUE)
se[2] <- sd(fg)

x <- rrayleigh(m , scale = 1.4) #using f2
fg <- (g(x)/(drayleigh(x,scale = 1.4)))*(x>1)
theta.hat[3] <- mean(fg , na.rm = TRUE)
se[3] <- sd(fg)

#5.15
M <- 1e4
k <- 5 #分成五个区间
r <- M/k 
N <- 50 #估计重复的次数
T2 <- numeric(k) #储存各区间抽样的结果
est <- matrix(0, N, 2)
#g(x)/f(x)=(1-exp(-1))/(1+x^2)
g<-function(t)(1-exp(-1))/(1+t^2)*(t>0)*(t<1) 
for (i in 1:N) 
{
  est[i, 1] <- mean(g(runif(M)))
  for(j in 1:k){T2[j]<-mean(g(runif(r,(j-1)/k,j/k)))} #分层
  est[i, 2] <- mean(T2)
}
apply(est,2,mean)
apply(est,2,var)
#example 5.10 result:6.504485e-08