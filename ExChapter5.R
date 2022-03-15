#Exercise5.3
m <- 10000
set.seed(0315)
##简单Monte Carlo估计
x1 <- runif(m,0,0.5)
g1 <- exp(-x1)
th.hat1 <- 0.5*mean(g1)   ### theta.hat估计
vg1 <- mean((g1-mean(g1))^2)/m   ### theta.hat方差
##从指数分布中抽样
x2 <- rexp(m)
g2 <- (x2 <= 0.5)
th.hat2 <- mean(g2)   ### theta.star估计
vg2 <- mean((g2-mean(g2))^2)/m   ### theta.star方差
##比较两种方法方差
vg1/vg2
###g1由均匀分布衍生,g2实际服从二项分布

#Exercise5.6
m <- 10000
##对偶方法计算
x <- runif(m/2)
u <- exp(x)
v <- exp(1-x)
result1 <- cov(u,v)
result2 <- var((u+v))/2
##简单蒙特卡洛方法
y <- runif(m)
result3 <- var(exp(y))
cat(100*(result3-result2)/result3,"%") ###计算减少的比例

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
