# Exercise 7.2
##initialize
library(bootstrap) ##get law data
n <- nrow(law)
B <- 2000
R <- numeric(B)
indices <- matrix(0, nrow = B, ncol = n)

##bootstrap for se(R)
for (b in 1:B) {
    #randomly select the indices
    i <- sample(1:n, size = n, replace = TRUE)
    LSAT <- law$LSAT[i] #i is a vector of indices
    GPA <- law$GPA[i]
    R[b] <- cor(LSAT, GPA)
    #save the indices for the jackknife
    indices[b,] <- i
}

#Jackknife-after-bootstrap
se.jack <- numeric(n)
for (i in 1:n) {
    #in i-th replicate omit all samples with x[i]
    keep <- (1:B)[apply(indices, MARGIN = 1,
                        FUN = function(k) {!any(k == i)})]
    se.jack[i] <- sd(R[keep])
}

print(sd(R))
print(sqrt((n-1) * mean((se.jack - mean(se.jack))^2)))

#Exercise 7.3
##bootstrap t interval method
boot.t.ci <-
    function(x, B = 500, R = 100, level = .95, statistic){
        #compute the bootstrap t CI
        x <- as.matrix(x)
        n <- nrow(x)
        stat <- numeric(B)
        se <- numeric(B)
        
        boot.se <- function(x, R, fun) {
            #local function to compute the bootstrap
            #estimate of standard error for statistic f(x)
            x <- as.matrix(x)
            m <- nrow(x)
            th <- replicate(R, expr = {
                i <- sample(1:m, size = m, replace = TRUE)
                fun(x[i, ])
            })
            return(sd(th))
        }
        
        ##re-sample for 
        for (b in 1:B) {
            j <- sample(1:n, size = n, replace = TRUE)
            y <- x[j, ]
            stat[b] <- statistic(y)
            se[b] <- boot.se(y, R = R, fun = statistic)
        }
        
        ##estimate the distribution for quantile
        stat0 <- statistic(x) #original statistic
        t.stats <- (stat - stat0) / se
        se0 <- sd(stat)
        alpha <- 1 - level
        
        Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
        names(Qt) <- rev(names(Qt))
        CI <- rev(stat0 - Qt * se0)
        return(CI)
    }

##use function boot.t.ci
library(bootstrap) # get law data
stat <- function(x) { cor(x[, 1], x[,2]) } # calculate the correlation
ci <- boot.t.ci(law, statistic = stat, B=2000, R=200)
print(ci) 

#Exercise 7.4
library(boot) # get aircondit data # get boot function
mle <- function(x,i){ 1/mean(x[i]) } 
aircondit.boot <- boot(data = as.matrix(aircondit), statistic = mle, R = 2000)

#Exercise 7.5
library(boot) # get aircondit data and function boot.ci
##a simple method
aircondit0 <- as.matrix(aircondit)
meantime <- function(x,i) mean(x[i])
##get boot class
aircondit.boot <- boot(data = aircondit0, statistic = meantime, R = 5000)
boot.ci(aircondit.boot, type = c("norm", "basic", "perc", "bca"))

#Exercise 7.7
library(bootstrap) # get scor data

# default: perform centralized and standardized process
theta <- function(data, center = TRUE, scale = TRUE) { 
    Sigma.hat <- 
        cov(scale(as.matrix(data),center = center, scale = scale))
    e.hat <- eigen(Sigma.hat) # calculate the eigenvalue and eigen-vector
    
    ## estimate the proportion of the first principal component
    theta.hat <- e.hat$values[1]/sum(e.hat$values)
    
    return(theta.hat)
}

## use all sample
theta.all <- theta(scor)

## use bootstrap method
 # initial
B <- 2000; n <- nrow(scor)
theta.B <- numeric(B)
 # bootstrap
theta.B <- replicate(B, expr = {
    i <- sample(1:n, size = n, replace = TRUE)
    scor.B <- as.matrix(scor)[i,]
    theta.B <- theta(scor.B)
})

## for bias and standard error
bias <- mean(theta.B-theta.all)
std <- sd(theta.B)
result <- c(theta.all,bias,std)
names(result) <- c("original","bias.B", "std.B")
print(result)
#Exercise 7.8

##use jackknife method
 # initial
n <- nrow(scor)
theta.jack <- numeric(n)
 # jackknife
for (i in 1:n) 
    theta.jack[i] <- theta(as.matrix(scor)[-i,])

## for bias and standard error
bias.jack <- (n-1)*mean(theta.jack-theta.all)
std.jack <- sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))

result.jack <- c(theta.all,bias.jack,std.jack)
names(result.jack) <- c("original","bias.jack", "std.jack")
print(result.jack)

#Exercise 7.10

library(DAAG); attach(ironslag) # get ironslag data

# for n-fold cross validation
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
# fit models on leave-one-out samples
for (k in 1:n) {
    y <- magnetic[-k]
    x <- chemical[-k]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
    e1[k] <- magnetic[k] - yhat1
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
        J2$coef[3] * chemical[k]^2
    e2[k] <- magnetic[k] - yhat2
    
    J3 <- lm(y ~ x + I(x^2) + I(x^3))
    yhat3 <- J3$coef[1] + J3$coef[2] * chemical[k] +
        J3$coef[3] * chemical[k]^2 + J3$coef[4] * chemical[k]^3
    e3[k] <- magnetic[k] - yhat3
    
    J4 <- lm(log(y) ~ x)
    logyhat4 <- J4$coef[1] + J4$coef[2] * chemical[k]
    yhat4 <- exp(logyhat4)
    e4[k] <- magnetic[k] - yhat4
}

crossv <- c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
best.cross.ind <- which(crossv == min(crossv), arr.ind = T)

# fit models on adjusted quadratic R
L1 <- summary(lm(magnetic ~ chemical))
L2 <- summary(lm(magnetic ~ chemical + I(chemical^2)))
L3 <- summary(lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3)))
L4 <- summary(lm(log(magnetic) ~ chemical))

adj.r <- 
    c(L1$adj.r.squared,L2$adj.r.squared,L3$adj.r.squared,L4$adj.r.squared)
best.adj.ind <- which(adj.r == max(adj.r),arr.ind = T)

best.ind <- c(best.cross.ind, best.adj.ind)
names(best.ind) <- c('best.cross.ind', 'best.adj.ind')
print(best.ind)