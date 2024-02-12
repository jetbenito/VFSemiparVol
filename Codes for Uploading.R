setwd(...)
#Load packages----
library(TSA)
library(fGarch)
library(mgcv)
library(beepr)
library(forecast)
#Create Evaluation functions:----
mse	<- function(x) {mean(x^2)*100}
rmse	<- function(x) {sqrt(mse(x))}
mad	<- function(x,y) {mean(abs(x-y))}
mdpe <- function(x,y) {100*median(abs((x-y)/x))}
ape <- function(x,y) {100*abs((x-y)/x)}
#Presetting----
x0	= 100			#Initial Price
mu	= 0.2			#Drift: Expected Annual Return
sigma	= 0.30	#Volatility: Expected Annual Volatility
psi	= 0.5			#Postulated effect of Change in Trading Volume
m.dly	= mu/252		#Daily Return
s.dly	= sigma/sqrt(252)	#Daily Volatility
#In between simulation settings, value of t is changed between 205, 510, and 1530
t = 0510 #Number of days
#Vector of Evaluation Statistics----
m1.mse<- c()	#Vector of the MSEs for Algo 1
m2.mse<- c()	#Vector of the MSEs for Algo 2
b1.mse<- c()	#Vector of the MSEs for the benchmark:GARCH
b2.mse<- c()	#Vector of the MSEs for the benchmark:GJR
m1.mad<- c()	#Vector of the MADs for Algo 1
m2.mad<- c()	#Vector of the MADs for Algo 2
b1.mad<- c()	#Vector of the MADs for the benchmark:GARCH
b2.mad<- c()	#Vector of the MADs for the benchmark:GJR
m1.rse<- c()	#Vector of the RMSEs for Algo 1
m2.rse<- c()	#Vector of the RMSEs for Algo 2
b1.rse<- c()	#Vector of the RMSEs for the benchmark:GARCH
b2.rse<- c()	#Vector of the RMSEs for the benchmark:GJR
f1.mpe<- c()  #Vector of the MdAPE for the out-of-sample forecasts from Algo 1
f2.mpe<- c()  #Vector of the MdAPE for the out-of-sample forecasts from Algo 2
b.mpe<- c() #Vector of the MdAPE for the out-of-sample forecasts from benchmark
J1 <- c() #Number of iterations in Algo 1
J2 <- c() #Number of iterations in Algo 2
fJ1 <- c() #Number of forecasts iterations in Algo 1
fJ2 <- c() #Number of forecasts iterations in Algo 2
#Replication Procedure----
start.time <- Sys.time()
  for(i in 1:100){
	  set.seed(i)
	  #DGP----
		  #Phase 1:
		  #Random Generation of s and u
      #In between simulations, s will be switched in between Independent and Autocorrelated settings
		  s	<- rnorm(t+20, m.dly, s.dly)		#Independent HF Data
		  # <- arima.sim(list(order=c(1,0,0),ar=0.5),rand.gen=rnorm,mean=m.dly,sd=s.dly,n=t+20) #Autoregressive HF Data
		  u	<- runif(t + 219,-1,1)			#For Volume
		  #Compute for daily volume and daily price
		  u1	<- c(0, cumsum(u))			#Percentage Change in Volume
		  w	<- abs(1000000*u1)			#Daily Trading Volume
		  lx	<- log(x0) + cumsum(s)			#Daily Return
		  #In between simulations, x will be switched in between the Linear and Exponential functional form settings
		  x	<- exp(lx) - (psi*u1[-c(1:200)])	#Linear Generation
		  # <- exp(lx - psi*u1[-c(1:200)])  #Exponential Generation
		  #Phase 2:
		  #Compile the Weekly Data
		  X	<- matrix(x,ncol=5,byrow=TRUE)
		  colnames(X) <- c("Mo","Tu","We","Th","Fr")
		  X	<- data.frame(X)				#Daily Price, compiled by week
		  #Volume
		  w.dly	<- w[-c(1:200)]
		  W	<- matrix(w.dly,ncol=5,byrow=TRUE)
		  V	<- rowSums(W)				#Weekly Trading Volume
		  #Compute the input and output variables
		  l <- t/5
		  Y	<- diff(log(X$Fr))
		  y <- Y[c(1:l-1)]			#C.C. Returns
		  v0	<- diff(log(V))				
		  v <- v0[c(1:l-1)] #Percentage Change in Weekly Volume
		  X.dly	<- X[c(1:l-1),]
		  Mo <- as.vector(X.dly$Mo)
		  Tu <- as.vector(X.dly$Tu)
		  We <- as.vector(X.dly$We)
		  Th <- as.vector(X.dly$Th)
		  Fr <- as.vector(X.dly$Fr)
		  #Out of sample Data
		  y.out <- tail(Y,n=4)  #Out of sample Returns
		  v.out <- tail(v0,n=4) #Out of sample % change in Weekly Volume
		  X.out0 <- tail(X,n=5)
		  X.out <- X.out0[-nrow(X.out0),] #Out of sample price
	  #VF-GARCH----
	  #ARMAX -> GJR -> GAM
		  #First Iteration
			  r.tfm <- 0
		    r.gjr <- 0
		    r.in1 <- 0
		    r.11 <- 0
		    r.21 <- 0
		    r.31 <- 0
		    y1 <- 0
		    res1 <- 0
		    #ARMA(1,1) + TFM:AR(1)
			  try({init1	<- arimax(y, order=c(1,0,1), xtransf=v, transfer=list(c(1,0)))
			  r.tfm	<- init1$residuals	
			  #GJR-GARCH(1,1)
			  in2.1	<- garchFit(~aparch(1,1), delta=2, data=r.tfm, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
			  r.gjr	<- in2.1@residuals
			  #GAM
			  in3.1	<- gam(r.gjr~s(X.dly$Mo,bs="cr")+s(X.dly$Tu,bs="cr")+s(X.dly$We,bs="cr")+s(X.dly$Th,bs="cr")+s(X.dly$Fr,bs="cr"))
			  r.in1	<- in3.1$residuals},silent=TRUE)
		  #Initial Values
			  d	<- mse(r.in1)	#test criteria
			  m1	<- 0			#mse
			  res1	<- r.in1		#residuals of final iteration
			  y1	<- y - res1		#updated response
			  j1 <- 1
		  #Iteration
		  	try(while(d > 0.005){
		  		mod11	<- arimax(y1, order=c(1,0,1), xtransf=v, transfer=list(c(1,0)))
		  		r.11	<- mod11$residuals
				
		  		mod21	<- garchFit(~aparch(1,1), delta=2, data=r.11, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
		  		r.21	<- mod21@residuals

		  		mod31	<- gam(r.21~s(X.dly$Mo,bs="cr")+s(X.dly$Tu,bs="cr")+s(X.dly$We,bs="cr")+s(X.dly$Th,bs="cr")+s(X.dly$Fr,bs="cr"))
		  		r.31	<- mod31$residuals

		  		res1	<- r.31
			  	y1	<- y - res1
		  		d	<- d - mse(res1)
		  		m1	<- mse(res1)
		  		j1 <- j1 + 1
		  	},silent=TRUE)
		  #Fitted Values of Algo 1
		  	f.11	<- y - r.11
		  	f.21	<- mod21@fitted
		  	f.31	<- mod31$fitted.values
		  	f1	<- f.11 + f.21 + f.31
	  	#Out of Sample Forecasting
		  	t.ar1 <- mod11$coef[4]
		  	t.ma1 <- mod11$coef[5]
		  	tf1 <- filter(v0,filter=t.ar1,method='recursive',sides=1)*t.ma1
		  	stock1 <- data.frame(cbind(Y, tf1))
		  	train1 <- stock1[1:l-1,]
		  	test1 <- tail(stock1,n=4)
		  	im.ar1 <- Arima(window(train1$Y),order=c(1,0,0))
		  	im.ma1 <- Arima(im.ar1$residuals,order=c(0,0,1))
		  	im.v1 <- glm(im.ma1$fitted~tf1,data=train1)
		  	it.f1 <- im.ar1$fitted + im.ma1$fitted + im.v1$fitted.values
		  	#Forecasting Iteration
		  	b <- mdpe(y,it.f1)
		  	past1 <- mdpe(y,it.f1)
		  	o1 <- 0
		  	y1 <- y - im.v1$residuals
		  	f1j <- 1
		  	while(b > 0.005){
		  	  m.ar1 <- Arima(y1,order=c(1,0,0))
		  	  p.ar1 <- predict(m.ar1, n.ahead=4)
		  	  m.ma1 <- Arima(m.ar1$residuals,order=c(0,0,1))
		  	  p.ma1 <- predict(m.ma1,n.ahead=4)
		  	  m.v1 <- glm(m.ma1$fitted~tf1, data=train1)
		  	  p.v1 <- predict.glm(m.v1, newdata=test1)
		  	  t.p1 <- p.ar1$pred + p.ma1$pred + as.vector(ts(p.v1))
		  	  curr1 <- mdpe(y.out, t.p1)
		  	  if(curr1>past1) o1 = past1 else o1 = curr1
		  	  
		  	  b <- b - mdpe(y.out, t.p1)
		  	  past1 <- curr1
		  	  y1 <- y - m.v1$residuals
		  	  f1j <- f1j + 1
		  	}
    #VF-ARMA----
	  #GAM -> ARMAX -> GJR
	  	#First Iteration
	  		e.gam <- 0
		  	y.upd <- 0
		  	e.tfm <- 0
		  	e.gjr <- 0
		  	r.12 <- 0
		  	r.22 <- 0
		  	r.32 <- 0
		  	res2 <- 0
		  	#GAM
		  	try({
	  	  init2	<- gam(y~s(X.dly$Mo,bs="cr")+s(X.dly$Tu,bs="cr")+s(X.dly$We,bs="cr")+s(X.dly$Th,bs="cr")+s(X.dly$Fr,bs="cr"))
	  		e.gam	<- init2$residuals
	  		y.upd	<- y - e.gam
			  #ARMAX
			  in2.2	<- arimax(y.upd, order=c(1,0,1), xtransf=v, transfer=list(c(1,0)))
			  e.tfm	<- in2.2$residuals
			  #GJR
			  in3.2	<- garchFit(~aparch(1,1), delta=2, data=e.tfm, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
			  e.gjr	<- in3.2@residuals
			  },silent = TRUE)
		  #Initial Values
		  	d	<- mse(e.gjr)	#test criteria
		  	m2	<- 0			#mse
		  	res2	<- e.gjr		#residuals of final iteration
		  	y2	<- y - res2		#updated response
		  	j2 <- 1
		  #Iteration
		  	try(while(d > 0.005){
		  		mod12	<- gam(y2~s(X.dly$Mo,bs="cr")+s(X.dly$Tu,bs="cr")+s(X.dly$We,bs="cr")+s(X.dly$Th,bs="cr")+s(X.dly$Fr,bs="cr"))
		  		r.12	<- mod12$residuals
		  		y.upd	<- y - r.12

		  		mod22	<- arimax(y.upd, order=c(1,0,1), xtransf=v, transfer=list(c(1,0)))
		  		r.22	<- mod22$residuals

		  		mod32	<- garchFit(~aparch(1,1), delta=2, data=r.22, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
		  		r.32	<- mod32@residuals

		  		res2	<- r.32
		  		y2	<- y - res2
		  		d	<- d - mse(res2)
		  		m2	<- mse(res2)
		  		j2 <- j2 + 1
		  	},silent=TRUE)
		  #Fitted Values of Algo 2
		  	f.12	<- mod12$fitted.values
		  	f.22	<- y - r.22
		  	f.32	<- mod32@fitted
		  	f2	<- f.12 + f.22 + f.32
		  #Out of Sample Forecasting
		  	o2 <- 0
		  	try({im.g <- gam(y~s(Mo,bs="cr")+s(Tu,bs="cr")+s(We,bs="cr")+s(Th,bs="cr")+s(Fr,bs="cr"))
		  	y.upd<- y - im.g$residuals
		  	t.ar2 <- mod22$coef[4]
		  	t.ma2 <- mod22$coef[5]
		  	tf2 <- filter(v0,filter=t.ar2,method='recursive',sides=1)*t.ma2
		  	stock2 <- data.frame(cbind(Y,tf2))
		  	train2 <- stock2[1:l-1,]
		  	test2 <- tail(stock2,n=4)
		  	im.ar2 <- Arima(y.upd,order=c(1,0,0))
		  	im.ma2 <- Arima(im.ar2$residuals,order=c(0,0,1))
		  	im.v2 <- glm(im.ma2$fitted~tf2,data=train2)
		  	it.f2 <- im.g$fitted.values + im.ar2$fitted + im.ma2$fitted + im.v2$fitted.values
		  #Forecasting Iteration
		  	b <- mdpe(y,it.f2)
		  	past2 <- mdpe(y,it.f2)
		  	y2 <- y - im.v2$residuals
		  	f2j <- 1
		  	while(b > 0.005){
		  	  m.g <- gam(y2~s(Mo,bs="cr")+s(Tu,bs="cr")+s(We,bs="cr")+s(Th,bs="cr")+s(Fr,bs="cr"))
		  	  p.g <- predict(m.g,newdata=X.out)
		  	  y.upd<- y - m.g$residuals
		  	  m.ar2 <- Arima(y.upd,order=c(1,0,0))
		  	  p.ar2 <- predict(m.ar2, n.ahead=4)
		  	  m.ma2 <- Arima(m.ar2$residuals,order=c(0,0,1))
		  	  p.ma2 <- predict(m.ma2, n.ahead=4)
		  	  m.v2 <- glm(m.ma2$fitted~tf2,data=train2)
		  	  p.v2 <- predict.glm(m.v2, newdata=test2)
		  	  t.p2 <- as.vector(ts(p.g)) + p.ar2$pred + p.ma2$pred + as.vector(ts(p.v2))
		  	  curr2 <- mdpe(y.out, t.p2)
		  	  if(curr2>past2) o2 = past2 else o2 = curr2
		  	  
		  	  b <- b - mdpe(y.out, t.p2)
		  	  past2 <- curr2
		  	  y2 <- y - m.v2$residuals
		  	  f2j <- f2j + 1
		  	}})
	  #Benchmark:	ARIMAX -> GARCH----
		  b.init	<- arimax(y, order=c(1,0,1), xtransf=v, transfer=list(c(1,0)))
		  r.tfm	<- b.init$residuals	
		  #GARCH and GJR
		  b.garch	<- garchFit(~garch(1,1), data=r.tfm, trace=FALSE,include.mean=FALSE)
		  r.b1	<- b.garch@residuals
		  
		  try({b.gjr	<- garchFit(~aparch(1,1), delta=2, data=r.tfm, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
		  r.b2	<- b.gjr@residuals})
		  
		  #Fitted Values of Bench 1
		  	f.tfm	<- y - r.tfm
		  	f.b21	<- b.garch@fitted
		  	f.b1	<- f.tfm + f.b21
		  #Out of Sample Forecasting
		  	bench.ar1 <- b.init$coef[4]
		  	bench.ma1 <- b.init$coef[5]
		  	tfb <- filter(v0,filter=bench.ar1,method='recursive',sides=1)*bench.ma1
		  	stockb <- data.frame(cbind(Y,tfb))
		  	trainb <- stockb[1:l-1,]
		  	testb <- tail(stockb,n=4)
		  	b.ar <- Arima(window(trainb$Y),order=c(1,0,0))
		  	p.arb <- predict(b.ar, n.ahead=4)
		  	b.ma <- Arima(b.ar$residuals,order=c(0,0,1))
		  	p.mab <- predict(b.ma, n.ahead=4)
		  	b.v <- glm(b.ma$fitted~tfb,data=trainb)
		  	p.vb <- predict.glm(b.v, newdata=testb)
		  	t.pb <- p.arb$pred + p.mab$pred + as.vector(ts(p.vb))
		  	
		  	ob <- mdpe(y.out,t.pb)
		  #Fitted Values of Bench 2
		  	f.b22	<- b.gjr@fitted
		  	f.b2	<- f.tfm + f.b22
	  #Compile the Evaluation Statistics of each model----
  		m1.mse<- c(m1.mse,m1)
  		m2.mse<- c(m2.mse,m2)
  		b1.mse<- c(b1.mse,mse(r.b1))
  		b2.mse<- c(b2.mse,mse(r.b2))	
  		m1.mad<- c(m1.mad,mad(y,f1))
  		m2.mad<- c(m2.mad,mad(y,f2))
  		b1.mad<- c(b1.mad,mad(y,f.b1))
  		b2.mad<- c(b2.mad,mad(y,f.b2))
  		m1.rse<- c(m1.rse,sqrt(m1))
  		m2.rse<- c(m2.rse,sqrt(m2))
  		b1.rse<- c(b1.rse,rmse(r.b1))
  		b2.rse<- c(b2.rse,rmse(r.b2))
  		f1.mpe<- c(f1.mpe,o1)
  		f2.mpe<- c(f2.mpe,o2)
  		b.mpe<- c(b.mpe,ob)
  		J1 <- c(J1,j1)
  		J2 <- c(J1,j2)
  		fJ1 <- c(fJ1,f1j)
  		fJ2 <- c(fJ2,f2j)
    }
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#Model Comparison----
  m1.rse[m1.rse==0] <- NA
  m1.stat <- cbind(m1.rse,m1.mad)
  m1.stat[which(is.na(m1.stat)),] <- NA
  m2.rse[m2.rse==0] <- NA
  m2.stat <- cbind(m2.rse,m2.mad)
  m2.stat[which(is.na(m2.stat)),] <- NA
  b1.stat <- cbind(b1.rse,b1.mad)
  b2.stat <- cbind(b2.rse,b2.mad)
  f2.mpe[f2.mpe==0] <- NA
  mods.sum <- c(mean(m2.stat[,1],na.rm=TRUE),mean(m1.stat[,1],na.rm=TRUE),
                mean(b1.stat[,1],na.rm=TRUE),mean(b2.stat[,1],na.rm=TRUE),
                mean(m2.stat[,2],na.rm=TRUE),mean(m1.stat[,2],na.rm=TRUE),
                mean(b1.stat[,2],na.rm=TRUE),mean(b2.stat[,2],na.rm=TRUE))
  sum.tab <- matrix(mods.sum, ncol=2)
  colnames(sum.tab) <- c("RMSE","MAD")
  rownames(sum.tab) <- c("VF-ARMA", "VF-GARCH", "GARCH", "GJR")
  sum.tab
#Import Output----
  write.csv(sum.tab,file="scen1.csv")
  summary(m1.stat)
  summary(m2.stat)
  summary(f2.mpe)
  c(length(m1.rse),length(m2.rse))
  c(mean(f2.mpe,na.rm=TRUE),mean(f1.mpe),mean(b.mpe))
beep("fanfare")

#Diebold-Mariano test----
#VF-ARMA vs GJR
dm.test(r.b1,res2, h = 4)

#VF-ARMA vs VF-GARCH
dm.test(res1,res2, h = 4)

#Plots----
x0	= 100			#Initial Price
mu	= 0.2			#Drift: Expected Annual Return
sigma	= 0.30	#Volatility: Expected Annual Volatility
psi	= 0.5			#Postulated effect of Change in Trading Volume
m.dly	= mu/252		#Daily Return
s.dly	= sigma/sqrt(252)	#Daily Volatility
#In between simulation settings, value of t is changed between 205, 510, and 1530
t = 1530 #Length of series

set.seed(1)
#In between simulations, s will be switched in between Independent and Autocorrelated settings
#	<- rnorm(t+20, m.dly, s.dly)		#Independent HF Data
s <- arima.sim(list(order=c(1,0,0),ar=0.5),rand.gen=rnorm,mean=m.dly,sd=s.dly,n=t+20) #Autoregressive HF Data
u	<- runif(t + 219,-1,1)			#For Volume
#Compute for daily volume and daily price
u1	<- c(0, cumsum(u))			#Percentage Change in Volume
w	<- abs(1000000*u1)			#Daily Trading Volume

#In between simulations, x will be switched in between the Linear and Exponential functional form settings
lx	<- log(x0) + cumsum(s)			#Daily Return
#	<- exp(lx) - (psi*u1[-c(1:200)])	#Linear Generation
x <- exp(lx - psi*u1[-c(1:200)])  #Exponential Generation

#Compile the Weekly Data
X	<- matrix(x,ncol=5,byrow=TRUE)
colnames(X) <- c("Mo","Tu","We","Th","Fr")
X	<- data.frame(X)				#Daily Price, compiled by week
#Volume
w.dly	<- w[-c(1:200)]
W	<- matrix(w.dly,ncol=5,byrow=TRUE)
V	<- rowSums(W)				#Weekly Trading Volume
#Compute the input and output variables
l <- t/5
Y	<- diff(log(X$Fr))
y <- Y[c(1:l-1)]			#C.C. Returns
v0	<- diff(log(V))				
v <- v0[c(1:l-1)] #Percentage Change in Weekly Volume
X.dly	<- X[c(1:l-1),]
Mo <- as.vector(X.dly$Mo)
Tu <- as.vector(X.dly$Tu)
We <- as.vector(X.dly$We)
Th <- as.vector(X.dly$Th)
Fr <- as.vector(X.dly$Fr)


par(mfrow=c(2,2))
plot(ts(x[1:t]),ylab="HF Data x",main="Uncorrelated and Linear")
plot(ts(x[1:t]),ylab="HF Data x",main="Uncorrelated and Exponential")
plot(ts(x[1:t]),ylab="HF Data x",main="Autocorrelated and Linear")
plot(ts(x[1:t]),ylab="HF Data x",main="Autocorrelated and Exponential")

plot(ts(Y),ylab="Response variable y",main="Uncorrelated and Linear")
plot(ts(Y),ylab="Response variable y",main="Uncorrelated and Exponential")
plot(ts(Y),ylab="Response variable y",main="Autocorrelated and Linear")
plot(ts(Y),ylab="Response variable y",main="Autocorrelated and Exponential")


#Empirical Study
#Stock Data----
#Import Data
ali <- read.csv(".../ali.csv")
xtra <- read.csv(".../extra data.csv")

#Converts daily data into a weekly format
X <- matrix(ali$Close,ncol=5,byrow=TRUE)
colnames(X) <- c("Mo","Tu","We","Th","Fr")
X	<- data.frame(X)

#Converts extra data into return

#PSEi
x1 <- diff(log(xtra$PSEi))
X1 <- matrix(x1,ncol=5,byrow=TRUE)
colnames(X1) <- c("Mo","Tu","We","Th","Fr")
X1	<- data.frame(X1)

#USDPHP
x2 <- xtra$USDPHP[-1]
X2 <- matrix(x2,ncol=5,byrow=TRUE)
colnames(X2) <- c("Mo","Tu","We","Th","Fr")
X2	<- data.frame(X2)

#Takes daily volume, and aggregates it into a weekly format
W <- matrix(as.numeric(ali$Volume),ncol=5,byrow=TRUE)
V <- rowSums(W)

#Converts Closing Price into Price Return and Weekly Volume into percentage change using log differences
y	<- diff(log(X$Fr))
v <- diff(log(V))

PS <- diff(log(X1$Fr))
UP <- diff(log(X2$Fr))

#Removes last row in the price data, simulating the unknown but upcoming prices of the week
X.dly <- X[-nrow(X),]
X1.dly <- X1[-nrow(X1),]
X2.dly <- X2[-nrow(X2),]

#Sets up data in terms of last week's price return and volume 
#Original stock data
stock <- cbind(y,v,X.dly)

#Stock data including Stock Index and Exchange rate
stock1 <- cbind(y,v,X1.dly)
stock2 <- cbind(y,v,X2.dly)

#Sets up training and test data set, train will consist of the past two years while test would be next month's data
train <- stock2[1:(length(y)-4),]
test <- tail(stock2,n=4)

Mo <- as.vector(train$Mo)
Tu <- as.vector(train$Tu)
We <- as.vector(train$We)
Th <- as.vector(train$Th)
Fr <- as.vector(train$Fr)


#Plots
par(mfrow=c(1,2))
plot(xts(ali$Close,as.Date(ali$Date, "%m/%d/%Y")),main="Daily Closing Price", ylab=NA)
plot(ts(X.dly$Fr),main="Weekly Closing Price",ylab=NA)

acf(ali$Close,main=NA)
pacf(ali$Close,main=NA)

plot(ts(V),type="l",main="Weekly Volume",ylab=NA)
plot(ts(v),main="Weekly % Change in Volume",ylab=NA)

#Modeling
#Can switch between Volume (v) or Stock Index (PS) or Exchange Rate (UP)
#VF-GARCH
init1	<- arimax(train$y, order=c(1,0,1), xtransf=train$v, transfer=list(c(1,0)))
r.tfm	<- init1$residuals	
#GJR-GARCH(1,1)
in2.1	<- garchFit(~aparch(1,1), delta=2, data=r.tfm, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
r.gjr	<- in2.1@residuals
#GAM
in3.1	<- gam(r.gjr~s(train$Mo,bs="cr")+s(train$Tu,bs="cr")+s(train$We,bs="cr")+s(train$Th,bs="cr")+s(train$Fr,bs="cr"))
r.in1	<- in3.1$residuals
d	<- mse(r.in1)	#test criteria
m1	<- 0			#mse
y1	<- train$y - r.in1		#updated response
#Iteration
try(while(d > 0.005){
  mod11	<- arimax(y1, order=c(1,0,1), xtransf=train$v, transfer=list(c(1,0)))
  r.11	<- mod11$residuals
  
  mod21	<- garchFit(~aparch(1,1), delta=2, data=r.11, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
  r.21	<- mod21@residuals
  
  mod31	<- gam(r.21~s(train$Mo,bs="cr")+s(train$Tu,bs="cr")+s(train$We,bs="cr")+s(train$Th,bs="cr")+s(train$Fr,bs="cr"))
  r.31	<- mod31$residuals
  
  res1	<- r.31
  y1	<- train$y - res1
  d	<- d - mse(res1)
  m1	<- mse(res1)
},silent=TRUE)
#Fitted Values of Algo 1
f.11	<- train$y - r.11
f.21	<- mod21@fitted
f.31	<- mod31$fitted.values
f1	<- f.11 + f.21 + f.31
c(rmse(res1),mad(train$y,f1))
#Out of Sample Forecasting
t.ar1 <- mod11$coef[4]
t.ma1 <- mod11$coef[5]
tf1 <- filter(v,filter=t.ar1,method='recursive',sides=1)*t.ma1
tf1 <- data.frame(tf1)
Stock1 <- cbind(stock2,tf1)
train1 <- data.frame(Stock1[1:(length(y)-4),])
test1 <- tail(Stock1,n=4)
im.ar1 <- Arima(train1$y,order=c(1,0,0))
im.ma1 <- Arima(im.ar1$residuals,order=c(0,0,1))
im.v1 <- glm(im.ma1$fitted~train1$tf1)
it.f1 <- im.ar1$fitted + im.ma1$fitted + im.v1$fitted.values
#Forecasting Iteration
b <- mdpe(train1$y,it.f1)
past1 <- mdpe(train1$y,it.f1)
o1 <- 0
y1 <- train1$y - im.v1$residuals
f1j <- 1
while(b > 0.005){
  m.ar1 <- Arima(y1,order=c(1,0,0))
  p.ar1 <- predict(m.ar1, n.ahead=4)
  m.ma1 <- Arima(m.ar1$residuals,order=c(0,0,1))
  p.ma1 <- predict(m.ma1,n.ahead=4)
  m.v1 <- glm(m.ma1$fitted~tf1, data=train1)
  p.v1 <- predict.glm(m.v1, newdata=test1)
  t.p1 <- p.ar1$pred + p.ma1$pred + as.vector(ts(p.v1))
  curr1 <- mdpe(test1$y, t.p1)
  if(curr1>past1) o1 = past1 else o1 = curr1
  
  b <- b - mdpe(test1$y, t.p1)
  past1 <- curr1
  y1 <- train1$y - m.v1$residuals
  f1j <- f1j + 1
}
c(rmse(res1),mad(train$y,f1),curr1)

#VF-ARMA
try({
  init2	<- gam(y~s(Mo,bs="cr")+s(Tu,bs="cr")+s(We,bs="cr")+s(Th,bs="cr")+s(Fr,bs="cr"),data=train)
  e.gam	<- init2$residuals
  y.upd	<- train$y - e.gam
  #ARMAX
  in2.2	<- arimax(y.upd, order=c(1,0,1), xtransf=train$v, transfer=list(c(1,0)))
  e.tfm	<- in2.2$residuals
  #GJR
  in3.2	<- garchFit(~aparch(1,1), delta=2, data=e.tfm, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
  e.gjr	<- in3.2@residuals
},silent = TRUE)
#Initial Values
d	<- mse(e.gjr)	#test criteria
m2	<- 0			#mse
res2	<- e.gjr		#residuals of final iteration
y2	<- train$y - res2		#updated response
#Iteration
try(while(d > 0.005){
  mod12	<- gam(y2~s(train$Mo,bs="cr")+s(train$Tu,bs="cr")+s(train$We,bs="cr")+s(train$Th,bs="cr")+s(train$Fr,bs="cr"))
  r.12	<- mod12$residuals
  y.upd	<- train$y - r.12
  
  mod22	<- arimax(y.upd, order=c(1,0,1), xtransf=train$v, transfer=list(c(1,0)))
  r.22	<- mod22$residuals
  
  mod32	<- garchFit(~aparch(1,1), delta=2, data=r.22, include.delta=FALSE, trace=FALSE, include.mean=FALSE)
  r.32	<- mod32@residuals
  
  res2	<- r.32
  y2	<- train$y - res2
  d	<- d - mse(res2)
  m2	<- mse(res2)
},silent=TRUE)
#Fitted Values of Algo 2
f.12	<- mod12$fitted.values
f.22	<- train$y - r.22
f.32	<- mod32@fitted
f2	<- f.12 + f.22 + f.32
c(rmse(res2),mad(train$y,f2))

#Out of Sample Forecasting
o2 <- 0
t.ar2 <- mod22$coef[4]
t.ma2 <- mod22$coef[5]
tf2 <- filter(v,filter=t.ar2,method='recursive',sides=1)*t.ma2
tf2 <- data.frame(tf2)
Stock2 <- data.frame(cbind(stock2,tf2))
train2 <- Stock2[1:(length(y)-4),]
test2 <- tail(Stock2,n=4)
try({im.g <- gam(train2$y~s(Mo,bs="cr")+s(Tu,bs="cr")+s(We,bs="cr")+s(Th,bs="cr")+s(Fr,bs="cr"))
y.upd<- train2$y - im.g$residuals
im.ar2 <- Arima(y.upd,order=c(1,0,0))
im.ma2 <- Arima(im.ar2$residuals,order=c(0,0,1))
im.v2 <- glm(im.ma2$fitted~tf2,data=train2)
it.f2 <- im.g$fitted.values + im.ar2$fitted + im.ma2$fitted + im.v2$fitted.values
#Forecasting Iteration
b <- mdpe(train2$y,it.f2)
past2 <- mdpe(train2$y,it.f2)
y2 <- train2$y - im.v2$residuals
f2j <- 1
while(b > 0.005){
  m.g <- gam(y2~s(Mo,bs="cr")+s(Tu,bs="cr")+s(We,bs="cr")+s(Th,bs="cr")+s(Fr,bs="cr"))
  p.g <- predict(m.g,newdata=test2)
  y.upd<- train2$y - m.g$residuals
  m.ar2 <- Arima(y.upd,order=c(1,0,0))
  p.ar2 <- predict(m.ar2, n.ahead=4)
  m.ma2 <- Arima(m.ar2$residuals,order=c(0,0,1))
  p.ma2 <- predict(m.ma2, n.ahead=4)
  m.v2 <- glm(m.ma2$fitted~tf2,data=train2)
  p.v2 <- predict.glm(m.v2, newdata=test2)
  t.p2 <- as.vector(ts(p.g)) + p.ar2$pred + p.ma2$pred + as.vector(ts(p.v2))
  curr2 <- mdpe(test2$y, t.p2)
  if(curr2>past2) o2 = past2 else o2 = curr2
  
  b <- b - mdpe(test2$y, t.p2)
  past2 <- curr2
  y2 <- train2$y - m.v2$residuals
  f2j <- f2j + 1
}})
c(rmse(res2),mad(train$y,f2),curr2)

#Benchmark:	ARIMAX -> GARCH----
b.init	<- arimax(train$y, order=c(1,0,1), xtransf=train$v, transfer=list(c(1,0)))
r.tfm	<- b.init$residuals	
#GARCH and GJR
bench	<- garchFit(~garch(1,1), data=r.tfm, trace=FALSE,include.mean=FALSE)
r.b1	<- bench@residuals

#Fitted Values of Bench 1
f.tfm	<- train$y - r.tfm
f.b21	<- bench@fitted
f.b1	<- f.tfm + f.b21
c(rmse(r.b1),mad(train$y,f.b1))
#Fitted Values of Bench 2
f.b22	<- in2.1@fitted
f.b2	<- f.tfm + f.b22
c(rmse(r.gjr),mad(train$y,f.b2))


#Out of Sample Forecasting
bench.ar1 <- b.init$coef[4]
bench.ma1 <- b.init$coef[5]
tfb <- filter(v,filter=bench.ar1,method='recursive',sides=1)*bench.ma1
tfb <- data.frame(tfb)
stockb <- cbind(stock,tfb)
trainb <- stockb[1:(length(y)-4),]
testb <- tail(stockb,n=4)
b.ar <- Arima(trainb$y,order=c(1,0,0))
p.arb <- predict(b.ar, n.ahead=4)
b.ma <- Arima(b.ar$residuals,order=c(0,0,1))
p.mab <- predict(b.ma, n.ahead=4)
b.v <- glm(b.ma$fitted~tfb,data=trainb)
p.vb <- predict.glm(b.v, newdata=testb)
t.pb <- p.arb$pred + p.mab$pred + as.vector(ts(p.vb))

ob <- mdpe(testb$y,t.pb)

