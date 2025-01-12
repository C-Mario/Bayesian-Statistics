library(invgamma)
options(scipen = 5)

datos <- c(495,541,1461,1555,1603,2201,2750,3468,3516,4319,6622,7728,13159,21194)
mean(datos)
sd(datos)
s <- sum(datos)
n <- length(datos)
ybar <- s/n

################## Punto 4. ##################

# hiperparámetros de la distribución previa:
a <- 8.25
b <- 32625

# parámetros de la distribución posterior:
ap <- a+n
bp <- b+s

# Distribución Gamma Inversa
GI <- function(lambda,alpha,beta){
  beta^(alpha)/gamma(alpha)*lambda^(-(alpha+1))*exp(-beta/lambda)}

# Gráfico de la distribución a priori y posterior
# png(file="punto4.png",width=13,height=10,units="cm",res=1000, pointsize=6)
par(bg='white')
par(cex=1.5)
curve(expr = dinvgamma(x,shape=a,rate=b), from = 0, to = 13000, col = "slateblue1", 
      xlab = expression(lambda), ylim=c(0,0.00042), 
      ylab = expression(Probabilidad), lwd=2, cex.axis=0.8,
      cex.lab=1.1,font.lab=3,col.lab="black", col.axis="gray25")
curve(expr = dinvgamma(x,shape=ap,rate=bp), col = "turquoise3",add = TRUE, lwd=2)
legend("topright", legend = c("Previa","Posterior"), col = c("slateblue1", "turquoise3"), 
       lty = 1, lwd = 2, bty = "o",pt.cex = cex, border = "black", bg = "ghostwhite")
# dev.off()

################# Punto 6. ##################

# IC bayesiano
meanp <- bp/(ap-1) # media posterior (estimación puntual)
varp <- bp^2/((ap-1)^2*(ap-2)) # varianza posterior
sdp <- sqrt(varp) # desviación estándar posterior
cvarp <- sdp/meanp # coeficiente de variacion posterior
meanp; round(cvarp*100,3)
# Intervalo de credibilidad bayesiano al 95%
LI_b <- qinvgamma(0.025, shape=ap, scale = 1/bp)
LS_b <- qinvgamma(0.975, shape=ap, scale = 1/bp)
LI_b; LS_b

# ICA frecuentista
mle <- ybar # estimación puntual
varmle <- ybar^2/n #Varianza del mle (I^-1)
sdmle <- sqrt(varmle) #desviación estándar del mle
cvarmle <- sdmle/mle #coeficiente de variación del mle
mle; round(cvarmle*100,3)
# Intervalo de confianza frecuentista al 95%
LI_f <- ybar-qnorm(0.975)*ybar/sqrt(n)
LS_f <- ybar+qnorm(0.975)*ybar/sqrt(n)
LI_f; LS_f

# IC Bootstrap
set.seed(9099)
boots<-function(x,n){
  MBootstraps <- replicate(50000,sample(x,n,replace=T))
  medias <- apply(MBootstraps,2,mean)}

medias <- boots(datos,14)
meanboos <- mean(medias) # estimación puntual bootstrap
cvarboos <- sd(medias)/meanboos # coeficiente de variación bootstrap
meanboos; round(cvarboos*100,3)
LI_boo <- as.numeric(quantile(x=medias,0.025))
LS_boo <- as.numeric(quantile(x=medias,0.975))
LI_boo; LS_boo

################### Punto 7. ###################

# probabilidad a posteriori de que lambda sea menor a 4000

round(pinvgamma(4000, shape=ap, scale=1/bp,lower.tail = TRUE),4)

# Método de monte carlo para aproximar la distribución predictiva posterior
# 1) simular 50000 valores del parámetro \lambda a partir de la distribución posterior

set.seed(9190)
B <- 50000
lambdasp <- rinvgamma(B, shape=ap, scale=1/bp)

# 2) generar B números aleatorios de una distribución exponencial condicionada en cada lambda

BN <- rexp(B,rate=1/lambdasp)
set.seed <- NULL

round(mean(BN<4000),4)

# BN <- matrix(nrow=B,ncol=14)
# for (i in 1:B){
# BN[i,] <- rexp(n,rate=1/lambdasp[i])
# }

################ Punto 8. #################
# Factor de bayes 10
lambda_0 <- 4000
b_10 <- exp(a*log(b)-ap*log(bp)+lgamma(ap)-lgamma(a)+n*log(lambda_0)+(s/lambda_0))
round(b_10,3)

################ Punto 9. #################

datos2 <- c(294,569,766,1576,1602,2015,2166,3885,8141,10285)
summary(datos2)
sd(datos2)
n2 <- length(datos2)
s2 <- sum(datos2)

# hiperparámetros de la distribución posterior (alambre tipo 2)
ap2 <- a+n2 
bp2 <- b+s2

b2_10 <- exp(a*log(b)+lgamma(ap)+lgamma(ap2)+(ap+n2)*log(bp+s2)-
               lgamma(a)-(ap)*log(bp)-(ap2)*log(bp2)-lgamma(ap+n2))
round(b2_10,3)


############### Punto 10. #################

B <- 50000
## Alambre Tipo 1
# Estadistico observado
t1_obs <- s/n
t1_obs
# Distribucion predictiva posterior
t1_mc <- NULL
set.seed(9091)
lambda1_mc <- rinvgamma(B,shape=ap,scale=1/bp)
set.seed(9091)
for(i in 1:B){
  y1_rep <- rexp(n,rate = 1/lambda1_mc[i])
  t1_mc[i] <- mean(y1_rep)
}
t1_mc

# Grafico
# png(file="punto10a.png",width=13,height=10,units="cm",res=1000, pointsize=6)
par(bg='white')
par(cex=1.5)
hist(x = t1_mc, freq = F, col = "gray60", border = "gray60", xlab = "t1", main="",
     ylab = expression(paste("p","(",t1," | ",y,")",sep="")),ylim=c(0,0.00030))
lines(density(t1_mc, adjust=2), col = "steelblue1", lwd = 2)
abline(v = t1_obs, col = "firebrick3", lwd = 2, lty = 1)
abline(v = quantile(x = t1_mc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = "darkolivegreen3")
legend("topright", legend = c("Posterior", "t1 obs", "IC 95%"), 
       col = c("steelblue1", "firebrick3", "darkolivegreen3"), lty = 1, lwd = 2,
       bty = "o",pt.cex = cex, border = "black", bg = "ghostwhite")
# dev.off()

#Valor predictivo posterior
round(mean(t1_mc > t1_obs),3)

## Alambre Tipo 2
# Estadistico observado
t2_obs <- s2/n2
t2_obs
# Distribucion predictiva posterior
t2_mc <- NULL
set.seed(9091)
lambda2_mc <- rinvgamma(B,shape=a+n2,scale=1/(b+s2))
set.seed(9091)
for(i in 1:B){
  y2_rep <- rexp(n,rate = 1/lambda2_mc[i])
  t2_mc[i] <- mean(y2_rep)
}
t2_mc

# Grafico
# png(file="punto10b.png",width=13,height=10,units="cm",res=1000, pointsize=6)
par(bg='white')
par(cex=1.5)
hist(x = t2_mc, freq = F, col = "gray90", border = "gray90", xlab = "t2", 
     ylab = expression(paste("p","(",t2," | ",y,")",sep="")), main = "", ylim = c(0,0.00035))
axis(1, at=seq(0, max(t2_mc), 2500))
lines(density(t2_mc, adjust=2), col = "coral3", lwd = 2)
abline(v = t2_obs, col = "darkslategray", lwd = 2, lty = 1)
abline(v = quantile(x = t2_mc, probs = c(0.025, 0.975)), lty = 2, lwd = 2, col = "mediumorchid1")
legend("topright", legend = c("Posterior", "t2 obs", "IC 95%"), 
       col = c("coral3", "darkslategray", "mediumorchid1"),
       lty = 1, lwd = 2, bty = "o",pt.cex = cex, border = "black", bg = "ghostwhite")
# dev.off()

#Valor predictivo posterior
round(mean(t2_mc > t2_obs),3)

