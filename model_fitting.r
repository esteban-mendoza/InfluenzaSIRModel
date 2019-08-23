# Librería 
library(tidyverse)
library(lubridate)

# Subestimación de beta e I_0
data = read.csv("infectados.csv", stringsAsFactors = FALSE) %>%
    mutate(fecha = dmy(fecha))

I = filter(data, fecha <= dmy("26-04-2009"))[,"I"]
t = 0:(length(I)-1)
Y = log(I)

model = lm(Y ~ t)

C = coef(model)

I_0 = exp(unname(C[1]))
beta_0 = unname(C[2])

h = function(t) {
    return(I_0 * exp(beta_0 * t))
}

plot(t, I, pch=19)
lines(t, h(t), col="red")

# Estimación de beta_0 e I_0
iters = 201000
alpha = 0.0000001

J_0 = function() {
    return(mean((I - h(t))^2)/2)
}

JI_0 = function() {
    return(mean((I - h(t))*(-exp(beta_0*t))))
}

Jbeta_0 = function() {
    return(mean((I - h(t))*(-I_0*t*exp(beta_0*t))))
}

tempI_0 = c(I_0)
tempBeta_0 = c(beta_0)

cost_0 = c(J_0())

for (i in 2:iters) {
    tempI_0 = c(tempI_0, I_0 - alpha*JI_0())
    tempBeta_0 = c(tempBeta_0, beta_0 - alpha*Jbeta_0())
    
    I_0 = tempI_0[i]
    beta_0 = tempBeta_0[i]
    
    cost_0 = c(cost_0, J_0())
}

I_0 = tempI_0[which.min(cost_0)]
beta_0 = tempBeta_0[which.min(cost_0)]

which.min(cost_0)

# I_0 = 0.2173231; beta_0 = 0.2516953
print(c(I_0, beta_0)) 
lines(t, h(t), col="blue")

# Propuesta de S_0
d05 = dmy("29-08-2005")
d09 = dmy("27-03-2009")
d10 = dmy("25-06-2010")

p05 = 103263388
p10 = 112336538

r = 100 * ( (p10 / p05)^(1 / as.integer(d10 - d05)) - 1 )
S_0 = p05 * (1 + r/100)^(as.integer(d09 - d05)) # S_0 = 109918558

# Estimación de mu y kappa
mu = 1 / 3

kappa_0 = (beta_0 + mu) / S_0 
R_0 = beta_0 / mu + 1

# kappa_0 = 5.322383e-09; R_0 = 1.755086e+00
print(c(kappa_0, R_0))

# Estimación lineal de beta_1 e I_1
I_ = filter(data, dmy("27-04-2009") <= fecha, fecha <= dmy("14-05-2009"))[,"I"]
t_ = 0:(length(I_)-1)
Y_ = log(I_)

model_ = lm(Y_ ~ t_)

C_ = coef(model_)

I_1 = exp(unname(C_[1]))
beta_1 = unname(C_[2])

# I_1 = 321.81661488; beta_1 = -0.06198531

g = function(t) {
    return(I_1 * exp(beta_1 * t))
}

plot(t_, I_, pch=19)
lines(t_, g(t_), col="red")

# Estimación de beta_1 e I_1
iters = 3e5
alpha_ = 0.0000001

J_1 = function() {
    return(mean((I_ - g(t_))^2)/2)
}

JI_1 = function() {
    return(mean((I_ - g(t_))*(-exp(beta_1*t_))))
}

Jbeta_1 = function() {
    return(mean((I_ - g(t_))*(-I_1*t_*exp(beta_1*t_))))
}

tempI_1 = c(I_1)
tempBeta_1 = c(beta_1)

cost_1 = c(J_1())

for (i in 2:iters) {
    tempI_1 = c(tempI_1, I_1 - alpha_*JI_1())
    tempBeta_1 = c(tempBeta_1, beta_1 - alpha_*Jbeta_1())
    
    I_1 = tempI_1[i]
    beta_1 = tempBeta_1[i]
    
    cost_1 = c(cost_1, J_1())
}

I_1 = tempI_1[which.min(cost_1)]
beta_1 = tempBeta_1[which.min(cost_1)]

# I_1 = 321.81661488; beta_1 = -0.06198531
print(c(I_1, beta_1))

lines(t_, g(t_), col="blue")

# Estimación de kappa_ y R_0_
kappa_1 = (beta_1 + mu) / S_0
R_1 = beta_1 / mu + 1

print(c(kappa_0, R_0))
print(c(kappa_1, R_1))

# Visualización de predicción de ambos modelos
I_v = filter(data, fecha <= dmy("14-05-2009"))[,"I"]
t_v = 0:(length(I_v)-1)

sir = read.csv("SIR.csv")
names(sir) = c("t", "S", "I", "R")

X = filter(sir, t <= 48)

plot(t_v, I_v, pch=19)
points(X[,"t"], X[,"I"], col="red")
points(X[,"t"], X[,"R"], col="green")
