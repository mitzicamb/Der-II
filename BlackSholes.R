#############################################################
### M\'etodos para modelo Black-Scholes y Black-Scholes76 ###
#############################################################

source("optimization.R")

# prima de una opci\'on americana por Black-Scholes -----------------------


BlackScholes<-function(strike,spot,t,vol,rate,call_put){
  moneyness<-spot/strike
  d1<-(log(moneyness)+(rate+vol^2/2)*t)/(vol*sqrt(t))
  d2<-d1-vol*sqrt(t)
  Discount<-exp(-rate*t)
  Forward<-exp(rate*t)*spot
  premium<-call_put*(Forward*pnorm(call_put*d1)-strike*pnorm(call_put*d2))*Discount
  return(premium)
}
#Ejemplo:
#BlackScholes(100,102,2,.25,.05,-1)



# Prima de una opci\'on americana por Black-Scholes76 ---------------------


BlackScholes76<-function(strike,forward,t,vol,rate,call_put){
  moneyness<-forward/strike
  d1<-(log(moneyness)+vol^2/2*t)/(vol*sqrt(t))
  d2<-d1-vol*sqrt(t)
  Discount<-exp(-rate*t)
  premium<-call_put*(forward*pnorm(call_put*d1)-strike*pnorm(call_put*d2))*Discount
  return(premium)
}
#Ejemplo:
#BlackScholes76(100,102,2,.25,.05,-1)



# Prima de una opci\'on americana por simulaci\'on -----------------------


BlackScholes_mc<-function(strike,spot,t,vol,rate,call_put,n=1000){
  S_t<-spot*exp(rnorm(n,(rate-vol^2/2)*t,vol*sqrt(t)))
  payoff<-pmax(call_put*(S_t-strike),0)
  Discount<-exp(-rate*t)
  return(mean(payoff)*Discount)
}
#Ejemplos:
#BlackScholes_mc(100,102,2,.25,.05,-1)
#BlackScholes_mc(100,102,2,.25,.05,-1,10000)



# M\'etodo para calcular Delta de Black-Scholes ---------------------------


BlackScholes_Delta<-function(strike,spot,t,vol,rate,call_put){
  moneyness<-spot/strike
  d1<-(log(moneyness)+(rate+vol^2/2)*t)/(vol*sqrt(t))
  return(call_put*pnorm(call_put*d1))
}
#Ejemplo:
#BlackScholes_Delta(100,102,2,.25,.05,-1)



# M\'etodo para calcular Delta de Black-Scholes76 ---------------------------


BlackScholes76_Delta<-function(strike,forward,t,vol,call_put){
  moneyness<-forward/strike
  d1<-(log(moneyness)+vol^2/2*t)/(vol*sqrt(t))
  return(call_put*pnorm(call_put*d1))
}
#Ejemplo:
#BlackScholes76_Delta(100,102,2,.25,.05,-1)



# M\'etodo para calcular Vega de Black-Scholes ----------------------------


BlackScholes_Vega<-function(strike,spot,t,vol,rate){
  moneyness<-spot/strike
  d2<-(log(moneyness)+(rate-vol^2/2)*t)/(vol*sqrt(t))
  return(strike*dnorm(d2)*sqrt(t))
}
#Ejemplo:
#BlackScholes_Vega(100,102,2,.25,.05,-1)



# M\'etodo para calcular Vega de Black-Scholes ----------------------------


BlackScholes76_Vega<-function(strike,forward,t,vol,rate){
  moneyness<-forward/strike
  d2<-(log(moneyness)-vol^2/2*t)/(vol*sqrt(t))
  return(strike*dnorm(d2)*sqrt(t))
}
#Ejemplo:
#BlackScholes76_Vega(100,102,2,.25,.05,-1)



# vol impl\'icita de Black-Scholes ----------------------------------------


BlackScholes_Impvol<-function(premium,strike,spot,t,rate,call_put,vol_init=0.2,max_iter=1e3,tolerance=1e-6,epsilon=1e-16,delta=1e-14){
  bs_mod<-function(vol){BlackScholes(strike,spot,t,vol,rate,call_put)}
  bs_vega_mod<-function(vol){BlackScholes_Vega(strike,spot,t,vol,rate)}
  if(premium<=0) stop("Premium must be positive.")
  vol_ret<-newton_raphson(vol_init,bs_mod,bs_vega_mod,premium,NULL,max_iter,tolerance,epsilon,delta)
  return(vol_ret)
}
#Ejemplo:
#BlackScholes_Impvol(8.505022,100,102,2,.05,-1)



# vol impl\'icita de Black-Scholes76 ----------------------------------------


BlackScholes76_Impvol<-function(premium,strike,forward,t,rate,call_put,vol_init=0.2,max_iter=1e3,tolerance=1e-6,epsilon=1e-16,delta=1e-14){
  bs_mod<-function(vol){BlackScholes76(strike,forward,t,vol,rate,call_put)}
  bs_vega_mod<-function(vol){BlackScholes76_Vega(strike,forward,t,vol,rate)}
  if(premium<=0) stop("Premium must be positive.")
  vol_ret<-newton_raphson(vol_init,bs_mod,bs_vega_mod,premium,NULL,max_iter,tolerance,epsilon,delta)
  return(vol_ret)
}
#Ejemplo:
#BlackScholes76_Impvol(11.93836,100,102,2,.05,-1)


