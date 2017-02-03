### M\'etodo para calcular la prima de una opci\'on Call o Put por Black-Scholes
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

### M\'etodo para calcular la prima de una opci\'on Call o Put por Black-Scholes
BlackScholes_mc<-function(strike,spot,t,vol,rate,call_put,n=1000){
  S_t<-spot*exp(rnorm(n,(rate-vol^2/2)*t,vol*sqrt(t)))
  payoff<-pmax(call_put*(S_t-strike),0)
  Discount<-exp(-rate*t)
  return(mean(payoff)*Discount)
}
#Ejemplos:
#BlackScholes_mc(100,102,2,.25,.05,-1)
#BlackScholes_mc(100,102,2,.25,.05,-1,10000)


