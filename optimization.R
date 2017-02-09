##############################################
### M\'etodos de optimizaci\'on num\'erica ###
##############################################

# Newton-Raphson ----------------------------------------------------------


newton_raphson<-function(x_init,f,f_der,zero_val=0,f_args=NULL,max_iter=1e3,tolerance=1e-6,epsilon=1e-16,delta=1e-14){
  n=0
  x_ret<-x_init
  x_0<-Inf
  f_x<-do.call(f,as.list(c(x_init,f_args)))
  f_der_x<-do.call(f_der,as.list(c(x_init,f_args)))
  while(n<max_iter && abs(f_x-zero_val)>tolerance && abs(f_der_x)>epsilon && abs(x_ret-x_0)>delta){
    x_0<-x_ret
    f_x<-do.call(f,as.list(c(x_0,f_args)))
    f_der_x<-do.call(f_der,as.list(c(x_0,f_args)))
    x_ret<-x_0-(f_x-zero_val)/f_der_x
    n<-n+1
  }
  if(abs(f_x-zero_val)<=tolerance) return(x_ret)
  stop("Method failed to converge.")
}


