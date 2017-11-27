costate <-function(t,p,parameters,input) {
	with (as.list(parameters), {

    a = parameters[1]
    omega = parameters[2]


		optW <- input$optW
		x <- sapply(input$interpX, mapply, t)
		y <- sapply(input$interpY, mapply, t)
		yhat <- sapply(input$interpyHat, mapply, t)
		q <- sapply(input$q, mapply, t)

				dp1 = optW[1]*(-p[1]*(0 )-p[2]*(1 )-p[3]*(0 )-p[4]*(0 ) - 2*( 1*q[1]*(y[1]-yhat[1])+0*q[2]*(y[2]-yhat[2])+0*q[3]*(y[3]-yhat[3])+0*q[4]*(y[4]-yhat[4]) ))
				dp2 = optW[2]*(-p[1]*(0 )-p[2]*(0 )-p[3]*(1 )-p[4]*(0 ) - 2*( 0*q[1]*(y[1]-yhat[1])+1*q[2]*(y[2]-yhat[2])+0*q[3]*(y[3]-yhat[3])+0*q[4]*(y[4]-yhat[4]) ))
				dp3 = optW[3]*(-p[1]*(0 )-p[2]*(0 )-p[3]*(0 )-p[4]*(1 ) - 2*( 0*q[1]*(y[1]-yhat[1])+0*q[2]*(y[2]-yhat[2])+1*q[3]*(y[3]-yhat[3])+0*q[4]*(y[4]-yhat[4]) ))
				dp4 = optW[4]*(-p[1]*(0 )-p[2]*(0 )-p[3]*(0 )-p[4]*(0 ) - 2*( 0*q[1]*(y[1]-yhat[1])+0*q[2]*(y[2]-yhat[2])+0*q[3]*(y[3]-yhat[3])+1*q[4]*(y[4]-yhat[4]) ))

				list(c(dp1 ,dp2 ,dp3 ,dp4 ))

  })
}
