costate <-function(t,p,parameters,input) {
	with (as.list(parameters), {

    alpha1 = parameters[1]
    alpha2 = parameters[2]
    alpha3 = parameters[3]
    beta1 = parameters[4]
    beta2 = parameters[5]


		x <- sapply(input$interpX, mapply, t)
		y <- sapply(input$interpY, mapply, t)
		yhat <- sapply(input$interpyHat, mapply, t)
		q <- sapply(input$q, mapply, t)

				dp1= -p[1]*(-alpha1 )-p[2]*(beta1 )-p[3]*(beta2 ) - 2*( 1*q[1]*(y[1]-yhat[1])+0*q[2]*(y[2]-yhat[2])+0*q[3]*(y[3]-yhat[3]) )
				dp2= -p[1]*(0 )-p[2]*(-alpha2 )-p[3]*(0 ) - 2*( 0*q[1]*(y[1]-yhat[1])+1*q[2]*(y[2]-yhat[2])+0*q[3]*(y[3]-yhat[3]) )
				dp3= -p[1]*(0 )-p[2]*(0 )-p[3]*(-alpha3 ) - 2*( 0*q[1]*(y[1]-yhat[1])+0*q[2]*(y[2]-yhat[2])+1*q[3]*(y[3]-yhat[3]) )

				list(c(dp1,dp2,dp3))

  })
}
