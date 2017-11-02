hiddenInputState  <- function(t, x, parameters, input) { 
 with (as.list(parameters),{ 

		optW <- input$optW
		w <- sapply(input$w, mapply, t)

    a = parameters[1]
    omega = parameters[2]


		dx1 = a * sin(omega * t) * cos(omega * t)+ optW[1] *w[1]
		dx2 = x[1]+ optW[2] *w[2]
		dx3 = x[2]+ optW[3] *w[3]
		dx4 = x[3]+ optW[4] *w[4]

	 	list(c(dx1 ,dx2 ,dx3 ,dx4 ))
  })
}
