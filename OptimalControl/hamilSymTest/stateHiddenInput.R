hiddenInputState  <- function(t, x, parameters, input) { 
 with (as.list(parameters),{ 

		w <- sapply(input, mapply, t)

    alpha1 = parameters[1]
    alpha2 = parameters[2]
    alpha3 = parameters[3]
    beta1 = parameters[4]
    beta2 = parameters[5]


		dx1= -alpha1 * x[1]+ t+ w[1]
		dx2= -alpha2 * x[2]+ beta1 * x[1]+ w[2]
		dx3= -alpha3 * x[3]+ beta2 * x[1]+ w[3]

	 	list(c(dx1,dx2,dx3))
  })
}
