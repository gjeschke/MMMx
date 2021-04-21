function prob = prob_Gaussian(rsim,meanr,sigr)

Garg = (meanr - rsim)/sigr;
prob = exp(-Garg^2);