# load the wmwpowR package
library(WebPower)

# Set the parameters for the power calculation
f <- 0.5     # Desired effect size (Cohen's f)
alpha <- 0.05 # Significance level
power <- 0.8  # Desired power
n <- 21       # Number of participants
ng <-1        # Number of groups
nm <- 3         # number of measurements same as: k <- 3        # Number of within-subjects levels
nscor	 <- 1      # Nonsphericity correction coefficient. The nonsphericity correction coefficient is a measure of the degree of sphericity in the population. A coefficient of 1 means sphericity is met, while a coefficient less than 1 means not met. The samller value of the coefficient means the further departure from sphericity. The lowest value of the coefficient is 1/(nm-1) where nm is the total number of measurements. Two viable approaches for computing the empirical nonsphericity correction coefficient are sggested. One is by Greenhouse and Geisser (1959), the other is by Huynh and Feldt (1976).
type <- 1    # The value "0" is for between-effect; "1" is for within-effect; and "2" is for interaction effect.

# Conduct power calculation
wp.rmanova(ng = ng, nm = nm, f = f, nscor = nscor,
           alpha = alpha, power=power, type = type)
