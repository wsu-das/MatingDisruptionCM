# This code emulates mating disruption effects of codling moth phenology models described in Pest Manag Sci 2021; 77: 1081–1093
# and implemented in the pesticide evaluator tool in DAS


# Functions to convert form C to F

source("./functions/FDD_CDD.R")
source("./functions/f_to_c.R")
source("./functions/calc_dd_vec.R")

require(ExtDist)
require(SuppDists)


# For illustration purposes, temperature data from the Pullman weather station during the growing season of 2023 will be used.

pullman <- read.csv("./data/pullman.csv")

# Convert F to C
pullman$maxC <- f_to_c(pullman$max)
pullman$minC <- f_to_c(pullman$min)
pullman$DDC <- FDD_CDD(pullman$DDF)

# Calculate daily degree-days

pullman$dds <- calc_dd_vec(tmax = pullman$maxC, tmin = pullman$minC, upper_threshold = 31, lower_threshold = 10, cutoff = "horizontal")

# Calculate daily delay in degree-days. The assumption is that there is a three-day delay in mating with MD.

a <- rep(NA, 100)
for(i in 1:100) {
  a[i] <- sum(pullman$dds[i:(i + 2)])
}

pullman$delD <- c(0, 0, a[1:98])

# Calculate proportion of larvae per calendar day using cumulative degree-days. Parameters of the JohnsonSB pdf from 
# Pest Manag Sci 2021; 77: 1081–1093 ("First summer egg hatch")

pullman$dlarva <- dJohnsonSB(pullman$DDC, 
                             params = list(gamma = 1.092, delta = 1.376, xi = 173.590, 
                                           lambda = 644.440))

# Function that delivers the proportion of no delay control. This was extracted from Fig 2 in Popul Ecol (2012) 54:421–429
# The original idea was to use the product of the estimated rm and lx, but these values do not match those in Fig 2.

mod_md1 <- function(x) {
  a = 0.372937
  b = 0.029193
  exp(a * (1 - exp(b*x)))
}

plot(seq(1, 100), pullman$dlarva, ylim = c(0, 0.0045), type = "l", lwd = 2, 
     ylab = "Proportion of larvae", xlab = "Calendar days", cex.lab  = 1.5,
     cex.axis = 1.5, lty = 2)
lines(seq(1, 100), pullman$dlarva * mod_md1(pullman$delD), col = "red", lwd = 2)

# Which produces a plot that is similar to that produced by the website for the first generation of larvae (Picture2.jpg attached). 
# For subsequent generations, just add the corresponding JohnsonSB.
# One thing that I would ask Vince is the reason for using the degree-day delays reported for the calendar days 
# when larvae hatched instead of using those reported when the adults emerged.  