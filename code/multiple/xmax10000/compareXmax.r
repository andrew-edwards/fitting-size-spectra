# compareXmax.r - comparing sample of 1,000 random PLB numbers
#  for xmax=1,000 and xmax=10,000 with same seed, to demonstrate
#  that have to be careful because the resulting samples are quite
#  similar (due to inverse method, and hardly any chance of being
#  >1,000 with xmax=10,000).
#  24th November 2015.

rm(list=ls())

source("../PLBfunctions.r")

n = 1000                  # sample size
b.known = -2              # known fixed value of b
xmin.known = 1            # known fixed value of xmin
# xmax.known = 10000         # known fixed value of xmax

set.seed(42)              # Usual seed.
x1000 = rPLB(n, b = b.known, xmin = xmin.known, xmax = 1000)

set.seed(42)
x10000 = rPLB(n, b = b.known, xmin = xmin.known, xmax = 10000)

ratio = x1000/x10000
print(ratio)

print(summary(ratio))

print(quantile(ratio, seq(0.1, 1, 0.1)))

print(sum(ratio > 0.99))

print(x1000[1])
print(x10000[1])
print(ratio[1])

print(pPLB(1000, b=-2, xmin=1, xmax=10000)) 
