
ND = 1000 # number of simulation replicates
N = 300 # number of subjects
K = 9 # number of biomarkers; excluding MMScore and thik
L = 2 # number of mixture component
Ni = 10
gap = 2 # one visit gap year
lambda1 = 0.35; lambda = c(lambda1, 1-lambda1)

rho=c(1.3, 0.8)

gamma1 = c(2,2,       1,1,1,2,            2,2,1)
gamma2 = c(3,3,       4,4,4,4,            2,3,4)
gamma3 = c(1,0.8,     -0.2,-0.2,0,0.2,    -2,-1,-1)
beta0 = c(0.2,0.1,    0,-0.3,-0.4,-0.3,   -1,-0.5,-0.1)
betat = c(-0.1,0.2,   0.2,0.1,0.2,0.4,    0.1,0.1,0.1)
betaX1 = c(-0.3,-0.4, 0,0.2,0.4,0.3,      0,0,0 ) # gender
betaX2 = c(-0.1,-0.2, 0.1,0,0,0,          0,0,0) # education
# var and sd
sig1 = c(0.8,0.7,     0.9,0.8,0.6,0.3,    1,1,1)
sig2 = c(0.8,0.7,     0.6,0.8,0.6,0.6,    0.04,0.04,0.5)
sqsig1 = sqrt(sig1)
sqsig2 = sqrt(sig2)
sqsig = rbind(sqsig1, sqsig2)

alphat = 0.4
alpha0 = 2