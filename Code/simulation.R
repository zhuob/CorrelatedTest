source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/simulationFunctions.R")
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/product.moments.R")
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/plot.rhoT.vs.rho.R")

mu1 <- c(0, 0);  sigma <- c(1, 1)
nreps <- 1000 # number of simulated data for each case
textsize = c(20, 20, 20, 20)
rhovec <- seq(-1, 1, by = 0.01)

wd <- 5
ht <- 7
FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Manuscript/Figures/"
n1 <- n2 <- 3 # small sample size
n3 <- n4 <- 10 # larger sample size


mu2 <- c(0, 0);
res1_1 <- compare_rho(mu1, mu2, rhovec, sigma, n1, n2, nreps)
res1_2 <- compare_rho(mu1, mu2, rhovec, sigma, n3, n4, nreps)
res1_3 <- sapply(rhovec, theore_rho,delta1 = mu2[1], delta2= mu2[2], sigma=sigma, n1=n1, n2=n2)
res1 <- comb1(res1_1, res1_2, res1_3)

mu2 <- c(0, 2)
res2_1 <- compare_rho(mu1, mu2, rhovec, sigma, n1, n2, nreps)
res2_2 <- compare_rho(mu1, mu2, rhovec, sigma, n1, n2, nreps)
res2_3 <- sapply(rhovec, theore_rho,delta1 = mu2[1], delta2= mu2[2], sigma=sigma, n1=n1, n2=n2)
res2 <- comb1(res2_1, res2_2, res2_3)

mu2 <- c(0.5, 2)
res3_1 <- compare_rho(mu1, mu2, rhovec, sigma, n1, n2, nreps)
res3_2 <- compare_rho(mu1, mu2, rhovec, sigma, n3, n4, nreps)
res3_3 <- sapply(rhovec, theore_rho,delta1 = mu2[1], delta2= mu2[2], sigma=sigma, n1=n1, n2=n2)
res3 <- comb1(res3_1, res3_2, res3_3)

mu2 <- c(1, 2)
res4_1 <- compare_rho(mu1, mu2, rhovec, sigma, n1, n2, nreps)
res4_2 <- compare_rho(mu1, mu2, rhovec, sigma, n3, n4, nreps)
res4_3 <- sapply(rhovec, theore_rho,delta1 = mu2[1], delta2= mu2[2], sigma=sigma, n1=n1, n2=n2)
res4 <- comb1(res4_1, res4_2, res4_3)



mu2 <- c(3, 2)
res5_1 <- compare_rho(mu1, mu2, rhovec, sigma, n1, n2, nreps)
res5_2 <- compare_rho(mu1, mu2, rhovec, sigma, n3, n4, nreps)
res5_3 <- sapply(rhovec, theore_rho,delta1 = mu2[1], delta2= mu2[2], sigma=sigma, n1=n1, n2=n2)
res5 <- comb1(res5_1, res5_2, res5_3)



mu2 <- c(-3, 2)
res6_1 <- compare_rho(mu1, mu2, rhovec, sigma, n1, n2, nreps)
res6_2 <- compare_rho(mu1, mu2, rhovec, sigma, n3, n4, nreps)
res6_3 <- sapply(rhovec, theore_rho,delta1 = mu2[1], delta2= mu2[2], sigma=sigma, n1=n1, n2=n2)
res6 <- comb1(res6_1, res6_2, res6_3)


smaller <- n1
larger <- n3

f1 <- plot_rho1(res1, textsize= textsize, include_invS = F, 
                case ="a", smaller = smaller, larger = larger)
f2 <- plot_rho1(res2, textsize= textsize, include_invS = F, 
                case ="b", smaller = smaller, larger = larger)
f3 <- plot_rho1(res3, textsize= textsize,
                case ="c", smaller = smaller, larger = larger)
f4 <- plot_rho1(res4, textsize= textsize,
                case ="d", smaller = smaller, larger = larger)
f5 <- plot_rho1(res5, textsize= textsize,
                case ="e", smaller = smaller, larger = larger)
f6 <- plot_rho1(res6, textsize= textsize,
                case ="f", smaller = smaller, larger = larger)

ggsave(paste(FigurePath,"case1.eps", sep =""), f1,
       width = wd, height = ht)
ggsave(paste(FigurePath,"case2.eps", sep =""), f2,
       width = wd, height = ht)
ggsave(paste(FigurePath,"case3.eps", sep =""), f3,
       width = wd, height = ht)
ggsave(paste(FigurePath,"case4.eps", sep =""), f4,
       width = wd, height = ht)
ggsave(paste(FigurePath,"case5.eps", sep =""), f5,
       width = wd, height = ht)
ggsave(paste(FigurePath,"case6.eps", sep =""), f6,
       width = wd, height = ht)



delta1_vector <- seq(-4, 4, 0.1)
sigma <- c(1, 1)
n1 <- n2 <- 1000
th1 <- theoretical_plot(delta1_vector, rho= 0.7, sigma = sigma, n1=n1, n2=n2)  
ggsave(paste(FigurePath,"th1.eps", sep =""), th1,
       width = wd, height = ht)
th2 <- theoretical_plot(delta1_vector, rho= 0.1, sigma = sigma, n1=n1, n2=n2)  
ggsave(paste(FigurePath,"th2.eps", sep =""), th2,
       width = wd, height = ht)
th3 <- theoretical_plot(delta1_vector, rho= 0.05, sigma = sigma, n1=n1, n2=n2)  
ggsave(paste(FigurePath,"th3.eps", sep =""), th3,
       width = wd, height = ht)
th4 <- theoretical_plot(delta1_vector, rho= -0.1, sigma = sigma, n1=n1, n2=n2)  
ggsave(paste(FigurePath,"th4.eps", sep =""), th4,
       width = wd, height = ht)


#FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/ThesisRelated/Dissertation/Figures/"
wd <- 5
ht <- 5
ct1 <- contour_plot(rho = 0.7, sigma= sigma, n1=n1, n2= n2, brks = seq(0, 0.7, by = 0.1))
ggsave(paste(FigurePath,"ct1.eps", sep =""), ct1,
       width = wd, height = ht)
ct2 <- contour_plot(rho = -0.7, sigma= sigma, n1=n1, n2= n2, brks = seq(-0.7, -0.1, by = 0.1))
ggsave(paste(FigurePath,"ct2.eps", sep =""), ct2,
       width = wd, height = ht)
ct3 <- contour_plot(rho = 0.1, sigma= sigma, n1=n1, n2= n2,brks = seq(0, 0.1, by = 0.01))
ggsave(paste(FigurePath,"ct3.eps", sep =""), ct3,
       width = wd, height = ht)
ct4 <- contour_plot(rho = -0.1, sigma= sigma, n1=n1, n2= n2, brks = seq(-0.1, 0, by = 0.01))
ggsave(paste(FigurePath,"ct4.eps", sep =""), ct4,
       width = wd, height = ht)





