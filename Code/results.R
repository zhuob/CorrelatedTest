
## the final result summary
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/simulationFunctions.R")
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/product.moments.R")
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Code/plot.rhoT.vs.rho.R")

mu1 <- c(0, 0);  
sigma <- c(1, 1)
nreps <- 5000 # number of simulated data for each case
textsize = c(20, 20, 20, 20)
rhovec <- seq(-1, 1, by = 0.01)
n1 <- n2 <- 3 # small sample size
n3 <- n4 <- 10 # larger sample size


case1 <- final_dat(mu1, mu2 = c(0, 0), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = DeltaY = 0
case2 <- final_dat(mu1, mu2 = c(0, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = 0, DeltaY = 2
case3 <- final_dat(mu1, mu2 = c(1, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = 1, DeltaY = 2
case4 <- final_dat(mu1, mu2 = c(3, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = 3, DeltaY = 2
case5 <- final_dat(mu1, mu2 = c(-1, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = -1, DeltaY = 2
case6 <- final_dat(mu1, mu2 = c(-3, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = -3, DeltaY = 2
case7 <- final_dat(mu1, mu2 = c(-5, 5), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = -5, DeltaY = 5
case8 <- final_dat(mu1, mu2 = c(5, 5), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = 5, DeltaY = 5

final_data <- bind_rows(case1, case2, case3, case4, case5, case6, case7, case8) %>% 
              mutate(setting = replace(setting, setting == "(0, 0)", "  (0, 0)"), 
                     setting = replace(setting, setting == "(0, 2)", "  (0, 2)"), 
                     setting = replace(setting, setting == "(1, 2)", "  (1, 2)"), 
                     setting = replace(setting, setting == "(3, 2)", "  (3, 2)"), 
                     setting = replace(setting, setting == "(-1, 2)", " (-1, 2)"), 
                     setting = replace(setting, setting == "(-3, 2)", " (-3, 2)"), 
                     setting = replace(setting, setting == "(-5, 5)", " (-5, 5)"), 
                     category = replace(category, category == "rhoSimu3", "  simulated, n=3"), 
                     category = replace(category, category == "rhoSimu10", " simulated, n=10"), 
                     category = replace(category, category == "rhoT3", "  theoretical, n=3"), 
                     category = replace(category, category == "rhoT10", " theoretical, n=10"), 
                     category = replace(category, category == "rhoTinf", "theoretical, asymptotic"))
saveRDS(object = final_data, file = paste("simuResults", nreps, ".rds", sep = ""))
final_data  <- readRDS(paste("simuResults", nreps, ".rds", sep = ""))


rho_figure <- ggplot(data = final_data, aes(x = rhoTrue, y = rho)) + 
          geom_line(aes(color = category, linetype = category)) + 
          facet_wrap(~setting, ncol= 4) + 
          theme(legend.position = "bottom") + 
          labs(x = "rho", y = "rhoT") + 
          guides(col = guide_legend(nrow = 2, byrow = T, keywidth = 2)) + 
          geom_abline(intercept = 0, slope=1, colour = "#000000") + ylim(-1, 1)

ggsave(paste(FigurePath,"sim.eps", sep =""), rho_figure,
       width = 10, height = 8)



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

