
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
FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/Manuscript/Figures/"


case1 <- final_dat(mu1, mu2 = c(0, 0), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = DeltaY = 0
case2 <- final_dat(mu1, mu2 = c(0, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = 0, DeltaY = 2
case3 <- final_dat(mu1, mu2 = c(1, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = 1, DeltaY = 2
case4 <- final_dat(mu1, mu2 = c(3, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = 3, DeltaY = 2
case5 <- final_dat(mu1, mu2 = c(-1, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = -1, DeltaY = 2
case6 <- final_dat(mu1, mu2 = c(-3, 2), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = -3, DeltaY = 2
case7 <- final_dat(mu1, mu2 = c(-5, 5), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = -5, DeltaY = 5
case8 <- final_dat(mu1, mu2 = c(5, 5), rhovec, sigma, n1, n2, n3, n4, nrep) ## DeltaX = 5, DeltaY = 5

final_data <- bind_rows(case1, case2, case3, case4, case5, case6, case7, case8) 
#saveRDS(object = final_data, file = paste("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/simuResults", nreps, ".rds", sep = ""))


final_data  <- readRDS(paste("/Users/Bin/Google Drive/Study/Thesis/Correlation/CorrelatedTest/simuResults", nreps, ".rds", sep = ""))


final_data <- final_data %>% 
                mutate(type = ifelse(grepl("Simu", category)>0, "simulated", " theoretical"), 
                       samplesize =  case_when(category %in% c("rhoSimu3", "rhoT3")~"n= 3", 
                                               category %in% c("rhoSimu10", "rhoT10")~"n=10", 
                                               category == "rhoTinf"~" asymptotic"),
                       setting = case_when(setting == "(0, 0)"~"  (0, 0)",
                                           setting == "(0, 2)"~"  (0, 2)",
                                           setting == "(1, 2)"~"  (1, 2)",
                                           setting == "(3, 2)"~"  (3, 2)",
                                           setting == "(-1, 2)"~" (-1, 2)",
                                           setting == "(-3, 2)"~" (-3, 2)",
                                           setting == "(-5, 5)"~" (-5, 5)",
                                           setting == "(5, 5)"~"(5, 5)")) %>% 
                arrange(type, samplesize, setting)
              


final_data$setting <- factor(final_data$setting, 
                             labels = c(substitute("("~delta[x]~","~delta[y]~")"~"="~"(0, 0)"),
                                        substitute("("~delta[x]~","~delta[y]~")"~"="~"(0, 2)"), 
                                        substitute("("~delta[x]~","~delta[y]~")"~"="~"(1, 2)"),
                                        substitute("("~delta[x]~","~delta[y]~")"~"="~"(3, 2)"),
                                        substitute("("~delta[x]~","~delta[y]~")"~"="~"(-1, 2)"),
                                        substitute("("~delta[x]~","~delta[y]~")"~"="~"(-3, 2)"),
                                        substitute("("~delta[x]~","~delta[y]~")"~"="~"(-5, 5)"),
                                        substitute("("~delta[x]~","~delta[y]~")"~"="~"(5, 5)")))

final_data$type <- factor(final_data$type, labels = c("theoretical", "simulated"))



rho_figure <- ggplot(data = final_data, aes(x = rhoTrue, y = rho)) + 
          geom_line(aes(color = samplesize, linetype = type)) + 
          facet_wrap(~setting, ncol= 2, labeller = label_parsed) + 
          theme(legend.position = "bottom") + 
          labs(x = expression(rho), y = expression(rho[T])) + 
          guides(colour = guide_legend(nrow = 1, byrow = T, title.position ="left", title ="" )) + 
          geom_abline(intercept = 0, slope=1, colour = "darkgrey") + 
          scale_linetype(guide = FALSE) +
          ylim(-1, 1) + 
          scale_color_manual(breaks = c("asymptotic", "n= 3", "n=10"),
                                      values=c("black", "magenta", "darkgreen"),
                             labels = expression(paste(n[1], "=", n[2], "=10"), 
                                                 paste(n[1], "=", n[2], "=3"),
                                                 paste("theoretical"))) 
     

rho_figure
ggsave(paste(FigurePath,"sim.eps", sep =""), rho_figure,
       width = 5, height = 9)



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

