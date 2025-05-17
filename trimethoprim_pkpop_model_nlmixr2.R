setwd('D:\\GitHub\\Farmakokinetika Sulfamethoxazole')

# load library
library(nlmixr2)
library(rxode2)
library(ggplot2)
library(readxl)
library(xpose)
library(xpose.nlmixr2)

# data preparation
inputnlmixr <- read_excel('tmp.xlsx')
tmpdat <- as.data.frame(inputnlmixr)
head(tmpdat)
ggplot(tmpdat, aes(TIME, DV)) + geom_point(aes(group=ID)) +
  scale_x_continuous("Time (h)") + scale_y_continuous("Du") +
  labs(title="Trimethoprim single-dose", 
       subtitle="Du vs. time by individual")

# pkpop function
onecomp.urine <- function(){
  ini({
    tka   <- log(1.2)   # log ka
    tkel  <- log(.3)    # log kel
    tknr  <- log(.2)    # log knr
    
    #Between subject variability .000002
    eta.ka   ~ .004
    eta.kel  + eta.knr  ~ c(.1, .005, .1)
    prop.err <- 0.001 # residual variability
  })
  
  model({
    #model kovariat
    ka   <- exp(tka  + eta.ka)
    kel   <- exp(tkel  + eta.kel) 
    knr   <- exp(tknr  + eta.knr)
    
    #model struktural
    d/dt (A1) = -ka*A1                # absorption
    d/dt (A2) = ka*A1 - kel*A2        #central
    d/dt (A3) = kel*A2 - knr*A2       #urine
    A3~ prop(prop.err)
  })
}

# fitting function to the data
fit <- nlmixr2(onecomp.urine, tmpdat, est="saem")

# the base model
fit

# covariate OBAT on kel (T = with lactoferrin-containing supplement, R = without supplement )
fit1 <- fit %>%
  model({kel <- exp(tkel + eff.TREATMENT*TREATMENT + eta.kel) }) %>%
  ini(eff.TREATMENT<- 0.1) %>%
  nlmixr2(tmpdat, est="saem")

# covariate OBAT on knr 
fit2 <- fit %>%
  model({knr <- exp(tknr + eff.TREATMENT*TREATMENT + eta.knr) }) %>%
  ini(eff.TREATMENT<- 0.1) %>%
  nlmixr2(tmpdat, est="saem")

# covariate OBAT on ka
fit3 <- fit %>%
  model({ka <- exp(tka + eff.TREATMENT*TREATMENT + eta.ka) }) %>%
  ini(eff.TREATMENT<- 0.1) %>%
  nlmixr2(tmpdat, est="saem")

# addCwres
addCwres(fit)
addCwres(fit1)
addCwres(fit2)
addCwres(fit3)

# calculate objective function value
fit$objDf[[1]]
fit1$objDf[[1]]
fit2$objDf[[1]]
fit3$objDf[[1]]

# log likelihood estimation
LL_base <- logLik(fit) # base model
LL_1 <- logLik(fit1)
LL_2 <- logLik(fit2)
LL_3 <- logLik(fit3)
rbind(LL_base, LL_1, LL_2, LL_3)

# chi-square test for asses the covariate effect
p_value_1 = 1 - pchisq( LL_base[1] - LL_1[1], df=1)
p_value_2 = 1 - pchisq( LL_base[1] - LL_2[1], df=1)
p_value_3 = 1 - pchisq( LL_base[1] - LL_3[1], df=1)
rbind(p_value_1, p_value_2, p_value_3)

# fit was chosen as the best model
# diagnostic plot

xpdb <- xpose_data_nlmixr(fit)

vpc_plot<- vpcPlot(fit, n=1000, show=list(obs_dv=T),
                   ylab = "The amount of unchanged trimethoprim excreted in urine (mg)", 
                   xlab = "Median of the interval sampling time  (h)",
                   log_y=T)
vpc_plot

ggsave(filename = "vpc_plot_high_res.png", plot = vpc_plot, dpi = 400, 
       width = 10, height = 8)

print(dv_vs_pred(xpdb) + 
        ylab("The observed amount of unchanged trimetoprim excreted in urine (mg)") +
        xlab("The population-based predicted amount of unchanged trimethoprim excreted in urine (mg)"))
print(dv_vs_ipred(xpdb) + 
        ylab("observed trimethoprim Excreted (mg)") +
        xlab("Population predicted trimethoprim Excreted (mg)"))
print(res_vs_pred(xpdb, res="CWRES") + 
        ylab("Conditional Weighted Residuals") +
        xlab("population Predicted Trimethoprim Excreted (mg)"))
print(res_vs_idv(xpdb, res="CWRES") + 
        ylab("Conditional Weighted Residuals") +
        xlab("T-mid (h)"))


#bootstrap
bt <- bootstrapFit(fit, nboot=200)
fit
