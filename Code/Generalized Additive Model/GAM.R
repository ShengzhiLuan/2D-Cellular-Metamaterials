#  GAM for 2D Cellular Metamaterials
#  Johns Hopkins University
#  Shengzhi Luan
#  02.08.2023

library(mgcv)
df = read.csv("D:/Johns Hopkins University/Research/2023/January-March/Machine Learning/GAM/Structure-Property-Refined-1P-Data.csv")
mod_lm <- gam((Relative.stiffness) ~ Relative.density + s(Mean.of.nodal.connectivity) + 
                                     s(Variance.of.strut.orientation) + s(Kurtosis.of.strut.orientation), 
                                     data = df, family = gaussian(link = "identity"))
summary(mod_lm)
plot(mod_lm, scheme = 3)
