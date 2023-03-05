# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: CALIBRACION DE ORIFICIO

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# //////////////////////////////////////////////////////////////////////////////////

# INFO:
# Cuantificacion_incertidumbre
# Propagacion_error
# Modelos_nls_(Nonlinear_Least_Squares)
# Analisis_grafico_ggplot2
# //////////////////////////////////////////////////////////////////////////////////

# Scientific notation is suppress
options(scipen = 0)

# Workspace is cleared
rm(list = ls())

# Working directory is selected
# setwd("/media/maikel/Trabajo/R_ITC/R_LABHYD/EXP_ORIFICIO")
setwd("C:/DATOS/R_ITC/R_LABHYD/EXP_ORIFICIO")

# CRAN libraries are loaded
# require(Agreement)
require(DescTools)
require(effects)
require(ggplot2)
require(ggalt)
require(MASS)
require(nls2)
require(nlstools)
require(pastecs)
require(reshape)
require(visreg)

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round data.frame to specif digits
# /////////////////////////////////////////////////////////////
round_df <- function(df, digits) {
  options(scipen = 0)
  options(scipen = -2)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# ////////////////////////////////////////////////////////
# BLOCK: Custom Function {as.lm.nls}
# ////////////////////////////////////////////////////////

as.lm.nls <- function(object, ...) {
  if (!inherits(object, "nls")) {
    w <- paste("expected object of class nls but got object of class:",
               paste(class(object), collapse = " "))
    warning(w)
  }
  
  gradient <- object$m$gradient()
  if (is.null(colnames(gradient))) {
    colnames(gradient) <- names(object$m$getPars())
  }
  
  response.name <- if (length(formula(object)) == 2) "0" else
    as.character(formula(object)[[2]])
  
  lhs <- object$m$lhs()
  L <- data.frame(lhs, gradient)
  names(L)[1] <- response.name
  
  fo <- sprintf("%s ~ %s - 1", response.name,
                paste(colnames(gradient), collapse = "+"))
  fo <- as.formula(fo, env = as.proto(   as.list(L)   ))
  
  do.call('lm', list(fo, offset = substitute(fitted(object))))
}
# ////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////
# BLOCK: CD Comparison
# ////////////////////////////////////////////////////////

# Input data is loaded and a data.frame is created
df.base <- read.table("test_orificio.txt", header = TRUE)

# Desc {DescTools} function is requested
Desc(df.base, plotit = TRUE)

# names {base} function is requested
names(df.base)

# A simple deviation (SD) column is created
df.base$SD <- (df.base$CD_TEO - df.base$CD_EXP)

# shapiro.test {stats} Normality Test is applied to df.base$SD
# if p-value > 0.05 then normality stands true, meaning that
# the variable is parametric
shapiro.test(df.base$SD)

# Descriptive statistics are requested and rounded to 5 decimals
df.base.desc <- round(stat.desc(df.base[, 2:8]),5)

# A correlation test is executed for X vs Y
cor.test(df.base$CD_TEO,df.base$CD_EXP)

# A ggplot2 object is created
fg01 <- ggplot(aes(x = CD_TEO,y = CD_EXP),data=df.base) +
  geom_abline(intercept = 0, slope = 1, colour = "gray", size = 1) +
  geom_point(aes(colour = GROUP, shape = GROUP),
             size = 5.5, alpha = 0.95) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Comparaci?n entre coeficientes teoricos y experimentales del orificio") +
  xlab("CD_TEO (-)") +
  ylab("CD_EXP (-)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg01

# melt {reshape} function is requested to convert data from "wide"
# to "long" format but extracting the first 5 columns
df.base.long <- melt(df.base[,6:8])

# Third column (reynolds) is repeated 3 times and a data.frame is created
df.reynolds <- as.data.frame(rep(df.base[,3],3))

# df.reynolds data.frame names are changed
names(df.reynolds) <- c("Reynolds")

# Selected columns are repeated "n" times 
rep1 <- as.character(rep(df.base[,1],1))
rep2 <- rep(c("no_aplica"),length(df.base[,1]))
rep3 <- as.character(rep(df.base[,1],1))

# Repetitions are c-binded
rep.chain <- c(rep1, rep2, rep3)

# A data.frame is created
df.group <- as.data.frame(rep.chain)

# df.group data.frame names are changed
names(df.group) <- c("GROUP")

# cbind {base} function is used to join data.frames
df.base.long <- cbind(df.group, df.reynolds, df.base.long)

# A ggplot object is created
fg02 <- ggplot(data = df.base.long, aes(x = Reynolds, y = value)) +
  geom_point(aes(colour = variable, shape = GROUP),size = 6.0,alpha = 0.95) +
  geom_line(aes(colour = variable,group = variable), size = 0.75) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  # geom_text(aes(label = value),parse = FALSE) +
  ggtitle("Comportamiento de Coef.CD en virtud del numero de Reynolds") +
  xlab("Numero de Reynolds (-)") +
  ylab("Coeficiente de Descarga CD (-)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg02          

# A new "CLASS" character column is created in data.frame data=df.base
df.base$CLASS <- c("Deviation")

# subsetting data.frame
df.base.boxplot <- subset(df.base.long, variable == "CD_EXP" | variable == "CD_TEO")

# A ggplot2 object is created
fg03 <- ggplot(aes(y = value,x = variable),data=df.base.boxplot) +
  geom_boxplot(aes(colour = variable),
               outlier.colour = '#ff0000',outlier.size = 4) +
  geom_point(aes(shape = GROUP),size = 5.5,
             position = position_jitter(width = 0.1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("Boxplot - CD Teorico + CD experimental") +
  xlab("Clase") +
  ylab("CD (-)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg03

# ////////////////////////////////////////////////////////
# BLOCK: Models
# ////////////////////////////////////////////////////////

# A ggplot2 object is created
fg04 <- ggplot(aes(x = DELTA_P_m,y = Q_LPS),data=df.base) +
  geom_smooth(size = 0.85, alpha = 0.35, method = lm, fullrange = TRUE) +
  geom_point(aes(colour = GROUP, shape = GROUP),size = 6.5, alpha = 0.90) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Modelo de calibracion del orificio (lm)") +
  xlab("Delta presion (m)") +
  ylab("Caudal Q (LPS)") +
  theme_bw(base_size = 22.0)

# A ggplot2 object is requested
fg04

# A lm {stats} Linear Model is created
mod.lm <- lm(Q_LPS ~ DELTA_P_m, data=df.base)

# A summary function is requested
summary (mod.lm)

# /////////////////////////////////////////////////
# lm SECTIONS:
# /////////////////////////////////////////////////
# 
# Residuals section:
# it provides a quick summary (min, 1Q, median, 3Q, max) of the distribution.
# 
# Coefficients section: each coefficient is a Gaussian random variable
# Estimate represents the mean distribution of the variable
# Std.Error displays the standard error of the variable
# the t value is Estimate divided by Std.Error
# the p value indicates the probability of getting a value larger than the t value
#
# Residual standard error outputs the standard deviation of residuals
#
# The degree of freedom indicates the differences between the observation 
# in training samples and the number used in the model
#
# Multiple R-squared is obtained by dividing the sum of squares.
#
# Adjusted R-squared uses an unbiased estimate, and will 
# be slightly less than multiple R-squared

# /////////////////////////////////////////////////

# model residuals are requested and its SD is calculated
mod.lm$residuals
sd(mod.lm$residuals)

# shapiro.test {stats} Normality Test is applied
# if p-value > 0.05 then normality stands true, meaning that
# the variable is parametric
shapiro.test(mod.lm$residuals)

# visreg {visreg} function is requested
visreg(mod.lm, partial = TRUE)

# A model plot is requested
# plot(mod.lm)

# A full lm-model diagnostic plot is requested
oldpar.lm <- par(oma=c(0,0,3,0), mfrow=c(2,2))
plot(mod.lm)
par(oldpar.lm)

# RESIDUALS VERSUS FITTED VALUES
# Residuals represent the vertical distance from a point to the regression line
# If all points fall exactly on the regression line, all residuals will fall exactly on the dotted gray line.
# The red line within the plot is a smooth curve with regard to residuals
# [SUGGESTS HOMOGENEITY,LINEARITY, CONSTANT VARIANCE (VARIABLES INDEPENDANCE)]
#
# QQ-PLOT (NORMAL OF RESIDUALS)
# This plot verifies the assumption that residuals were normally distributed.
# Thus, if the residuals were normally distributed, they should lie exactly on the gray dash line
# [SUGGESTS NORMALITY]
#
# SCALE-LOCATION PLOT
# It measures the square root of the standardized residuals against the fitted value
# Therefore, if all dots lie on the regression line, the value of y should be close to zero.
# Since it is assumed that the variance of residuals does not change the distribution substantially,
# if the assumption is correct, the red line should be relatively flat.
# [SUGGESTS THAT VARIANCE HOMOGENEITY HOLDS TRUE]
#
# STANDARDIZED RESIDUALS VERSUS LEVERAGE
# The leverage is a measurement of how each data point influences the regression.
# It is a measurement of the distance from the centroid of regression and level of isolation
# [SUGGESTS THAT VARIANCE IS SENSITIVE TO OUTLIERS]
# 
# COOK'S DISTANCE
# It measures how regression would change if a single point is deleted.
# Cook's distance is affected by HIGH leverage and LARGE residuals.
# For a perfect fit regression, the red line should be close to the dashed
# line with no points over 0.5 in Cook's distance contour.

# A rlm {MASS} Linear Model is created (Robust Fitting of Linear Models)
mod.rlm <- rlm(Q_LPS ~ DELTA_P_m, data=df.base)

# A summary function is requested
summary(mod.rlm)

# model residuals are requested
mod.rlm$residuals

# shapiro.test {stats} Normality Test is applied
# if p-value > 0.05 then normality stands true, meaning that
# the variable is parametric
shapiro.test(mod.rlm$residuals)

# A full rlm-model diagnostic plot is requested
oldpar.rlm <- par(oma=c(0,0,3,0), mfrow=c(2,2))
plot(mod.rlm)
par(oldpar.rlm)

# A nls {stats} Nonlinear Least Squares model is fitted
# f(x) = a + b*Ln(x); Remember, in R Ln() = log !!!
# parameters a and b initial values must be defined
mod.nls <- nls(Q_LPS ~ a + b*(log(DELTA_P_m)),
               data=df.base,
               start=list(a=1,b=1))

# A summary function is requested
summary(mod.nls)

# A model plot is requested
plot(mod.nls)

#------------------------------------------
# R2 Comparison
#------------------------------------------

# Residual sum of squares is calculated for each model
RSS.p1 <- sum(residuals(mod.lm)^2)
RSS.p2 <- sum(residuals(mod.rlm)^2)
RSS.p3 <- sum(residuals(mod.nls)^2)

# Total sum of squares is calculated regardless of the model
TSS <- sum((df.base$Q_LPS - mean(df.base$Q_LPS))^2)


# R-squared is calculated for all models
R2.lm <- (1 - (RSS.p1/TSS))
R2.rlm <- (1 - (RSS.p2/TSS))
R2.nls <- (1 - (RSS.p3/TSS))

# R2 is requested for all models
print(R2.lm)
print(R2.rlm)
print(R2.nls)

# Based on R2 results, nls is selected

# CIs are calculated using custom function {as.lm.nls}
df.predCI <- as.data.frame(predict(as.lm.nls(mod.nls),
                                   interval = 'confidence',
                                   level = 0.95))

# df.predCI data.frame is rounded to 3 decimals
df.predCI <- round(df.predCI,3)

# cbind {base} function is used to join data.frames
df.output <- cbind(df.base,df.predCI)

# confint2 {nlstools} function is requested to 
# confidence intervals in nonlinear regression coefficients
confint2(mod.nls)

# A ggplot object is created
fg05 <- ggplot() +
  geom_point(aes(x = DELTA_P_m,y = Q_LPS, colour = GROUP, shape = GROUP),
             data=df.output, size = 6.5) +
  geom_line(aes(x = DELTA_P_m,y = fit),
            data=df.output,colour = '#0000ff',size = 0.95,alpha = 0.9) +
  geom_line(aes(x = DELTA_P_m,y = lwr),
            data=df.output,colour = '#666666',size = 0.95,linetype = 2,alpha = 0.9) +
  geom_line(aes(x = DELTA_P_m,y = upr),
            data=df.output,colour = '#666666',size = 0.95,linetype = 2,alpha = 0.9) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Modelo de calibracion del orificio (nls)") +
  xlab("Delta presion (m)") +
  ylab("Caudal Q (LPS)") +
  theme_bw(base_size = 22.0)

# A ggplot object is requested
fg05

# round_df function is applied to relevant data.frames
df.base.desc <- round_df(df=df.base.desc, digits=3)
df.output <- round_df(df=df.output, digits=3)

# Objects to export
# Desc(2REYNOLDS, 4DELTA_P_m, 5CD_EXP)
# shapiro.test(df.base$SD), df.base.desc, cor.test(df.base$CD_TEO,df.base$CD_EXP),
# fg01, fg02, fg03, fg04, fg05 
# summary (mod.lm)
# full lm-model diagnostic plot
# summary(mod.rlm)
# full rlm-model diagnostic plot
# summary(mod.nls)
# R2.lm, R2.rlm, R2.nls
# df.output, df.base.desc
# shapiro.test(mod.lm$residuals)
# shapiro.test(mod.rlm$residuals)

write.csv(df.base.desc, file = "df.base.desc.csv")
write.csv(df.output, file = "df.output.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////
