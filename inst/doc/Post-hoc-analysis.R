## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning = FALSE, message=FALSE-------------------------------
library(postHoc)

## ------------------------------------------------------------------------
str(DeIdentifiedExample)

## ------------------------------------------------------------------------
table(DeIdentifiedExample$Treatment)

## ------------------------------------------------------------------------
FULLmodel <- lm(Y ~ Treatment + 0, data = DeIdentifiedExample)
coef(FULLmodel)

## ------------------------------------------------------------------------
NULLmodel <- lm(Y ~ 1, data = DeIdentifiedExample)
anova(NULLmodel, FULLmodel) 

## ------------------------------------------------------------------------
TT <- posthoc(Model = FULLmodel, EffectLabels = LETTERS[1:7], digits = 1)
summary(TT)

## ------------------------------------------------------------------------
print(TT)

## ---- fig.height = 5, fig.width = 6, fig.align ="center"-----------------
barplot(TT, ylim = c(5,21))
abline(h = 5)

## ---- fig.height = 5, fig.width = 6, fig.align ="center"-----------------
lines(TT, ylim = c(5,21))

## ------------------------------------------------------------------------
ZZ <- posthoc(Model = FULLmodel, EffectIndices = c(2,3,4,5,7), digits = 2)

## ------------------------------------------------------------------------
ZZ$PvaluesMatrix

## ---- fig.height = 5, fig.width = 6, fig.align ="center"-----------------
set.seed(143)
plot(ZZ)

## ------------------------------------------------------------------------
summary(ZZ)

## ---- fig.height = 5, fig.width = 6, fig.align ="center"-----------------
set.seed(1243)
plot(TT)

