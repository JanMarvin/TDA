
library(tdaR)

## 8.3.2
tda_evalfi("sin(x)", x = c(0, 1))

## 8.4.1
tda_gmin("sin(x)", start = list(c(1, 0, 10)), use_derivatives = TRUE)$run

## 8.4.2
tda_range("sin(x)", start = list(c(1, 0, 10)), use_derivatives = TRUE)

## 8.5.1
id1dat <- read.table("../examples/exam/id1.dat")
tda_sddf(id1dat, self_consistent = TRUE)

id1adat <- read.table("../examples/exam/id1a.dat")
tda_iddf(~iv(V1, V2), data = id1adat, self_consistent = TRUE)

## 8.5.2
id2dat <- read.table("../examples/exam/id2.dat")
id2dat$XL <- id2dat$V1 - 30L
id2dat$XH <- id2dat$V1 + 30L
tda_idf(~iv(XL, XH), data = id2dat)

## 8.6.1
tda_imean(~iv(XL, XH), data = id2dat)

## 8.6.2
tda_ivariance(~iv(V1, V2), data = id1adat)
