library(tdaR)

## 6.1.2
# tda_cwt() missing

## 6.2.1
ds1dat <- read.table("../examples/exam/ds1.dat")
tda_dstat(ds1dat)

## 6.2.2
gdf1dat <- read.table("../examples/exam/gdf1.dat")
names(gdf1dat) <- c("SEL", "ID", "L1", "YL", "CEN")

gdf1 <- tda_gdf(~YL, data = gdf1dat[gdf1dat$SEL == 0, ], censor = "CEN")
gdf2 <- tda_gdf(~YL, data = gdf1dat[gdf1dat$SEL == 0, ], censor = "CEN", what = "survivor")
gdf3 <- tda_gdf(~YL, data = gdf1dat[gdf1dat$SEL == 0, ], censor = "CEN", what = "expected_values")


## 6.2.3
tda_quant(cars)

## 6.2.4
tda_freq(ds1dat)
tab <- tda_freq2(ds1dat$V2, ds1dat$V1)

## 6.2.5
atab <- tda_atab(ds1dat[c("V3", "V4")], breaks = seq(-10, 10, 0.5))
atab$table

## 6.2.6
tda_corr(ds1dat)

## 6.2.7
tab <- tda_freq2(ds1dat$V2, ds1dat$V1, contingency = TRUE)
tab

## 6.2.8
gen <- tda_rng()
n <- 200
X <- Y <- numeric(n)
for (i in seq_len(n)) {
    X[i] <- gen$rd(a = 0, b = 3)

    Y[i] <- sin(X[i]) + gen$rd()
}
df <- data.frame(x = X, y = Y)

### scatter plot
p <- tda_ps(data = df, xlim = c(0,3), ylim = c(-1,2.5))
p <- tda_pl_scatter(p, "x", "y", symbol = 5, size = 1)
p <- tda_pl_axes(p, sc = c(1), ic = c(10, 0))
plot(p)

# # ## FIXME: throws error
# p <- tda_ps(data = df, xlim = c(0,3), ylim = c(-1,2.5))
# # # Error: TDA could not draw this: Syntax error or undefined variables.
# p <- tda_pl_scatter(p, x = df$X, y = df$Y, symbol = 5, size = 1)
# p <- tda_pl(p, "plxa", sc = 1, ic = 10)
# p <- tda_pl(p, "plya", sc = 1, ic = 0)
# plot(p)

### Sunflower plot
## random dot at axis intersection with sunflower plot
p <- tda_ps(data = df, xlim = c(0,3), ylim = c(-1,2.5))
p <- tda_pl_scatter(p, "x", "y", type = "sunflower", grid = c(15, 10), size = 3)
p <- tda_pl_axes(p, sc = c(1), ic = c(10, 1))
plot(p)

### lowess regression
p <- tda_ps(data = df, xlim = c(0,3), ylim = c(-1,2.5))
# scatter does not care about xlim, ylim and the no clip command does not work?
p <- tda_pl_scatter(p, "x", "y", symbol = 5,
                    size = 1, type = "lowess")
p <- tda_pl(p, "plxa", sc = 1, ic = 10)
p <- tda_pl(p, "plya", sc = 1, ic = 0)
plot(p)


## 6.2.9
dh1dat <- read.table("../examples/exam/dh1.dat")
p <- tda_ps(data = dh1dat, xlim = c(0,8), ylim = c(0,0.4))
p <- tda_pl(p, "plxa", sc = 1)
p <- tda_pl(p, "plya", sc = .1)
p <- tda_pl_histogram(p, x = "V1", breaks = c(0,2,4,6,8))
plot(p)

## 6.2.10
p <- tda_ps(data = dh1dat, xlim = c(0,8), ylim = c(0,0.6))
p <- tda_pl(p, "plxa", sc = 1)
p <- tda_pl(p, "plya", sc = .1)
p <- tda_pl_density(p, x = "V1", kernel = "triangle", at = seq(0, 8, 0.1))
plot(p)


## 6.4.2
tda_ineq(ds1dat)

## 6.4.3
ds2dat <- read.table("../examples/exam/ds2.dat")
names(ds2dat) <- c("G", paste0("X", 1:4))
tda_segr(ds2dat, group = "G")

## 6.5.1
rrdat1 <- tda_rrdat()
rrdat1$DES <- ifelse(rrdat1$TFin == rrdat1$TI, 0, 1)
rrdat1$DUR <- rrdat1$TFin - rrdat1$TStart + 1

ltb <- tda_ltb(Surv(DUR, DES) ~1, data = rrdat1, tp = seq(0, 500, 30))
# box 3
ltb$summary
# box 4
ltb$table
tda_survivor(ltb)

rrdat2 <- read.table("../examples/exam/rrdat.2")
names(rrdat2) <- c("DUR", "CEN", "WT", "GRP")
fit <- tda_ple(Surv(DUR, CEN) ~ as.factor(GRP), data = rrdat2, weights = "WT")
fit$blocks

ple <- tda_ple(Surv(TFP, DES) ~ 1, data = tda_rrdat(states = 4))
ple$summary


srv <- tda_survivor(ple)
p <- tda_ps(srv, width = 80, height = 40, xlim = c(0, 300), ylim = c(0, 1))
p <- tda_pl(p, "plxa", sc = 60, ic = 5)
p <- tda_pl(p, "plya", sc = 1, ic = 10)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "0,1", lty = 1)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "0,2", lty = 5)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "0,3", lty = 8)
plot(p)

ple <- tda_ple(Surv(TFP, DES) ~ 1, data = tda_rrdat(states = 4), quantiles = seq(0.1, .9, by = 0.1))
ple$quantiles

## 6.5.4
pl2 <- tda_ple(Surv(TFP, DES) ~ as.factor(SEX), tda_rrdat(), compare = TRUE)
sb <- tda_survivor(pl2)
sb <- sb[sb$time < 290, ]
p <- tda_ps(sb, width = 90, height = 50, xlim = c(0,300), ylim = c(0,1))
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "1",
                  band = c("upper", "lower"), lty = 1)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "2",
                            band = c("upper", "lower"), lty = 5)
p <- tda_pl(p, "plxa", sc = 60, ic = 1)
p <- tda_pl(p, "plya", sc = 1, ic = 10)
p <- tda_pl_frame(p)
plot(p)

## 6.6.1
seqd1 <- read.table("../examples/exam/seq.d1")
names(seqd1) <- c("ID", paste0("Y", seq_len(ncol(seqd1)-1L)))

tda_seqlg(seqd1, id = "ID")

## 6.6.2
tda_seqgc(seqd1, id = "ID")

## 6.6.3
tda_seqsd(seqd1[, -1])
tda_seqen(seqd1[,-1])

## 6.6.4
tda_seqpm(seqd1[,-1L], patterns = list("-", c(3, 3), c(3, "*", 3)))

## 6.7.2.3
seqmd1 <- read.table("../examples/exam/seqm.d1")
seqmd1a <- read.table("../examples/exam/seqm.d1a")
tda_seqm(seqmd1)
s <- tda_seqm(seqmd1, print = "sequential")
attributes(s)
s <- tda_seqm(seqmd1, print = "lcs")
attributes(s)


seqmd3 <- read.table("../examples/exam/seqm.d3")

scost <- matrix(0, 10, 10)
for (i in 1:10) for (j in 1:10) scost[i, j] <- 0.1 * abs(i - j)

r <- tda_seqm(seqmd3, indel = 1, subcost = scost, dp_matrix = TRUE)
attr(r, "dp_matrix")[["2_1"]]$D

## 6.7.2.4
seqmd4 <- read.table("../examples/exam/seqm.d4")
r <- tda_seqm(seqmd4, dp_matrix = TRUE)
r <- tda_seqm(seqmd4, indel = c(1, 0.5), dp_matrix = TRUE)
r <- tda_seqm(seqmd4, indel = c(1, 0), dp_matrix = TRUE)

## 6.7.2.5
r <- tda_seqm(seqmd4, subcost = 2, dp_matrix = TRUE)


## 6.9.1
lsreg1dat <- read.table("../examples/exam/lsreg1.dat")
names(lsreg1dat) <- c("Height", "Weight")
reg <- tda_lsreg(Weight ~Height, data = lsreg1dat)
summary(reg)

## 6.9.1.2
lsreg2dat <- read.table("../examples/exam/lsreg2.dat")
names(lsreg2dat) <- c("X1", "X2", "Y")

reg <- tda_lsreg(Y ~ 0 + X1 + X2, data = lsreg2dat,
                 equality   = "b1 + b2 = 6",
                 inequality = c("b1 = 2", "b2 = 2"))

## 6.9.1.3
lsreg3dat <- read.table("../examples/exam/lsreg3.dat")
names(lsreg3dat) <- c("X1", "G", "Y", "NE", "NC", "SO", "WE")
lsreg3dat$Y2 <- lsreg3dat$Y * lsreg3dat$Y

region <- c("NE", "NC", "SO", "WE")
reg <- tda_lsreg(G ~ Y + Y2, data = lsreg3dat, dgroup = region)
summary(reg)


## 6.10.1
l1reg4dat <- read.table("../examples/exam/l1reg4.dat")
names(l1reg4dat) <- c("Y", "X1", "X2")

tda_l1reg(Y ~X1 + X2, data = l1reg4dat)

## 6.10.2
gen <- tda_rng()
n <- 256
X <- E <- numeric(n)
for (i in seq_len(n)) {
    X[i] <- gen$rd()
    E[i] <- 0.32 * gen$rdn()
}
Z <- sin(2 * pi * X * X * X)
Y <- Z * Z * Z + E
npregdat <- data.frame(X = X, Y = Y)

np <- tda_npreg(Y ~ X, data = npregdat, method = "mean", kernel = "quartic",
                bandwidth = 0.1, x = seq(0, 1, 0.01))

p <- tda_ps(npregdat, width = 90, height = 50, xlim = c(0, 1),
            ylim = c(-1.2, 1.8))
p <- tda_pl_axes(p, sc = c(0.1, 0.5), ic = 0)
p <- tda_pl_points(p, "X", "Y", symbol = 5, size = 0.5, lty = 0)
p <- tda_pl_function(p, "sin(2*pi*x*x*x)*sin(2*pi*x*x*x)*sin(2*pi*x*x*x)",
                     range = c(0, 1), step = 0.01, lty = "dashed")
p <- tda_pl_lines(p, x = np$table$X, y = np$table$YM, lty = "solid", lw = 1)

plot(p)



## 6.11.2
dfile <- tda_rrdat()
dfile$DUR <- dfile$TFin - dfile$TStart + 1L
dfile$W <- as.integer(dfile$SEX == 2)

ml <- tda_fml({
    rate = exp(a0 + COHO2 * a1 + COHO3 * a2 + W * a3)
    l1 = ifelse(DES, log(rate), 0)
    fn = l1 - rate * DUR
}, data = dfile, start = c(-4, 0, 0, 0))
ml


## 6.12.2
qr1dat <- read.table("../examples/exam/qr1.dat")
names(qr1dat) <- c("Dose", "Weight", "Response")
qr1dat$Log10Dose <- log(qr1dat$Dose) / log(10)

# logit
logit <- tda_qreg(Response ~ Log10Dose, data = qr1dat, weights = "Weight", predictions = TRUE, standardized = TRUE)
logit
logit$standardized
logit$predictions

# probit
probit <- tda_qreg(Response ~ Log10Dose, data = qr1dat, weights = "Weight", model = 2)


## 6.12.2.1
qr5dat <- read.table("../examples/exam/qr5.dat")
names(qr5dat) <- c("Y", "N", "S", "X")
qr5dat$SX = qr5dat$S * qr5dat$X

## TODO: Signif has just 5 digits precision, can we tweak this too?
# get start values for random effects
fml <- tda_fml(
    {
        xb = beta0 + S * beta1 + X * beta2 + SX * beta3
        ee = bc(N,Y) * (exp(xb)^Y) / ((1 + exp(xb))^N)
        fn = log(ee)
    },
    data = qr5dat,
    options = list(mfmt=24.16, tfmt=24.16, pfmt=24.16)
)

re <- tda_fml(
    {
        xb = beta0 + S * beta1 + X * beta2 + SX * beta3
        ee = intn(7,bc(N,Y) * (exp(xb+t*gam)^Y)/((1 + exp(xb+t*gam))^N) )
        fn = log(ee)
    },
    data = qr5dat,
    start = c(coef(fml), gam = .1),
    options=list(mina = 4,mfmt=24.16, tfmt=24.16)
)

## 6.12.3
qr2dat <- read.table("../examples/exam/qr2.dat",
                     col.names = c("Y", "X1", "X2"))

# ordinal logit
qr_ol <- tda_qreg(Y ~ X1 + X2, data = qr2dat, model = 3, predictions = TRUE)
qr_ol
qr_ol$predictions

# ordinal probit
qr_pb <- tda_qreg(Y ~ X1 + X2, data = qr2dat, model = 4)
qr_pb

fml_ol <- tda_fml(
    {
        xb = X1 * beta1 + X2 * beta2
        theta0 = 1
        theta1 = eexp(alpha1 + xb)
        theta2 = eexp(alpha2 + xb)
        theta3 = eexp(alpha3 + xb)
        ll = ifelse(Y == 0, theta0 - theta1,
                    ifelse(Y == 1, theta1 - theta2,
                           ifelse(Y == 2, theta2 - theta3,
                                  theta3
                           )))
        fn = log(ll)
    },
    data = qr2dat,
    start = c(-2, -3, -4, 0, 0),
    options = list(mfmt=24.16, tfmt=24.16)
)
fml_ol

## 6.12.4
qr3dat <- read.table("../examples/exam/qr3.dat",
                     col.names = c("X", "Weight", "Y"))

tda_qreg(Y ~ X, data = qr3dat, weights = "Weight", model = 5, nq = 5)

ml <- tda_fml(
    {
        xb2 = exp(beta20 + X * beta21)
        xb3 = exp(beta30 + X * beta31)
        xb4 = exp(beta40 + X * beta41)
        xb5 = exp(beta50 + X * beta51)
        xbb = 1 + xb2 + xb3 + xb4 + xb5
        ll = ifelse(Y == 1, 1 / xbb,
                    ifelse(Y == 2, xb2 / xbb,
                           ifelse(Y == 3,xb3 / xbb,
                                  ifelse(Y == 4, xb4 / xbb, xb5 / xbb))))
        fn = Weight * log(ll)
    },
    data = qr3dat
)
ml

qr4dat <- read.table("../examples/exam/qr4.dat",
                     col.names = c("Z1", "Z2", "Z3", "Y"))

# cross section
cml <- tda_qreg(Y ~ 0 + lvl(Z1, Z2, Z3), data = qr4dat, model = 5, nq = 3, predictions = TRUE)
cml
head(cml$predictions, 8)

## 6.12.5
qr3dat <- read.table("../examples/exam/qr3.dat",
                     col.names = c("X", "Weight", "Y"))

tda_qreg(Y ~ X, data = qr3dat, model = 6, nq = 5,
         constraints = c("b9 = 0", "b10 = 0", "b11 = 0", "b12 = 0", "b13 = 0",
                         "b14 = 0", "b15 = 0", "b16 = 0", "b17 = 0", "b18 = 0"),
         weights = "Weight")


qr4dat <- read.table("../examples/exam/qr4.dat",
                     col.names = c(paste0("Z", 1:3), "Y"))
fit <- tda_qreg(Y ~ 0 + lvl(Z1, Z2, Z3), data = qr4dat, model = 6, nq = 3,
                start = c(-0.1716, 0, 0, 0), constraints = c("b3 = 0", "b4 = 0"))
fit


## 6.14.1
cd1dat <- read.table("../examples/exam/cd1.dat",
                     col.names = c("NDI", "Service", "B", "C", "D", "E", "C60", "C65", "C70", "P75"))
cd1dat$LOGS <- log(cd1dat$Service)
cd1dat <- cd1dat[cd1dat$NDI >= 0, ]

tda_fml({
    xb = b0 + B * bb + C * bc + D * bd + E * be +
        C60 * bc60 + C65 * bc65 + C70 * bc70 + P75 * bp75 + LOGS * blogs
    exb = exp(xb)
    fn = poisson(exb,NDI)
}, data = cd1dat)

tda_fml({
    xb = b0 + B * bb + C * bc + D * bd + E * be +
        C60 * bc60 + C65 * bc65 + C70 * bc70 + P75 * bp75 + LOGS * blogs
    exb = exp(xb)
    fn = poisson(exb,NDI)
}, data = cd1dat, constraints = "b9 = 1")


## 6.15.2.1
lsreg1dat <- read.table("../examples/exam/lsreg1.dat",
                        col.names = c("Height", "Weight"))

tda_glm(Weight ~ Height, data = lsreg1dat)

gen <- tda_rng()
n <- 100
X <- Y <- numeric(n)
for (i in seq_len(n)) {
    X[i] <- gen$rd(a = -3, b = 3)

    Y[i] <- exp(3 + X[i]) + gen$rd()
}

glm2dat <- data.frame(X, Y)

tda_glm(Y ~ X, data = glm2dat, start = c(5, 5), link = 2)


## 6.15.2.2
qr1dat <- read.table("../examples/exam/qr1.dat")
names(qr1dat) <- c("Dose", "Weight", "Response")
qr1dat$Log10Dose <- log(qr1dat$Dose) / log(10)

tda_glm(Response ~ Log10Dose, data = qr1dat, weights = "Weight", family = 2)

glm3dat <- read.table("../examples/exam/glm3.dat", col.names = c("D", "N", "Y"))
glm3dat$D <- log(glm3dat$D)

## TODO trials could be just "N"?
fit <- tda_glm(Y ~ D, data = glm3dat, family = 2, link = 3, trials = glm3dat$N)
fit

fit <- tda_glm(Y ~ D, data = glm3dat, trials = glm3dat$N, family = 2, custom_link = "log(mue / (1 - mue))")
fit


## 6.15.2.3
glm5dat <- read.table("../examples/exam/glm5.dat",
                      col.names = c("X", "Y"))

tda_glm(Y ~ X, data = glm5dat, family = "poisson")


cd1dat <- read.table("../examples/exam/cd1.dat",
                     col.names = c("NDI", "Service", "B", "C", "D", "E", "C60", "C65", "C70", "P75"))
cd1dat$LOGS <- log(cd1dat$Service)
cd1dat <- cd1dat[cd1dat$NDI >= 0, ]

fit <- tda_glm(NDI ~ B + C + D + E + C60 + C65 + C70 + P75 + LOGS, data = cd1dat, family = "poisson",
               equality = c("b9 = 1"), control = tda_control(maxit = 100))

fit


## 6.15.2.4
gen <- tda_rng()
n <- 1000
X <- Y <- numeric(n)
for (i in seq_len(n)) {
    X[i] <- gen$rd()

    Y[i] <- -exp(1 + X[i]) * (log(gen$rd()) + log(gen$rd()))/2
}

glm7dat <- data.frame(X, Y)

## Box 1
fit <- tda_glm(Y ~ X, data = glm7dat, family = "Gamma",
               link = 2, start = c(1, 1), control = tda_control(maxit = 50))
fit

X <- Y <- numeric(n)
for (i in seq_len(n)) {
    X[i] <- gen$rd()

    Y[i] <- -(log(gen$rd()) + log(gen$rd())) / ((1 + X[i]) * 2)
}

glm7datb <- data.frame(X, Y)
fit <- tda_glm(Y ~ X, data = glm7datb, family = "Gamma",
               start = c(1, 1), control = tda_control(maxit = 50))
fit

## Box 4
glm8adat <- read.table("../examples/exam/glm8a.dat",
                       col.names = c("U", "Y", "L"))
glm8adat$X <- log(glm8adat$U)
glm8adat$LX <- glm8adat$L * glm8adat$X

fit <- tda_glm(Y ~ X, data = glm8adat[glm8adat$L == 0,], family = "Gamma",
               start = c(-.02, .02))
fit

fit <- tda_glm(Y ~ X, data = glm8adat[glm8adat$L == 1,], family = "Gamma",
               start = c(-.02, .02))
fit

## 6.15.2.5
gen <- tda_rng()
n <- 1000
X <- Mu <- Z <- X1 <- Y <- numeric(n)
for (i in seq_len(n)) {
    X[i] <- gen$rd()
    Mu[i] <- 1 + X[i]
    Z[i] <- gen$rdn() ^ 2
    X1[i] <- Mu[i] + Mu[i] ^ 2 * Z[i] - Mu[i] * sqrt(2 * Mu[i] * Z[i] + Mu[i] ^ 2 * Z[i] ^ 2)
    Y[i] <- ifelse(gen$rd() <= (Mu[i] / (Mu[i] + X1[i])), X1[i], Mu[i] ^ 2 / X1[i])
}

glm9dat <- data.frame(
    X, Mu, Z, X1, Y
)

fit <- tda_glm(Y ~ X, data = glm9dat, family = "inverse.gaussian", link = 1, start = c(.5, .5))
fit

summary(fit)

## 6.16.1
source("../R_checks/freg1dat.R", echo = TRUE)
fit <- tda_freg(
    {
        ax = a1 * X1 + a2 * X2 + a4 * exp(a3 * X3)
        fn = (Y - ax) ^ 2
    },
    data = freg1dat,
    start = c(-.1, 1, -1, -1),
    control = tda_control(algorithm = 5)
)

fit

## 6.17.2

rrdat1 <- tda_rrdat()
rrdat1$W <- as.integer(rrdat1$SEX == 2)

## 6.17.2.1
fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 2, residuals = TRUE,
                prate = list(tp = seq(0, 100, 5), COHO3 = 1, W = 1))
head(fit$residuals, 4)
fit$rates

# m
rrdat1_4 <- tda_rrdat(states = 4)
rrdat1_4$W <- as.integer(rrdat1_4$SEX == 2)
rrdat1_4$DES <- 2L
# is TDA using rounded values in nvar()?
rrdat1_4$DES[rrdat1_4$PRESN / rrdat1_4$PRES - 1L >= 0.19999999] <- 1L
rrdat1_4$DES[rrdat1_4$PRESN / rrdat1_4$PRES - 1L < 0.0] <- 3L
rrdat1_4$DES[rrdat1_4$TFin == rrdat1_4$TI] <- 0L

fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1_4, model = 2, relative_risk = TRUE)
fit$relative_risk

## 6.17.2.2
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 3, tp = seq(0, 96, 12))

# c
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 3, tp = seq(0, 96, 12),
         constraints = c("b1 - b2 = 0", "b2 - b3 = 0", "b3 - b4 = 0",
                         "b4 - b5 = 0", "b6 - b7 = 0", "b8 - b9 = 0"))

## 6.17.2.3
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 16, tp = seq(0, 96, 12))

## 6.17.3.1
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 4, degree = 2)
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 5, degree = 2)

## 6.17.3.2 (I cannot say that I understand this, is this some kind of SUR?)
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 6, on = c("xb", "xc"))

## 6.17.3.3
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 7)

## 6.17.3.4
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 8)

## 6.17.3.5
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 9)

## 6.17.3.6
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 12)

## 6.17.3.7
tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 13)

## 6.17.3.8
fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1, model = 14)


## 6.17.5.1
rrdat1 <- tda_rrdat()
rrdat1$W <- as.integer(rrdat1$SEX == 2)

fit <- tda_frml(Surv(TFP, DES) ~ COHO2 + COHO3 + W, {
    rate = exp(a0 + COHO2 * a1 + COHO3 * a2 + W * a3)
    l1 = ifelse(DES,log(rate),0)
    fn = l1 - rate * TFP
},
data = rrdat1,
start = c(-4, 0, 0, 0)
)
fit

rrdat1 <- tda_rrdat()
rrdat1$W <- as.integer(rrdat1$SEX == 2)

## box 3
fit <- tda_frml(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    {
        rate1 = exp(a10 + COHO2 * a11 + COHO3 * a12 + W * a13)
        rate2 = exp(a20 + COHO2 * a21 + COHO3 * a22 + W * a23)
        rate3 = exp(a30 + COHO2 * a31 + COHO3 * a32 + W * a33)
        rate = rate1 + rate2 + rate3
        l1 = ifelse(DES == 1, log(rate1),
                    ifelse(DES == 2, log(rate2),
                           ifelse(DES == 3, log(rate3),0)))
        fn = l1 - rate * TFP
    },
    data = rrdat1,
    start = c(-4, 0, 0, 0,
              -4, 0, 0, 0,
              -4, 0, 0, 0),
    control = tda_control(maxit = 100)
)
fit

## box 4
fit <- tda_frml(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    {
        bb = COHO2 * b1 + COHO3 * b2 + W * b3
        p1 = tf >= 0 & tf < 12
        p2 = tf >= 12 & tf < 24
        p3 = tf >= 24 & tf < 36
        p4 = tf >= 36 & tf < 48
        p5 = tf >= 48 & tf < 60
        p6 = tf >= 60 & tf < 72
        p7 = tf >= 72 & tf < 84
        p8 = tf >= 84 & tf < 96
        p9 = tf >= 96
        rr = ifelse(p1, a1,
                    ifelse(p2, a2,
                           ifelse(p3, a3,
                                  ifelse(p4, a4,
                                         ifelse(p5, a5,
                                                ifelse(p6, a6,
                                                       ifelse(p7, a7,
                                                              ifelse(p8, a8,
                                                                     ifelse(p9, a9, 0)))))))))
        s1 = 12 * exp(a1 + bb)
        s2 = s1 + 12 * exp(a2 + bb)
        s3 = s2 + 12 * exp(a3 + bb)
        s4 = s3 + 12 * exp(a4 + bb)
        s5 = s4 + 12 * exp(a5 + bb)
        s6 = s5 + 12 * exp(a6 + bb)
        s7 = s6 + 12 * exp(a7 + bb)
        s8 = s7 + 12 * exp(a8 + bb)
        lsurv = ifelse(p1, tf * exp(a1 + bb),
                       ifelse(p2, s1 + (tf - 12) * exp(a2 + bb),
                              ifelse(p3, s2 + (tf - 24) * exp(a3 + bb),
                                     ifelse(p4, s3 + (tf - 36) * exp(a4 + bb),
                                            ifelse(p5, s4 + (tf - 48) * exp(a5 + bb),
                                                   ifelse(p6, s5 + (tf - 60) * exp(a6 + bb),
                                                          ifelse(p7, s6 + (tf - 72) * exp(a7 + bb),
                                                                 ifelse(p8, s7 + (tf - 84) * exp(a8 + bb),
                                                                        ifelse(p9, s8 + (tf - 96) * exp(a9 + bb), 0)))))))))
        l1 = ifelse(DES, rr + bb, 0)
        fn = l1 - lsurv
    },
    data = rrdat1,
    start = c(rep(-4, 9), 0, 0, 0),
    control = tda_control(maxit = 50)
)
fit

## box 5
fit <- tda_frml(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    {
        aa = exp( a0 + COHO2 * a1 + COHO3 * a2 + W * a3 )
        bb = exp( b0 )
        l1 = ifelse(DES,log(bb * aa^bb * tf^(bb - 1)),0)
        fn = l1 - (aa * tf)^bb
    },
    data = rrdat1,
    start = c(-4, 0, 0, 0, 0),
    control = tda_control(maxit = 100)
)
fit

## box 6
fit <- tda_frml(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    {
        aa = a0 + COHO2 * a1 + COHO3 * a2 + W * a3
        bb = exp( b0 )
        zz = (log(tf) - aa) / bb
        qz = 1 - nd(zz)
        vz = ndf(zz) / qz
        l1 = ifelse(DES,log(vz / (bb * tf)),0)
        fn = l1 + log(qz)
    },
    data = rrdat1,
    start = c(4, 0, 0, 0, 0),
    control = tda_control(maxit = 100)
)
fit


## box 7
fit <- tda_frml(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    {
        kk = 1
        aa = a0 + COHO2 * a1 + COHO3 * a2 + W * a3
        bb = exp(b0)
        zz = (log(tf) - aa) / bb
        qt = kk * exp(zz / sqrt(kk))
        rr = (kk - 0.5) * log(kk) + (sqrt(kk) * zz - qt) - (b0 + log(tf) + lgam(kk))
        icc = icg(qt,kk)
        fn = ifelse(DES,rr,log(1 - icc))
    },
    data = rrdat1,
    start = c(5, 0, 0, 0, 0),
    control = tda_control(maxit = 100, algorithm = 4)
)
fit

## box 8
fit <- tda_frml(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    {
        rate = exp(a0 + COHO2 * a1 + COHO3 * a2 + W * a3)
        l1 = ifelse(des,log(rate),0)
        fn = l1 - rate * (tf - ts)
    },
    data = rrdat1,
    start = c(-4, 0, 0, 0),
    control = tda_control(maxit = 100)
)
fit

## 6.17.5.2
fit <- tda_frml(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    {
        sigma = exp(s)
        lambda = exp(l)/(1 + exp(l))
        beta = exp(b)
        ff = sigma * exp(beta * lambda^tf)
        ff1 = ff + 1
        dens = -(ff * beta * lambda^tf * log(lambda)) / (ff1 * ff1)
        surv = 1 - 1 / ff1
        fn = ifelse(des,log(dens),log(surv))
    },
    data = rrdat1,
    control = tda_control(maxit = 100, algorithm = 4)
)
fit

## 6.17.5.3
fit <- tda_frml(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    {
        bg = 0.604
        sigma = exp(s) / (1 + exp(s))
        lambda = exp(l)
        mue = exp(m)
        beta = bg * lambda
        tl = exp(-lambda * (tf - mue))
        ldens = log(sigma) + l - beta * (tf - mue) - tl - lgam(bg)
        surv = 1 - sigma * (1 - icg(tl,bg))
        fn = ifelse(DES,ldens,log(surv))
    },
    data = rrdat1, start = c(-3, 2),
    control = tda_control(maxit = 100, algorithm = 4)
)
fit

## 6.17.6.2
tda_rate(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    data = rrdat1, model = 20
)

## 6.17.6.3
tda_rate(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    data = rrdat1, model = 21
)

## 6.17.7.1
tda_rate(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    data = rrdat1, model = 1
)

## box 4
tda_rate(
    Surv(TFP, DES) ~ COHO2 + COHO3 + W,
    data = rrdat1_4, model = 1
)

## 6.17.7.2
rrdat1 <- tda_rrdat()
rrdat1$W <- as.integer(rrdat1$SEX == 2)
# MDATE is static (doesn't reference `time`), so it can be pre-computed
# directly in R, matching the manual's own nvar()-level formula exactly:
rrdat1$MDATE <- ifelse(rrdat1$TMAR <= 0, 10000, rrdat1$TMAR - rrdat1$TStart)

fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W + MARR, data = rrdat1,
                model = 1,
                define = list(MARR = "gt(time,MDATE)"),
                helpers = "MDATE")
fit

## 6.17.7.3 - we dont have this

## 6.17.7.4
rrdat1 <- tda_rrdat()
rrdat1$W <- as.integer(rrdat1$SEX == 2)

fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W + WTEST, data = rrdat1,
                model = 1,
                define = list(WTEST = "W * (log(time) - 4.22)"))

## box 3

rrdat5 <- read.table("../examples/exam/rrdat.5", col.names = c("TF", "DES", "G"))
rrdat5$G1 <- as.integer(rrdat5$G == 1)
rrdat5$G2 <- as.integer(rrdat5$G == 2)

tda_rate(Surv(TF, DES) ~ G1, data = rrdat5, tp = c(0, 170.1, 354.1, 535.1), model = 1)

## box 5
fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1,
                model = 1, tp = seq(0, 96, 12))
fit$gof

fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rrdat1_4,
                model = 1, tp = seq(0, 96, 12))
fit$gof

## 6.17.7.5
rrdat1 <- tda_rrdat()
fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + strata(SEX),
                data = rrdat1, model = 1)

## 6.17.7.6
rrdat1 <- tda_rrdat()
rrdat1$W <- as.integer(rrdat1$SEX == 2)
## prate not supported?
fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W,
                data = rrdat1, model = 1, prate = list(COHO3 = 1))
fit <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W,
                data = rrdat1, model = 1, prate = list(COHO3 = 1, W = 1))
fit$rates

## 6.18.2
seqd4 <- read.table("../examples/exam/seq.d4")
names(seqd4) <- c("ID", paste0("Y", 0:5), paste0("S", 0:5), "V1", "V2")

multi <- list(paste0("Y", 0:5), paste0("S", 0:5))
# box 3
tda_seq_info(multi, data = seqd4)
tda_seqev(multi, data = seqd4)

# box 4
tda_seqevd(multi, data = seqd4, sn = 2)

## 6.18.3

# Box 1
tda_seqmd(paste0("Y", 0:5), event = c(1, 2), data = seqd4, summary = TRUE)

# Box 2
tda_seqmd(paste0("Y", 0:5), event = c(1, 2), data = seqd4)

# Box 3
tda_seqmd(paste0("Y", 0:5), event = c(1, 2), data = seqd4,
          covariates = c("V1", "V2"))

## 6.19.1
ll1dat <- read.table("../examples/exam/ll1.dat", col.names = c("X1", "X2", "H"))
# TODO why 5 coefficients?
fit <- tda_loglin(~ X1 + X2, data = ll1dat, weights = "H")
fit$table

## 6.19.2
ll2dat <- read.table("../examples/exam/ll2.dat", col.names = c("X1", "X2", "F", "X4"))
fit <- tda_loglin(~ X1 + X2, data = ll2dat, weights = "F", residuals = TRUE)
fit
coef(fit)
residuals(fit)

