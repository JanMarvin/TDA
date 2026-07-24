# Independence probe: how many readers no longer need TDA's printed
# output or its written files at all.
#
# The console output is emptied and every written file deleted AT
# CONSTRUCTION TIME, before the reader sees the run.  A reader that
# still produces its result was BUILT from exports; one that errors or
# returns NULL still routes through a parse, even if every value in its
# result came from an export afterwards.
#
# This is a different measure from tools/phase3_scan.R.  That one asks
# whether the NUMBERS are the doubles TDA computed (they are, 94/94).
# This one asks whether the READER can work without the text.  Both are
# needed: precision and provenance are not the same claim.
#
#     TDA_EXAMPLES=... Rscript tools/independence_probe.R
options(digits = 17, warn = -1)
suppressMessages(library(tdaR))
EX <- Sys.getenv("TDA_EXAMPLES", "examples")
ns <- asNamespace("tdaR"); orig <- get("tda_run", ns)
strip <- function(res) { res$output <- character(0)
                         unlink(list.files(res$dir, full.names = TRUE)); res }
unlockBinding("tda_run", ns); assign("tda_run", function(...) strip(orig(...)), ns)

ok <- err <- character(0)
t_ <- function(nm, f) {
  r <- try(suppressWarnings(f()), silent = TRUE)
  if (inherits(r, "try-error") || is.null(r)) err <<- c(err, nm) else ok <<- c(ok, nm)
}
d <- read.table(file.path(EX,"exam","lsreg1.dat")); names(d) <- c("Height","Weight")
ds1 <- read.table(file.path(EX,"exam","ds1.dat"))
rr <- tda_rrdat(); rr$W <- as.integer(rr$SEX == 2)
qr1 <- read.table(file.path(EX,"exam","qr1.dat")); names(qr1) <- c("Dose","Weight","Response")
qr1$L <- log(qr1$Dose)/log(10)
ll2 <- read.table(file.path(EX,"exam","ll2.dat"), col.names=c("X1","X2","F","X4"))
sq <- data.frame(Y0=c(1,1,2),Y1=c(1,2,1),Y2=c(2,1,1),Y3=c(1,1,2),Y4=c(1,1,2),Y5=c(3,3,3))

t_("dstat", function() tda_dstat(d));      t_("quant", function() tda_quant(d))
t_("corr", function() tda_corr(ds1));      t_("cov", function() tda_cov(ds1))
t_("freq", function() tda_freq(ds1));      t_("freq2", function() tda_freq2(ds1$V2, ds1$V1))
t_("atab", function() tda_atab(ds1["V3"], breaks=seq(-10,10,0.5)))
t_("ineq", function() tda_ineq(ds1))
t_("lsreg", function() tda_lsreg(Weight~Height, data=d))
t_("glm", function() tda_glm(Weight~Height, data=d))
t_("l1reg", function() tda_l1reg(Weight~Height, data=d))
t_("rate", function() tda_rate(Surv(TFP,DES)~COHO2, data=rr, model=2))
t_("qreg", function() tda_qreg(Response~L, data=qr1, weights="Weight"))
t_("loglin", function() tda_loglin(~X1+X2, data=ll2, weights="F"))
t_("ple", function() tda_ple(Surv(TFP,DES)~1, data=rr))
t_("ltb", function() tda_ltb(Surv(TFP,DES)~1, data=rr, tp=seq(0,500,30)))
t_("dple", function() tda_dple(Surv(TFP,DES)~1, data=rr))
t_("dltb", function() tda_dltb(Surv(TFP,DES)~1, data=rr))
t_("seqgc", function() tda_seqgc(sq)); t_("seqlg", function() tda_seqlg(sq))
t_("seqsd", function() tda_seqsd(sq)); t_("seqen", function() tda_seqen(sq))
t_("seqsi", function() tda_seqsi(sq, tp="0(1)5"))
t_("seq_info", function() tda_seq_info(list(names(sq)), data=sq))
t_("seqm", function() tda_seqm(sq))
t_("npreg", function() tda_npreg(Y~X, data=data.frame(X=seq(0,1,length.out=40), Y=sin(seq(0,1,length.out=40))), method="mean", bandwidth=0.2, x=seq(0,1,0.2)))
t_("segr", function() tda_segr(data.frame(G=c(0,0,1,1),X1=c(1,2,1,2)), group="G"))
t_("mds", function() tda_mds(dist(matrix(c(1,2,3,8,9,10), ncol=1))))

cat(sprintf("independent of output+files: %d of %d\n", length(ok), length(ok)+length(err)))
cat("  OK :", paste(ok, collapse=", "), "\n")
cat("  DEP:", paste(err, collapse=", "), "\n")
