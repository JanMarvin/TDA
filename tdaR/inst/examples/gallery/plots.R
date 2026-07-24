dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plots.cf using only documented tda_pl_* API.
# TDA's special-character text: the "@nnn"/"\nnn" escape syntax
# plots.cf itself uses directly, matching exactly what the shipped
# .cf file contains -- TDA writes these through two different fonts
# (FT, Times-Roman re-encoded via a table embedded in every TDA file's
# own prolog; FS, the standard Adobe Symbol font), and this package's
# own PostScript reader now decodes both correctly rather than
# falling back to the raw byte value regardless of which font drew
# it. Checked against the real, unmodified plots.ps: the
# same "@141 = 2.5 \245 @142" text this uses decodes to "alpha = 2.5
# (bullet) beta" through the real TDA binary, matching this script's
# own output.
library(tdaR)

p <- tda_ps(xlim = c(0, 10), ylim = c(0, 5), width = 120, height = 60)
p <- tda_pl_axes(p, sc = 1, ic = c(0, 2))
p <- tda_pl_frame(p)
p <- tda_pl_text(p, "@141 = 2.5 \\245 @142", at = c(1, 4))
p <- tda_pl_text(p, "pltext (xy=1,4) = \\324\\100141 = 2.5 \\\\245 \\100142\\324",
                 at = c(4, 4))
p <- tda_pl_text(p, "@123 = @173 @161 @174 @161 @316 @121 @175", at = c(1, 2))
p <- tda_pl_text(p, paste("pltext (xy=1,2) = \\324\\100123 = \\100173 \\100161",
                          "\\100174 \\100161 \\100316 \\100121 \\100175\\324"),
                 at = c(2.2, 0.9))
p <- tda_pl_text(p, "Sm\\277rebr\\277d", at = c(1, 3))
p <- tda_pl_text(p, "pltext (xy=1,3) = \\324Sm\\\\277rebr\\\\277d\\324", at = c(4, 3))
p <- tda_pl_lines(p, x = c(3.5, 2.7), y = c(4.05, 4.05), a = "1,1")
p <- tda_pl_lines(p, x = c(3.5, 2.7), y = c(3.05, 3.05), a = "1,1")
p <- tda_pl_arc(p, x = c(1.5, 3.2, 3.2), y = c(1.5, 1.5, 2.1), curvature = -90)
p <- tda_pl_lines(p, x = c(3.2, 2.9), y = c(2.1, 2.1), a = "1,1")
p <- tda_pl_arc(p, x = c(1.5, 1.5), y = c(1.5, 0.95), curvature = -90)
p <- tda_pl_lines(p, x = c(1.5, 2), y = c(0.95, 0.95))
png("out/plots_r.png", width = 700, height = 420)
plot(p)
dev.off()
