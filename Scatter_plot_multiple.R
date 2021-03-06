library(gplots)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(magrittr)
Base_respi_ordenada$LGMD <- factor(Base_respi_ordenada$LGMD, levels = c("0", "1", "2", "3", "4"), labels = c("CAPN3", "ANO5", "DYSF", "FKRP", "SGC"))
p <- ggplot(data = Base_respi_ordenada, aes(x = Diseasetime, y = sitting_FVC, color = LGMD)) + geom_point() + scale_y_continuous(name="FVC% sitting") + geom_hline(yintercept=80, col="grey") + scale_x_continuous(name="Time of disease progression") + geom_hline(yintercept=80, col="grey") + geom_smooth(method = "lm", se=TRUE)
g <- p + facet_wrap(~LGMD, nrow = 1) + theme_light()
a <- ggplot(data = Base_respi_ordenada, aes(x = Diseasetime, y = sitting_FEV1, color = LGMD)) + geom_point() + scale_y_continuous(name="FEV1% sitting") + geom_hline(yintercept=80, col="grey") + scale_x_continuous(name="Time of disease progression") + geom_hline(yintercept=80, col="grey") + geom_smooth(method = "lm", se=TRUE)
b <- a + facet_wrap(~LGMD, nrow = 1) + theme_light()
b
c <- ggplot(data = Base_respi_ordenada, aes(x = Diseasetime, y = Lying_FVC, color = LGMD)) + geom_point() + scale_y_continuous(name="FVC% lying") + geom_hline(yintercept=80, col="grey") + scale_x_continuous(name="Time of disease progression") + geom_hline(yintercept=80, col="grey") + geom_smooth(method = "lm", se=TRUE)
d <- c + facet_wrap(~LGMD, nrow = 1) + theme_light()
d
e <- ggplot(data = Base_respi_ordenada, aes(x = Diseasetime, y = Lying_FEV1, color = LGMD)) + geom_point() + scale_y_continuous(name="FEV1% lying") + geom_hline(yintercept=80, col="grey") + scale_x_continuous(name="Time of disease progression") + geom_hline(yintercept=80, col="grey") + geom_smooth(method = "lm", se=TRUE)
f <- e + facet_wrap(~LGMD, nrow = 1) + theme_light()
f
gggarrange(g, b, d, f, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
