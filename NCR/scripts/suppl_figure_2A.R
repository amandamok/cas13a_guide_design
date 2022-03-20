rm(list=ls())

library(here)
library(ggplot2)

project_dir <- file.path(here(), "NCR")
data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

# load regression data ----------------------------------------------------

regression_data <- read.csv(file.path(data_dir, "regression_data.csv"),
                            stringsAsFactors=T)

# one-hot encoding
regression_data_onehot <- as.matrix(mltools::one_hot(data.table::data.table(regression_data)))
regression_data_onehot <- subset(regression_data_onehot,
                                 select=!grepl("_structured", colnames(regression_data_onehot)))

# compute random forest regression ----------------------------------------

rf_regression <- randomForest::randomForest(x=regression_data_onehot[, -1],
                                            y=regression_data_onehot[, 1],
                                            importance=T)
rf_regression_importance <- data.frame(rf_regression$importance, 
                                       rf_regression$importanceSD,
                                       matrix(unlist(strsplit(rownames(rf_regression$importance), 
                                                              split="_")), 
                                              ncol=3, byrow=T))
colnames(rf_regression_importance) <- c("increaseMSE", "increaseNodePurity", 
                                        "increaseMSE_stderr", "type", "position", "coef")
rf_regression_importance$position <- with(rf_regression_importance,
                                          as.numeric(sub("\\.", "-", position)))
rf_regression_r2 <- 1-sum((regression_data_onehot[,1]-rf_regression$predicted)^2) / 
  sum((regression_data_onehot[,1]-mean(regression_data_onehot[,1]))^2)

# generate plot -----------------------------------------------------------

suppl_figure_2A <- ggplot(rf_regression_importance, aes(x=position, y=increaseMSE, fill=coef)) + 
  geom_col(position="stack") + geom_hline(yintercept=0) +
  xlab("position relative to spacer") + ylab("MSE increase") + 
  theme_classic(base_size=10) + labs(fill="") +
  scale_x_continuous(breaks=c(1, 10, 20)) +
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(4, "Set1"), "grey")) +
  theme(legend.position="bottom") +
  annotate("rect", xmin=0.5, xmax=19.5, ymin=-18, ymax=-10, fill="black") + 
  annotate("text", x=10.5, y=-14, label="protospacer", col="white", size=2.5) +
  annotate("rect", xmin=-5.5, xmax=0.5, ymin=-18, ymax=-10, fill="darkgrey") +
  annotate("text", x=-2.5, y=-14, label="anti-tag", col="white", size=2.5) 

ggsave(filename=file.path(figure_dir, "suppl_figure_2A.pdf"),
       plot=suppl_figure_2A,
       device="pdf", width=6.5, height=2, units="in")
