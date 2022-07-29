rm(list=ls())

library(here)
library(ggplot2)
library(patchwork)
library(glmnet)

project_dir <- file.path(here(), "NCR")
figure_dir <- file.path(project_dir, "figures")

ucb_colors <- as.list(read.csv(file.path(project_dir, "palette.csv"), header=F)[,2])
names(ucb_colors) <- read.csv(file.path(project_dir, "palette.csv"), header=F)[,1]

# functions ---------------------------------------------------------------

optimize_glmnet_hyperparameters <- function(input_data, fold_ids, 
                                            metric, family="gaussian") {
  # input_data: matrix; first column is response variable
  # family: character; "family" argument for glmnet()
  # metric: character; one of "r2" or "mse"
  # fold_ids: numeric vector; "foldid" argument for cv.glmnet()
  num_folds <- length(unique(fold_ids))
  # 1. fit all models
  alphas <- seq(0, 1, by=0.1)
  all_models <- lapply(alphas,
                       function(alpha) {
                         glmnet::cv.glmnet(input_data[, -1], input_data[, 1], 
                                           foldid=fold_ids, alpha=alpha)
                       })
  names(all_models) <- alphas
  # 2. calculate r^2 and pearson correlation per split across all models
  all_models_stats <- lapply(seq_along(all_models),
                             function(alpha_id) {
                               lambdas <- all_models[[alpha_id]]$lambda
                               do.call(rbind,
                                       lapply(seq(num_folds),
                                              function(fold_id) {
                                                tmp_training <- subset(input_data,
                                                                       fold_ids != fold_id)
                                                tmp_testing <- subset(input_data,
                                                                      fold_ids == fold_id)
                                                tmp_fit <- glmnet::glmnet(tmp_training[, -1],
                                                                          tmp_training[, 1],
                                                                          alpha=alphas[alpha_id],
                                                                          lambda=lambdas)
                                                tmp_training_predict <- lapply(lambdas,
                                                                               function(lambda) {
                                                                                 predict(tmp_fit,
                                                                                         newx=tmp_training[,-1],
                                                                                         s=lambda)
                                                                               })
                                                tmp_testing_predict <- lapply(lambdas,
                                                                              function(lambda) {
                                                                                predict(tmp_fit,
                                                                                        newx=tmp_testing[,-1],
                                                                                        s=lambda)
                                                                              })
                                                tmp_testing_null_deviance <- sum((tmp_testing[,1] -
                                                                                    mean(tmp_testing[,1]))^2)
                                                tmp_testing_r2 <- sapply(tmp_testing_predict,
                                                                         function(x) {
                                                                           1 - sum((tmp_testing[,1]-x)^2) / 
                                                                             tmp_testing_null_deviance
                                                                         })
                                                tmp_training_cor <- sapply(tmp_training_predict,
                                                                           function(x) {
                                                                             cor(x, tmp_training[,1])
                                                                           })
                                                tmp_testing_cor <- sapply(tmp_testing_predict,
                                                                          function(x) {
                                                                            cor(x, tmp_testing[,1])
                                                                          })
                                                return(data.frame(alpha=alphas[alpha_id],
                                                                  fold_id=fold_id,
                                                                  lambda=rep(lambdas, times=2),
                                                                  r2=c(tmp_fit$dev.ratio, 
                                                                       tmp_testing_r2),
                                                                  pearson_r=c(tmp_training_cor,
                                                                              tmp_testing_cor),
                                                                  split=rep(c("training", "testing"),
                                                                            each=length(lambdas))))
                                              }))
                             })
  all_models_stats <- do.call(rbind, all_models_stats)
  all_models_stats$alpha <- factor(as.character(all_models_stats$alpha),
                                   levels=as.character(alphas))
  all_models_stats$split <- factor(all_models_stats$split,
                                   levels=c("training", "testing"))
  all_models_r2 <- aggregate(r2 ~ alpha + lambda + split, 
                             data=all_models_stats, FUN=mean)
  r2_plot <- ggplot(all_models_r2, aes(x=lambda, y=r2, col=alpha)) + 
    geom_line() + theme_bw() + scale_x_log10() + 
    facet_grid(split ~ ., scales="free_y") + 
    xlab(expr(lambda)) + ylab(bquote("mean"~r^2)) + labs(col=expr(alpha)) + 
    ggtitle(paste0("Elastic net regression, ", num_folds, "-fold cross-validation"))
  all_models_pearson_r <- aggregate(pearson_r ~ alpha + lambda + split, 
                                    data=all_models_stats, FUN=mean)
  pearson_r_plot <- ggplot(all_models_pearson_r, aes(x=lambda, y=pearson_r, col=alpha)) + 
    geom_line() + theme_bw() + scale_x_log10() + 
    facet_grid(split ~ ., scales="free_y") + 
    xlab(expr(lambda)) + ylab(bquote("pearson"~r)) + labs(col=expr(alpha)) + 
    ggtitle(paste0("Elastic net regression, ", num_folds, "-fold cross-validation"))
  # 3. calculate mean squared error across all models
  all_models_MSE <- lapply(seq_along(all_models),
                           function(x) {
                             tmp_fit <- all_models[[x]]
                             data.frame(lambda=tmp_fit$lambda,
                                        cv_error=tmp_fit$cvm,
                                        cv_upper=tmp_fit$cvup,
                                        cv_lower=tmp_fit$cvlo,
                                        num_nonzero=tmp_fit$nzero,
                                        alpha=alphas[x],
                                        dev_ratio=tmp_fit$glmnet.fit$dev.ratio)
                           })
  all_models_MSE <- do.call(rbind, all_models_MSE)
  all_models_MSE$alpha <- factor(all_models_MSE$alpha, 
                                 levels=as.character(alphas))
  mse_plot <- ggplot(all_models_MSE, aes(x=lambda, y=cv_error, col=alpha)) + 
    geom_line() + theme_bw() + scale_x_log10() + scale_y_log10() + 
    xlab(expr(lambda)) + labs(col=bquote(alpha)) + 
    ylab(paste0("cross-validated error (k=", num_folds, ")"))
  # 4. calculate best r^2 and pearson correlation per split per alpha
  # best r^2
  all_models_best_r2 <- do.call(rbind,
                                lapply(alphas,
                                       function(a) {
                                         tmp_testing <- subset(all_models_r2, 
                                                               split=="testing" & alpha==a)
                                         tmp_lambda <- tmp_testing$lambda[which.max(tmp_testing$r2)]
                                         subset(all_models_r2, alpha==a & lambda==tmp_lambda)
                                       }))
  best_r2_plot <- ggplot(all_models_best_r2, 
                         aes(x=as.numeric(as.character(alpha)), 
                             y=r2, col=split)) + 
    geom_point() + geom_line(alpha=0.5) + theme_bw() +
    xlab(expr(alpha)) + ylab(bquote("mean"~r^2)) + 
    scale_color_manual(values=c("#E69F00", "#56B4E9"))
  # best pearson r
  all_models_best_pearson_r <- do.call(rbind,
                                       lapply(alphas,
                                              function(a) {
                                                tmp_testing <- subset(all_models_pearson_r, 
                                                                      split=="testing" & alpha==a)
                                                tmp_lambda <- tmp_testing$lambda[which.max(tmp_testing$pearson_r)]
                                                subset(all_models_pearson_r, alpha==a & lambda==tmp_lambda)
                                              }))
  best_pearson_r_plot <- ggplot(all_models_best_pearson_r, 
                                aes(x=as.numeric(as.character(alpha)),
                                    y=pearson_r, col=split)) + 
    geom_point() + geom_line(alpha=0.5) + theme_bw() + 
    xlab(expr(alpha)) + ylab(bquote("mean pearson"~r)) + 
    scale_color_manual(values=c("#E69F00", "#56B4E9"))
  # pearson r, optimizing r^2
  match_var <- c("alpha", "lambda", "split")
  all_models_best_r2_pearson_r <- all_models_pearson_r[prodlim::row.match(all_models_best_r2[, match_var],
                                                                          all_models_pearson_r[, match_var]),]
  best_r2_pearson_r_plot <- ggplot(all_models_best_r2_pearson_r, 
                                   aes(x=as.numeric(as.character(alpha)),
                                       y=pearson_r, col=split)) + 
    geom_point() + geom_line(alpha=0.5) + theme_bw() + 
    xlab(expr(alpha)) + ylab(bquote("mean pearson"~r)) + 
    scale_color_manual(values=c("#E69F00", "#56B4E9"))
  # r^2, optimizing pearson r
  all_models_best_pearson_r_r2 <- all_models_r2[prodlim::row.match(all_models_best_pearson_r[, match_var],
                                                                   all_models_r2[, match_var]),]
  best_pearson_r_r2_plot <- ggplot(all_models_best_pearson_r_r2, aes(x=as.numeric(as.character(alpha)),
                                                                     y=r2, col=split)) + 
    geom_point() + geom_line(alpha=0.5) + theme_bw() + 
    xlab(expr(alpha)) + ylab(bquote("mean"~r^2)) + 
    scale_color_manual(values=c("#E69F00", "#56B4E9"))
  # 5. choose optimal hyperparameters
  if(metric=="mse") {
    best_param <- do.call(rbind, 
                          lapply(alphas,
                                 function(a) {
                                   tmp_fit <- subset(all_models_MSE, alpha==a)
                                   tmp_fit[which.min(tmp_fit$cv_error),]
                                 }))
    best_param <- best_param[which.min(best_param$cv_error),]
  } else {
    min_r2 <- with(subset(all_models_best_r2, split=="testing"), max(r2) - sd(r2))
    best_param <- subset(all_models_best_r2, split=="testing" & r2 >= min_r2)$alpha
    best_param <- subset(all_models_best_r2, split=="training" & alpha %in% best_param)
    best_param <- best_param[which.max(best_param$r2), ]
  }
  best_lambda <- best_param$lambda
  best_alpha <- as.numeric(as.character(best_param$alpha))
  best_r2_plot <- best_r2_plot + 
    geom_point(data=subset(all_models_best_r2, alpha==best_alpha & lambda==best_lambda),
               aes(x=as.numeric(as.character(alpha)), y=r2), col="red") + 
    ggtitle("optimal hyperparameters:", 
            subtitle=bquote(alpha~"="~.(best_alpha)~";"~lambda~"="~.(best_lambda)))
  best_pearson_r_plot <- best_pearson_r_plot + 
    geom_point(data=subset(all_models_best_pearson_r, alpha==best_alpha & lambda==best_lambda),
               aes(x=as.numeric(as.character(alpha)), y=pearson_r), col="red")
  r2_plot <- r2_plot + 
    geom_point(data=subset(all_models_best_r2, alpha==best_alpha & lambda==best_lambda),
               aes(x=lambda, y=r2), col="red")
  mse_plot <- mse_plot + 
    geom_point(data=subset(all_models_MSE, alpha==best_alpha & lambda==best_lambda),
               aes(x=lambda, y=cv_error), col="red")
  
  return(list(all_models=all_models, all_models_stats=all_models_stats, MSE=all_models_MSE,
              best_alpha=best_alpha, best_lambda=best_lambda, metric=metric,
              r2_plot=r2_plot, pearson_r_plot=pearson_r_plot, mse_plot=mse_plot, 
              best_r2_plot=best_r2_plot, best_pearson_r_plot=best_pearson_r_plot,
              best_r2_pearson_r_plot=best_r2_pearson_r_plot,
              best_pearson_r_r2_plot=best_pearson_r_r2_plot))
}

plot_coefs <- function(fit_obj, best_alpha, best_lambda) {
  # fit_obj: cv.glmnet object
  # best_alpha: numeric
  # best_lambda: numeric
  coefs <- coef(fit_obj, s=best_lambda)[, 1]
  coefs <- data.frame(matrix(unlist(strsplit(names(coefs)[-1],
                                             split="_")),
                             ncol=3, byrow=T),
                      coefs[-1],
                      row.names=NULL)
  colnames(coefs) <- c("type", "position", "coef", "value")
  coefs$position <- as.numeric(sub("\\.", "-", coefs$position))
  coef_plot <- ggplot(coefs, aes(x=position, y=value, fill=coef)) + 
    geom_col(position="stack") + geom_hline(yintercept=0) + 
    xlab("position relative to spacer") + ylab("regression coefficient") + 
    theme_classic(base_size=8) + labs(fill="") + theme(legend.position="bottom") + 
    scale_fill_manual(values=c(RColorBrewer::brewer.pal(4, "Set1"), "grey")) +
    annotate("rect", xmin=0.5, xmax=19.5, 
             ymin=max(coefs$value)+1, ymax=max(coefs$value+2), fill="black") + 
    annotate("text", x=10, y=max(coefs$value)+1.5, 
             label="protospacer", col="white") +
    annotate("rect", xmin=-3.5, xmax=0.5, 
             ymin=max(coefs$value)+1, ymax=max(coefs$value+2), fill="darkgrey") + 
    annotate("text", x=-1.5, y=max(coefs$value)+1.5, 
             label="anti-tag", col="white")
  if(max(coefs$position) >= 50) {
    coef_plot <- coef_plot + 
      annotate("rect", xmin=30.5, xmax=49.5, 
               ymin=max(coefs$value)+1, ymax=max(coefs$value)+2, 
               fill="darkgrey") + 
      annotate("text", x=40.5, y=max(coefs$value)+1.5, 
               label="cis-cleavage neighborhood", col="white")
  }
  return(coef_plot)
}

plot_fold_fits <- function(input_data, fold_ids, best_alpha, best_lambda) {
  # input_data: matrix; first column is response variable
  # fold_ids: numeric vector; folds assigned to each row of input_data
  # best_alpha: numeric
  # best_lambda: numeric
  num_folds <- length(unique(fold_ids))
  fold_plots <- lapply(seq(num_folds),
                       function(x) {
                         # split data into training and testing
                         tmp_data <- split(seq(nrow(input_data)), fold_ids==x)
                         tmp_data <- lapply(tmp_data,
                                            function(indices) {
                                              data.frame(input_data[indices,])
                                            })
                         names(tmp_data) <- c("training", "testing")
                         # compute fit
                         tmp_fit <- glmnet(as.matrix(tmp_data$training)[, -1],
                                           tmp_data$training[, 1],
                                           alpha=best_alpha, lambda=best_lambda)
                         # predict outcome
                         tmp_data$training$predicted <- predict(tmp_fit,
                                                                newx=as.matrix(tmp_data$training)[, -1],
                                                                s=best_lambda)
                         tmp_data$testing$predicted <- predict(tmp_fit,
                                                               newx=as.matrix(tmp_data$testing)[, -1],
                                                               s=best_lambda)
                         # generate plot
                         tmp_plot_data <- rbind(data.frame(tmp_data$training, split="training"),
                                                data.frame(tmp_data$testing, split="testing"))
                         tmp_r2 <- sapply(tmp_data,
                                          function(x) {
                                            with(x, 1-sum((rate-predicted)^2)/sum((rate-mean(rate))^2))
                                          })
                         tmp_cor <- sapply(tmp_data,
                                           function(x) {
                                             with(x, cor(rate, predicted))
                                           })
                         tmp_plot <- ggplot(tmp_plot_data, aes(x=rate, y=predicted, col=split)) + 
                           geom_point() + geom_abline(slope=1, intercept=0) + theme_bw() + 
                           ggtitle(paste("fold:", x), subtitle=paste("# non-zero:", tmp_fit$df)) + 
                           xlab("true rate (RFU/min)") + ylab("predicted rate (RFU/min)") + 
                           scale_color_manual(values=c(ucb_colors$`California Gold`, ucb_colors$Lawrence)) + 
                           annotate(geom="text", x=Inf, y=-Inf, col=ucb_colors$Lawrence,
                                    label=bquote(r^2~"="~.(round(tmp_r2["training"], digits=3))),
                                    hjust="inward", vjust=-2) + 
                           annotate(geom="text", x=Inf, y=-Inf, col=ucb_colors$`California Gold`,
                                    label=bquote(r^2~"="~.(round(tmp_r2["testing"], digits=3))),
                                    hjust="inward", vjust=-1) + 
                           annotate(geom="text", x=Inf, y=Inf, col="black", label="y=x",
                                    hjust="inward", vjust="inward") + 
                           annotate(geom="text", x=-Inf, y=Inf, col=ucb_colors$Lawrence,
                                    label=bquote(rho~"="~.(round(tmp_cor["training"], digits=2))),
                                    hjust="inward", vjust=1.5) +
                           annotate(geom="text", x=-Inf, y=Inf, col=ucb_colors$`California Gold`,
                                    label=bquote(rho~"="~.(round(tmp_cor["testing"], digits=2))),
                                    hjust="inward", vjust=2.5)
                         return(tmp_plot)
                       })
  fold_plots <- lapply(seq_along(fold_plots),
                       function(x) {
                         if(x != length(fold_plots)) {
                           return(fold_plots[[x]] + guides(col="none"))
                         } else {
                           return(fold_plots[[x]])
                         }
                       })
  return(wrap_plots(fold_plots, nrow=1))
}

# load guide rates from primary vRNA screen -------------------------------

guide_rate <- read.csv(file.path(project_dir, "data", "guide_rate.csv"))

# format data for regression ----------------------------------------------

window_start <- -5
window_stop <- 5

# load viral genome sequence
viral_genome <- read.csv(file.path(project_dir, "data", "viral_genome_features.csv"))

# annotate sequence features
guide_sequence <- t(sapply(guide_rate$start,
                           function(x) {
                             tmp_positions <- (x+19+window_stop):(x+window_start)
                             viral_genome$base[tmp_positions]
                           }))
guide_sequence <- data.frame(lapply(seq(ncol(guide_sequence)),
                                    function(x)  {
                                      factor(guide_sequence[, x], 
                                             levels=c("A", "C", "G", "U")) 
                                    }))
colnames(guide_sequence) <- paste0("seq_", (1+window_start):(20+window_stop))

# annotate structure features
guide_structure <- t(sapply(guide_rate$start,
                            function(x) {
                              tmp_positions <- (x+19+window_stop):(x+window_start)
                              viral_genome$invivo[tmp_positions]
                            }))
guide_structure <- data.frame(lapply(seq(ncol(guide_structure)),
                                     function(x) {
                                       factor(guide_structure[, x],
                                              levels=c("structured", "unstructured"))
                                     }))
colnames(guide_structure) <- paste0("str_", (1+window_start):(20+window_stop))

# wrangle data for model
regression_data <- data.frame(rate=guide_rate$Estimate, guide_sequence, guide_structure)
regression_data <- subset(regression_data, nchar(guide_rate$spacer)==20)
regression_data_onehot <- as.matrix(mltools::one_hot(data.table::data.table(regression_data)))
regression_data_onehot <- subset(regression_data_onehot,
                                 select=!grepl("_structured", colnames(regression_data_onehot)))

write.csv(regression_data, 
          file=file.path(project_dir, "data", "regression_data.csv"), 
          row.names=F)

# compute elastic net regression ------------------------------------------

# set folds
set.seed(21)
num_folds <- 5
fold_ids <- sample(cut(seq(sum(nchar(guide_rate$spacer)==20)), breaks=num_folds, labels=F))

# optimize hyperparameters
model_fit <- list(r2=optimize_glmnet_hyperparameters(regression_data_onehot,
                                                     fold_ids,
                                                     metric="r2"),
                  MSE=optimize_glmnet_hyperparameters(regression_data_onehot,
                                                      fold_ids,
                                                      metric="mse"))
model_fit_hyperparameter_plot <- ((model_fit$r2$best_r2_plot + ggtitle(bquote("optimize"~r^2))) + 
                                    (model_fit$MSE$best_r2_plot + ggtitle("optimize MSE"))) /
  (model_fit$r2$mse_plot + model_fit$MSE$mse_plot)
model_fit_optimize_plot <- ((model_fit$MSE$best_r2_plot + ggtitle("optimizing testing r^2")) / 
                              model_fit$MSE$best_r2_pearson_r_plot) | 
  ((model_fit$MSE$best_pearson_r_r2_plot + ggtitle("optimize testing pearson r")) /
     model_fit$MSE$best_pearson_r_plot)
# plot coefficients
model_fit_bestfit <- model_fit$MSE$all_models[[which(names(model_fit$MSE$all_models) == model_fit$MSE$best_alpha)]]
model_fit_coef_plot <- plot_coefs(model_fit_bestfit, model_fit$MSE$best_alpha, model_fit$MSE$best_lambda)
# plot fit across folds
model_fit_fold_plot <- plot_fold_fits(regression_data_onehot, fold_ids, 
                                      model_fit$MSE$best_alpha, model_fit$MSE$best_lambda)
model_fit_plot <- ((model_fit$MSE$mse_plot + 
                      ggtitle("Model 1", subtitle="rate ~ sequence + structure (+/- 5nt)")) +
                     (model_fit$MSE$best_r2_plot) +
                     (model_fit$MSE$best_pearson_r_plot)) / 
  (model_fit_coef_plot + ggtitle("", subtitle="") + 
     labs(col="sequence", shape="structure")) /
  model_fit_fold_plot

# generate plot -----------------------------------------------------------

fill_colors <- c(RColorBrewer::brewer.pal(4, "Set1"), "grey")

model_fit_coefs <- coef(model_fit_bestfit, s=model_fit$MSE$best_lambda)[, 1]
model_fit_coefs <- data.frame(matrix(unlist(strsplit(names(model_fit_coefs)[-1],
                                                     split="_")),
                                     ncol=3, byrow=T),
                              model_fit_coefs[-1],
                              row.names=NULL)
colnames(model_fit_coefs) <- c("type", "position", "coef", "value")
model_fit_coefs$position <- as.numeric(sub("\\.", "-", model_fit_coefs$position))

figure_4A <- ggplot(model_fit_coefs, aes(x=position, y=value, fill=coef)) + 
  geom_col(position="stack") + geom_hline(yintercept=0) +
  xlab("position relative to spacer") + ylab("regression\ncoefficient") + 
  theme_classic(base_size=8) + labs(fill="") +
  scale_x_continuous(breaks=c(1, 10, 20)) +
  scale_fill_manual(values=c("G"=fill_colors[which(unique(model_fit_coefs$coef)=="G")],
                             "A"=fill_colors[which(unique(model_fit_coefs$coef)=="A")],
                             "C"=fill_colors[which(unique(model_fit_coefs$coef)=="C")],
                             "U"=fill_colors[which(unique(model_fit_coefs$coef)=="U")],
                             "unstructured"=fill_colors[which(unique(model_fit_coefs$coef)=="unstructured")])) +
  theme(legend.position="bottom")
  # annotate("rect", xmin=0.5, xmax=19.5, ymin=2, ymax=4, fill="black") + 
  # annotate("text", x=10.5, y=3, label="protospacer", col="white", size=2.5) +
  # annotate("rect", xmin=-5.5, xmax=0.5, ymin=2, ymax=4, fill="darkgrey") +
  # annotate("text", x=-2.5, y=3, label="anti-tag", col="white", size=2.5) 

ggsave(filename=file.path(figure_dir, "figure_4A.pdf"),
       plot=figure_4A,
       device="pdf", width=3.25, height=2, units="in")
