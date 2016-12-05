# Copyright (C) 2016 Genome Research Ltd.

#' @import dplyr
#' @import ggplot2
NULL

l3_model <- function (x, xmid, scal)
{
  .expr1 <- x - xmid
  .expr4 <- exp(-.expr1/(scal))
  .expr5 <- 1 + .expr4
  .expr9 <- .expr5^2
  .value <- 1/.expr5
  return(.value)
}

l3_model2 <- function(lx, maxc, xmid, scal){
  x <- getXfromConc(exp(lx), maxc)
  yhat <- 1 - l3_model(x, xmid, scal)
  return(yhat)
}

getConcFromX <- function(x, maxc) {
  xc <- maxc * 2 ^ (x - 9)
  return(xc)
}

getXfromConc <- function(xc, maxc) {
  x <- (log(xc / maxc)/log(2))+ 9
  return(x)
}

#' plotResponse
#'
#' Plot dose reponse curve.
#' @param model_stats dataframe of fitted values from the GDSC nlme logistic
#' dose response model
#' @param cl as identified from the \code{CL} column of the model stats data
#' frame.
#' @param drug character string as identified by the \code{drug} column of the
#' model stats data frame.
#'
#' @return ggplot of dose response curve.
#' @export
#'
#' @examples
#' plotResponse()
plotResponse <- function(model_stats, cell_line, drug_identifier) {
  plot_data <- model_stats %>%
    filter_(~CL == cell_line, ~drug == drug_identifier)
  stopifnot(nrow(plot_data) > 0)

  IC50 <- unique(plot_data$IC50)
  stopifnot(length(IC50) == 1)

  auc <- unique(plot_data$auc)
  stopifnot(length(auc) == 1)

  rmse <- unique(plot_data$RMSE)
  stopifnot(length(rmse) == 1)

  drug_id <- unique(plot_data$DRUG_ID)
  stopifnot(length(drug_id) == 1)

  cell_line_name <- unique(plot_data$CELL_LINE_NAME)
  cell_line_name <- ifelse(is.null(cell_line_name),  "", cell_line_name)
  stopifnot(length(cell_line_name) == 1)

  max_conc <- unique(plot_data$maxc)
  stopifnot(length(max_conc) == 1)

  plot_data <- plot_data %>%
    mutate_(lx = ~log(getConcFromX(x, maxc)),
            lxmid = ~log(getConcFromX(xmid, maxc))
            )

  plot_xmid <- plot_data %>% select_(~xmid) %>% distinct()
  plot_scal <- plot_data %>% select_(~scal) %>% distinct()
  plot_maxc <- plot_data %>% select_(~maxc) %>% distinct()

  plot_low_x <- 1 - plot_scal$scal * log((1 - 1e-3) / 1e-3) + plot_xmid$xmid
  plot_low_x <- log(getConcFromX(plot_low_x, +plot_maxc$maxc))
  plot_low_x <- min(c(plot_data$lx, plot_low_x))

  plot_high_x <- 1 - plot_scal$scal * log(1e-3 / (1 - 1e-3)) + plot_xmid$xmid
  plot_high_x <- log(getConcFromX(plot_high_x, plot_maxc$maxc))
  plot_high_x <- max(c(plot_data$lx, plot_high_x))

  p <- ggplot(plot_data) + aes(x = lx, y = 1 - yhat) + geom_point(shape = 3) +
  scale_x_continuous(limits = c(plot_low_x, plot_high_x))

  p <- p + stat_function(aes(x = lx),
                         fun = l3_model2,
                         args = list(maxc = plot_maxc$maxc,
                                     xmid = plot_xmid$xmid,
                                     scal = plot_scal$scal))
  p <- p + annotate("rect",
                     xmin = min(plot_data$lx),
                     xmax = max(plot_data$lx),
                     ymin = 0,
                     ymax =  max(c(1, 1 - plot_data$yhat)),
                     alpha = 0.2
  )

  p <- p + aes(x = lx, y = 1 - yhat) + geom_point(shape = 4)

  p <- p + geom_point(aes(x = lx, y = 1 - y), shape = 1)
  p <- p + geom_point(aes(x = lxmid, y = 0.5), colour = "red") + theme(legend.position = "none")
  p <- p + geom_point(aes(x = lxmid, y = 0.5), colour = "green", shape = 5, size = 3) + theme(legend.position = "none")
  p <- p + geom_linerange(aes(x = IC50, ymax = 0.5, ymin = 0, colour = "red"), linetype = "dashed")
  p <- p + annotate("label",
                    x = IC50 + 1,
                    y = 0.5,
                    hjust = "left",
                    label = sprintf("IC50==%.3f~log[e]~mu*M", IC50), parse = T)
  p <- p + annotate("label",
                    x = min(plot_data$lx) + 0.1,
                    y = 0.25,
                    hjust = "left",
                    label = sprintf("auc = %.3f", auc))
  p  <- p + ylab("Response: normalized intensity") +
    theme(axis.title.y = element_text(size=14))
  p <- p + xlab(expression(Dose/log[e]~mu*M)) +
    theme(axis.title.x = element_text(size=14))

  p <- p + annotate("label",
                     # plot label at 2/3 x width
                     x = plot_low_x + ((plot_high_x - plot_low_x) * 2 / 3),
                     y =  max(c(1, 1 - plot_data$yhat)),
                     hjust = "left",
                     vjust = "inward",
                     # label = sprintf("rmse = %.3f", rmse)
                     label = sprintf("Cell-line: %s\nDrug id: %d\nMax dose = %.3f uM\nrmse = %.3f",
                                    cell_line_name, drug_id, max_conc, rmse))
  # label = paste(sprintf("Cell-line:~%s", cell_line_name), "\n", sprintf("Drug~id:~%d", drug_id)), parse = T)

  p <- p + ggtitle(paste("Dose response: cell line ", cell_line, "; drug ", drug_identifier, "."))

  return(p)
}
