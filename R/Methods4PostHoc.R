# ============================================================================ #
#' Print methods for objects of class PostHoc
#'
#' @author Rodrigo Labouriau
#' @param x an object of class PostHoc to be printed
#' @param digits number of digits in the output (default = 4)
#' @param ... further arguments passed to or from other methods.
#' @return a dataframe with two variables, Levels a factor containing the levels
#'  of the effects and ParaeterCI which is a factor with the effects and the
#'  corresponding confidence intervals and the grouping combined in a character
#'  constructed in such a way that when printing this dataframe yields a table
#'  arranged in a suitable format.
#' @examples MM <- glm(Y ~ Treatment + 0,  data = DeIdentifiedExample)
#' @examples GG <- posthoc(MM)
#' @examples print(GG)
#' @export
print.PostHoc <- function(x, digits = 4, ...){
  object <- x

  CI <- object$CI

  if(is.null(CI)){
    stop("The confidence intervals were not available.")
  }

  CI <- round(CI , digits = digits)

  Out <- data.frame(CI, object$Grouping)
  Low <- object$SignificanceLevel/2*100
  Up  <- 100-object$SignificanceLevel/2*100
  names(Out) <- c("Est.",
                  paste(Low,"%", sep=""),
                  paste(Up,"%", sep=""),
                  "Group"
  )

  if(is.null(row.names(CI))) row.names(CI) <- paste("Level", 1:(dim(CI)[1]))
  Out <-  data.frame(Levels=row.names(CI),
                    ParameterCI=paste( CI[ ,1], "(", CI[ ,2], "-",
                                       CI[ ,3], ")",object$Grouping , sep="")
    )

    return(Out)
  }

# ============================================================================ #

# ============================================================================ #
#' Plot method for objects of class PostHoc
#'
#' @author Rodrigo Labouriau
#' @param x an object of class PostHoc to be printed.
#' @param y an object of class PostHoc to be printed.
#' @param ... further arguments passed to or from other methods.
#' @return none
#' @examples MM <- glm(Y ~ Treatment + 0,  data = DeIdentifiedExample)
#' @examples GG <- posthoc(MM)
#' @examples plot(GG)
#' @export
plot.PostHoc <- function(x, y, ...)
{
  object <- x

  Cliques <- object$Cliques
  plot(object$Gr, mark.groups = Cliques, ...)
}
# ============================================================================ #

# ============================================================================ #
#' Summary method for objects of class PostHoc
#'
#' @author Rodrigo Labouriau
#' @param object an object of class PostHoc to be printed.
#' @param ... further arguments passed to or from other methods.
#' @return a dataframe constructed in such a way that when printing this
#'  dataframe yields a table arranged in a suitable format. The summary,
#'  differently than the print method displays also the matrix of p-values
#'  of all the pairwise comparisons.
#' @examples MM <- glm(Y ~ Treatment + 0,  data = DeIdentifiedExample)
#' @examples GG <- posthoc(MM)
#' @examples summary(GG)
#' @export
summary.PostHoc <- function(object, ...){

  if (class(object) != "PostHoc"){
    stop("Sorry, but the object should be of class PostHoc")
  }

  if(is.null(object$CI)){
    stop("The confidence intervals were not available.")
  }

  Pvalues <- as.character(round(object$PvaluesMatrix, digits=object$digits))
  Pvalues[is.na(Pvalues)] <- ""
  dim(Pvalues) <- dim(object$PvaluesMatrix)
  Pvalues <- Pvalues[ , -dim(object$PvaluesMatrix)[2] ]
  Out <- data.frame(object$CI, object$Grouping, Pvalues )
  Low <- object$SignificanceLevel/2*100
  Up  <- 100-object$SignificanceLevel/2*100
  names(Out) <- c("Est.",
                  paste(Low,"%", sep=""),
                  paste(Up,"%", sep=""),
                  "Group",
                  row.names(object$CI)[-length(row.names(object$CI))]
  )
  return(Out)
}
# ==============================================================================

# ============================================================================ #
#' Barplot method for objects of class PostHoc
#'
#' @author Rodrigo Labouriau
#' @param height an object of class PostHoc to be printed
#' @param col the colour of the bars (default = "lightblue")
#' @param labelsCol the colour of the bars (default = "black")
#' @param space2max space between the upper limit of the confidence interval and the label¨
#' @param UseGrouping should the grouping be added to the plots (default = TRUE)
#' @param ylim range of the vertical axis
#' @param main character with the title of the plot (default = '')
#' @param ylab label of the vertical axis
#' @param xlab label of the horizontal axis
#' @param lty type of line
#' @param drawAxis should the axis be drawn (default = TRUE)
#' @param rightshift a number specifying a (small) right shift of the line
#' @param additionalTextGrouping character vector with additional text to the grouping
#' @param superpose should the graph be superposed to an existing graph (default = FALSE)
#' @param cex.grouping size of the labels of groups
#' @param cex.ticks size of the thicks defining the the limits of the confidence intervals
#' @param cex.lab size of the labels
#' @param ylog should the vertical axis be expressed in a logarithmic scale (default = FALSE)
#' @param ... parameters to be passed to the barlot function
#' @return none
#' @examples MM <- glm(Y ~ Treatment+0,  data=DeIdentifiedExample)
#' @examples GG <- posthoc(MM)
#' @examples barplot(GG)
#' @import graphics
#' @export
barplot.PostHoc <- function(height,
                            col = "lightblue",
                            labelsCol = "black",
                            space2max = 0,
                            UseGrouping = TRUE,
                            ylim = NULL,
                            main = "",
                            ylab = "",
                            xlab = "",
                            lty = 1,
                            drawAxis = TRUE,
                            rightshift = 0,
                            additionalTextGrouping = "",
                            superpose = FALSE,
                            cex.grouping = 1.0,
                            cex.ticks = 0.1,
                            cex.lab = 1,
                            ylog = FALSE,
                            ...){

  object <- height
  if (class(object) != "PostHoc"){
    stop("Sorry, but the object should be of class PostHoc")
  }

  if(is.null(object$CI)){
    stop("The confidence intervals were not available.")
  }

  Groups <- object$Grouping
  CI <- object$CI
  Effects <- CI[ ,1]

  if(is.null(ylim)){
    DeltaSpace <- (abs(max(CI - min(CI))) * 0.10) +  space2max
    ylim <- c(min(0, min(CI)), max(CI)+DeltaSpace)
  }

  bp <- barplot(Effects, ylim = ylim, col = col, xpd = FALSE,
                ylab = ylab, xlab = xlab, main = main)

  arrows(x0 = bp, x1 = bp, y0 = CI[ ,2], y1 = CI[ ,3], lwd = 1,
         angle = 90, col = labelsCol, code = 3, length = cex.ticks)

  text(x = bp, y = CI[ ,3] + space2max/2,
       labels = paste(Groups, additionalTextGrouping, sep =""),
       pos = 3, col = labelsCol, cex = cex.grouping)

}
# ============================================================================ #

# ============================================================================ #
# ============================================================================ #
#' Lines method for objects of class PostHoc
#'
#' @author Rodrigo Labouriau
#' @param x an object of class PostHoc to be printed
#' @param col the colour of the lines (default = "black")
#' @param labelsCol the colour of the bars (default = "black")
#' @param space2max space between the upper limit of the confidence interval
#'   and the label¨
#' @param UseGrouping should the grouping be added to the plots (default = TRUE)
#' @param ylim range of the vertical axis
#' @param main character with the title of the plot (default = '')
#' @param ylab label of the vertical axis
#' @param xlab label of the horizontal axis
#' @param lty type of line
#' @param drawAxis should the axis be drawn (default = TRUE)
#' @param rightshift a number specifying a (small) right shift of the line
#' @param additionalTextGrouping character vector with additional text to the
#'  grouping
#' @param superpose should the graph be superposed to an existing graph
#'  (default = FALSE)
#' @param cex.grouping size of the labels of groups
#' @param cex.ticks size of the thicks defining the the limits of the
#'  confidence intervals
#' @param cex.lab size of the labels
#' @param ylog should the vertical axis be expressed in a logarithmic scale
#'  (default = FALSE)
#' @param ... parameters to be passed to the function lines
#' @return none
#' @examples MM <- glm(Y ~ Treatment+0,  data = DeIdentifiedExample)
#' @examples GG <- posthoc(MM)
#' @examples lines(GG)
#' @import graphics
#' @export
lines.PostHoc <- function(x,
                          col = "black",
                          labelsCol = "black",
                          space2max = 0,
                          UseGrouping = TRUE,
                          ylim = NULL,
                          main = "",
                          ylab = "",
                          xlab = "",
                          lty = 1,
                          drawAxis = TRUE,
                          rightshift = 0,
                          additionalTextGrouping = "",
                          superpose = FALSE,
                          cex.grouping = 1.0,
                          cex.ticks = 0.1,
                          cex.lab = 1,
                          ylog = FALSE,
                          ...)
{
  object <- x
  if (class(object) != "PostHoc"){
    stop("Sorry, but the object should be of class GroupClusterEff")
  }

  Groups <- object$Grouping

  CI <- object$CI
  Effects <- CI[ ,1]

  Nlevels <- length(Effects)
  # bp <- barplot(Effects, plot = FALSE)
  bp <- 1:Nlevels
  # if(rightshift != 0){
  bp <- bp + rightshift
  # }

  if(is.null( ylim)){
    DeltaSpace <- (abs(max(CI - min(CI))) * 0.10) +   space2max
    ylim <- c(min(0, min(CI)), max(CI)+DeltaSpace)
  }

  if(is.null( ylim)){
    ylim <- c(min( CI), max( CI))
  } else{
    ylim <-  ylim
  }

  if( superpose){ ################
    points(bp, Effects, type = "l", lwd = 3, lty =  lty,
           col =  labelsCol)
  } else{
    plot(bp, Effects, ylim = ylim, type = "l", lwd = 3, axes = FALSE,
         ylab =  ylab, xlab =  xlab, main =  main,
         lty =  lty,
         col =  labelsCol, ylog =  ylog)
  }

    if( drawAxis){
      Axis(side = 1, at = bp, labels = names(Effects),
                               cex.lab =  cex.lab)
      Axis(side = 2, labels = TRUE)
    }

    arrows(x0 = bp, x1 = bp, y0 = CI[ ,2], y1 = CI[ ,3], lwd = 1,
           angle = 90, col =  labelsCol, code = 3, length =  cex.ticks)

    text(x = bp, y = CI[ ,3] +  space2max/2,
         labels = paste(Groups,  additionalTextGrouping, sep =""),
         pos = 3, col =  labelsCol, cex =  cex.grouping)

}
# ============================================================================ #
