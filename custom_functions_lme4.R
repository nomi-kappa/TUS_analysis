#### All functions for analysis code ####

# ============================================================================ #
#### Recode character variables to proper factors ####

recode_char2fac <- function(data){
  #' Turn string variables into proper factors again
  #' @param data    data frame with trial-level data
  #' @return data   same data frame with all character variables turned into proper factors
  
  cat("Recode any character variable to factor\n")
  
  for (iCol in 1:ncol(data)){
    if (is.character(data[,iCol])){ # if character
      cat(paste0("Recode variable ",names(data)[iCol],"\n"))
      data[,iCol] <- factor(data[,iCol]) # recode to proper factor
    }
  }
  
  cat("Finished :-)\n")
  return(data)
}

# ============================================================================ #
#### Recode all numerical variables to factors ####

recode_num2fac <- function(data,variables = NULL){
  #' Recode selected numerical variable to factor
  #' @param data    data frame with trial-level data
  #' @variables     vector of strings, numerical variables to turn into factors (if not provided: all numerical variables in data frame)
  #' @return data   same data frame with all numerical variables also available as factor ("_f")
  
  cat("Recode selected numerical variables to factors\n")
  
  ### Determine if input variables given, otherwise take all numerical variables: 
  if(is.null(variables)){
    cat("No variables given, take all numerical variables\n")
    variables <- c()
    for (iCol in 1:ncol(data)){
      if (is.numeric(data[,iCol])){ # if numeric
        variables <- c(variables,names(data)[iCol]) # add variable name
      }
    }
  }
  
  ### Number variables:
  nVar <- length(variables)
  
  ### Loop through variables:  
  for (iVar in 1:nVar){ # loop through variables
    varName <- variables[iVar]
    cat(paste0("Recode variable ",varName,"\n"))
    newVarName <- paste0(varName,"_f") # new variable name
    data[,newVarName] <- factor(data[,varName]) # recode to proper factor
  }
  cat("Finished :-)\n")
  return(data)
}

# ============================================================================ #
#### Recode formula in Wilkinson notation to handle that can be used for saving models: ####

formula2handle <- function(formula){
  #' Create name based on formula with underscores instead of spaces without random-effects part to be used in plot names.
  #' @param formula   string, formula in Wilkinson notation.
  #' @return handle   string, converted to be without spaces but with underscores, with random-effects part removed.
  require(stringr)
  
  cat(paste0("Input: ", formula, "\n"))

  handle <- formula # copy over
  
  # ------------------------------------------------------------------------ #
  ## Extract until random effects parentheses:
  
  handle  <- sub("\\(.*", "", handle) # delete everything after (
  
  # ------------------------------------------------------------------------ #
  ## Delete only very last (!) plus before random-effects part:

  if(grepl( "+", str_sub(handle, -2, -1), fixed = T)){
    handle <- substring(handle, 1, nchar(handle) - 3)
    
  }
  
  ## Replace every * by x:
  handle <- gsub("\\*", "x", handle) # extract until last plus
  
  ## Substitute every space with underscore:
  handle <- gsub(" ", "_", handle)
  
  cat(paste0("Output: ", handle, "\n"))
  return(handle)
  
}

# ============================================================================ #
#### Print probability of direction for all effects in brms model: ####

print_Pd <- function(mod, nPad = 40){
  #' Print probability of effect (Pd), 
  #' i.e., percentage of posterior that is above/below zero (whatever is smaller),
  #' for each fixed effect in brms model.
  #' @param mod       model object, model fitted in brms.
  #' @param nPad      numeric scalar, number of characters to pad effect name with spaces (default: 40).
  #' @return none, just printing to console.
  
  # Manually:
  # mean(c(as_draws_matrix(brms_mod, variable = "b_reqAction_f1")) < 0) # Probability of direction
  
  print(mod$formula)
  
  ## Extract names of all available fixed effects:
  draws_fit <- as_draws_array(mod) # convert model to large array.
  varVec <- dimnames(draws_fit)$variable # all variables in model.
  effVec <- varVec[grep("b_", varVec, fixed = T)] # retrieve names of fixed effects terms.
  
  ## Loop over effects:
  for (iEff in effVec){ # iEff <- effVec[5]
    
    ## Name of effect, padded:
    iEff_str <- str_pad(iEff, nPad, side = "right", pad = " ")
    
    ## Compute both possible Pds:
    posPd <- mean(c(as_draws_matrix(mod, variable = iEff)) < 0) # Pd for positive effects.
    negPd <- mean(c(as_draws_matrix(mod, variable = iEff)) > 0) # Pd for negative effects.
    
    ## Compare:
    isPos <- posPd < negPd # whether positive or negative effect.
    effDir <- ifelse(isPos, "positive", "negative") # label if positive or negative.
    selPd <- ifelse(isPos, posPd, negPd) # Pd that is smaller.
    
    ## Print:
    cat(paste0(iEff_str, ": ", effDir, " effect, Pd = ", selPd, "\n"))
    
  } # end of iEff for loop
  
} # end of function

# ============================================================================ #
#### Determine y-axis limits ####

determine_ylim_data_y <- function(data){
  #' Determine optimal y-axis limits based on some input heuristics
  #' @param data data frame, aggregated per subject, long format, with variable \code{y}
  #' @return yLim vector with to elements: minimal and maximal y-axis limit
  
  require(plyr)
  
  # Determine minimum and maximum:
  yMin <- min(data$y, na.rm = T)
  yMax <- max(data$y, na.rm = T)
  nRound <- 3
  cat(paste0("Detected yMin = ", round(yMin, nRound), ", yMax = ", round(yMax, nRound), "\n"))
  
  ## if likely probability
  if(yMin >= 0 & yMax <= 1){ 
    cat("Set automatic yLim to [0, 1]\n")
    yLim <- c(0, 1)
    
    ## If positive number, but not huge:
  } else if(yMin >= 0 & yMax <= 10){ # if rather small positive number: plot from 0 onwards to 110% of yMax
    cat("Set automatic yLim to [0, round(yMax*1.1,1)]\n")
    yLim <- c(0, round(yMax*1.1,1))
    
    ## If positive number and huge:
  } else if ((yMin > 0) & (yMax > 100)) { # if very big number: plot from 0 onwards to next hundred
    cat("Set automatic yLim to [0, hundres]\n")
    yMin <- 0 # from zero onwards
    yMax <- round_any(yMax, 100, f = ceiling) # round up to next hundred
    yLim <- c(yMin, yMax)
    
    ## Else (if yMin negative): enforce symmetric yLim using the bigger of yMin and yMax
  } else { # take the numerically bigger one, symmetric
    cat("Set automatic yLim to be symmetric, use bigger of yMinAbs and yMaxAbs\n")
    yMaxAbs <- ceiling(c(abs(yMin), yMax)) # round to next integer
    yMaxAbs <- yMaxAbs[1] # only first entry
    yLim <- c(-1*yMaxAbs, 1*yMaxAbs)
  }
  
  ## Check if only 2 output elements:
  if (length(yLim) < 2){stop("yLim output has less than 2 elements, please check input\n")}
  if (length(yLim) > 2){stop("yLim output has more than 2 elements, please check input\n")}
  
  return(yLim)
} 

# ============================================================================ #
#### Plot single bar plot with individual data points: ####

custom_singlebar <- function(d, selVar, yLim = c(0,1), color = "grey80",
                             xLab = "x", yLab = "y", main = "Plot"){
  
  ## Load packages:
  require(ggplot2)
  require(ggbeeswarm)
  
  ## Aggregate again with Rmisc:
  # library(Rmisc)
  # summary_d <- summarySEwithin(d,measurevar="ACC", idvar = "sID", na.rm = T)
  
  if(!(selVar %in% names(data))){stop(paste0("Variable ", selVar, " not found in data"))}
  
  lineWidth <- 1.3
  fontSize <- 30
  colAlpha <- .70
  lineWidth <- 1.5
  
  ## Add grouping variable (containing all subjects):
  nSub <- nrow(d)
  d$x <- rep(1,nSub)
  d$j <- jitter(rep(0,nrow(d)), amount=.06)
  d$xj <- d$x + d$j
  
  ## Repeat selected variable:
  d$y <- d[,selVar]
  
  p <- ggplot(data = d, aes(x = x, y = y)) # initialize
  
  ## Bars:
  p <- p + stat_summary(fun = mean, geom = "bar", fill = color,
                        color = 'black', width = 0.3, lwd = lineWidth)
  
  ## Confidence intervals: 
  p <- p + stat_summary(fun.data = mean_cl_normal, geom =
                          "errorbar", width = 0.05, lwd = lineWidth)
  
  ## Individual data points:
  p <- p + geom_beeswarm(color = "black", fill = color, alpha = colAlpha)
  # p <- p + geom_point(data = d, aes(x = xj), color = "black", fill = "grey40", 
  #                     shape = 21, size = 4,
  #                     alpha = colAlpha)
  
  ## Other settings:
  
  if (mean(yLim) == 0.5){ # if conditional probabilities:
    p <- p + geom_hline(yintercept=0.5, linetype=2, color="black")
  }
  
  ## X-axis:
  p <- p + scale_x_continuous(limits = c(0.5,1.5), breaks = c(0,1,2), labels = c("","",""))
  
  # Y-axis:
  p <- p + scale_y_continuous(limits = yLim, breaks = seq(yLim[1],yLim[-1],0.25)) 
  
  # Axis labels:
  p <- p + xlab(xLab) + ylab(yLab) +
    
    ## Add title:
    ggtitle(main) +
    
    ## Add theme:
    theme_classic()
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size=fontSize),
                 axis.title = element_text(size=fontSize), 
                 plot.title = element_text(size=fontSize, hjust = 0.5), # center title 
                 legend.text = element_text(size=fontSize))
  
  # Print plot in the end:
  print(p)
  return(p)
}

# =============================================================================================== #
#### Plot flexible raincloud plot: #####

plot_raincloud <- function(data, xVar, yVar, subVar = "subject_n",
                           isBar = F, isPoint = T, isViolin = F, isBoxplot = F, isMean = T, 
                           isMirror = F, useCond = NULL, jitterNum = 0.09,
                           yLim = NULL, color, Labels = NULL, xLab = "x", yLab = "y", main = NULL, fontSize = NULL){
  #' Make raincloud plot
  #' @param data data frame, with variables \code{variables}
  #' @param xVar string, name of variable that goes on x-axis. Factor or numeric.
  #' @param yVar string, name of variable that goes on y-axis. Variable needs to be numeric.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param isBar Boolean, plot means per condition as bar (default: FALSE)
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: TRUE)
  #' @param isViolin Boolean, plot distribution per condition as violin plot (default: TRUE)
  #' @param isBoxplot Boolean, plot median and interquartile per condition with boxplot (default: FALSE)
  #' @param isMean Boolean, plot means per condition as thick connected points with error bars (default: TRUE) 
  #' @param isMirror Determine position of boxplots and violin plots: always on the left (FALSE, default) or mirred to the middle (TRUE)
  #' @param useCond vector of numbers, use only subset of input conditions (default: all conditions)
  #' @param jitterNum numeric, amount of jitter to use for points (default: 0.9).
  #' @param isMirror Determine position of boxplots and violin plots: always on the left (FALSE, default) or mirred to the middle (TRUE)
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data).
  #' @param color vector of strings (HEX colors), colors for input conditions (d$x)
  #' @param Labels vector of strings, labels for input conditions (d$x) on the x-axis
  #' @param xLab string, label for x-axis (default: "x")
  #' @param yLab string, label for y-axis (default: "y")
  #' @param main string, title of plot (default: "NULL")
  #' @param fontSize integer, font size for axes ticks and labels and title.
  #' @return makes raincloud plot
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.Rmd
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.pdf
  
  require(ggplot2)
  require(ggstatsplot) # for geom_flat_violin
  require(gghalves) # for half plots
  
  # General settings:
  colAlpha <- .50
  lineWidth <- 1.5
  
  # ---------------------------------------------------------------------------- #
  ### Check for NAs in xVar:
  
  if (sum(is.na(data[,xVar]))>0){
    cat("xVar contains NAs, remove\n")
    data <- data[!(is.na(data[,xVar])),]
  }
  
  # ---------------------------------------------------------------------------- #
  ### Aggregate data per condition per subject into long format:
  
  d <- prepare_data_plot_long(data, x = xVar, y = yVar, subVar = subVar, jitterNum = jitterNum)
  
  # ---------------------------------------------------------------------------- #
  ### Determine defaults:
  
  cat("Set defaults\n")
  
  # Determine y limits if not given:
  if(is.null(yLim)){
    cat("Automatically determine y-axis\n")
    yLim <- determine_ylim_data_y(d)
  }
  
  # Determine number of conditions to use:
  if(is.null(useCond)){ # if not defined: use all conditions
    nCond <- length(unique(d$x))
    useCond <- 1:nCond #  sort(unique(d$x))
    d$condition <- as.numeric(as.factor(d$x)) # from 1 to nCond
  } else {
    cat(paste0("Use only conditions ",toString(useCond)),"\n")
    nCond <- length(useCond)
    # Loop through selected conditions, recompute x so it is consecutive:
    ix <- 0 # initialize
    for (iCond in useCond){
      ix <- ix + 1 # increment
      d[which(d$condition==iCond),"x"] <- ix
    } 
    d$xj <- d$x + d$j # recompute xj
  }
  
  if (is.null(Labels)){
    cat("Automatically extract x-axis labels\n")
    Labels <- sort(unique(d$x))
  }
  
  if(length(color)!=nCond){
    cat("Less colors provided than conditions, repeat first color for all conditions\n")
    color <- rep(color[1],nCond)
  }
  if(length(Labels)!=nCond){
    cat("Less labels provided than conditions, use just numbers from 1 to number of conditions\n")
    Labels <- 1:nCond
  }
  
  if (is.null(fontSize)){
    ## Font sizes for ordinary viewing: 15
    # fontSize <- 15
    ## Font sizes for saving: 30
    fontSize <- 30
    cat(paste0("No font size provided, use font size ",fontSize),"\n")
  }
  
  # --------------------------------------------------------------------- #
  ## Compute summary per condition:
  require(Rmisc)
  dsel = subset(d, condition %in% useCond) # select used conditions:
  summary_d <- summarySEwithin(dsel, measurevar="y", withinvar = "x", idvar = "sID", na.rm = T)
  # automatically rescales x to 1:nCond
  # summary_d$x <- as.numeric(summary_d$x)
  if (is.numeric(d$x)){summary_d$x <- sort(unique(d$x))} # back to original numerical values
  summary_d$condition <- as.numeric(as.factor(summary_d$x)) # consecutive count (for labels)
  
  # --------------------------------------------------------------------- #
  ## Determine position of boxplots and violin plots: always on the left (FALSE, default) or mirred to the middle (TRUE)
  if(isMirror == T){ # if mirror images for box plots and violin plots
    sideVec <- c(rep("l",floor(nCond/2)),rep("r",ceiling(nCond/2)))
  } else { # otherwise default
    sideVec <- rep("l",nCond)
  }
  
  # --------------------------------------------------------------------- #
  ## Start ggplot: 
  p <- ggplot(data = d, aes(y = y)) # initialize
  
  # --------------------------------------------------------------------- #
  ## Bars:
  if (isBar == T){
    cat("Make bar plot\n")
    
    for(iCond in useCond){
      p <- p + geom_bar(data = subset(summary_d, condition == iCond), aes(x = x, y = y), stat = "identity",
                        fill = color[iCond], col = "black", width = 0.4, lwd = lineWidth, alpha = 0.3)
      p <- p + geom_errorbar(data = subset(summary_d, condition==iCond), 
                             aes(x = x, y = y, ymin = y-ci, ymax = y+ci),
                             color = "black", width = 0.10, lwd = lineWidth, alpha = 0.6) 
    }
  }
  
  # --------------------------------------------------------------------- #
  ## Individual data points:
  if(isPoint == T){
    cat("Make individual data points\n")
    
    for(iCond in useCond){
      p <- p + geom_point(data = subset(d, condition == iCond), aes(x = xj), color = color[iCond], size = 1.5, 
                          alpha = .35) # alpha = .50
    }
    # Add lines to combine points:
    p <- p + geom_line(data = subset(d, condition %in% useCond), aes(x = xj, group = subject), 
                       size = 1.0, color = 'grey40', alpha = 0.35) # lightgray
  }
  # Till 2021-11-22: color = grey70, size = 1
  
  # --------------------------------------------------------------------- #
  ## Box plots:
  if (isBoxplot == T){
    
    cat("Make box plot\n")
    
    for(iCond in useCond){
      p <- p + geom_half_boxplot(
        data = subset(d, condition == iCond), aes(x=x, y = y), position = position_nudge(x = .15), 
        side = sideVec[iCond], outlier.shape = NA, center = TRUE, errorbar.draw = TRUE, width = .1, 
        fill = color[iCond], alpha = colAlpha)
    }
  }
  
  # Violin plots:
  if (isViolin == T){
    
    cat("Make violin plot\n")
    
    for(iCond in useCond){
      p <- p + geom_half_violin(
        data = subset(d, condition == iCond),aes(x = x, y = y), position = position_nudge(x = -.12),
        side = sideVec[iCond], fill = color[iCond], alpha = colAlpha, trim = FALSE)
    }
  }
  
  # --------------------------------------------------------------------- #
  ## Mean and CI per condition:
  posNudge = 0 # .13
  if (isMean == T){
    
    cat("Make point for mean\n")
    
    ## Thick line connecting means (plot line first and points on top):
    p <- p + geom_line(data = summary_d, aes(x = x, y = y), color = color[iCond],
                       position = position_nudge(x = 1*posNudge), size = 1.5)
    
    ## Error shades:
    p <- p + geom_ribbon(data = summary_d, aes(x = x, y = y, ymin = y-ci, ymax = y+ci),
                         color = color[iCond], alpha=0.15, linetype = 0)
    
    ## Points for means and error bars:
    for(iCond in useCond){
      p <- p + geom_point(data = subset(summary_d, condition==iCond), aes(x = x, y = y), # point
                          position = position_nudge(x = -1*posNudge), color = color[useCond[iCond]], alpha = .6, size = 3) # size = 2
      
      p <- p + geom_errorbar(data = subset(summary_d, condition==iCond), aes(x = x, y = y, # error bar
                                                                             ymin = y-ci, ymax = y+ci),
                             position = position_nudge(-1*posNudge), color = color[useCond[iCond]], width = 0.05, size = 1.5, alpha = .6)
    }
  }
  
  if (mean(yLim) == 0.5){ # if conditional probabilities:
    p <- p + geom_hline(yintercept=0.5, linetype=2, color="black", size = 1) # dotted line in the middle
  }
  if (yLim[1] == 0 & yLim[2] == 1){
    # p <- p + scale_y_break(c(0, 0.5, 1))
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  
  # --------------------------------------------------------------------- #
  cat("Adjust axes, labels\n")
  
  # Y-axis:
  # p <- p + scale_y_continuous(limits = yLim) 
  p <- p + coord_cartesian(ylim=yLim) 
  
  # Labels:
  p <- p + scale_x_continuous(breaks=sort(unique(d$x)), labels=Labels[useCond]) +
    xlab(xLab) + ylab(yLab) + theme_classic()
  
  # Other settings:
  if (!(is.null(main))){
    cat("Add title\n")
    p <- p + ggtitle(main) # title off for printing for poster
  }  
  
  # p + theme_classic()
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size=fontSize),
                 axis.title = element_text(size=fontSize), 
                 plot.title = element_text(size=fontSize, hjust = 0.5), # center title 
                 legend.text = element_text(size=fontSize))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
  return(p)
}

# =============================================================================================== #
#### Plot beeswarm plot: #####

plot_beeswarm <- function(data, xVar, yVar, subVar = "subject_n",
                          nRound = NULL, yLim = NULL, color, Labels, 
                          xLab = "x", yLab = "y", main = "Plot", fontSize = NULL){
  #' Make raincloud plot
  #' @param data data frame, with variables \code{variables}
  #' @param xVar string, name of variable that goes on x-axis. Factor or numeric.
  #' @param yVar string, name of variable that goes on y-axis. Variable needs to be numeric.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param nRound numeric, digits used for rounding dependent variable y.
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data)
  #' @param color vector of strings (HEX colors), colors for input conditions (d$x)
  #' @param Labels vector of strings, labels for input conditions (d$x) on the x-axis
  #' @param xLab string, label for x-axis (default: "x")
  #' @param yLab string, label for y-axis (default: "y")
  #' @param main string, title of plot (default: "Plot")
  #' @param fontSize integer, font size for axes ticks and labels and title.
  #' @return makes beeswarm and box plot for each condition.
  
  require(Rmisc)
  require(ggplot2)
  require(ggbeeswarm)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  ## General settings:
  colAlpha <- .20
  lineWidth <- 1.5
  
  ## Font sizes for ordinary viewing and saving:
  if (is.null(fontSize)){
    ## Font sizes for ordinary viewing: 15
    # fontSize <- 15
    ## Font sizes for saving: 30
    fontSize <- 30
    cat(paste0("No font size provided, use font size ",fontSize),"\n")
  }
  
  # ---------------------------------------------------------------------------- #
  ### Aggregate data per condition per subject into long format:
  
  d <- prepare_data_plot_long(data, x = xVar, y = yVar, subVar = subVar)
  
  # ---------------------------------------------------------------------------- #
  ### Determine defaults:
  
  cat("Set defaults\n")
  
  # Determine y limits if not given:
  if(is.null(yLim)){
    yLim <- determine_ylim_data_y(d)
  }
  
  nCond <- length(unique(d$condition))
  condNames <- sort(unique(d$condition))
  
  # --------------------------------------------------------------------- #
  ## Compute summary per condition:
  
  summary_d <- summarySEwithin(d, measurevar = "y", withinvar = "x", idvar = "sID", na.rm = T)
  # automatically rescales x to 1:nCond
  # summary_d$x <- as.numeric(summary_d$x)
  if (is.numeric(d$x)){summary_d$x <- sort(unique(d$x))} # back to original numerical values
  summary_d$condition <- as.numeric(as.factor(summary_d$x)) # consecutive count (for labels)
  
  if (!is.null(nRound)){
    cat(paste0("Round DV on ",nRound," digits"))
    d$y <- round(d$y, nRound)
    # summary_d$y <- round(summary_d$y, nRound)
  }
  
  # --------------------------------------------------------------------- #
  ## Start ggplot: 
  
  p <- ggplot(data = d, aes(x = x, y = y)) # initialize
  
  # --------------------------------------------------------------------- #
  ## Make beeswarm:
  
  cat("Make beeswarm plot\n")
  
  p <- p + geom_beeswarm()
  
  # --------------------------------------------------------------------- #
  ## Box plots:
  
  cat("Make box plot\n")
  
  for(iCond in 1:nCond){
    p <- p + geom_boxplot(
      data = subset(d, condition == condNames[iCond]), aes(x = x, y = y), 
      width = 0.3,
      color = color[iCond], fill = color[iCond], alpha = colAlpha)
    
    # position = position_nudge(x = .15), 
    # side = sideVec[iCond], outlier.shape = NA, center = TRUE, errorbar.draw = TRUE, width = .1, 
    # fill = color[iCond], alpha = colAlpha)
  }
  
  # --------------------------------------------------------------------- #
  cat("Adjust axes, labels\n")
  
  # Y-axis:
  # p <- p + scale_y_continuous(limits = yLim) 
  p <- p + coord_cartesian(ylim=yLim) 
  
  # Labels:
  p <- p + scale_x_continuous(breaks=sort(unique(d$x)), labels=Labels) +
    xlab(xLab) + ylab(yLab) +
    
    # Other settings:
    ggtitle(main) + # title off for printing for poster
    theme_classic()
  
  ## Font sizes:
  p <- p + theme(axis.text=element_text(size=fontSize),
                 axis.title=element_text(size=fontSize), 
                 title = element_text(size=fontSize), 
                 legend.text = element_text(size=fontSize))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
}

# ============================================================================ #
#### Barplot 1 IV: Aggregate per condition per subject, plot (1 IV on x-axis): ####

custom_barplot1 <- function(data, xVar = NULL, yVar = NULL, subVar = "subject_n", 
                            xLab = "Condition", yLab = "p(Go)", main = NULL, selCol = "red",
                            isPoint = F, isBeeswarm = F, yLim = NULL, fontSize = NULL,
                            savePNG = F, saveEPS = F, prefix = NULL, suffix = NULL){
  #' Make bar plot with error bars and individual-subject data points.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. If numeric, it will be converted to an (ordered) factor.
  #' @param yVar string, name of variable that goes on y-axis. Needs to be numeric.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param main string, overall plot label (optional).
  #' @param selCol vector of strings (HEX colors), colors for bars (default: "red").
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: FALSE).
  #' @param isBeewswarm Boolean, plot individual data points per condition as beeswarm densities (default: FALSE).
  #' @param yLim vector of two numbers, y-axis (default: automatically determined by ggplot).
  #' @param fontSize scalar integer, font size to use.
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @param prefix string, string to add at the beginning of plot name (optional).
  #' @param suffix string, string to add at the end of plot name (optional).
  #' @return creates (and saves) plot.
  
  # -------------------------------------------------------------------------- #
  ## Load required packages:
  
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Check inputs:
  
  if(!(yVar %in% names(data))){stop("yVar not found in data")}
  if(!(xVar %in% names(data))){stop("xVar not found in data")}
  if(!(subVar %in% names(data))){stop("subVar not found in data")}
  
  # -------------------------------------------------------------------------- #
  ## Fixed plotting settings:
  
  SEweight <- 1
  lineWidth <- 1.5
  dodgeVal <- 0.6
  colAlpha <- 0.5 # 1
  
  if (is.null(fontSize)){
    fontSize <- 30 # 30 or 15?
  }
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  
  data$x <- data[,xVar]
  data$y <- data[,yVar]
  data$subject <- data[,subVar]
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data per subject per condition:
  
  aggrData <- ddply(data, .(subject, x), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  cat(paste0("Min = ", round(min(aggrData$y), 3), "; Max = ", round(max(aggrData$y), 3)), "\n")
  
  # -------------------------------------------------------------------------- #
  ## Add jittered x-axis variable for points:
  
  aggrData$xpos <- as.numeric(aggrData$x) # to numeric
  aggrData$xpos <- aggrData$xpos - min(aggrData$xpos) + 1 # plot starts at lowest level of x-variable
  aggrData$j <- jitter(rep(0, nrow(aggrData)), amount = .09) # jitter 0.09
  aggrData$xj <- aggrData$xpos + aggrData$j # add jitter
  
  # -------------------------------------------------------------------------- #
  ## Determine y limits if not given:
  
  if(is.null(yLim)){
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate across subjects with Rmisc:
  
  summary_d <- summarySEwithin(aggrData, measurevar = "y", idvar = "subject", na.rm = T,
                               withinvars = c("x"))
  
  # -------------------------------------------------------------------------- #
  ## Control settings for saving:
  
  ## Additions:
  if(is.null(prefix)){prefix <- ""} else {prefix <- paste0(prefix, "_")}
  if(is.null(suffix)){suffix <- ""} else {suffix <- paste0("_", suffix)}
  
  ## Name:
  plotName <- paste0("custombarplot1_", prefix, yVar, "_", xVar)
  if (isPoint){plotName <- paste0(plotName, "_points")} 
  plotName <- paste0(plotName, suffix)
  cat(paste0("Start plot ", plotName, "\n"))
  
  # Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)}
  
  # -------------------------------------------------------------------------- #
  ## Start ggplot:
  p <- ggplot(summary_d,aes(x, y))
  
  ## Bars of means:
  cat("Add group-level bars \n")
  p <- p + stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.6,
                        lwd = lineWidth, fill = selCol, color = "black") + 
    
    ## Error bars:
    cat("Add error bars \n")
  p <- p + geom_errorbar(data = summary_d,
                         aes(x = x, y = y, ymin = y - se * SEweight, ymax = y + se * SEweight),
                         position = position_dodge(width = dodgeVal), width = 0.1,
                         lwd = lineWidth, color = "black", alpha = 1)
  
  # -------------------------------------------------------------------------- #
  ## Individual data points:
  
  cat("Add per-subject data points\n")
  if (isPoint){
    ## Colored dots:
    p <- p + geom_point(data = aggrData, aes(x = xj, fill = x), shape = 21, size = 2, stroke = 1.2, # size = 0.6, 
                        color = "black", alpha = colAlpha)
    p <- p + scale_fill_manual(values = selCol)
    ## Grey dots:
    # p <- p + geom_point(data = aggrData, aes(x = xj), shape = 21, size = 2, fill = NA, stroke = 1.5, # size = 0.6, 
    #                     color = "grey", alpha = colAlpha) # color = black, grey60,
  }
  
  if (isBeeswarm){
    p <- p + geom_beeswarm(data = aggrData, aes(x = xpos), shape = 1, size = 2, stroke = 1, # size = 0.6, 
                           color = "black", alpha = colAlpha)
  }
  
  # -------------------------------------------------------------------------- #
  ## Settings:
  
  if (yLim[1] == 0 & yLim[2] == 1){
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  if(!(is.null(yLim))){p <- p + coord_cartesian(ylim=yLim)}
  
  ## Axis labels:
  p <- p + labs(x=xLab, y = yLab)
  
  # Add title:
  if (!(is.null(main))){
    cat("Add title\n")
    p <- p + ggtitle(main)  
  }
  
  ## Theme:
  p <- p + theme_classic()
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = fontSize),
                 axis.title = element_text(size = fontSize), 
                 title = element_text(size = fontSize),
                 legend.position = "none",
                 axis.line = element_line(colour = 'black', size = lineWidth)) # fixed font sizes
  print(p)
  if(savePNG | saveEPS){
    dev.off(); 
    print(p)
  }
  return(p)
  cat("Finished :-)\n")
}

# ============================================================================ #
#### Barplot 2 IVs: Aggregate per condition per subject, plot (2 IVs, 1 on x-axis, 1 as color): ####

custom_barplot2 <- function(data, xVar, yVar, zVar, subVar = "subject_n", 
                            xLab = "Condition", yLab = "p(Go)", zLab = "Action", main = NULL,
                            selCol = c("blue","red"), isPoint = F, isBeeswarm = F, yLim = NULL, 
                            savePNG = F, saveEPS = F, prefix = NULL, suffix = NULL){
  #' Make bar plots with 2 IVs: x-axis and color.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. If numeric, it will be converted to an (ordered) factor.
  #' @param yVar string, name of variable that goes on y-axis. Needs to be numeric.
  #' @param zVar string, name of variable that determines bar coloring. Needs to be a factor.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param zLab string, label for color legend (default: "z").
  #' @param main string, overall plot label (optional).
  #' @param selCol vector of strings (HEX colors), colors for input levels of zVar (default: c("blue","red")).
  #' @param yLim vector of two numbers, y-axis (default: automatically determined by ggplot).
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: TRUE).
  #' @param isBeeswarm Boolean, plot individual data points per condition as beeswarm density (default: FALSE).
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @param prefix string, string to add at the beginning of plot name (optional).
  #' @param suffix string, string to add at the end of plot name (optional).
  #' @return creates (and saves) plot.
  
  # -------------------------------------------------------------------------- #
  ## Load packages:
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  require(ggbeeswarm) # for ggbeeswarm
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Check inputs:
  
  if(!(yVar %in% names(data))){stop("yVar not found in data")}
  if(!(xVar %in% names(data))){stop("xVar not found in data")}
  if(!(zVar %in% names(data))){stop("zVar not found in data")}
  if(!(subVar %in% names(data))){stop("subVar not found in data")}
  if(length(selCol) != length(unique(data[,zVar]))){stop("Length selCol and number levels zVar do not match")}
  
  if(is.numeric(data[, xVar])){data[, xVar] <- factor(data[, xVar]); cat("Convert xVar to factor\n")}
  
  # -------------------------------------------------------------------------- #
  ## Fixed plotting settings:
  
  lineWidth <- 1.5 # 1.3
  fontSize <- 30 # 30
  if (length(unique(data[, xVar])) > 10){fontSize <- 15; cat("Lower fontSize from 30 to 15 because many x levels\n")}
  dodgeVal <- 0.6
  colAlpha <- 1
  condCol <- rep(selCol, length(unique(data[,xVar])))
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  
  cat("Create new variables x, y, z, subject based on inputs\n")
  data$x <- data[, xVar]
  data$y <- data[, yVar]
  data$z <- data[, zVar]
  data$subject <- data[, subVar]
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data per subject per condition:
  cat("Aggregate data per subject\n")
  aggrData <- ddply(data, .(subject, x, z), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  # Wide format: each subject/condition1/condition2 in one line, variables subject, x, y, z
  
  ## Add condition variable:
  cat("Create condition variable\n")
  nZlevel <- length(unique(data$z))
  posScale <- 0.05 * (nZlevel + 1)
  aggrData$cond <- as.numeric(aggrData$x)*nZlevel - nZlevel + as.numeric(aggrData$z) # necessary for proper axis positions
  # aggrData$cond <- 1 + as.numeric(aggrData$x)*2 + as.numeric(aggrData$z)
  nCond <- length(unique(aggrData$cond))
  if (length(condCol) < nCond){condCol <- rep(condCol, length.out = nCond)}
  
  ## Add jittered x-axis for points:
  cat("Add jitter for points\n")
  aggrData$j <- jitter(rep(0, nrow(aggrData)), amount = .05) # pure jitter .05
  if (nZlevel == 2){
    aggrData$xpos <- as.numeric(aggrData$x) + (as.numeric(aggrData$z) - 1.5) * 2 * posScale # convert to [1 2], to [-0.5 0.5], * 2 so [-1 1], scale by 0.15
  } else if (nZlevel == 3) {
    zMid <- round(mean(as.numeric(data$z)))
    aggrData$xpos <- as.numeric(aggrData$x) + ((as.numeric(aggrData$z) - zMid)) * posScale  # demean, scale by 0.20
  } else {
    zMid <- round(mean(as.numeric(data$z)))
    zScale <- ceiling(max(as.numeric(data$z))) - zMid
    aggrData$xpos <- as.numeric(aggrData$x) + ((as.numeric(aggrData$z) - zMid)) * zScale * posScale  # demean, bring min/max to 1, scale by 0.20
    warning(paste0("Not yet implement for z variable with ", nZlevel, " levels\n"))
  }
  aggrData$xj <- aggrData$xpos + aggrData$j # add jitter to xpos
  
  ## Determine y limits if not given:
  if(is.null(yLim)){
    cat("Automatically y-axis limits based on per-subject-per-condition means\n")
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate across subjects with Rmisc:
  cat("Aggregate data across subjects\n")
  summary_d <- summarySEwithin(aggrData, measurevar = "y", idvar = "subject", na.rm = T,
                               withinvars = c("x", "z"))
  # Aggregated over subjects, one row per condition, variables x, z, N, y, sd, se, ci
  
  # -------------------------------------------------------------------------- #
  ### Start plot:
  
  ## Additions:
  if(is.null(prefix)){prefix <- ""} else {prefix <- paste0(prefix, "_")}
  if(is.null(suffix)){suffix <- ""} else {suffix <- paste0("_", suffix)}
  
  ## Name:
  plotName <- paste0("custombarplot2_", prefix, yVar, "_", xVar, "_", zVar)
  if (isPoint){plotName <- paste0(plotName, "_points")} 
  if (isBeeswarm){plotName <- paste0(plotName, "_beeswarm")} 
  plotName <- paste0(plotName, suffix)
  cat(paste0("Start plot ", plotName, "\n"))
  
  ## Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)}
  
  ## Start plot:
  p <- ggplot(summary_d, aes(x, y, fill = z))
  
  ## Bars of means:
  cat("Add bars\n")
  p <- p + stat_summary(fun = mean, geom = "bar", position = "dodge", width = dodgeVal,
                        lwd = lineWidth, color = "black")
  
  ## Error bars:
  cat("Add error bars\n")
  p <- p + geom_errorbar(data = summary_d,
                         aes(x = x, y = y, ymin = y - se, ymax = y + se),
                         position = position_dodge(width = dodgeVal), width = 0.2,
                         lwd = lineWidth, color = "black", alpha = 1)
  
  # Individual data points:
  if (isPoint){
    cat("Start adding per-subject points \n")
    for(iCond in 1:nCond){ # add separately per condition
      p <- p + geom_point(data = aggrData[aggrData$cond == iCond, ],
                          aes(x = xj), # position = "dodge",
                          shape = 21, size = 2, stroke = 1.2, color = "black", fill = condCol[iCond],
                          alpha = 0.5) # colAlpha)
    }
  }
  
  # Beeswarm style plots:
  if (isBeeswarm){
    cat("Start adding beeswarm \n")
    for(iCond in 1:nCond){ # add separately per condition
      p <- p + geom_beeswarm(data = aggrData[aggrData$cond == iCond, ],
                             aes(x = xpos), # position = "dodge",
                             # priority = "ascending",
                             shape = 21, size = 2, stroke = 1.2, color = "black", fill = condCol[iCond],
                             alpha = 0.5) # colAlpha)
    }
  }
  
  # Add title:
  if (!(is.null(main))){
    cat("Add title\n")
    p <- p + ggtitle(main)  
  }
  
  # Settings:
  if (yLim[1] == 0 & yLim[2] == 1){
    cat("Add y-axis ticks for 0, 0.5, 1\n")
    # p <- p + scale_y_break(c(0, 0.5, 1))
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  if(!(is.null(yLim))){p <- p + coord_cartesian(ylim = yLim)}
  # if(!(is.null(yLim))){p <- p + scale_y_continuous(limits = yLim, breaks = seq(yLim[1], yLim[-1], (yLim[-1] - yLim[1])/2))}
  
  # Add theme, fontsizes:
  cat("Add axis labels, colors, theme, font sizes\n")
  require(ggthemes)
  p <- p + labs(x = xLab, y = yLab, fill = zLab) +
    scale_fill_manual(values = selCol) + 
    theme_classic() + 
    theme(axis.text = element_text(size = fontSize),
          axis.title = element_text(size = fontSize), 
          plot.title = element_text(size = fontSize, hjust = 0.5), 
          legend.text = element_text(size = fontSize),
          legend.title = element_text(size = fontSize),
          axis.line = element_line(colour = 'black', size = lineWidth)) # fixed font sizes
  #          legend.title = element_blank(), legend.position = "none")
  
  print(p)
  if(savePNG | saveEPS){dev.off(); print(p)}
  cat("Finished :-)\n")
  return(p)
}

# ============================================================================ #
#### Barplot 3 IVs: Aggregate per condition per subject, plot (3 IVs, 1 on x-axis, 1 as color, 1 as facets): ####

custom_barplot3 <- function(data, yVar, xVar, zVar, splitVar = NULL, subVar = "subject_n", 
                            yLab = "DV", xLab = "IV1", zLab = "IV2", main = NULL,
                            xLevels = c("Go", "NoGo"), zLevels = c("Win", "Avoid"), splitLevels = c("On", "Off"),
                            selCol = greenred, isPoint = F, yLim = NULL, savePNG = F, saveEPS = F){
  #' Make bar plots with 3 IVs: x-axis and color and facetwrap.
  #' Can add points with geom_point, 
  #' but not beeswarm plots because no position argument (hence no dodge) and manual x-position not compatible with facet_wrap.
  #' @param data data frame, trial-by-trial data.
  #' @param yVar string, name of variable that goes on y-axis. Needs to be numeric.
  #' @param xVar string, name of variable that goes on x-axis. If numeric, it will be converted to an (ordered) factor.
  #' @param zVar string, name of variable that determines bar coloring. Needs to be a factor.
  #' @param splitVar string, name of variable by which to split plot (facetwrap, optional).
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param yLab string, label for y-axis (default: "y").
  #' @param xLab string, label for x-axis (default: "x").
  #' @param zLab string, label for color legend (default: "z").
  #' @param main string, title of plot (optional).
  #' @param selCol vector of strings (HEX colors), colors for input levels of zVar (default: greenred).
  #' @param yLim vector of two numbers, y-axis (default: automatically determined by ggplot).
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: TRUE).
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @return creates (and saves) plot.  #' Make bar plot per subject
  #' @param data data frame, with variables \code{variables}
  
  # -------------------------------------------------------------------------- #
  ## Load packages:
  
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  require(ggbeeswarm) # for ggbeeswarm
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Check inputs:
  
  if(!(yVar %in% names(data))){stop("yVar not found in data")}
  if(!(xVar %in% names(data))){stop("xVar not found in data")}
  if(!(zVar %in% names(data))){stop("zVar not found in data")}
  if(!(subVar %in% names(data))){stop("subVar not found in data")}
  if(!is.null(splitVar) & !(subVar %in% names(data))){stop("splitVar not found in data")}
  if(length(selCol) != length(unique(data[,zVar]))){stop("Length selCol and number levels zVar do not match")}
  
  if (is.null(xLab)){xLab <- str_to_title(xVar)}
  if (is.null(zLab)){zLab <- str_to_title(zVar)}
  
  # -------------------------------------------------------------------------- #
  ## General settings:
  
  lineWidth <- 1.3 # 1.5
  if (savePNG | saveEPS){fontSize <- 24} else {fontSize <- 15} # 30
  dodgeVal <- 0.6
  barWidth <- 0.15
  colAlpha <- 0.6
  
  SEweight <- 1.96
  
  condCol <- rep(selCol, length(unique(data[,xVar])))
  
  ## Set default y-axis limits:
  if (is.null(yLim)){
    yLim <- c(0, 1)
  }
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  cat("Create new variables x, y, z, subject based on inputs\n")
  
  data$y <- data[, yVar]
  data$x <- as.numeric(data[, xVar])
  data$z <- as.numeric(data[, zVar])
  if (!is.null(splitVar)){data$split <- data[, splitVar]}
  data$subject <- data[, subVar]
  
  if(sum(is.na(data$x)) > 0){stop("Found NAs in x variable\n")}
  if(sum(is.na(data$z)) > 0){stop("Found NAs in z variable\n")}
  if(!is.null(splitVar) & sum(is.na(data$split)) > 0){stop("Found NAs in split variable\n")}
  
  # table(data$FRR_f, data$x) # poor becomes 2, mid becomes 1, rich becomes 0
  # table(data$BRR_f, data$z) # poor becomes 1, rich becomes 0
  # table(data$valence_f, data$split) # gain becomes 1, loss becomes 0
  
  ## Recode to 0, 1 if necessary (for combining into condition variable):
  
  xLevelVec <- sort(unique(data$x))
  nXlevels <- length(xLevelVec)
  if (nXlevels > 3){stop("Plot not implemented for x variable with > 3 levels")}
  
  zLevelVec <- sort(unique(data$z))
  nZlevels <- length(zLevelVec)
  if (nZlevels > 2){stop("Plot not implemented for z variable with > 2 levels")}
  
  ## Recode variables to go from (nLevels - 1) to zero:
  if(nXlevels == 2){if(all(xLevelVec == c(1, 2))){data$x <- nXlevels - data$x}} # to 0 or 1
  if(nXlevels == 3){if(all(xLevelVec == c(1, 2, 3))){data$x <- nXlevels - data$x}} # to 0 or 1 or 2
  if(all(zLevelVec == c(1,2))){data$z <- nZlevels - data$z} # to 0 or 1,
  xLevelVec <- rev(sort(unique(data$x))) # update, descending order
  zLevelVec <- rev(sort(unique(data$z))) # update, descending order
  
  if (!is.null(splitVar)){
    sLevelVec <- sort(unique(data$split))
    nSlevels <- length(sLevelVec)
    if (nSlevels > 3){stop("Plot not implemented for split variable with > 2 levels")}
    if(all(sLevelVec == c(1,2))){data$split <- nSlevels - data$split}
    sLevelVec <- rev(sort(unique(data$split))) # update, descending order
  } # to 0 - 1
  
  # sort(unique(data$x)) # should be 0, 1 (, 2)
  # sort(unique(data$z)) # should be 0, 1
  # sort(unique(data$split)) # should be 0, 1
  
  # -------------------------------------------------------------------------- #
  ## Combine into single condition variable:
  
  if (!is.null(splitVar)){ # if splitVar: 8 conditions
    data$condition <- 1 + data$split*nXlevels*nZlevels + data$x*nZlevels + data$z
  } else { # No splitVar: 4 conditions
    data$condition <- 1 + data$x*nZlevels + data$z
  }
  # sort(unique(data$condition)) # from 1 to nCond
  nCond <- length(unique(data$condition))
  stopifnot(min(data$condition)  == 1)
  stopifnot(max(data$condition)  == nCond)
  
  # table(data$condition, data$x) # always nZlevels after each other
  # table(data$condition, data$z) # odd and even
  # table(data$condition, data$split) # 1 for upper half, 0 for lower half
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data per subject per condition:
  
  cat("Aggregate data per subject\n")
  aggrData <- ddply(data, .(subject, condition), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  
  ## Recover original variables:
  aggrData$x <- (ceiling(aggrData$condition/nZlevels) - 1) %% nXlevels
  aggrData$z <- aggrData$condition %% nZlevels
  if (!is.null(splitVar)){aggrData$split <- 1 - ceiling(aggrData$condition/(nXlevels*nZlevels)) %% nSlevels}
  head(aggrData, n = nCond)
  
  ## Create factors:
  aggrData$x_f <- factor(aggrData$x, levels = xLevelVec, labels = xLevels)
  aggrData$z_f <- factor(aggrData$z, levels = zLevelVec, labels = zLevels)
  if (!is.null(splitVar)){aggrData$split_f <- factor(aggrData$split, levels = sLevelVec, labels = splitLevels)}
  
  cat(paste0("Assume ", paste0(xLevelVec, collapse = ", "), " corresponds to ", paste0(xLevels, collapse = ", "), "\n"))
  cat(paste0("Assume ", paste0(zLevelVec, collapse = ", "), " corresponds to ", paste0(zLevels, collapse = ", "), "\n"))
  if (!is.null(splitVar)){cat(paste0("Assume ", paste0(sLevelVec, collapse = ", "), " corresponds to ", paste0(splitLevels, collapse = ", "), "\n"))}
  
  aggrData
  # table(data$FRR_f, data$x) # poor becomes 2, mid becomes 1, rich becomes 0
  # table(data$BRR_f, data$z) # poor becomes 1, rich becomes 0
  # table(data$valence_f, data$split) # gain becomes 1, loss becomes 0
  
  # -------------------------------------------------------------------------- #
  ## Determine y limits if not given:
  
  if(is.null(yLim)){
    cat("Automatically y-axis limits based on per-subject-per-condition means\n")
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate across subjects with Rmisc:
  cat("Aggregate data across subjects\n")
  
  d <- summarySEwithin(aggrData, measurevar = "y", withinvar = "condition", idvar = "subject", na.rm = T)
  d$condition <- as.numeric(as.factor(d$condition)) # condition back to numeric
  
  ## Recover independent variables from condition:
  d$x <- ceiling(d$condition/nZlevels - 1) %% nXlevels
  d$z <- d$condition %% nZlevels
  if (!is.null(splitVar)){d$split <- 1 - ceiling(d$condition/(nXlevels*nZlevels)) %% nSlevels}
  
  ## Create factors:
  d$x_f <- factor(d$x, levels = xLevelVec, labels = xLevels)
  d$z_f <- factor(d$z, levels = zLevelVec, labels = zLevels)
  if (!is.null(splitVar)){d$split_f <- factor(d$split, levels = sLevelVec, labels = splitLevels)}
  
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  ## Start plot:
  
  ## Name:
  plotName <- paste0("custombarplot3_",  yVar, "~", xVar, "_", zVar, "_")
  if (!(is.null(splitVar))){plotName <- paste0(plotName, splitVar)}
  if (isPoint){plotName <- paste0(plotName, "_points")} 
  cat(paste0("Start plot ", plotName, "\n"))
  
  ## Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)}
  
  ## Initialize ggplot object:
  p <- ggplot(d, aes(x = x_f, y = y, fill = z_f))
  
  ## Bars of means:
  cat("Add bars\n")
  p <- p + geom_bar(position = "dodge", stat = "summary", fun = "identity",
                    color = "black", width = dodgeVal, lwd = lineWidth)
  
  ## Error bars:
  cat("Add error bars\n")
  p <- p + geom_errorbar(data = d,
                         aes(x = x_f, y = y, ymin = y - se, ymax = y + se),
                         position = position_dodge(width = dodgeVal), width = barWidth,
                         lwd = lineWidth, color = "black", alpha = 1)
  
  ## Individual data points:
  if (isPoint){
    cat("Start adding per-subject points \n")
    p <- p + geom_point(data = aggrData,
                        position = position_dodge(width = dodgeVal),
                        shape = 21, size = 2, stroke = 1.2, color = "black",
                        alpha = 0.5) 
  }
  
  ## Facet wrap:
  if (!is.null(splitVar)){
    cat("Start adding facet_wrap\n")
    p <- p + facet_wrap(vars(split_f))
  }
  
  ## Y axis limits:
  cat("Start adding y-axis ticks\n")
  p <- p + scale_y_continuous(limits = yLim, breaks = seq(yLim[1], yLim[-1], (yLim[-1] - yLim[1])/2)) 
  
  ## Add labels:
  cat("Start labels\n")
  p <- p + labs(x = xLab, fill = zLab,
                y = yLab)
  
  ## Add color:
  cat("Start colors for fill\n")
  p <- p + scale_fill_manual(values = rep(selCol, 4))
  
  ## Add theme:
  cat("Start theme\n")
  p <- p + theme_classic()
  
  # Add title:
  if (!(is.null(main))){
    cat("Add title\n")
    p <- p + ggtitle(main)  
  }
  
  ## Font sizes:
  cat("Start line width and font size \n")
  p <- p + theme(axis.line = element_line(colour = 'black', size = lineWidth),
                 axis.text = element_text(size = fontSize),
                 axis.title = element_text(size = fontSize), 
                 plot.title = element_text(size = fontSize, hjust = 0.5), # center title 
                 legend.title = element_text(size = fontSize),
                 strip.text.x = element_text(size = fontSize), # facetwrap fontsize
                 legend.text = element_text(size = fontSize))
  print(p)
  if (savePNG | saveEPS){ dev.off(); print(p); Sys.sleep(1)}
  cat("Finished :-)\n")
  
  return(p)
  
} # end of function

# ================================================================================================================================================ #
#### Lineplot 1 IV: Aggregate per time point per condition per subject, plot (1 IV for any condition, time on x-axis): ####

custom_lineplot <- function(data, xVar="counter", yVar="response_cleaned", zVar="condition_f", subVar="sID", 
                            xLab = "Time (trial number)", yLab = "p(Go)", main = "",
                            selCol = c("#009933","#CC0000","#009933","#CC0000"), selLineType = c(1,1,2,2),
                            SEweight = 1, yLim = NULL, savePNG = F, saveEPS = F){
  #' Make line plot with group-level and individual lines.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. Variable needs to be numeric.
  #' @param yVar string, name of variable that goes on y-axis. Variable needs to be numeric.
  #' @param zVar string, name of variable that determines bar coloring. Variable needs to be a factor.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param main string, title of plot (optional).
  #' @param selCol vector of strings (HEX colors), colors for input levels of zVar (default: c("#009933","#CC0000","#009933","#CC0000")).
  #' @param selLineType vector of numerics, line types to use (default: c(1,1,2,2))
  #' @param SEweight scalar, weight to use for error shades (how many times SE; default: 1).
  #' @param yLim vector of two numbers, y-axis (default: NULL).
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @return creates (and saves) plot.
  
  # -------------------------------------------------------------------------- #
  ## Load packages:
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  
  # -------------------------------------------------------------------------- #
  ## Fixed plotting settings:
  LWD <- 3 # axes of plot
  CEX <- 1.5 # axes ticks and labels
  lineWidth <- 3 # line
  fontSize <- 20
  dodgeVal <- 0.6
  colAlpha <- 1
  
  # -------------------------------------------------------------------------- #
  ## Create variables under standardized names:
  data$x <- data[,xVar]
  data$y <- data[,yVar]
  data$z <- data[,zVar]
  data$subject <- data[,subVar]
  
  # -------------------------------------------------------------------------- #
  ## Aggregate data:
  aggrData <- ddply(data, .(subject, x, z), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  # Wide format: each subject/x condition/z condition in one line, variables subject, x, y, z
  
  ## Determine y limits if not given:
  if(is.null(yLim)){
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  # -------------------------------------------------------------------------- #
  ## Aggregate across subjects with Rmisc:
  summary_d <- summarySEwithin(aggrData, measurevar="y", idvar = "sID", na.rm = T,
                               withinvars = c("x","z"))
  # Aggregated over subjects, one row per condition, variables x, z, N, y, sd, se, ci
  
  # -------------------------------------------------------------------------- #
  ## Data dimensions:
  xVec <- unique(sort(as.numeric(summary_d$x)))
  xMax <- max(xVec)
  condNames <- unique(summary_d$z)
  nCond <- length(unique(summary_d$z))
  
  # -------------------------------------------------------------------------- #
  ## Plot name:
  
  plotName <- paste0("lineplot_",yVar,"_",xVar,"_",zVar)
  
  # -------------------------------------------------------------------------- #
  # Saving:
  
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)}
  
  # -------------------------------------------------------------------------- #
  ### Start plot:
  
  par(mar = c(5.1, 5.1, 4.1, 2.1)) # bottom, left, top, right
  
  # dev.off()
  for (iCond in 1:nCond){ # iCond <- 1
    condName <- condNames[iCond] # name of condition
    yVec <- summary_d$y[summary_d$z == condName] # y-variable
    seVec <- summary_d$se[summary_d$z == condName] # se variable
    plot(xVec, yVec, type = "l", 
         col = selCol[iCond], lty = selLineType[iCond], axes = F,
         lwd = LWD, cex.lab=CEX, cex.axis=CEX, cex.main=CEX,
         xlab = xLab, ylab = yLab, main = main,
         xlim = c(0, xMax), ylim = yLim)
    axis(side = 1, lwd = LWD, cex.axis = CEX, at = seq(0,xMax,5), line = 0)
    axis(side = 2, lwd = LWD, cex.axis = CEX, at = c(0, 0.5, 1))
    polygon(c(xVec,rev(xVec)),
            c(yVec-SEweight*seVec,rev(yVec+SEweight*seVec)),col = alpha(selCol[iCond],0.2), border = F)
    par(new = TRUE)
    
  }
  
  # Add legend:
  legend("top", legend=condNames,
         col=selCol, lty = selLineType, border = 0, lwd = LWD, cex = CEX, horiz=TRUE, bty = "n")
  
  if(savePNG | saveEPS){dev.off()}
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # bottom, left, top, right
  
}

# ============================================================================ #
#### Lineplot 1 IV with ggplot: Aggregate per time point per condition per subject, plot (1 IV for any condition, time on x-axis): ####

custom_lineplot_gg <- function(data, xVar="counter", yVar="response_cleaned", zVar="condition_f", subVar="sID", 
                               xLab = "Time (trial number)", yLab = "p(Go)", main = NULL,
                               selCol = c("#009933","#CC0000","#009933","#CC0000"), selLineType = c(1,1,2,2),
                               SEweight = 1, yLim = NULL, savePNG = F, saveEPS = F){
  #' Make line plot with group-level lines plus shades in ggplot.
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. Variable needs to be numeric.
  #' @param yVar string, name of variable that goes on y-axis. Variable needs to be numeric.
  #' @param zVar string, name of variable that determines bar coloring. Variable needs to be a factor.
  #' @param subVar string, name of variable containing subject identifier (default: subject).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param main string, title of plot (optional).
  #' @param selCol vector of strings (HEX colors), colors for input levels of zVar (default: c("#009933","#CC0000","#009933","#CC0000")).
  #' @param selLineType vector of numerics, line types to use (default: c(1,1,2,2))
  #' @param SEweight scalar, weight to use for error shades (how many times SE; default: 1).
  #' @param yLim vector of two numbers, y-axis (default: NULL).
  #' @param savePNG Boolean, save as .png file.
  #' @param saveEPS Boolean, save as .eps file.
  #' @return creates (and saves) plot.
  
  ## Load packages:
  require(plyr) # for ddply
  require(Rmisc) # for summarySEwithin
  
  ## Fixed plotting settings:
  lineWidth <- 1.3
  fontSize <- 30
  colAlpha <- 1
  
  ## Create variables under standardized names:
  data$x <- data[,xVar]
  data$y <- data[,yVar]
  data$z <- data[,zVar]
  data$subject <- data[,subVar]
  
  ## Aggregate data:
  aggrData <- ddply(data, .(subject, x, z), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  # Wide format: each subject/x condition/z condition in one line, variables subject, x, y, z
  
  ## Determine y limits if not given:
  if(is.null(yLim)){
    yLim <- determine_ylim_data_y(aggrData)
  }
  
  ## Aggregate across subjects with Rmisc:
  summary_d <- summarySEwithin(aggrData, measurevar="y", idvar = "sID", na.rm = T,
                               withinvars = c("x","z"))
  summary_d$x <- as.numeric(summary_d$x) # back to numeric to get continuous x-axis
  # Aggregated over subjects, one row per condition, variables x, z, N, y, sd, se, ci
  
  # Data dimensions:
  xVec <- unique(sort(as.numeric(summary_d$x)))
  xMax <- max(xVec)
  condNames <- unique(summary_d$z)
  nCond <- length(unique(summary_d$z))
  
  ## 3) ggplot:
  # Name:
  plotName <- paste0("lineplot_gg_",yVar,"_",xVar,"_",zVar)
  
  # Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(dirs$plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)}
  
  # Start plot:
  # par(mar = c(5.1, 5.1, 4.1, 2.1)) # bottom, left, top, right
  p <- ggplot(data = summary_d, aes(x = x, y = y))
  
  for (iCond in 1:nCond){ # iCond <- 1
    condData <- subset(summary_d, z == condNames[iCond]) # select data for this condition
    condData$ymin <- condData$y - SEweight * condData$se # lower edge of shade
    condData$ymax <- condData$y + SEweight * condData$se # upper edge of shade
    ## Shade:
    p <- p + geom_ribbon(data = condData, aes(x = x, y = y, ymin = ymin, ymax = ymax, group = 1),
                         fill = alpha(selCol[iCond],0.2))
    ## Line:
    p <- p + geom_path(data = condData, aes(x = x, y = y, group = 1), 
                       col = selCol[iCond], linetype = selLineType[iCond], size = lineWidth) 
  }
  
  # Add title:
  if (!(is.null(main))){
    p <- p + ggtitle(main)  
  }
  
  ## X-axis:
  p <- p + scale_x_continuous(limits = c(0, xMax), breaks = seq(0,xMax,5))
  
  ## Y-axis:
  if (yLim[1] == 0 & yLim[2] == 1){
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  if(!(is.null(yLim))){p <- p + coord_cartesian(ylim=yLim)}
  
  # Add theme, fontsizes:
  require(ggthemes)
  p <- p + labs(x=xLab, y = yLab, fill = condNames) +
    scale_fill_manual(values=selCol) + 
    theme_classic() + 
    theme(axis.text=element_text(size=fontSize),
          axis.title=element_text(size=fontSize), 
          plot.title = element_text(size=fontSize, hjust = 0.5), 
          legend.text = element_text(size=fontSize),
          legend.title=element_blank(), legend.position = "none")
  
  print(p)
  if(savePNG | saveEPS){dev.off()}
  # par(mar = c(5.1, 4.1, 4.1, 2.1)) # bottom, left, top, right
  
}

# =============================================================================================== #
#### REGRESSION LINES 1 IV: Plot regression line per group and per subject based on model output: #####

custom_regressionline1 <- function(mod, selEff, xLim = NULL, yLim = NULL, useEffect = TRUE, xVec = NULL,
                                   selCol = "red", margin = NULL, xLab = "x", yLab = "y", main = NULL, fontSize = NULL){
  #' Plot group-level regression line and subject-level regression lines based on 1 continuous predictor.
  #' @param mod model fitted with lme4
  #' @param selEff string, name of predictor to plot
  #' @param xLim vector of two numbers for y-axis limits
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data)
  #' @param subVar string, name of grouping variable (default: subject)
  #' @param useEffect boolean, extract upper and lower bound of confidence interval from effects(mod) (TRUE) or compute yourself based on vcov(mod) (FALSE)
  #' @param xVec vector of two numbers for x-axis ticks, to be used only if useEffect = FALSE.
  #' @param selCol strings (HEX colors), colors for line and error shade (default: red)
  #' @param margin vector of 4 numbers, margin of plot (default: NULL)
  #' @param xLab string, label for x-axis (default: "x")
  #' @param yLab string, label for y-axis (default: "y")
  #' @param main string, title of plot (default: "Plot")
  #' @param fontSize integer, font size for axes ticks and labels and title.
  #' @return makes regression line plot
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.Rmd
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.pdf
  
  # --------------------------------------------------------------------- #
  ## Load packages:
  require(ggplot2)
  require(lme4)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## General settings:
  
  colAlpha <- .95 # transparency
  lineWidth <- 1.5 # line width
  
  if (is.null(fontSize)){
    ## Font sizes for ordinary viewing: 15
    # fontSize <- 15
    ## Font sizes for saving: 30
    fontSize <- 30
    cat(paste0("No font size provided, use font size ",fontSize),"\n")
  }
  
  # --------------------------------------------------------------------- #
  ## Extract relevant data from model:
  
  ## Extract group-level and subject-level data sets:
  groupCoefs <- fixef(mod) # fixed effects
  if(length(coef(mod)) > 1){warning("Model features > 1 grouping variable, extract random effects for 1st one")}
  subCoefs <- coef(mod)[[1]] # fixed + random effects (only random effects for first structure)
  
  ## Extract effect names, check:
  effNames <- names(subCoefs)
  if (sum(grepl(selEff, effNames, fixed = T)) == 0){stop("selEff not a predictor in mod")}
  selEff1 <- effNames[grep(selEff, effNames, fixed = T)]
  selEff1 <- selEff1[1] # first effect
  
  ## Check presence of interactions:
  if(sum(grepl(":", effNames, fixed = T)) > 0){
    useEffect = FALSE
    warning("Model contains interactions; set useEffect to FALSE, check if predictions make sense")
  }
  
  ## Locate effect, count subjects:
  iCol <- which(names(subCoefs) == selEff1) # localize where in subCoefs effect of interest is
  allCols <- 1:length(subCoefs)
  otherCols <- allCols[allCols != iCol]
  nSub <- nrow(subCoefs)
  
  # --------------------------------------------------------------------- #
  ## Generate x-axis values for which to generate plots:
  
  tmp <- effect(selEff, mod) # retrieve objects from effects
  
  if (useEffect){ # if x samples automatically generated by effects() to be used
    
    xVecEval <- as.numeric(tmp$model.matrix[, iCol])
    
  } else { # if use by hand
    
    if (is.null(xVec)){ # if no xVec provided
      
      cat("No x-axis samples provided, extract from effects(mod)\n")
      xVecEval <- as.numeric(tmp$model.matrix[, iCol])
      
    } else {
      
      xVecEval <- seq(xVec[1], xVec[2], (xVec[2] - xVec[1]) / 100) # 100 samples between min and max
      
    }
    
  }
  
  xLen <- length(xVecEval) # number of x-axis samples
  
  if (is.null(xLim)){ # if no x-axis limits provided
    xLim <- c(min(xVecEval), max(xVecEval))
  }
  
  # --------------------------------------------------------------------- #
  ## Start ggplot: 
  
  d <- data.frame(x = xVecEval, y = xVecEval) #  just to initialize ggplot; also assign xVecEval to y
  p <- ggplot(data = d) # initialize
  
  # --------------------------------------------------------------------- #
  ## Loop over subjects, create single subject lines + ribbons:
  cat("Draw random-effects lines\n")
  
  inputData <- mod@frame # extract input data
  subIdx <- inputData[, length(inputData)] # extract subject indices to compute marginal means
  
  for (iSub in 1:nSub){ # iSub <- 1
    
    margMean <- colMeans(model.matrix(mod)[subIdx == iSub, ]) # mean of each regressor in design matrix
    iIntercept <- as.numeric(as.numeric(subCoefs[iSub, otherCols]) %*% as.numeric(margMean[otherCols])) # intercept is estimate of y at mean of all other regressors (marginal effect)
    # iIntercept <- subCoefs[iSub, 1] # extract intercept --> works only if other variables centered or very small
    iSlope <- subCoefs[iSub, iCol] # extract slope
    yVecEval <- iIntercept + xVecEval * iSlope
    if (isGLMM(mod)){yVecEval <-mod@resp$family$linkinv(yVecEval)} # bring to response scale
    
    ## Create single data frame (2 points should be enough to draw line, i.e. xmin and xmax):
    d <- data.frame(x = xVecEval, y = yVecEval)
    
    ## Thick line connecting means (plot line first and points on top):
    p <- p + geom_path(data = d, aes(x = x, y = y), color = 'grey40', # color = 'grey70'
                       alpha = 0.35, size = 1.0)
    # Used till 2021-11-21: grey50, size = 1
    
  }
  
  # --------------------------------------------------------------------- #
  ## Overall group-level line:
  cat("Draw fixed-effects line\n")
  
  ## Compute group-level predicted variables:
  
  # iIntercept <- as.numeric(groupCoefs[1]) # extract intercept
  margMean <- colMeans(model.matrix(mod)) # mean of each regressor in design matrix
  iIntercept <- as.numeric(groupCoefs[otherCols] %*% margMean[otherCols]) # intercept is estimate of y at mean of all other regressors (marginal effect)
  # iIntercept <- groupCoefs[1] # extract intercept --> works only if other variables centered or very small
  iSlope <- as.numeric(groupCoefs[iCol]) # extract slope
  yVecEval <- iIntercept + xVecEval * iSlope # recompute before transform
  if (isGLMM(mod)){yVecEval <- mod@resp$family$linkinv(yVecEval)} # bring to response scale
  
  ## Based on effects output:
  # yVecEval2 <- as.numeric(tmp$fit) # y-axis coordinates (untransformed)
  # if (isGLMM(mod)){yVecEval2 <- mod@resp$family$linkinv(yVecEval2)} # bring to response scale
  
  
  ## Create single data frame (2 points should be enough to draw line, i.e. xmin and xmax):
  d <- data.frame(x = xVecEval, y = yVecEval)
  
  ## Thick line connecting means (plot line first and points on top):
  p <- p + geom_path(data = d, aes(x = x, y= y), color = selCol, size = 1.5)
  
  # ----------------------------------------------------- #
  ## Error shades:
  if (useEffect){
    
    # ------------------------------- #
    ## Option A: Extract from effect() object:
    
    tmp <- effect(selEff, mod)
    
    # Mean/ line itself:
    yVecEval <- as.numeric(tmp$fit) # y-axis coordinates (untransformed)
    if (isGLMM(mod)){yVecEval <-mod@resp$family$linkinv(yVecEval)} # transform to response scale
    
    # Lower and upper limit of CI interval:
    ymin <- t(tmp$lower) # lower bound of CI interval
    ymax <- t(tmp$upper) # upper bound of CI interval
    if (isGLMM(mod)){ymin <-mod@resp$family$linkinv(ymin)} # transform to response scale
    if (isGLMM(mod)){ymax <-mod@resp$family$linkinv(ymax)} # transform to response scale
    
  } else {
    
    # ------------------------------- #
    ## Option B: Compute SE yourself:
    # https://github.com/cran/effects/blob/master/R/Effect.R
    
    mmat <- matrix(0, xLen, ncol(subCoefs)) # initialize design matrix: by default all regressors 0
    mmat[,1] <- rep(1, xLen) # add intercept
    mmat[, iCol] <- xVecEval # add regressor of interest
    V <- vcov(mod, complete=FALSE) # covariance matrix from model
    vcov <- mmat %*% V %*% t(mmat) # multiply with design matrix
    var <- diag(vcov) # variance
    se <- sqrt(var) # se # see tmp$se
    conf <- 1.96
    
    # iIntercept and iSlope computed above:
    yVecEval <- iIntercept + xVecEval * iSlope # recompute before transform
    ymin <- yVecEval - se*conf # lower bound of CI interval
    ymax <- yVecEval + se*conf # upper bound of CI interval
    if (isGLMM(mod)){ymin <-mod@resp$family$linkinv(ymin)} # bring to response scale
    if (isGLMM(mod)){ymax <-mod@resp$family$linkinv(ymax)} # bring to response scale
    # ymin;ymax
  }
  
  # ------------------------------- #
  ## Plot error bars/ shades:
  
  d <- data.frame(x = xVecEval, y = yVecEval, ymin = ymin, ymax = ymax)
  p <- p + geom_ribbon(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                       fill = selCol, alpha = 0.20, linetype = 0)
  
  # --------------------------------------------------------------------- #
  cat("Adjust axes, labels\n")
  
  # Y-axis:
  p <- p + coord_cartesian(xlim = xLim) 
  if (!is.null(yLim)){p <- p + coord_cartesian(ylim = yLim)} 
  
  # X-axis:
  xTickVec <- round(seq(xLim[1], xLim[2],(xLim[2] - xLim[1]) / 2), 2)
  p <- p + scale_x_continuous(breaks=xTickVec, labels=xTickVec)
  
  # Labels:
  p <- p + xlab(xLab) + ylab(yLab)
  
  # Other settings:
  if (!is.null(main)){
    cat("Add title\n")
    p <- p + ggtitle(main)
  }
  
  p <- p + theme_classic()
  if (!is.null(margin)){
    cat("Adjust margin\n")
    p <- p + theme(plot.margin = unit(margin, "cm"))
  }
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = fontSize),
                 axis.title = element_text(size = fontSize), 
                 plot.title = element_text(size = fontSize, hjust = 0.5), 
                 legend.text = element_text(size = fontSize))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
  return(p)
}  

# =============================================================================================== #
#### REGRESSION LINES 2 IV: Plot regression line per condition per group and per subject based on model output: #####

custom_regressionline2 <- function(mod, xVar, zVar, 
                                   xLim = NULL, yLim = NULL, xVec = NULL,
                                   selCol = c("red", "blue"), margin = NULL, 
                                   xLab = "x", yLab = "y", main = NULL, fontSize = NULL){
  #' Plot group-level regression line and subject-level regression lines based on 1 continuous predictor and 1 binary predictor.
  #' @param mod model fitted with lme4.
  #' @param xVar string, name of continuous predictor to plot on x-axis.
  #' @param zVar string, name of binary predictor to plot with different colors.
  #' @param subVar string, name of grouping variable (default: subject).
  #' @param xLim vector of two numbers for y-axis limits.
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data).
  #' @param xVec vector of two numbers for x-axis ticks.
  #' @param selCol vector of strings (HEX colors), colors for line and error shade (default: blue and red).
  #' @param margin vector of 4 numbers, margin of plot (default: NULL).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param main string, title of plot (default: "Plot").
  #' @param fontSize integer, font size for axes ticks and labels and title.
  #' @return makes regression line plot.
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.Rmd
  #'  See: https://github.com/jorvlan/open-visualizations/blob/master/R/repmes_tutorial_R.pdf
  
  # --------------------------------------------------------------------- #
  ## Load packages:
  
  require(ggplot2)
  require(lme4)
  
  # --------------------------------------------------------------------- #
  ## General settings:
  
  alphaSub <- 0.10
  alphaShade <- 0.20
  sizeGroup <- 1.5
  sizeSub <- 1
  lineWidth <- 1.5
  
  if (is.null(fontSize)){
    ## Font sizes for ordinary viewing: 15
    # fontSize <- 15
    ## Font sizes for saving: 30
    fontSize <- 30
    cat(paste0("No font size provided, use font size ",fontSize),"\n")
  }
  
  # --------------------------------------------------------------------- #
  ## Extract relevant data from model:
  
  ## Extract group-level and subject-level data sets:
  groupCoefs <- fixef(mod) # fixed effects
  
  if(length(coef(mod)) > 1){warning("Model features > 1 grouping variable, extract random effects for 1st one")}
  subCoefs <- coef(mod)[[1]] # fixed + random effects (only random effects for first structure)
  
  ## Extract effect names, check:
  effNames <- names(subCoefs)
  
  if (sum(grepl(xVar, effNames, fixed = T)) == 0){stop("xVar not a predictor in mod")}
  xVar1 <- effNames[grep(xVar, effNames, fixed = T)] # check version in effect names
  xVar1 <- xVar1[1] # only first match
  
  if (sum(grepl(zVar, effNames, fixed = T)) == 0){stop("zVar not a predictor in mod")}
  zVar1 <- effNames[grep(zVar, effNames, fixed = T)] # check version in effect names
  zVar1 <- zVar1[1] # only first match
  
  intVar <- paste0(xVar, ":", zVar) # combine as name should look like
  intVar1 <- effNames[grep(":", effNames, fixed = T)] # check version in effect names
  intVar1 <- paste0(intVar,"1") # only first match
  
  ## Indices of effects:
  xCol <- which(names(subCoefs) == xVar) # localize where in subCoefs effect of interest is
  zCol <- which(names(subCoefs) == zVar1) # localize where in subCoefs effect of interest is
  intCol <- which(names(subCoefs) == intVar1) # localize where in subCoefs effect of interest is
  
  ## Number of subjects:
  nSub <- nrow(subCoefs)
  
  ## X-axis for which to generate plots:
  tmp <- effect(intVar, mod) # retrieve objects from effects
  if (is.null(xVec)){
    xVecEval <- sort(unique(as.numeric(tmp$model.matrix[, xCol])))
  } else {
    xVecEval <- xVec # copy over input
  }
  
  xLen <- length(xVecEval) # number of x-axis samples
  
  if (is.null(xLim)){ # if no x-axis limits provided
    xLim <- c(xVecEval[1], xVecEval[xLen])
  }
  
  ## Z-axis with how levels of factor are coded:
  zVecEval <- unique(as.numeric(tmp$model.matrix[, zCol]))
  if(length(zVecEval) > 2){stop("zVec is variable with > 2 levels, must have only 2 levels")}
  
  # --------------------------------------------------------------------- #
  ## Initialize empty ggplot: 
  
  d <- data.frame(x = xVecEval, y = xVecEval) #  just to initialize ggplot; assign xVecEval also to y
  p <- ggplot(data = d) # initialize
  
  # --------------------------------------------------------------------- #
  ## Loop over subjects, create single subject lines:
  
  cat("Draw random-effects lines\n")
  
  for (iSub in 1:nSub){ # iSub <- 1
    
    ## Intercepts:
    iInter1 <- subCoefs[iSub, 1] + zVecEval[1] * subCoefs[iSub, zCol] # extract intercept condition 1
    iInter2 <- subCoefs[iSub, 1] + zVecEval[2] * subCoefs[iSub, zCol] # extract intercept condition 2
    
    ## Slopes:
    iSlope1 <- subCoefs[iSub, xCol] + zVecEval[1] * subCoefs[iSub, intCol] # extract slope condition 1
    iSlope2 <- subCoefs[iSub, xCol] + zVecEval[2] * subCoefs[iSub, intCol] # extract slope condition 2
    
    # Simulated y-data per subject:
    yVecEval1 <- iInter1 + xVecEval * iSlope1
    yVecEval2 <- iInter2 + xVecEval * iSlope2
    
    if (isGLMM(mod)){yVecEval1 <- mod@resp$family$linkinv(yVecEval1)} # bring to response scale
    if (isGLMM(mod)){yVecEval2 <- mod@resp$family$linkinv(yVecEval2)} # bring to response scale
    
    ## Create single data points (2 should be enough, i.e. xmin and xmax)
    d <- data.frame(x = xVecEval, y1 = yVecEval1, y2 = yVecEval2)
    
    ## Thick line connecting means (plot line first and points on top):
    p <- p + 
      geom_path(data = d, aes(x = x, y = y1), color = selCol[1],
                alpha = alphaSub, size = sizeSub) + 
      geom_path(data = d, aes(x = x, y = y2), color = selCol[2],
                alpha = alphaSub, size = sizeSub)
    
  }
  
  # --------------------------------------------------------------------- #
  ## Overall group-level line:
  
  cat("Draw fixed-effects line\n")
  
  ## Extract from effect() object:
  tmp <- effect(intVar, mod)
  idx1 <- 1:xLen # first half of predicted values
  idx2 <- (xLen+1):(2*xLen) # second half of predicted values
  
  # Line itself (means): 
  yVecEval1 <- as.numeric(tmp$fit)[idx1] # y-axis coordinates (untransformed)
  yVecEval2 <- as.numeric(tmp$fit)[idx2] # y-axis coordinates (untransformed)
  if (isGLMM(mod)){yVecEval1 <- mod@resp$family$linkinv(yVecEval1)} # transform to response scale
  if (isGLMM(mod)){yVecEval2 <- mod@resp$family$linkinv(yVecEval2)} # transform to response scale
  
  # Lower and upper limit of CI interval:
  ymin1 <- t(tmp$lower)[idx1] # lower bound of CI interval condition 1
  ymin2 <- t(tmp$lower)[idx2] # lower bound of CI interval condition 2
  ymax1 <- t(tmp$upper)[idx1] # upper bound of CI interval condition 1
  ymax2 <- t(tmp$upper)[idx2] # upper bound of CI interval condition 2
  if (isGLMM(mod)){ymin1 <- mod@resp$family$linkinv(ymin1)} # transform to response scale
  if (isGLMM(mod)){ymin2 <- mod@resp$family$linkinv(ymin2)} # transform to response scale
  if (isGLMM(mod)){ymax1 <- mod@resp$family$linkinv(ymax1)} # transform to response scale
  if (isGLMM(mod)){ymax2 <- mod@resp$family$linkinv(ymax2)} # transform to response scale
  
  # --------------------------------------------------------------------- #
  ## Thick line connecting group-level means:
  
  d <- data.frame(x = xVecEval, y1 = yVecEval1, y2 = yVecEval2)
  p <- p + 
    geom_path(data = d, aes(x = x, y = y1), color = selCol[1], size = sizeGroup) + 
    geom_path(data = d, aes(x = x, y = y2), color = selCol[2], size = sizeGroup)
  
  # --------------------------------------------------------------------- #
  ## Error shades:
  
  ## Plot error bars/ shades:
  d <- data.frame(x = xVecEval, y = yVecEval1, ymin = ymin1, ymax = ymax1)
  p <- p + geom_ribbon(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                       fill = selCol[1], alpha = alphaShade, linetype = 0)
  d <- data.frame(x = xVecEval, y = yVecEval2, ymin = ymin2, ymax = ymax2)
  p <- p + geom_ribbon(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                       fill = selCol[2], alpha = alphaShade, linetype = 0)
  
  # --------------------------------------------------------------------- #
  ## Further plot settings:
  
  cat("Adjust axes, labels\n")
  
  ## Y-axis:
  p <- p + coord_cartesian(xlim = xLim, ylim = yLim) 
  
  ## X-axis:
  # xTickVec <- round(seq(xLim[1], xLim[2], (xLim[2] - xLim[1]) / 2), 2) # extremes, middle
  xTickVec <- round(xVecEval, 2) # extremes, middle
  p <- p + scale_x_continuous(breaks=xTickVec, labels=xTickVec)
  
  ## Labels:
  p <- p + xlab(xLab) + ylab(yLab)
  
  ## Title:
  if (!is.null(main)){
    cat("Add title\n")
    p <- p + ggtitle(main)
  }
  
  ## Theme:  
  p <- p + theme_classic()
  
  ## Margin:
  if (!is.null(margin)){
    cat("Adjust margin\n")
    p <- p + theme(plot.margin = unit(margin, "cm"))
  }
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = fontSize),
                 axis.title = element_text(size = fontSize), 
                 plot.title = element_text(size = fontSize, hjust = 0.5), 
                 legend.text = element_text(size = fontSize))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
  return(p)
}  

# =============================================================================================== #
#### REGRESSION BARS 1 IV: Plot regression bars per group and per subject based on model output: #####

custom_regressionbar1 <- function(mod, selEff, selCol = "red", 
                                  xLab = "x", yLab = "y", main = NULL, xLabels = NULL,
                                  margin = NULL, fontSize = NULL, yLim = NULL){
  #' Plot group-level regression bar and subject-level regression lines based on 1 binary predictor.
  #' @param mod model fitted with lme4.
  #' @param selEff string, name of predictor to plot.
  #' @param subVar string, name of grouping variable (default: subject).
  #' @param selCol strings (HEX colors), colors for line and error shade (default: "red" for all).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param main string, title of plot (default: NULL).
  #' @param xLabels vector of strings, x-axis ticks (optional; otherwise retrieved from model).
  #' @param margin vector of 4 numbers, margin of plot (default: NULL).
  #' @param fontSize integer, font size for axes ticks and labels and title.
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data).
  #' @return makes regression line plot
  
  # --------------------------------------------------------------------- #
  ## Load packages:
  
  require(ggplot2)
  require(lme4)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  # --------------------------------------------------------------------- #
  ## General settings:
  
  colAlpha <- .95
  lineWidth <- 1.5
  
  if (is.null(fontSize)){
    ## Font sizes for ordinary viewing: 15
    # fontSize <- 15
    ## Font sizes for saving: 30
    fontSize <- 30
    cat(paste0("No font size provided, use font size ",fontSize),"\n")
  }
  
  # --------------------------------------------------------------------- #
  ## Extract relevant data from model:
  groupCoefs <- fixef(mod) # fixed effects
  if(length(coef(mod)) > 1){warning("Model features > 1 grouping variable, extract random effects for 1st one")}
  subCoefs <- coef(mod)[[1]] # fixed + random effects
  
  ## Extract effect names, check:
  effNames <- names(subCoefs)
  if (sum(grep(selEff, effNames, fixed = T)) == 0){stop("selEff not a predictor in mod")}
  selEff1 <- effNames[grep(selEff, effNames, fixed = T)] # check version in effect names
  selEff1 <- selEff1[1] # only first match
  
  ## Check presence of interactions:
  if(sum(grepl(":", effNames, fixed = T)) > 0){message("Model contains interactions; predictions might be off")}
  
  ## Locate effect:
  iCol <- which(names(subCoefs) == selEff1) # localize where in subCoefs effect of interest is
  
  ## Count subjects:
  nSub <- nrow(subCoefs)
  
  ## x-coordinates for which to plot:
  tmp <- effect(selEff, mod) # retrieve objects from effects
  xVecEval <- as.numeric(tmp$model.matrix[, iCol])
  xLen <- length(xVecEval)
  
  ## In case of sum-to-zero coding: Flip everything
  sum2zero <- F
  if(all(xVecEval == c(1, -1))){sum2zero <- T; message("Sum-to-zero coding detected, flip y values, colors, x-axis labels")}
  
  ## Complete color vector:
  selCol <- rep(selCol, length.out = xLen)
  if (sum2zero){selCol <- rev(selCol)} # flip of sum-to-zero coding
  
  ## X-axis tick labels:
  if (is.null(xLabels)){
    xLabels <- tmp$variables[[1]]$levels
  }
  if (sum2zero){xLabels <- rev(xLabels)} # flip of sum-to-zero coding
  
  # --------------------------------------------------------------------- #
  ## Start ggplot: 
  
  d <- data.frame(x = xVecEval, y = xVecEval) #  just to initialize ggplot
  p <- ggplot(data = d) # initialize
  
  # --------------------------------------------------------------------- #
  ## Loop over subjects, create single subject lines + ribbons:
  cat("Draw random-effects lines\n")
  
  for (iSub in 1:nSub){ # iSub <- 1
    
    iInter <- subCoefs[iSub, 1] # extract intercept
    iSlope <- subCoefs[iSub, iCol] # extract slope
    yVecEval <- iInter + xVecEval * iSlope # compute y-axis value
    if (isGLMM(mod)){yVecEval <- mod@resp$family$linkinv(yVecEval)} # bring to response scale
    if (sum2zero){yVecEval <- rev(yVecEval)} # flip of sum-to-zero coding
    
    ## Create single data frame (2 data points should be enough, i.e. xMin and xMax):
    d <- data.frame(x = xVecEval, y = yVecEval)
    
    ## Thick line connecting means (plot line first and points on top):
    p <- p + geom_path(data = d, aes(x = x, y= y), color = 'grey40', # color = 'grey70'
                       alpha = 0.35, size = 1)
    
  }
  
  # --------------------------------------------------------------------- #
  ## Overall group-level mean and error-bar:
  cat("Draw fixed-effects line\n")
  
  # groupSE <- VarCorr(mod)
  # iCol <- which(colnames(groupCoefs)==selEff) # localize where in groupCoefs effect of interest is
  iInter <- as.numeric(groupCoefs[1]) # extract intercept
  iSlope <- as.numeric(groupCoefs[iCol]) # extract slope
  yVecEval <- iInter + xVecEval * iSlope
  if (isGLMM(mod)){yVecEval <-mod@resp$family$linkinv(yVecEval)} # bring to response scale
  if (sum2zero){yVecEval <- rev(yVecEval)} # flip of sum-to-zero coding
  
  ## Create single data frame (2 data points should be enough, i.e. xMin and xMax):
  d <- data.frame(x = xVecEval, y = yVecEval)
  
  # ------------------------------- #
  ## Thick line connecting means (plot line first and points on top):
  p <- p + geom_path(data = d, aes(x = x, y = y), color = "black", size = 1.5)
  
  # ------------------------------- #
  ## Point for mean:
  p <- p + geom_point(data = d, aes(x = x, y = y), # point
                      color = selCol, alpha = 1, size = 5) # size = 2
  
  # ------------------------------- #
  ## Error shades:
  
  ymin <- tmp$lower 
  ymax <- tmp$upper 
  if (isGLMM(mod)){ymin <- mod@resp$family$linkinv(ymin)} # bring to response scale
  if (isGLMM(mod)){ymax <- mod@resp$family$linkinv(ymax)} # bring to response scale
  if (all(xVecEval == c(1, -1))){yVecEval <- rev(yVecEval)} # flip of sum-to-zero coding
  if (all(xVecEval == c(1, -1))){ymin <- rev(ymin)} # flip of sum-to-zero coding
  if (all(xVecEval == c(1, -1))){ymax <- rev(ymax)} # flip of sum-to-zero coding
  d <- data.frame(x = xVecEval, y = yVecEval, ymin = ymin, ymax = ymax)
  p <- p + geom_errorbar(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                         color = selCol, width = 0.15, size = 1.5, alpha = .6)
  
  # --------------------------------------------------------------------- #
  ## Further plot settings:
  
  cat("Adjust axes, labels\n")
  
  ## Y-axis:
  if (!is.null(yLim)){
    p <- p + coord_cartesian(ylim = yLim) 
    if (yLim[1] == 0 & yLim[2] == 1){
      p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
    }
  }
  
  ## X-axis:
  p <- p + coord_cartesian(xlim = c(min(xVecEval) - 0.5, max(xVecEval) + 0.5)) 
  p <- p + scale_x_continuous(breaks = xVecEval, labels = xLabels)
  
  ## Labels:
  p <- p +  xlab(xLab) + ylab(yLab)
  
  ## Title:
  if (!is.null(main)){
    p <- p + ggtitle(main) # title off for printing for poster
  }
  
  ## Theme:
  p <- p + theme_classic() # theme
  
  ## Margin:
  if (!is.null(margin)){
    p <- p + theme(plot.margin = unit(margin, "cm"))
  }
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = fontSize),
                 axis.title = element_text(size = fontSize), 
                 plot.title = element_text(size = fontSize, hjust = 0.5), # center title 
                 legend.text = element_text(size = fontSize))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
  return(p)
}  

# =============================================================================================== #
#### REGRESSION BARS 1 IV flat: Plot regression bars per group based on model output: #####

custom_regressionbar_flat1 <- function(mod, selEff, selCol = "red", 
                                       xLab = "x", yLab = "y", main = NULL, xLabels = NULL,
                                       margin = NULL, fontSize = NULL, yLim = NULL){
  #' Plot values predicted by lm/glm model for given categorical conditions based on 1 binary predictor.
  #' @param mod model fitted with base (lm or glm).
  #' @param selEff string, name of predictor to plot.
  #' @param selCol strings (HEX colors), colors for line and error shade (default: "red" for all).
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param main string, title of plot (default: NULL).
  #' @param xLabels vector of strings, x-axis ticks (default: numbered through).
  #' @param margin vector of 4 numbers, margin of plot (default: NULL).
  #' @param fontSize integer, font size for axes ticks and labels and title.
  #' @param yLim vector of two numbers for y-axis limits (default: determine based on min and max of input data).
  #' @return makes regression line plot.
  
  require(ggplot2)
  require(effects)
  
  ## General settings:
  colAlpha <- .95
  lineWidth <- 1.5
  
  if (is.null(fontSize)){
    ## Font sizes for ordinary viewing: 15
    # fontSize <- 15
    ## Font sizes for saving: 30
    fontSize <- 30
    cat(paste0("No font size provided, use font size ",fontSize),"\n")
  }
  
  # --------------------------------------------------------------------- #
  ## Rely on output of effect(mod) object:
  cat("Extract y-axis values from effect object\n")
  tmp <- effect(selEff, mod) # extract CI end points from effect object
  
  ## Retrieve fitted y-values:
  yVecEval <- tmp$fit
  if(class(mod)[1] == "glm"){yVecEval <-mod$family$linkinv(yVecEval)}
  nVal <- length(yVecEval) # number of fitted values
  
  ## Count number of y-axis values, create corresponding number of x-axis values:
  xVecEval <- seq(1, nVal, 1)
  selCol <- rep(selCol, length.out = nVal) # repeated selCol until enough colors for conditions
  
  ## Compute end points of whiskers:
  ymin <- tmp$lower # lower CI end
  ymax <- tmp$upper  # upper CI end
  if(class(mod)[1] == "glm"){ymin <-mod$family$linkinv(ymin)} # bring to response scale
  if(class(mod)[1] == "glm"){ymax <-mod$family$linkinv(ymax)} # bring to response scale
  
  ## Create single data points (2 should be enough, i.e. xmin and xmax)
  d <- data.frame(x = xVecEval, y = yVecEval, ymin = ymin, ymax = ymax) # concatenate to data
  
  # --------------------------------------------------------------------- #
  ## Start ggplot: 
  p <- ggplot(data = d) # initialize
  
  # ------------------------------- #
  ## Thick line connecting means (plot line first and points on top):
  cat("Add line connecting condition means\n")
  p <- p + geom_path(data = d, aes(x = x, y = y), color = "black", size = 1.5)
  
  # ------------------------------- #
  ## Point for mean:
  cat("Add points for condition means\n")
  p <- p + geom_point(data = d, aes(x = x, y = y), # point
                      color = selCol, alpha = 1, size = 5) # size = 2
  
  # ------------------------------- #
  ## Error shades:
  cat("Add error shades\n")
  p <- p + geom_errorbar(data = d, aes(x = x, y = y, ymin = ymin, ymax = ymax),
                         color = selCol, width = 0.15, size = 1.5, alpha = .6)
  
  # --------------------------------------------------------------------- #
  cat("Adjust axes, labels\n")
  
  ## Y-axis:
  if (yLim[1] == 0 & yLim[2] == 1){
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  p <- p + coord_cartesian(xlim = c(xVecEval[1] - 0.1, xVecEval[length(xVecEval)] + 0.1), ylim = yLim) 
  
  ## X-axis:
  if (is.null(xLabels)){xLabels <- xVecEval} # if no names provided: use position indices
  p <- p + scale_x_continuous(breaks = xVecEval, labels = xLabels)
  
  ## Labels:
  p <- p +  xlab(xLab) + ylab(yLab)
  
  ## Title:
  if (!is.null(main)){
    p <- p + ggtitle(main) # title off for printing for poster
  }
  
  ## Theme:  
  p <- p + theme_classic() # theme
  
  ## Margin:
  if (!is.null(margin)){
    p <- p + theme(plot.margin = unit(margin, "cm"))
  }
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(size = fontSize),
                 axis.title = element_text(size = fontSize), 
                 plot.title = element_text(size = fontSize, hjust = 0.5), # center title 
                 legend.text = element_text(size = fontSize))
  
  # Print plot in the end:
  print(p)
  cat("Finished :-)\n")
  return(p)
}  

# =============================================================================================== #
#### Plot group-level and subject-level coefficients as dots in horizontal dot-plot: #####

custom_coefplot <- function(mod, plotSub = TRUE, plotText = FALSE, dropIntercept = FALSE, revOrder = FALSE,
                            xLab = "Regression weight", yLab = "Predictor", main = NULL,
                            selCol = "blue", yLabels = NULL, xLim = NULL){
  #' Plot group-level  and subject-level coefficients as dots in horizontal dot-plot.
  #' @param mod model fitted with lme4.
  #' @param plotSub Boolean, whether to plot per-subject effects (TRUE) or not (FALSE; default: TRUE).
  #' @param plotText Boolean, whether to print value of group-level coefficient next to dot (TRUE) or not (FALSE; default: FALSE).
  #' @param dropIntercept Boolean, do not plot intercept (TRUE; default: false).
  #' @param revOrder Boolean, revert order of predictors (first one on top, TRUE) or not (last one on top FALSE) (default: false).
  #' @param xLab string, label for x-axis (default: "Regression weight").
  #' @param yLab string, label for y-axis (default: "Predictor").
  #' @param main string, title of plot (default: NULL).
  #' @param selCol strings (HEX colors), colors for group-level dots and error lines (default: "blue" for all).
  #' @param yLabels vector of strings, y-axis ticks (default: terms extracted from mod).
  #' @param xLim vector of two numerics, x-axis limits (optional).
  #' @return coefplot created with ggplot.
  
  # -------------------------------------------------------------------------- #
  ## Required packages:
  
  require(ggplot2)
  require(lme4)
  require(arm) # for se.fixef
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Fixed settings:
  
  SEweight <- 1.96 # width of whiskers
  lineWidth <- 1.5 # linewidth of axes
  fontSize <- 20 # font size for all text: 30 or 15 
  colAlpha <- 0.6 # transparency of per-subject dots
  subDisplacement <- 0 # systematic vertical displacement of per-subject dots
  nJitter <- 0.07 # amount of jitter added to per-subject dots
  nRound <- 3 # how much to round plotted text.
  
  # -------------------------------------------------------------------------- #
  ## Extract group-level information from input:
  
  meanVec <- as.numeric(fixef(mod))
  seVec <- se.fixef(mod)
  
  if (is.null(yLabels)){ # if not provided
    labelsVec <- colnames(mod@pp$X)
    # labels <- c("Intercept", attr(terms(mod), "term.labels"))
    # labels <- names(coef(mod)[[1]])
    labelsVec <- substitute_label(labelsVec) # translate labels
  } else { # if provided
    labelsVec <- yLabels
  }
  
  ## Concatenate group-level values to data frame:
  coefData <- data.frame(labelsVec, meanVec, seVec)
  names(coefData) <- c("label", "mean", "se")
  
  ## Text labels with significance stars based on z-values:
  # (1 - pnorm(1.64))*2
  # (1 - pt(1.96, df = 1000))*2
  coefData$z <- abs(coefData$mean / coefData$se) # z-value to evaluate significance
  coefData$zLabel <- as.character(round(coefData$mean, nRound)) # copy over, to string
  coefData$zLabel <- ifelse(coefData$z > 1.64 & coefData$z < 1.96, paste0(coefData$zLabel, "\U207A"), coefData$zLabel) # latin cross
  coefData$zLabel <- ifelse(coefData$z > 1.96, paste0(coefData$zLabel, "*"), coefData$zLabel)
  coefData$zLabel <- ifelse(coefData$z > 3, paste0(coefData$zLabel, "*"), coefData$zLabel)
  
  # Alternatives to cross for marginally significant effects:
  # https://unicode-table.com/en/sets/crosses/
  # for (iLabel in 1:nrow(coefData)){
  #   if (coefData$z[iLabel] > 1.64 & coefData$z[iLabel] < 1.96){
  #     coefData$zLabel[iLabel] <- expression(paste(eval(coefData$zLabel[iLabel]), "^+"))
  #     # coefData$zLabel[iLabel] <- substitute(expression(n "^+"), list(n = coefData$zLabel[iLabel]))
  #     # coefData$zLabel[iLabel] <- bquote(.(coefData$zLabel[iLabel])^+)
  #   }
  # } 
  
  # substitute(expression(a + b), list(a = coefData$zLabel[iLabel]))
  # substitute(expression(a ^+), list(a = coefData$zLabel[iLabel]))
  
  # coefData$zLabel <- ifelse(coefData$z > 1.64 & coefData$z < 1.96, expression(paste0(coefData$zLabel, "^+")), coefData$zLabel) # latin cross
  # coefData$zLabel <- ifelse(coefData$z > 1.64 & coefData$z < 1.96, paste0(coefData$zLabel, "\U2670"), coefData$zLabel) # latin cross
  # coefData$zLabel <- ifelse(coefData$z > 1.64 & coefData$z < 1.96, paste0(coefData$zLabel, "\U207A"), coefData$zLabel) # superscript + (rather small)
  
  ## Drop intercept or not:
  if(dropIntercept){
    coefData <- coefData[2:nrow(coefData), ] 
    # labels <- labels[2:length(labels)]
  }
  
  ## Determine final number of effects:
  nEff <- nrow(coefData)
  
  ## Adjust number of colors:
  if (length(selCol) > nEff){selCol <-  selCol[1:nEff]} # if too many colors: only use first
  selCol <- rep(selCol, length.out = nEff) # if too few colors: repeat until enough
  # selCol <- rev(COL2("RdBu", nEff))
  
  ## Reverse order (first regressor will be plotted on top):
  if (revOrder){
    coefData <- coefData[nrow(coefData):1,]
    selCol <- rev(selCol)
  }
  
  ## Compute index and lower/upper confidence bound:
  coefData$idx <- seq(1, nrow(coefData), 1) # numerical index of each effect to loop through (in correct order)
  coefData$lower <- coefData$mean - coefData$se * SEweight
  coefData$upper <- coefData$mean + coefData$se * SEweight
  
  ## Determine x axis limits:
  xMin <- min(coefData$lower)
  xMax <- max(coefData$upper)
  # --> overwrite later if individual subject points plotted
  
  # -------------------------------------------------------------------------- #
  ## Group-level plot:
  
  p <- ggplot(coefData, aes(x = mean, y = label)) # define ggplot and axed
  
  ## Add error bar lines:
  cat("Plot error bar whiskers\n")
  for (iEff in 1:nEff){ # iEff <- 1
    ## For this effect: extract index, upper and lower end of whisker
    effData <- data.frame(x = c(coefData$lower[iEff], coefData$upper[iEff]),
                          y = rep(coefData$idx[iEff], 2) + subDisplacement)
    p <- p + geom_line(data = effData, aes(x = x, y = y), lineWidth = 1.2, color = selCol[iEff])
  }
  
  # -------------------------------------------------------------------------- #
  ## Add fixed-effects points:
  
  cat("Plot fixed-effect coefficients\n")
  p <- p + geom_point(aes(x = mean, y = idx, color = factor(idx)), size = 5) +  # points for point estimates; size = 5
    scale_color_manual(values = selCol)
  
  # -------------------------------------------------------------------------- #
  ## Add subject-level points:
  
  if (plotSub){
    
    ## Extract subject-level information from input:
    subCoefs <- coef(mod)[[1]]
    nSub <- nrow(subCoefs)
    
    ## Drop intercept or not:
    if(dropIntercept){
      subCoefs <- subCoefs[, 2:ncol(subCoefs)] 
      # labels <- labels[2:length(labels)]
    }
    
    ## Reverse order (first regressor will be plotted on top):
    if (revOrder){
      subCoefs <- subCoefs[, ncol(subCoefs):1]
    }
    if (is.vector(subCoefs)){subCoefs <- data.frame(subCoefs)} # ensure it has a column dimension
    
    ## New axis limits:
    xMin <- min(subCoefs)
    xMax <- max(subCoefs)
    
    cat("Plot effect per subject\n")
    for (iEff in 1:nEff){ # iEff <- 1
      ## For this effect: extract and plot effects per subject
      effData <- data.frame(x = subCoefs[, iEff], # per-subject effect
                            y = rep(coefData$idx[iEff], nSub) + subDisplacement) # y-axis position
      effData$y <- jitter(effData$y, amount = nJitter) # add jitter to distinuish subjects
      p <- p + geom_point(data = effData, aes(x = x, y = y), size = 2, 
                          # shape = 16, color = "gray30", # all grey
                          # shape = 21, color = "black", fill = "gray70", stroke = 1.2, # black edge, white fill
                          shape = 1, color = "black", # or color = selCol[iEff],
                          alpha = colAlpha)
    }  
  }
  
  # -------------------------------------------------------------------------- #
  ## Determine x-axis dimensions for scaling:
  
  xRange <- xMax - xMin
  iMag <- log10(xRange)
  if (iMag < 0){
    iMag <- floor(iMag); 
    iMag <- 10 ^ (iMag - 1)
    xStep <- iMag * 5
  } else if (iMag > 0){
    iMag <- ceiling(iMag)
    iMag <- 10 ^ (iMag - 1)
    xStep <- iMag
  } # round order of magnitude
  
  # -------------------------------------------------------------------------- #
  ## Add group-level coefficient as text:
  
  if (plotText){
    textDisplacement <- (xMax - xMin) * 0.10 # 10% of x-axis width
    cat("Print values of group-level effects as text\n")
    p <- p + geom_text(data = coefData, aes(x = mean, y = idx, label = zLabel),  
                       nudge_x = textDisplacement, nudge_y = 0.3, na.rm = T, check_overlap = T, size = fontSize/3) # nudge_x = 0.20
  }
  
  # -------------------------------------------------------------------------- #
  ### Other settings:
  
  ## Horizontal line at x = 0:  
  p <- p + geom_vline(xintercept = 0, 
                      linetype = "dashed", colour = "#949494", lwd = 1) # line at zero
  
  ## Erode axis limits:
  cat(paste0("xMin = ", round(xMin, 3), ", xMax = ", round(xMax, 3), "\n"))
  if(xMin > 0){xMin <- xMin * 0.90} else (xMin <- xMin * 1.10)
  if(xMax > 0){xMax <- xMax * 1.10} else (xMax <- xMax * 0.90)
  
  ## Determine breaks:
  if(is.null(xLim)){
    
    ## X-axis ticks:
    # ------------------------------------------------------------------------ #
    ## Optional A: round x limits, put 5 points in between (might be quirky digits)
    
    xBreakMin <- floor(xMin/iMag)*iMag # copy over, round
    xBreakMax <- ceiling(xMax/iMag)*iMag # copy over, round
    
    ## X-axis limits:
    xLim <- c(xBreakMin, xBreakMax)
    
    ## Distance between x-axis ticks:
    expVec <- seq(-10, 10, 1) # exponents for candidate ticks 
    xTickVec <- sort(c(10^expVec, 5 * 10^expVec)) # candidate distances
    nTick <- 5 # number desired ticks on x-axis
    xDist <- ceiling((xBreakMax - xBreakMin)) / nTick
    xStepIdx <- which(abs(xTickVec - xDist) == min(abs(xTickVec - xDist)))
    xStepIdx <- xStepIdx[1] # first in case of multiple minima
    xStep <- xTickVec[xStepIdx]
    
    ## New version:
    # xBreakRange <- xBreakMax - xBreakMin # range
    # xStep <- xBreakRange / 5 # 5 labels
    # breakVec <- seq(xBreakMin, xBreakMax, xStep)
    
    ## Old version:
    # iMag <- log10(xBreakMax) # order of magnitude
    # if (iMag < 0){iMag <- floor(iMag)} else if (iMag > 0){iMag <- floor(iMag)} # round order of magnitude
    # if (iMag < 0){xBreakMin <- round(xBreakMin, -1*iMag); xBreakMax <- round(xBreakMax, -1*iMag); xStep <- 10^iMag * 5} # * 2?
    # if (iMag >= 0){xBreakMin <- ceiling(xBreakMin) - 1; xBreakMax <- ceiling(xBreakMax); xStep <- 10^iMag}
    # breakVec <- seq(xBreakMin, xBreakMax, xStep)
    
    # ------------------------------------------------------------------------ #
    ## Option B: regular break points, selected between xMin and xMax (regular digits)
    
    breakVec <- seq(floor(xMin), ceiling(xMax), xStep) # just very broad break points, aligned to magnitude
    # breakVec <- breakVec[breakVec > xMin & breakVec < xMax]
    # if (iMag < 0){breakVec <- round(breakVec, -1*iMag)}
    
    # xLim <- c(xBreakMin - iMag, xBreakMax + iMag)
    
  } else { # of xLim provided:
    
    breakVec <- sort(c(xLim[1], mean(xLim), xLim[2]))
    
  }
  # cat(paste0("xLim = ", xLim[1], ", ", xLim[2], "\n"))
  # cat(paste0("breakVec = ", paste0(breakVec, collapse = ", "), "\n"))
  
  ## Axes:
  p <- p + coord_cartesian(xlim = xLim) # set x limits
  p <- p + scale_x_continuous(breaks = breakVec)
  p <- p + scale_y_continuous(breaks = 1:nEff, labels = coefData$label)
  
  ## Labels:
  p <- p + labs(x = xLab,
                y = yLab)
  
  ## Add title:
  if (!(is.null(main))){
    p <- p + ggtitle(main)  
  }
  
  ## Theme:
  p <- p + theme_classic() # base_size = 14
  
  ## Font sizes:
  p <- p + theme(axis.text = element_text(colour = "black", size = fontSize),
                 axis.line = element_line(colour = "black", size = lineWidth), # fixed font sizes
                 axis.title = element_text(colour = "black", size = fontSize), 
                 plot.title = element_text(colour = "black", size = fontSize, hjust = 0.5), # center title 
                 legend.position = "none")
  
  print(p)
  cat("Finished :-)\n")
  return(p)
  
}

# ============================================================================ #
#### Plot intercorrelations of regressors in design matrix: #####

corrplot_regressors <- function(mod, perSub = F, varNames = NULL, savePNG = TRUE){
  #' Plot intercorrelations between regressors in design matrix using coefplot.
  #' @param mod model object fitted with lme4 or another package.
  #' @param perSub compute correlation between regressors separately per subject, than average (TRUE) (default: FALSE).
  #' @param varNames vector of strings, names of variables to use for rows/ columns of design matrix.
  #' @param savePNG Boolean, save as PNG to dirs.plotDir (TRUE, default) or not (FALSE).
  #' @return nothing, plots to console and saves in dirs$plotDir.
  
  require(psych) # for fisherz and fisherz2r
  require(corrplot)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ## Extract design matrix:
  
  DM <- model.matrix(mod) # extract design matrix
  DM <- DM[, 2:ncol(DM)] # drop intercept
  
  # -------------------------------------------------------------------------- #
  ## Compute correlation:
  
  if (perSub){
    
    cat("Separate DM per subject, average\n")
    subIdx <- mod@frame[, ncol(mod@frame)] # extract subject indices
    stopifnot(nrow(DM) == length(subIdx)) # assure dimensions match
    subVec <- unique(subIdx)
    nSub <- length(subVec)
    
    Mlist <- vector(mode = "list", length = nSub) # initialize
    
    for (iSub in 1:nSub){ # iSub <- 1
      
      subID <- subVec[iSub] # subject ID
      DMsub <- DM[subIdx == subID, ] # select part of design matrix for this subject
      M <- cor(DMsub) # correlation
      MF <- fisherz(M) # Fisher-z transform
      diagIdx <- which(MF == Inf) # detect diagonal elements (infinite)
      MF[diagIdx] <- 1 # temporarily overwrite to 1, correct later after transforming back
      
      Mlist[[iSub]] <- MF # store
    }
    
    M <- Reduce("+", Mlist)/nSub # mean across subjects
    M <- fisherz2r(M) # transform back
    M[diagIdx] <- 1 # set diagonal back to 1 
    
  } else {
    M <- cor(DM)
  }
  
  ## Print range to console:
  Mvec <- as.vector(M) # to vector
  diagVec <- seq(1, length(Mvec), nrow(M) + 1) # identify indices of diagonal
  Mvec[diagVec] <- NA # set diagonal to NA
  nRound <- 2
  cat(paste0("All correlations between r = ", round(min(Mvec, na.rm = T), nRound), " and r = ", round(max(Mvec, na.rm = T), nRound), "\n"))
  
  # -------------------------------------------------------------------------- #
  ## Overwrite variables names:
  
  if (is.null(varNames)){
    rownames(M) <- substitute_label(rownames(M)) # substitute for known variable names
  } else {
    stopifnot(nrow(M) == length(varNames)) # check if same length
    rowNames(M) <- varNames # overwrite
  }
  colnames(M) <- rownames(M)
  
  # -------------------------------------------------------------------------- #
  ## Title and name for saving:
  
  # https://stackoverflow.com/questions/14671172/how-to-convert-r-formula-to-text
  if (class(mod) %in% c("glmerMod")){
    
    formulaStr <- mod@call$formula
    
  } else {
    
    formulaStr <- attr(mod@frame, "formula")
    
  }
  formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
  
  # https://stackoverflow.com/questions/40509217/how-to-have-r-corrplot-title-position-correct
  # titleStr <- paste0("Intercorrelation regressors for \n", deparse(formulaStr), width.cutoff = 20))
  
  plotName <- paste0("corr_reg_", formula2handle(formulaStr))
  if (perSub){plotName <- paste0(plotName, "_perSub")}
  plotName <- paste0(plotName, ".png")
  
  # -------------------------------------------------------------------------- #
  ## Make corrplot:
  
  # mar = c(0,0,1,0)  
  if(savePNG) {png(paste0(dirs$plotDir, "ic_reg/", plotName), width = 480, height = 480)}
  
  # https://stackoverflow.com/questions/40352503/change-text-color-in-corrplot-mixed
  
  # corrplot(M, method = "circle", col = rev(COL2('RdBu', 200))) # colored dots of different size
  # corrplot(M, method = "number", col = rev(COL2('RdBu', 200))) # numerals of different color
  
  ## Uper half colored dots of different size, lower half black numerals, variable names in diagonal:
  # corrplot.mixed(M, lower = "number", upper = "circle", lower.col = "black", upper.col = rev(COL2('RdBu', 200)), tl.col = "black")
  
  ## Colors dots of different size with black numerials in them:
  corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
  # corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.offset = 1, tl.srt = 0) # column labels higher, not rotated
  
  ## Also numerals in color, different color scale (uniform), variable names in diagonal:
  # corrplot.mixed(M, lower = "number", upper = "circle", lower.col = COL1('YlOrRd', 200), upper.col = rev(COL2('RdBu')))
  # corrplot.mixed(M, lower = "number", upper = "circle", lower.col = COL1('YlOrRd', 200), upper.col = COL1('YlOrRd', 200))
  
  # print(p)
  if(savePNG){
    dev.off(); 
    cat(paste0("Saved under \n", plotName, " :-)\n"))
    corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
  }
  
}

# ============================================================================ #
#### Plot intercorrelations of coefficients from model: #####

corrplot_coefficients <- function(input, varNames = NULL, savePNG = TRUE){ 
  #' Plot intercorrelations between regressors in design matrix using coefplot.
  #' @param mod model object fitted with lme4 or another package.
  #' @param varNames vector of strings, names of variables to use for rows/ columns of design matrix.
  #' @param savePNG Boolean, save as PNG to dirs.plotDir (TRUE, default) or not (FALSE).
  #' @return nothing, plots to console and saves in dirs$plotDir.
  
  require(corrplot)
  
  # -------------------------------------------------------------------------- #
  ## Close any open plots:
  
  # if (length(dev.list()!=0)){dev.off()}
  
  # -------------------------------------------------------------------------- #
  ### Detect class and extract coefficients:
  
  ## Detect model class:  
  modClass <- class(input)
  if (modClass == "list"){modClass <- class(input$modList[[1]]); modClass <- modClass[1]}
  cat(paste0("Input model of class ", modClass, "\n"))
  
  ## Extract coefficients, parameter names, formula:
  
  if(modClass %in% c("glmerMod", "lmerMod", "lmerTest", "lmerModLmerTest")){ # from lme4
    
    coefMat <- coef(input)[[1]]
    
    parNamesVec <- colnames(coefMat)
    
    if (modClass %in% c("glmerMod")){
      
      formulaStr <- input@call$formula
      
    } else {
      
      formulaStr <- attr(input@frame, "formula")
      
    }
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
    # lmer
  } else if (modClass %in% c("mixed")){ # from afex
    
    coefMat <- coef(input$full_model)[[1]]
    parNamesVec <- rownames(coefMat)
    
    formulaStr <- input@call$formula
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
  } else if (modClass %in% c("lm", "glm")){ # from lm
    
    coefMat <- input$bMat
    parNamesVec <- names(coef(input$modList[[1]]))
    if (modClass == "glm"){
      formulaStr <- input$modList[[1]]$formula
    } else if (modClass == "lm"){
      formulaStr <- eval(input$modList[[1]]$call[[2]]) 
    } else {
      stop("Unknown model class")
    }
    
  } else if (modClass %in% c("brmsfit")){ # from brms
    
    ## Parameter names and formula:
    parNamesVec <- row.names(fixef(input)) # names of all predictors
    nParam <- length(parNamesVec)
    
    formulaStr <- input$formula
    formulaStr <- formulaStr[[1]] # extract only first object
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
    ## Exttract correlation estimates:
    brmsVarCorr <- VarCorr(input)[[1]]$cor 
    brmsVarCorr <- VarCorr(input)$subject_f$cor 
    
    coefMat <- matrix(NA, nParam, nParam) # initialize
    rownames(coefMat) <- parNamesVec
    colnames(coefMat) <- parNamesVec
    for (iParam1 in 1:nParam){ # iParam1 <- 1
      for (iParam2 in 1:nParam){ # iParam2 <- 2
        coefMat[iParam1, iParam2] <- brmsVarCorr[iParam1, 1, iParam2]
      }
    }
    
  } else {
    
    stop("Unknown model class")
    
  }
  
  # -------------------------------------------------------------------------- #
  ## Compute correlation:
  
  M <- cor(coefMat)
  
  ## Print range to console:
  Mvec <- as.vector(M) # to vector
  diagVec <- seq(1, length(Mvec), nrow(M) + 1) # identify indices of diagonal
  Mvec[diagVec] <- NA # set diagonal to NA
  nRound <- 2
  cat(paste0("All correlations between r = ", round(min(Mvec, na.rm = T), nRound), " and r = ", round(max(Mvec, na.rm = T), nRound), "\n"))
  
  # -------------------------------------------------------------------------- #
  ## Overwrite variables names:
  
  if (is.null(varNames)){
    
    rownames(M) <- parNamesVec
    rownames(M) <- substitute_label(rownames(M)) # substitute for known variable names
    
  } else { # overwrite with inputs
    
    stopifnot(nrow(M) == length(varNames)) # check if same length
    rowNames(M) <- varNames # overwrite
  }
  
  colnames(M) <- rownames(M) # copy over
  
  # -------------------------------------------------------------------------- #
  ## Title and name for saving:
  
  # https://stackoverflow.com/questions/14671172/how-to-convert-r-formula-to-text
  
  plotName <- paste0("corr_coef_", modClass, "_", formula2handle(formulaStr), ".png")
  plotNameFull <- paste0(dirs$plotDir, "ic_coef/", plotName)
  cat(paste0("File path has ", nchar(plotNameFull), " characters\n"))
  if (nchar(plotNameFull) > 260){
    warning("File path too long, shorten\n")
    plotNameFull <- paste0(substr(plotNameFull, 1, 255)[1], ".png")
  }
  
  # -------------------------------------------------------------------------- #
  ## Visualize correlation matrix:
  # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
  
  if(savePNG) {png(plotNameFull, width = 480, height = 480)}
  
  # corrplot(M, method = "circle", col = rev(COL2('RdBu', 200)))
  # corrplot(M, method = "number", col = rev(COL2('RdBu', 200)))
  corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
  
  if(savePNG){
    dev.off(); 
    cat(paste0("Saved under \n", plotName, " :-)\n"))
    ## Plot again:
    corrplot::corrplot(M, addCoef.col = 'black', col = rev(COL2('RdBu')), tl.col = "black", tl.pos = "lt")
  }
  
}

# ============================================================================ #
#### Save coefficients from model: #####

save_coefficients <- function(input){ 
  #' Save coefficients from model output as .csv-file.
  #' @param mod model object fitted with lme4 or another package.
  #' @return nothing, plots to console and saves in dirs$plotDir.
  
  # -------------------------------------------------------------------------- #
  ### Detect class and extract coefficients:
  
  ## Detect model class:  
  modClass <- class(input)
  if (modClass == "list"){modClass <- class(input$modList[[1]]); modClass <- modClass[1]}
  cat(paste0("Input model of class ", modClass, "\n"))
  
  ## Extract coefficients and formula:
  
  if(modClass %in% c("glmerMod", "lmerMod", "lmerTest", "lmerModLmerTest")){ # from lme4
    
    coefMat <- coef(input)[[1]]
    
    if (modClass %in% c("glmerMod")){
      
      formulaStr <- input@call$formula
      
    } else {
      
      formulaStr <- attr(input@frame, "formula")
      
    }
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
    # lmer
  } else if (modClass %in% c("mixed")){ # from afex
    
    coefMat <- coef(input$full_model)[[1]]
    
    formulaStr <- input@call$formula
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
  } else if (modClass %in% c("lm", "glm")){ # from lm
    
    coefMat <- input$bMat
    colnames(coefMat) <- names(coef(input$modList[[1]])) # add variable names
    if (modClass == "glm"){
      formulaStr <- input$modList[[1]]$formula
    } else if (modClass == "lm"){
      formulaStr <- eval(input$modList[[1]]$call[[2]]) 
    } else {
      stop("Unknown model class")
    }
    
  } else if (modClass %in% c("brmsfit")){ # from brms
    
    ## Parameter names and formula:
    
    formulaStr <- input$formula
    formulaStr <- formulaStr[[1]] # extract only first object
    formulaStr <- Reduce(paste, deparse(formulaStr, width.cutoff = 500)) # convert from formula to string object
    
    ## Exttract mean of per-subject posterior:
    coefArray <- coef(input)[[1]] # extract array of dimensions nSub x 4 x nParam
    nSub <- dim(coefArray)[1]
    nParam <- dim(coefArray)[3]
    
    coefMat <- matrix(NA, nSub, nParam) # initialize
    colnames(coefMat) <- row.names(fixef(input))
    for (iSub in 1:nSub){ # iParam1 <- 1
      for (iParam in 1:nParam){ # iParam2 <- 2
        coefMat[iSub, iParam] <- coefArray[iSub, 1, iParam]
      }
    }
    
  } else {
    
    stop("Unknown model class")
    
  }
  
  # -------------------------------------------------------------------------- #
  ## Name for saving:
  
  fileName <- paste0("coefMat_", modClass, "_", formula2handle(formulaStr), ".csv")
  cat(paste0("Save data under ", fileName, "\n"))
  fileNameFull <- paste0(dirs$outputDir, fileName)
  cat(paste0("File path has ", nchar(fileNameFull), " characters\n"))
  if (nchar(fileNameFull) > 260){
    warning("File path too long, shorten\n")
    fileNameFull <- paste0(substr(fileNameFull, 1, 255)[1], ".csv")
  }
  
  # -------------------------------------------------------------------------- #
  ## Save:
  
  write.csv(coefMat, fileNameFull, row.names = F)
  
  cat("Finished :-)\n")
  
}

# ============================================================================ #
#### Quick CIs based on SEs: ####

quickCI <- function(mod, selEff = NULL, level = 0.95, nRound = 2){
  #' Compute CIs for given lme4 model given SEs from model.
  #' @param data mod model objected fitted with lme4.
  #' @param selEff vector of integers, index of effect in model for which to compute effect size (default: 2).
  #' @param savePNG numeric, 0-1, level of CIs (default: 0.95).
  #' @param nRound integer, number of digits after comma to round to (default: 2).
  #' @return print to console.
  require(arm) # for se.fixef
  
  twoSideLevel <- 1 - (1 - level) / 2 # correct for two-sided test
  zVal <- qnorm(twoSideLevel) # respective z-level threshold
  
  ## If no effect selected: print all effects in model (skipping intercept)
  if (is.null(selEff)){selEff <- 2:length(fixef(mod))}
  
  ## Loop through effects:
  for (iEff in selEff){
    
    # print(round(c(fixef(mod)[iEff] - zVal*se.fixef(mod)[iEff], fixef(mod)[iEff] + zVal*se.fixef(mod)[iEff]), nRound))
    tmp <- round(c(fixef(mod)[iEff] - zVal*se.fixef(mod)[iEff], 
                   fixef(mod)[iEff] + zVal*se.fixef(mod)[iEff]), 
                 nRound)  
    cat(paste0(level*100, "%-CIs for ", colnames(model.matrix(mod))[iEff], ": ", 
               paste(tmp, collapse = " "), "\n"))
  }
}

# =============================================================================================== #
#### Print effect from lme4 model: #####

print_effect <- function(mod, eff, nDigit = 3){
  #' Print selected effect from lme4 model
  #' @param mod fitted model
  #' @param eff string, name of effect for which to print effect
  #' @param nDigit integer, number of digits to round after comma, default 2
  #' @return nothing returned, but printed
  require(stringr)
  
  nPad <- nDigit + 2
  
  if (str_sub(eff,-1) == "f"){eff <- paste0(eff,"1")} # add 1 at the end if eff is factor
  
  # Extract output of fixed effects:
  coefs <- summary(mod)$coefficients # extract coefficients
  idx <- which(rownames(coefs)==eff) # find effect back
  if (length(idx)==0){stop(paste0("Effect ", eff, " not found"))} 
  
  ## Retrieve coefficients:
  if (summary(mod)$objClass=="glmerMod"){ # glmer
    
    # Extract relevant info:
    b <- coefs[idx,1]
    se <- coefs[idx,2]
    zScore <- coefs[idx,3]
    pVal <- coefs[idx,4]
    
  } else if (summary(mod)$objClass=="lmerModLmerTest"){ # lmer
    
    # Extract relevant info:
    b <- coefs[idx,1]
    se <- coefs[idx,2]
    dfs <- coefs[idx,3]
    zScore <- coefs[idx,4]
    pVal <- coefs[idx,5]
    
  } else {
    stop("Unknown model type")
  }
  
  # Variable padding of b based on sign:
  bPad <- ifelse(b > 0,nPad,nPad+1) # pad to 5 digits if negative
  zPad <- ifelse(zScore > 0,nPad,nPad+1) # pad to 5 digits if negative
  
  ## Handle b:
  if (round(b,nDigit) == 0){
    bText <- "0"
  } else {
    bText <- str_pad(round(b,nDigit), bPad, side="right", pad="0")
  }
  
  ## Handle se:
  if (round(se,nDigit) == 0){
    seText <- "0"
  } else {
    seText <- str_pad(round(se,nDigit), nPad, side="right", pad="0")
  }
  
  ## Handle statistic for given object:
  if (summary(mod)$objClass=="glmerMod"){
    zStat <- ", z = "
  } else {
    zStat <- paste0(", t(",round(dfs,nDigit),") = ")
  }
  
  if (round(zScore,nDigit) == 0){
    zText <- "0"
  } else {
    zText <- str_pad(round(zScore,nDigit), zPad, side="right", pad="0")
  }
  
  ## Handle very small p-values:
  if (pVal < 0.001){
    pText <- "p < .001"
  } else {
    pText <- paste0("p = ", str_pad(round(pVal,(nDigit+1)), 5, side="right", pad="0")) # p-value: always 5 digits
  }
  
  # Print to console:
  cat(paste0("b = ", bText,
             ", se = ", seText,
             zStat, zText,
             ", ", pText,"\n"))
}

# ============================================================================ #
#### Fit lm per subject: #####

loop_lm_subject <- function(data, formula, isBinomial = F, family = "binomial"){
  #' Perform lm separately for each subject, store coefficients and models, 
  #' one-sample t-test across subjects for each effect, return.
  #' @param data data frame with variable subject and DVs and IVs.
  #' @param formula string with formula to fit in Wilkinson notation.
  #' @param isBinomial boolean, fit generalized lm with binomial link function (T) or not (F).
  #' @param family distribution of DV to use (default: binomial).
  #' @return output list with elements:
  #' output$bVec: vector of size nSub x nEff, b weights for each effect for each subject.
  #' output$modList: list of models of each subject.
  #' Prints results from t-test across subjects for each effect.
  
  # ----------------------------------------------- #
  ## Fixed variables:
  subVar <- "subject_n"
  if (!(subVar %in% names(data))){stop(paste0("subVar ", subVar, "not contained in data"))}
  nDigit <- 3
  
  # ----------------------------------------------- #
  ## Determine number of subjects:
  subVec <- unique(data[, subVar])
  nSub <- length(subVec)
  cat(paste0("Found data from ", nSub, " subjects\n"))
  
  # ----------------------------------------------- #
  ## Fit model for first subject available to determine number of coefficients:
  subIdx <- which(data[, subVar] == subVec[1])
  subData <- data[subIdx, ]
  mod <- lm(formula = formula, data = subData) # fit model
  nEff <- length(mod$coefficients)
  
  # ----------------------------------------------- #
  ## Initialize matrix to store coefficients:
  bMat <- matrix(NA, nrow = nSub, ncol = nEff) # initialize
  modList <- list()
  
  ## Loop through all subjects:
  for (iSub in 1:nSub){ # iSub <- 1
    
    ## Select subject and data:    
    subID <- subVec[iSub]
    cat(paste0("Start subject ", subID, "\n"))
    subIdx <- which(data[, subVar] == subID)
    subData <- data[subIdx, ]
    
    ## Fit model:
    if (isBinomial) {
      if (family == "binomial"){
        mod <- glm(formula = formula, data = subData, family = binomial())
      }
      else if (family == "poisson"){
        mod <- glm(formula = formula, data = subData, family = poisson())
      } else {
        stop("Unknown family for DV")
      }
    } else {
      mod <- lm(formula = formula, data = subData)
    }
    
    bMat[iSub, ] <- mod$coefficients # store data
    modList[[iSub]] <- mod
    
  } # end iSub
  
  # ----------------------------------------------- #
  ## Perform 1-sample t-test across subjects for each effect:
  
  if(nEff > 0) {
    for (iEff in 1:nEff){ # iEff <- 1
      out <- t.test(bMat[, iEff]) # one-sample t-test
      cat(paste0("Effect for ", names(mod$coefficients)[iEff], 
                 ": t(", out$parameter, ") = ", round(out$statistic, nDigit), 
                 ", p = ", out$p.value, "\n"))
    }
    
  } else {
    cat("Intercept-only model, no effects printed\n")
  }
  
  cat("Finished! :-)\n")
  
  # ----------------------------------------------- #
  ## Output:
  output <- list()
  output$bMat <- bMat
  output$modList <- modList
  
  return(output)
}

# =============================================================================================== #
#### Fit lm per subject: #####

loop_glm_subject <- function(data, formula, family = "binomial"){
  #' Wrapper on loop_lm_subject to perform glm separately for each subject.
  #' @param data data frame with variable subject and DVs and IVs.
  #' @param formula string with formula to fit in Wilkinson notation.
  #' @param family distribution of DV to use (default: binomial).
  #' @return modList list with elements:
  #' modList$bVec: vector of size nSub x nEff, b weights for each effect for each subject.
  #' modList$modList: list of models of each subject.
  #' Prints results from t-test across subjects for each effect.
  
  ## Call loop_lm_subject with isGLM = T:
  out <- loop_lm_subject(data, formula, isBinomial = T, family = family)
  
  return(out)
}

# ============================================================================ #
#### Perform likelihood ratio test on two lists of lm models: #####

LRT_on_modList <- function(modList1, modList2){
  #' Loop over both lists of models, compute log-likelihoods, perform LRT. 
  #' @param modList1 list with lm() models as elements.
  #' @param modList2 list with lm() models as elements.
  #' @return modList list with fields:
  
  cat("Perform likelihood ratio test\n")
  
  # modList1 <- modList1$modList
  # modList2 <- modList2$modList
  
  # ----------------------------------------------- #
  ## Check if same number of elements in input lists:
  nSub1 <- length(modList1)
  nSub2 <- length(modList2)
  if (nSub1 != nSub2){stop("List have different numbers of elements")}
  
  # ----------------------------------------------- #
  ## Loop over both lists, extract log-likelihoods:
  
  logLikMat <- matrix(NA, nrow = nSub1, ncol = 2) # initialize
  
  for (iSub in 1:nSub1){
    logLikMat[iSub, 1] <- as.numeric(logLik(modList1[[iSub]]))
    logLikMat[iSub, 2] <- as.numeric(logLik(modList2[[iSub]]))
  }
  
  # ----------------------------------------------- #
  ## Perform one-sided LRT, print, return:
  
  # sum log-likelihood over trials, why not also over subjects?
  # https://stats.stackexchange.com/questions/526936/how-to-aggregate-log-likelihood-score-of-many-modelsmod
  
  ## x^2-value:
  cat(paste0("Sum log-likelihood across subjects\n"))
  logLikVec <- colSums(logLikMat) # aggregate across subjects
  # logLikVec <- colMeans(logLikMat) # aggregate across subjects
  x2val <- -2*(logLikVec[1] - logLikVec[2]) # -2*log likelihood difference is x^2-distributed
  
  ## degrees of freedom:
  df1 <- modList1[[1]]$iter
  if (is.null(df1)){df1 <- 0}
  df2 <- modList2[[1]]$iter
  dfDif <- df2 - df1
  
  ## p-value:
  p <- pchisq(x2val, df = dfDif, lower.tail = FALSE)
  cat(paste0("LRT: x^2(", dfDif, ") = ", x2val, ", p = ", p, "\n"))
  
  cat("Finished :-)\n")
  return(p)
  
}

# ============================================================================ #
#### Perform Bayesian model selection on two lists of lm models: #####

BMS_on_modList <- function(modList1, modList2, n_samples = 1e6){
  #' Loop over both lists of models, compute BICs, perform random-effects Bayesian model selection
  #' using VB_bms() from mattelisi's bmsR package.
  #' @param modList1 list with lm() models as elements.
  #' @param modList2 list with lm() models as elements.
  #' @param n_samples integer, number of iterations in VB_bms (default: 1e6).
  #' @return bms_mod with alpha (Dirichlet parameters), r (expected model frequencies), xp (exceedance probability),
  #' bor (Bayesian omnibus risk), pxp (protected exceedance probabilities).
  
  ## For details on BMS implementation, see:
  # https://github.com/mattelisi/bmsR
  # remotes::install_github("mattelisi/mlisi")
  # require(mlisi)
  # remotes::install_github("mattelisi/bmsR")
  require(bmsR)
  cat("Perform Bayesian model selection using VB_BMS from mattelisi's bmsR package\n")
  
  # ----------------------------------------------- #
  ## Check if same number of elements in input lists:
  nSub1 <- length(modList1)
  nSub2 <- length(modList2)
  if (nSub1 != nSub2){stop("List have different numbers of elements")}
  
  # ----------------------------------------------- #
  ## Loop over both lists, extract BIC per model:
  
  BICmat <- matrix(NA, nrow = nSub1, ncol = 2)
  
  for (iSub in 1:nSub1){
    BICmat[iSub, 1] <- BIC(modList1[[iSub]])
    BICmat[iSub, 2] <- BIC(modList2[[iSub]])
  }
  
  cat(paste0("Mean BICs per model are ", paste0(colMeans(BICmat), collapse = ", "), "\n"))
  
  ## Perform BMS:
  bms_mod <- VB_bms(BICmat, n_samples = n_samples)
  
  cat("Finished :-)\n")
  return(bms_mod)
  
}

# ============================================================================ #
#### Permutation test: #####
permutation_test <- function(x, y, n = 10000){
  #' Permutation test of 2-sided paired samples t-test
  #' Inputs:
  #' x,y = two vectors of equal length that will be used for permutation. Hypothesis x != y is tested
  #' n = number of permutations (default: 10000)
  #' Outputs: 
  #' p = p-value (number of samples in permutation distribution more extreme than critical value)
  
  ## Check inputs:
  if (length(x) != length(y)){
    stop("Error: x and y of different length")
  } 
  
  ## Preparation:
  nS <-  length(x) # number of samples in vectors 
  diff <- x-y # difference vector of x and y 
  testval <- mean(diff, na.rm = T)/sd(diff, na.rm = T) # empirical test statistics (here: Cohen's d)
  permdist <-  rep(NA, n) # initialize
  
  ## Loop over permutations: 
  for(i in 1:n){
    sign <- sample(c(1,-1), nS, replace = T) # sample vectors of signs with -1 or 1
    testVec <- (diff) * sign # multiply differences between x and y with randomly sampled signs
    permdist[i] <- mean(testVec, na.rm = T)/sd(testVec, na.rm = T) # Take the mean, divide by std (ie Cohen's d), store in permutation distribution
  }
  
  ## Compute p-value: 
  p1 <- sum(permdist > testval) / length(permdist) # number samples in permutation distribution that are larger than empirical test statistic. 
  p2 <- sum(permdist < testval) / length(permdist) # number samples in permutation distribution that are smaller than empirical test statistic
  p <- min(p1,p2)*2 # times 2 because 2-sided test
  return(p)
}

# END