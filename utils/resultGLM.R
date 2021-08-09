compute = function(resultT){
    res = Reduce("+", resultT)/length(resultT)
    resultS2 = lapply(resultT, function(tmp) tmp^2)
    res2 = Reduce("+", resultS2)/length(resultT)
    resSD = sqrt(res2 - res^2)/sqrt(length(resultT))
    return(list(res = res, resSD = resSD))
}


computeAll = function(resultT){
    resultMean = lapply(resultT, function(tmp) mean(tmp))
    res = Reduce("+", resultMean)/length(resultMean) 
    
    resultS2 = lapply(resultMean, function(tmp) tmp^2)
    res2 = Reduce("+", resultS2)/length(resultT)
    resSD = sqrt(res2 - res^2)/sqrt(length(resultT))
    return(list(res = res, resSD = resSD))
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}



