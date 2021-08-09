rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
nsample = as.numeric(args[1])
nu = as.numeric(args[2])
sigma = as.numeric(args[3])
method = args[4]

source("./utils/resultGLM.R")

savepath = paste0("nsample-", nsample, "nu-", nu, "sigma-", sigma, "method-",method, ".Rdata")

if(file.exists(savepath)) stop("File exists")

try({
    save(result, file = savepath)
})

compute = function(resultT){
    
}


mergeResult = function(fileList){
    resultT = list()
    for(f in fileList){
        load(f)
        resultT = c(resultT, result)
    }
    return(resultT)
}

outputNew = function(finalMean, finalSD, nn, pp){
    rk = 2
    cat("\\hline\n")
    cat("\\multicolumn{7}{l}{$n =", nn, "\\quad p = ", pp ,
        "$}  \\\\\n\\hline\n", sep = "")
    
    for(rowI in 1:5){
        if(rowI==1) cat("$\\mathrm{SME}_1$")
        if(rowI==2) cat("$\\mathrm{SME}_2$")
        if(rowI==3) cat("TPR")
        if(rowI==4) cat("TNR")
        if(rowI==5) cat("Rank")
        
        for(colJ in 1:6){
            if(rowI==5 && colJ < 5){
                cat( " &  --  ")
            }else{
                cat(" & ")
                cat(format(round(finalMean[rowI,colJ],rk), nsmall = rk))
                cat("(",format(round(finalSD[rowI,colJ],rk), nsmall = rk),")", sep="")
            }
        }
        cat("\\\\", "\n")
        
    }# rowI
    
    cat("\\hline\n")
}


collectResults = function(settingsID, simu_p, repID, allMethod){
    finalMean = matrix(0, 5, 6)
    finalSD = matrix(0,5,6)
    colSel = 1
    for(m in 1:length(allMethod)){
        method = allMethod[m]
        fileList = paste0("./data/setting-",settingsID,
                          "-method-", method,
                          "-p-", simu_p, 
                          "-", repID ,".RData")
        #print(fileList)
        load(fileList)
        resultT = mergeResult(fileList)
        resultT = compute(resultT)
        res = resultT$res
        resSD = resultT$resSD
        # the last two method has two columns
        if(m == length(allMethod) - 1) colSel = c(length(allMethod) - 1, length(allMethod) )
        if(m == length(allMethod)) colSel = c(length(allMethod) + 1, length(allMethod) + 2)
        if(length(res)==4) res = c(res, 0)
        if(length(resSD)==4) resSD = c(resSD, 0)
        finalMean[,colSel] = t(res)
        finalSD[,colSel] = t(resSD)
        colSel = colSel + 1
    }
    return(list(finalMean = finalMean, finalSD = finalSD))
}


finalMean = matrix(0, 5, 6)
finalSD = matrix(0,5,6)
outputNew(finalMean, finalSD, 200, 50)

 
##### Output table #####

M = 6
settingsID = 2
simu_n = 200
#method = "GLMLowRank"
method = "GLMFPCA"
simu_knots = 20
repID = 2
DataType = "GLMhybrid"
sigma_e = 0.5

collectResults_Mixed  = function(M, settingsID, simu_n, method, simu_knots, DataType){
    colSel = 1
    finalMean_Est = matrix(0, M, 1)
    finalSD_Est = matrix(0, M, 1)
    finalMean_Pred = matrix(0, M, 1)
    finalSD_Pred = matrix(0, M, 1)
    fileList = paste0("./data/setting-", settingsID,"-",simu_n, "-", method, "-", simu_knots, "-",repID, "-",
                      DataType, "-", sigma_e, ".RData")
    load(fileList)
    result_est = compute(Final$Est)
    result_pred = compute(Final$Pred)
    resEst = result_est$res
    resSDEst = result_est$resSD
    resPred = result_pred$res
    resSDPred = result_pred$resSD
    
    finalMean_Est[, colSel] = resEst
    finalSD_Est[, colSel] = resSDEst
    finalMean_Pred[, colSel] = resPred
    finalSD_Pred[, colSel] = resSDPred
    return(list(finalMean_Est = finalMean_Est, finalSD_Est = finalSD_Est, finalMean_Pred = finalMean_Pred,
                finalSD_Pred = finalSD_Pred))
}

outputNewMixed = function(finalMean_Est, finalSD_Est, finalMean_Pred, finalSD_Pred, simu_n){
    rk = 3
    cat("\\hline\n")
    # nrow(finalMean_Est)
    cat("\\multicolumn{1}{c}{$n =", simu_n, "\\\\\n\\hline\n", sep = "")
    for (rowI in 1:nrow(finalMean_Est)){
        cat(paste0("$ \\| \\hat{\\beta}_{", rowI, "} - \\beta_{0" ,rowI, "} \\|_{L_2}^2$" ))
        for (colJ in 1:ncol(finalMean_Est)){
            cat(" & ")
            cat(format(round( finalMean_Est[rowI, colJ],rk), nsmall = rk))
            cat("(",format(round(finalSD_Est[rowI, colJ],rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
    cat("\\hline\n")
    for (rowI in 1:nrow(finalMean_Pred)){
        cat(paste0("$ \\| \\hat{\\beta}_{", rowI, "} - \\beta_{0" ,rowI, "} \\|_{X_2}^2$" ))
        for (colJ in 1:ncol(finalMean_Pred)){
            cat(" & ")
            cat(format(round( finalMean_Pred[rowI, colJ],rk), nsmall = rk))
            cat("(",format(round(finalSD_Pred[rowI, colJ],rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
    cat("\\hline\n")
    
}

###### Grpahical data simulation ###

M = 20
settingsID = 3
simu_n = 200
method = "Graph"
simu_knots = 20
repID = 2
DataType = "Graph"

Final_Graph = collectResults_Mixed(M, settingsID, simu_n, method, simu_knots, DataType)
finalMean_Est = Final_Graph$finalMean_Est
finalSD_Est = Final_Graph$finalSD_Est
finalMean_Pred = Final_Graph$finalMean_Pred
finalSD_Pred = Final_Graph$finalSD_Pred

outputNewMixed(finalMean_Est, finalSD_Est, finalMean_Pred, finalSD_Pred, simu_n)


#### Multitask regression simulation ###
settingID = 1
simu_n = 256
simu_knots = 20
allMethod = c("Tikhonov", "FPCA", "RKHS", "Reduce")
repID = 2
DataType = "Reg"
sigma_e = 0.5
collectResults_Reg = function(){
    colSel = 1
    finalMean_Est = matrix(0, 5, 4)
    finalSD_Est = matrix(0, 5, 4)
    finalMean_Pred = matrix(0, 5, 4)
    finalSD_Pred = matrix(0, 5, 4)
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        fileList = paste0("dataServer/setting-", settingID, "-", simu_n, "-" , method, "-", simu_knots, "-",
         repID, "-", "knots", simu_knots, "-", repID, "-", DataType, "-", sigma_e, ".RData")
        load(fileList)
        result_est = compute(Final$Est)
        result_pred = compute(Final$Pred)
        resEst = result_est$res
        resSDEst = result_est$resSD
        resPred = result_pred$res
        resSDPred = result_pred$resSD
        
        finalMean_Est[, colSel] = resEst
        finalSD_Est[, colSel] = resSDEst
        finalMean_Pred[, colSel] = resPred
        finalSD_Pred[, colSel] = resSDPred
        colSel = colSel + 1
    }
    return(list(finalMean_Est = finalMean_Est, finalSD_Est = finalSD_Est, finalMean_Pred = finalMean_Pred,
                finalSD_Pred = finalSD_Pred))
}

outputNewReg = function(){
    rk = 3
    cat("\\hline\n")
    # nrow(finalMean_Est)
    cat("\\multicolumn{1}{c}{$n =", simu_n, "\\\\\n\\hline\n", sep = "")
    for (rowI in 1:nrow(finalMean_Est)){
        cat(paste0("$ \\| \\hat{\\beta}_{", rowI, "} - \\beta_{0" ,rowI, "} \\|_{L_2}^2$" ))
        for (colJ in 1:ncol(finalMean_Est)){
            cat(" & ")
            cat(format(round( finalMean_Est[rowI, colJ],rk), nsmall = rk))
            cat("(",format(round(finalSD_Est[rowI, colJ],rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
    cat("\\hline\n")
    for (rowI in 1:nrow(finalMean_Pred)){
        cat(paste0("$ \\| \\hat{\\beta}_{", rowI, "} - \\beta_{0" ,rowI, "} \\|_{X_2}^2$" ))
        for (colJ in 1:ncol(finalMean_Pred)){
            cat(" & ")
            cat(format(round( finalMean_Pred[rowI, colJ],rk), nsmall = rk))
            cat("(",format(round(finalSD_Pred[rowI, colJ],rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
    cat("\\hline\n")
}


#### Output graph ###
settingID = 3
simu_n = 200
simu_knots = 20
#allMethod = c("Tikhonov", "FPCA", "RKHS", "Reduce")
allMethod = c("Tikhonov", "FPCA","RKHS", "GraphMulti")
repID = 2
#DataType = "Reg"
sigma_e = 0.5

collectResults_Graph = function(){
    colSel = 1
    finalMean_Est = matrix(0, 1, length(allMethod))
    finalSD_Est = matrix(0, 1, length(allMethod))
    finalMean_Pred = matrix(0, 1, length(allMethod))
    finalSD_Pred = matrix(0, 1, length(allMethod))
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        fileList = paste0("dataServer/setting-", settingID, "-", simu_n, "-" , method, "-", simu_knots, 
                          "-", repID, "-", sigma_e, ".RData")
        load(fileList)
        result_est = computeAll(Final$Est)
        result_pred = computeAll(Final$Pred)
        resEst = result_est$res
        resSDEst = result_est$resSD
        resPred = result_pred$res
        resSDPred = result_pred$resSD
        
        finalMean_Est[, colSel] = resEst
        finalSD_Est[, colSel] = resSDEst
        finalMean_Pred[, colSel] = resPred
        finalSD_Pred[, colSel] = resSDPred
        colSel = colSel + 1
    }
}

outputNew_Graph = function(){
    rk = 4
    cat("\\hline\n")
    cat("\\multicolumn{1}{c}{$n =", simu_n, "\\\\\n\\hline\n", sep = "")
    for (rowI in 1:nrow(finalMean_Est)){
        cat(paste0("$ \\| \\hat{\\beta}_{", rowI, "} - \\beta_{0" ,rowI, "} \\|_{L_2}^2$" ))
        for (colJ in 1:ncol(finalMean_Est)){
            cat(" & ")
            cat(format(round( finalMean_Est[rowI, colJ],rk), nsmall = rk))
            cat("(",format(round(finalSD_Est[rowI, colJ],rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
    cat("\\hline\n")
    for (rowI in 1:nrow(finalMean_Pred)){
        cat(paste0("$ \\| \\hat{\\beta}_{", rowI, "} - \\beta_{0" ,rowI, "} \\|_{X_2}^2$" ))
        for (colJ in 1:ncol(finalMean_Pred)){
            cat(" & ")
            cat(format(round( finalMean_Pred[rowI, colJ],rk), nsmall = rk))
            cat("(",format(round(finalSD_Pred[rowI, colJ],rk),nsmall = rk),")", sep="" )
        }
        cat("\\\\", "\n")
    }
    cat("\\hline\n")
}


#### Spectra Data ####
allMethod = c("Tikhonov", "FPCA","RKHS", "Reduce")
nKnots = 30
collectResults_Spectra = function(){
    colSel = 1
    final = matrix(0, 3, length(allMethod))
    for (m in 1:length(allMethod)){
        method = allMethod[m]
        fileList = paste0("./SpectraData/method-", method, "-", nKnots, ".RData")
        load(fileList)
        
        final[, m] = Final$MSP
    }
    return(final)
}

output_Spectra = function(final){
    rk = 3
    cat("\\hline\n")
    for (rowI in 1:nrow(final)){
        cat("\\text{MSP}")
        for (colJ in 1:ncol(final)){
            cat(" & ")
            cat(format(round( final[rowI, colJ],rk), nsmall = rk))
        }
    cat("\\\\", "\n")
    }
}

