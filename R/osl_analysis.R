#' Convert Paleodose and Sample Dose Rate to Age (before present)
#'
#' @param DeGy numeric or vector of numeric values of sample paleodose
#' @param DoseRate numeric or vector of numeric values of sample dose rate
#'
#' @return Age (numeric or vector of numeric values of sample age in years before sampled day)
#' @export
#'
#' @examples
#' doseToAge(6.44, 2040)
doseToAge <- function(DeGy,DoseRate){
  Age <- 1000000 * DeGy / DoseRate
  return(Age)
}


#' Convert sample ages to an age-density data frame
#'
#' @param df data frame with sample results
#' @param AgeField string, name of the Sample Age field
#' @param SDField string, name of the Sample Age Standard Deviation field
#' @param IDField string, name of the Sample ID field
#' @param grainSD numeric greater than zero, value in standard deviations of to calculate the age-density grid
#'
#' @return data frame of age-densities
#' @export
#'
#' @examples
#' myDF <- data.frame(Age = c(1000, 3000, 5000), SD = c(50, 40, 100), Sample = c('x1','x2','x3'))
#' ageDensities <- sampleDensity(myDF)
sampleDensity <- function(df,
                           AgeField = 'Age',
                           SDField='SD',
                           IDField='Sample',
                           grainSD = 2){

  #progress bar setup
  pb<-utils::txtProgressBar(min = 0, max = 1, initial = 0)
  #tpb<-length(levels(df[,IDField]))
  tpb <- length(df[,IDField])
  ipb <- 0

  #extra column names
  cnam <- colnames(df)
  cnam <- cnam[which(cnam != AgeField & cnam != SDField & cnam != IDField)]

  #process data
  outDF <- data.frame(sample = character(),
                      ageGrid = numeric(),
                      densities = numeric(),
                      dy = numeric())
  for(i in levels(df[,IDField])){
    subOSL <- subset(df,df[,IDField] == i)

    compDens <- data.frame(ageGrid=numeric(),densities=numeric())
    for(j in 1:length(subOSL[,IDField])){
      ageGrid <- seq(floor(subOSL[j,AgeField]-grainSD*subOSL[j,SDField]),
                     ceiling(subOSL[j,AgeField]+grainSD*subOSL[j,SDField]),
                     1)
      densities <- stats::dnorm(ageGrid,
                         subOSL[j,AgeField],
                         subOSL[j,SDField])
      tempDF <- data.frame(ageGrid=ageGrid,densities=densities)
      compDens <- rbind(compDens,tempDF)

      ipb <- ipb + 1
      utils::setTxtProgressBar(pb, ipb/tpb)
    }
    sampDens <- stats::aggregate(densities ~ ageGrid, data = compDens, sum)
    sampDens$densities <- sampDens$densities / sum(sampDens$densities)

    dy <- sampDens$densities / max(sampDens$densities)

    subDF <- data.frame(sample = rep(subOSL[1,IDField],length(sampDens$ageGrid)),
                        ageGrid = sampDens$ageGrid,
                        densities = sampDens$densities,
                        dy = dy)
    #insert extra columns
    if(length(cnam)>0){
      for(k in cnam){
        if(is.factor(df[,k])){
          tempvar <- rep(as.character(subOSL[1,k]),length(sampDens$ageGrid))
        }else{tempvar <- rep(subOSL[1,k],length(sampDens$ageGrid))}
        subDF <- cbind(subDF,as.data.frame(tempvar))
      }
    }

    outDF <- rbind(outDF,subDF)
  }

  close(pb) #close progress bar
  if(length(cnam>0)){colnames(outDF) <- c('sample','ageGrid','densities','dy',cnam)}
  return(outDF)
}


#' Calculate the percent of a sample that is younger than a given date
#'
#' @param df data frame of age densities
#' @param AgeField string, name of the Age field
#' @param DensField string, name of the Density field
#' @param IDField string, name of the Sample ID field
#' @param yearBP integer, value of the date to calculate percent younger than
#'
#' @return data frame containing percent younger than value for each sample
#' @export
#'
#' @examples
pct.younger <- function(df,
                        AgeField = 'ageGrid',
                        DensField = 'densities',
                        IDField = 'sample',
                        yearBP = 1000){
  cnam <- colnames(df)
  cnam <- cnam[which(cnam != AgeField & cnam != DensField & cnam != IDField & cnam != 'dy')]

  outDF <- data.frame(sample = factor(),
                      younger = numeric())
  for(i in levels(df[,IDField])){
    subDF <- subset(df,df[,IDField] == i)
    sumDens <- sum(subDF[,DensField])
    youngDens <- subDF[which(subDF[,AgeField] <= yearBP),]
    sumYoung <- sum(youngDens[,DensField])

    tempDF <- data.frame(sample = as.factor(subDF[1,IDField]),
                         younger = sumYoung / sumDens)

    if(length(cnam)>0){
      for(k in cnam){
        if(is.factor(df[,k])){
          tempvar <- as.character(subDF[1,k])
        }else{tempvar <- subDF[1,k]}
        tempDF <- cbind(tempDF,as.data.frame(tempvar))
      }
    }

    outDF <- rbind(outDF,tempDF)
  }
  if(length(cnam>0)){colnames(outDF) <- c('sample','younger',cnam)}
  return(outDF)
}


#' Fast sampling a normal distribution for plotting individual grain ages
#'
#' @param df data frame
#' @param IDField string, name of the Sample ID field
#' @param ID string, name of the sample ID to be analyzed
#' @param AgeField string, name of the Sample Age field
#' @param SDField string, name of the Sample Standard Deviation field
#' @param nSamples integer, number of samples to draw from a normal distribution
#' @param checkAge integer, age of sample to return T or F if the sample distribution is below or above a value
#'
#' @return
#' @export
#'
#' @examples
single.osl <- function(df,
                       IDField,
                       ID,
                       AgeField = 'Age',
                       SDField='SD',
                       nSamples = 5000,
                       checkAge = 0){
  subDF <- subset(df,df[,IDField] == ID)
  subDF <- subDF[order(subDF[,AgeField]),]
  outDF <- data.frame(Age = numeric(), Grain = numeric(),ZeroDose = logical(),CrossZero = logical(),OlderThan = logical())
  for(i in 1:length(subDF[,IDField])){
    grainDist <- stats::rnorm(nSamples,subDF[i,AgeField],subDF[i,SDField])
    if(subDF[i,AgeField] <= 0){ZD <- T}else{ZD <- F}
    if(subDF[i,AgeField]-2*subDF[i,SDField] <= 0){CZ <- T}else{CZ <- F}
    if(subDF[i,AgeField]-2*subDF[i,SDField] <= checkAge){OT <- F}else{OT <- T}
    tempDF <- data.frame(Age = grainDist, Grain = rep(i,length(grainDist)),ZeroDose = rep(ZD,length(grainDist)),CrossZero = rep(CZ,length(grainDist)), OlderThan = rep(OT,length(grainDist)))
    outDF <- rbind(outDF,tempDF)
  }
  return(outDF)
}


# this is going to be pushed to check the git upload.
