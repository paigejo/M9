setwd("~/git/M9/")

# clean leonardDatMod
cleanDR1 = function() {
  dat = read.csv("leonardDatMod_DR1.csv", header=TRUE)
  
  # Notes on data just loaded:
  # plusMinus is ±
  # lots of ?'s.  Maybe just remove those data points for now
  # bp: number of calendar years before 1950
  # rcybp: radiocarbon years before 1950
  
  # make everything strings instead of factor variables
  for(c in 1:ncol(dat)) {
    dat[,c] = as.character(dat[,c])
  }
  
  # fix index 801 before removing data or 375 after, which has unnecessary semicolons at end
  dat[801,7] = "2465-2972"
  dat[801,8] = "2650plusMinus90"
  
  # remove 802nd row
  dat = dat[-802,]
  
  # set soilNo for empty rows to 1 (107-111)
  dat[107:111,6] = "1"
  
  # change the weird "1?" soilNo values to 1
  dat[dat[,6] == "1?",6] = "1"
  
  attach(dat)
  
  # remove rows where event is unknown
  
  n = nrow(dat)
  
  ##### convert rcybp to 2 columns, a high and a low
  
  # split strings around "pm"
  rcybpRoot=strsplit(as.character(rcybp), "plusMinus")
  bpRoot=strsplit(as.character(bp), "-")
  
  # get left part of strings (central estimate) and right part (margin of error)
  rcybpCntr = sapply(rcybpRoot, "[", 1)
  rcybpMOE = sapply(rcybpRoot, "[", 2)
  bpL = as.numeric(sapply(bpRoot, "[", 1))
  bpR = as.numeric(sapply(bpRoot, "[", 2))
  bpCntr = apply(cbind(bpL, bpR), 1, mean)
  bpMOE = (bpR - bpL)*0.5
  
  # convert from string to number (don't worry about warnings)
  rcybpCntr = as.numeric(rcybpCntr)
  rcybpMOE = as.numeric(rcybpMOE)
  
  # these datapoints have max and min instead of center and MOE estimated for some reason
  #NOTE: assume max minus central estimate is like a 95% MOE and the central estimate is a good estimate
  # maxMinI = c(240, 242) # these are the indices after removing data
  maxMinI = c(340, 342) # these are the indices before removing data
  rcybpCntr[maxMinI] = c(mean(c(763, 855)), mean(c(1530, 1700)))
  rcybpMOE[maxMinI] = c(855-mean(c(763, 855)), 1700-mean(c(1530, 1700)))
  
  ##### get ready to fill in blanks of table using this function
  
  # function for filling in blanks in table
  fillBlanks = function(varTable, varInds, changeInds) {
    # loop to fill in blanks in dataset
    for(i in 1:length(varInds)) {
      ind = varInds[i]
      
      ##### get the set of indices with the same dates as the data at the given index
      # first data point
      if(i == 1) {
        inds = 1:changeInds[1]
      }
      # last data point
      else if(ind > max(changeInds)) {
        inds = (changeInds[length(changeInds)]+1):n
      }
      # general case
      else {
        indR = which(changeInds >= ind)[1]
        indL = indR-1
        inds = (changeInds[indL]+1):changeInds[indR]
      }
      
      ##### fill in the blanks
      tmp = t(varTable[inds,])
      tmp[,] = varTable[ind,]
      varTable[inds,] = t(tmp)
    }
    
    return(varTable)
  }
  
  ##### convert certain numbers to numbers if R assumes they are factor variable
  soilNo = as.numeric(as.character(soilNo))
  Lat = as.numeric(as.character(Lat))
  Uncertainty = as.numeric(as.character(Uncertainty))
  
  ##### fill in blanks for lon and lat data using soilNo variable
  #for testing purposes
  orig = cbind(Lon, Lat)
  
  #compute data indices and when groups change
  newLoc = soilNo == 1
  locInds = (1:n)[newLoc]
  changeInds = locInds[2:length(locInds)]-1
  
  # fill in the blanks
  toFillIn = cbind(Lon, Lat)
  filledIn = fillBlanks(toFillIn, locInds, changeInds)
  Lon = filledIn[,1]
  Lat = filledIn[,2]
  
  # check results:
  ids = cbind(Region, Site, soilNo, event)
  test = cbind(Lon, Lat)
  print(cbind(ids[1:30,], test[1:30,], orig[1:30,]))
  print(cbind(ids[(n-20):n,], test[(n-20):n,], orig[(n-20):n,]))
  
  ##### make new data table
  newDat = data.frame(cbind(Region, Site, Lat, Lon, coreNo, soilNo, bpCntr, bpMOE, 
                            rcybpCntr, rcybpMOE, event, subsidence, Uncertainty, 
                            source, quality, method))
  
  ##### remove bad data
  
  # remove the "?" event data points (NOTE may change this in the future)
  isQu = newDat$event == "?"
  newDat = newDat[!isQu,]
  
  # remove the ? subsidence data
  isQu = is.na(as.numeric(as.character(newDat$subsidence)))
  newDat = newDat[!isQu,]
  
  # remove the NA Uncertainty/subsidence/event data
  cols = c(11:13)
  isNA = is.na(newDat[,cols])
  isNA = apply(isNA, 1, any)
  newDat = newDat[!isNA,]
  
  # test newDat:
  print(dat[1:30,])
  print(newDat[1:30,])
  print(dat[(nrow(dat)-30):nrow(dat),])
  print(newDat[(nrow(newDat)-30):nrow(newDat),])
  
  ##### write file
  detach(dat)
  rownames(newDat)=NULL
  
  #write as excel file
  library(xlsx)
  write.xlsx(newDat, file="DR1.xlsx", row.names = FALSE, showNA = FALSE)
  
  #make sure all variable types are right
  numericCols = c(3:10, 12:13)
  factorCols = c(1:2, 11, 14, 15)
  for(c in numericCols) { newDat[,c] = as.numeric(as.character(newDat[,c])) }
  for(c in factorCols) { newDat[,c] = as.factor(as.character(newDat[,c])) }
  
  # Lon values are negative of what they should be
  newDat[,4] = -newDat[,4]
  
  # Don't remove this event!  It's a full or nearly full rupture that is agreed upon existing
#   #remove T10R1 event.  Not sure what it is
#   newDat = newDat[newDat$event != "T10R1",]
  
  # save as rdata file
  dr1 = newDat
  save(dr1, file="DR1.RData")
  
  invisible(NULL)
}
cleanDR1()

