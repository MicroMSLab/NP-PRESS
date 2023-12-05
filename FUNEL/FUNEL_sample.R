# CODE UTF-8
# LICENSE Apache-2.0
# AUTHOR

library(xcms)
library(CluMSID)
library(CAMERA)
library(dplyr)
library(stringr)
library(rjson)
library(RSQLite)

#update on 20211213, context-based adduct/fragment calculation
#update on 20211218, fix bug in adduct calculation threshold usage
#update on 20220103, fix bug in DBsearch() function
#update on 20230908, Added mult tabel support. untested.

configs <- as.data.frame(fromJSON(file = ""))

DB <- dbConnect(SQLite(), configs$DatabaseFile)
setwd(configs$wd)
MergeControlPeaklistFile <- configs$MergeControlPeaklistFile
MergeControlPeaklist2File <- configs$MergeControlPeaklist2File
sfolder <- configs$SampleFolder
sfileExtension = configs$SampleFileExtension
DBtype <- configs$dbtype
rulesfile <- configs$rulesfile
isotopet <- as.logical(configs$IsotopeMatch)
sThreshold <- as.double(configs$SampleThreshold)
method <- configs$method
ppm = as.double(configs$ppm)
DeltaMZ <- as.double(configs$DeltaMZ)
peakwidth = c(as.double(configs$lowPeakwidth), as.double(configs$highPeakwidth))
rtCheck = as.double(configs$rtCheck)
snoise = as.double(configs$SampleNoise)
calAdductS = as.logical(configs$calAdductS)
snthresh = as.double(configs$snthresh)
perfwhm = as.double(configs$perfwhm)
consecMissedLimit = as.double(configs$consecMissedLimit)
integrate = as.double(configs$integrate)
looseCtrRemove = as.logical(configs$looseCtrRemove)
looseIsoFilter = as.logical(configs$looseIsoFilter)
rtThrFilter = as.double(configs$rtThrFilter)
foldchangeFilter = as.double(configs$foldchangeFilter)
Cutoff = as.double(configs$cutoff)
mzdiffCal = as.logical(configs$mzdiffCal)
polarity = as.character(configs$polarity)

if (foldchangeFilter == 0) {
  foldchangeFilter <- FALSE
}

if (DeltaMZ == 0) {
  Delta <- FALSE
} else {
  Delta <- TRUE
}

if (Cutoff == 0) {
  remopeak <- FALSE
} else {
  remopeak <- TRUE
}

tabel_list = unlist(strsplit(DBtype, "&"))

calMZdiff <- function(pl, DeltaMZ = 0.01) {
  # compute mz differences with high frequencies within the same pcgroups (similar retentiontime)
  # can be used for context based adduct/fragment calculation
  pcgroups <- levels(factor(pl[, 'pcgroup']))
  pcgroups <- as.numeric(pcgroups)
  pcgroups <- pcgroups[!is.na(pcgroups)]
  pcgroups <- sort(pcgroups)

  tryCatch({
    mzdiffs <- c()
    for (pcgroup in pcgroups) {
      wgroup <- pl[pl[, 'pcgroup'] == pcgroup, 'mz']
      difflist <-
        lapply(1:max((length(wgroup) - 1), 1), function(x) {
          diff(wgroup, lag = x)
        })
      mzdiffs <- c(mzdiffs, unlist(difflist))
    }

    if (length(mzdiffs) < 1) {
      return(FALSE)
    }
    mzdiffs <- abs(mzdiffs)
    # bin mzdiffs
    binMzdiffs <-
      as.data.frame(table(cut(mzdiffs, ceiling((max(mzdiffs) - min(mzdiffs)) /
                                                 (DeltaMZ / 2)
      ))))
    colnames(binMzdiffs) <- c('levels', 'count')
    binMzdiffs <-
      binMzdiffs[order(binMzdiffs$count, decreasing = TRUE),]
    # calculate the mzdifference (mean of bin boundries)
    lower <- as.numeric(sub("\\((.+),.*", "\\1", binMzdiffs$levels))
    upper <-
      as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", binMzdiffs$levels))
    meandiff <- (lower + upper) / 2
    binMzdiffs <-
      data.frame(lower, upper, meandiff, count = binMzdiffs$count)
    binMzdiffs <-
      binMzdiffs[binMzdiffs[, 'count'] > 0 &
                   binMzdiffs[, 'meandiff'] > 3.1,]
    count <- binMzdiffs$count
    # to determine which mzdiff is meaningful, most of the mzdiff will have similar small
    # counts, indicating random values. meaniful mzdiff will have count numbers that are
    # outliers. determine outliers using IQR method.
    binCount <-
      data.frame(table(cut(count, breaks = (min(
        count
      ) - 1):max(count))))
    binCount$countValue <- min(count):max(count)
    # if the counts of many different mzdiffs are the same, such mzdiffs are noise. so find outlier
    # of countvalues that are extremely large
    binCountFreq <- binCount$Freq[binCount$Freq > 0]
    thrCountFreq <-
      (quantile(binCountFreq, 0.75) - quantile(binCountFreq, 0.25)) * 1.5 + quantile(binCountFreq, 0.75)

    thrCountValue <-
      max(binCount$countValue[binCount$Freq > thrCountFreq])
    binMzdiffs <- binMzdiffs[binMzdiffs[, 'count'] > thrCountValue,]
    count <- binMzdiffs$count
    thrCountValue <-
      (quantile(count, 0.75) - quantile(count, 0.25)) * 1.5 + quantile(count, 0.75)
    binMzdiffs <-
      binMzdiffs[binMzdiffs[, 'count'] >= thrCountValue,]
    # at maximum, take the top 20 mzdiff
    if (nrow(binMzdiffs) > 20) {
      binMzdiffs <- binMzdiffs[1:20,]
    }
    if (nrow(binMzdiffs) < 1) {
      return(FALSE)
    } else {
      return(binMzdiffs)
    }
  },
  warning = function(w) {
    print(paste("MY WARNING: ", w))
    return(FALSE)
  },
  error = function(e) {
    print(paste("MY ERROR: ", e))
    return(FALSE)
  })
}

extractPeaklist <- function(filelist, ppm = 20, Delta = FALSE, DeltaMZ = 0.0, perfwhm = 1, peakwidth = c(5, 20),
                            noise = 10000, snthresh = 10, rtCheck = 30, rules = NULL,
                            calAdduct = TRUE, filter = TRUE, method = 'centWave',
                            consecMissedLimit = 1, integrate = 2, remopeak = FALSE, Cutoff = 1000000) {
  if (Delta == FALSE) {
    if (method == 'centWave') {
      xs.s <- xcmsSet(
        filelist,
        method = "centWave",
        ppm = ppm,
        peakwidth = peakwidth,
        noise = noise,
        snthresh = snthresh
      )
    }
    if (method == 'massifquant') {
      xs.s <- xcmsSet(
        filelist,
        method = "massifquant",
        ppm = ppm,
        peakwidth = peakwidth,
        noise = noise,
        snthresh = snthresh,
        consecMissedLimit = consecMissedLimit,
        integrate = integrate
      )
    }
  }
  else {
    if (method == 'centWave') {
      xs.s <- xcmsSet(
        filelist,
        method = "centWave",
        ppm = ppm,
        mzdiff = DeltaMZ,
        peakwidth = peakwidth,
        noise = noise,
        snthresh = snthresh
      )
    }
    if (method == 'massifquant') {
      xs.s <- xcmsSet(
        filelist,
        method = "massifquant",
        ppm = ppm,
        mzdiff = DeltaMZ,
        peakwidth = peakwidth,
        noise = noise,
        snthresh = snthresh,
        consecMissedLimit = consecMissedLimit,
        integrate = integrate
      )
    }
  }

  if (length(filelist) > 1) {
    xs.s <- group.nearest(xs.s, mzVsRTbalance = 1, mzCheck = 0.01, rtCheck = rtCheck)
  }
  xs.s <- xsAnnotate(xs.s, polarity = polarity)
  #for shorter gradient, set perfwhm parameter low. less retention time variance,
  #then lower perfwhm needed. in general, higher perfwhm gives better results
  xs.s <- groupFWHM(xs.s, perfwhm = perfwhm)
  # filter based on theoratical C12/C13 ratio, if M+1 peak has intensity higher
  # then possible value with max possible carbon numbers, or lower than value with
  # one carbon atom, the corresponding M peak will not be considered monoisotopic
  # peak

  #xs.s <- findIsotopes(xs.s, mzabs = 0.008, minfrac = 0.1, filter = filter)
  xs.s <- findIsotopesWithValidation(xs.s, mzabs = 0.008)
  # calculate rules-allowed adduct
  if (calAdduct) {
    xs.s <- groupCorr(xs.s, cor_eic_th = 0.75)
    if (remopeak == TRUE) {
      rindex = which(xs.s@groupInfo[, 'maxo'] < Cutoff)
      if (length(rindex) > 0) {
        nuinde <- c()
        for (i in 1:length(xs.s@pspectra)) {
          reindex <- which(xs.s@pspectra[[i]] %in% rindex)
          xs.s@pspectra[[i]] <- xs.s@pspectra[[i]][-reindex]
          if (length(xs.s@pspectra[[i]]) > 0) {
            nuinde <- c(nuinde, i)
          }
        }
        xs.s@pspectra <- xs.s@pspectra[nuinde]
      }
    }
    xs.s <- findAdducts(xs.s, polarity = polarity, rules = rules)

  }
  return(xs.s)
}

removeControl <- function(c, s, ppm = 20, Delta = FALSE, DeltaMZ = 0.0, rt = 0, isotope = FALSE,
                          foldchange = FALSE) {
  flag <- list()
  if (isotope) {
    # exatract isotope information from 'isotope' column, [65][M+H]2+ to [M+H]2+
    isotopeList <- str_extract(c[, 'isotopes'], '\\[M.*\\].*')
    isotopeList[is.na(isotopeList)] <- '0'
  }
  for (i in 1:nrow(s)) {
    # if the intensity of peaks in sample is lower than ctrThreshold(threshold
    # for control peak extraction), definitely no match in control peaks. 

    # if the peak has no isotopes related, not need to compare, as all control
    # peaks extracted have isotope related
    # if (s[i,'isotopes'] == ''){
    #   next
    # }
    # if the mz and rt both fit within the range, and the intensity in samples is at
    # least foldchange (to avoid removal because of carry-over) stronger than in the
    # control, then remove the signal from samples.
    if (Delta == FALSE) {
      match <- (abs(c[, 'mz'] - s[i, 'mz']) < (ppm * 1e-06 * c[, 'mz']))
    }
    else {
      match <- (abs(c[, 'mz'] - s[i, 'mz']) < (ppm * 1e-06 * c[, 'mz']))
      match <- match & (abs(c[, 'mz'] - s[i, 'mz']) < (DeltaMZ))
    }
    if (rt > 0) {
      match <- match & (abs(c[, 'rt'] - s[i, 'rt']) < rt)
    }
    # if no hit after mz and rt filter, then skip the rest
    if (!(TRUE %in% match)) {
      next
    }
    if (foldchange != FALSE) {
      match <- match &
        ((c[, 'maxo'] / s[i, 'maxo']) > (1 / foldchange)) &
        (c[, 'maxo'] > 0)
    }
    if (!(TRUE %in% match)) {
      next
    }
    if (isotope) {
      sisotope <- str_extract(s[i, 'isotopes'], '\\[M.*\\].*')
      match <- match & (isotopeList == sisotope)
    }
    if (TRUE %in% match) {
      flag <- c(flag, i)
    }
  }
  if (length(flag) > 0) {
    return(s[-unlist(flag),])
  }
  else {
    print('no control signal removed this time')
    return(s)
  }

}

addExactMass <- function(pl,
                         xsa,
                         mzdiffs = FALSE,
                         DeltaMZ = 0.01) {
  pl[, 'exactmass'] <- '-1'
  pl[, 'charge'] <- '-1'
  iso <- xsa@isotopes # isotope information
  dions <- xsa@derivativeIons # adduct information

  #calculate exactmass for the features with adducts assigned
  for (i in 1:nrow(pl)) {
    if (length(dions) == 0) {
      next
    }
    else {
      if (is.null(dions[[as.numeric(rownames(pl[i, ]))]])) {
        next
      } else {
        # if adducts calculated, combine all possible exactmass
        d <- unlist(dions[[as.numeric(rownames(pl[i,]))]])
        d <- unlist(d[grepl('mass.mass', names(d))])
        d <- as.vector(d) #remove names attribute
        d <- round(as.numeric(d), 4)
        pl[i, 'exactmass'] <- paste(d, collapse = ',')
      }
    }
  }

  # if require context-based adduct calculation, extract the entries with exactmass
  # calculated

  # now start to deal with the entries without exactmass calculated
  rowIndex <- c()
  for (i in 1:nrow(pl)) {
    if (pl[i, 'exactmass'] == '-1') {
      rowIndex <- c(rowIndex, i)
    }
  }
  # only entries with exactmass calculated will be anchors for the context-based
  # calculation, otherwise, the exactmass generated from context-based method will be
  # used for context-based exactmass calculation for following entries
  rowEM <- pl[, 'exactmass'] != '-1'

  for (i in rowIndex) {
    # the entries with pcgroup assigned, with pcgroup same as the target mz, with exactmass
    # calculated from the first round will be used for context-based adduct calculation
    index <-
      (pl$pcgroup > 0) & rowEM &
      (pl$pcgroup == pl[i, 'pcgroup']) & (pl$exactmass != '-1')
    if ((typeof(mzdiffs) != 'logical') & (TRUE %in% index)) {
      mzWithEM <- pl[index, 'mz']
      EMList <- pl[index, 'exactmass']
    }

    if (length(iso) == 0 |
        is.null(iso[[as.numeric(rownames(pl[i, ]))]])) {
      # no isotopes are calculated or the entry has no isotope assigned,
      # give exactmass the value of -1
      if ((typeof(mzdiffs) != 'logical') & (TRUE %in% index)) {
        for (mzdiff in mzdiffs) {
          # calculate the difference between the mz and the mz entries with exactmass calculated,
          # if the difference matches the context-based mzdiffs, then the two mz are from the same
          # species
          if (TRUE %in% (abs(abs(mzWithEM - pl[i, 'mz']) - mzdiff) < (DeltaMZ /
                                                                      2))) {
            pl[i, 'exactmass'] <-
              paste(EMList[(abs(abs(mzWithEM - pl[i, 'mz']) - mzdiff) < (DeltaMZ / 2))], collapse = ",")
            pl[i, 'adduct'] <- mzdiff
            # mzdiff list is ordered by likelihood. if the most possible context-based mzdiff 
            # can be matched, then stop searching for other less likely ones
            break
          }
        }
      }
    }
    else {
      charge <- iso[[as.numeric(rownames(pl[i,]))]]$charge
      pl[i, 'charge'] <- charge
      if (length(dions) == 0) {
        # length(dions)==0 means adducts are not calculated at all, then calculate with
        # context-based adduct rules or calculate exactmass assuming it's [M+nH]n+ ion
        flag <- FALSE
        if ((typeof(mzdiffs) != 'logical') & (TRUE %in% index)) {
          for (mzdiff in mzdiffs) {
            # calculate the difference between the mz and the mz entries with exactmass calculated,
            # if the difference matches the context-based mzdiffs, then the two mz are from the same
            # species
            if (TRUE %in% (abs(abs(mzWithEM - pl[i, 'mz']) - mzdiff) < (DeltaMZ /
                                                                        2))) {
              pl[i, 'exactmass'] <-
                paste(EMList[(abs(abs(mzWithEM - pl[i, 'mz']) - mzdiff) < (DeltaMZ / 2))], collapse = ",")
              pl[i, 'adduct'] <- mzdiff
              # mzdiff list is ordered by likelihood. if the most possible context-based mzdiff 
              # can be matched, then stop searching for other less likely ones
              break
            }
          }
        }
        if (!flag) {
          pl[i, 'exactmass'] <-
            as.character(round((pl[i, 'mz'] - 1.0072) * charge, 4))
        }
      }
      else {
        if (is.null(dions[[as.numeric(rownames(pl[i, ]))]])) {
          # if no adducts calculated, calculated exactmass based on context-based rules or charge
          flag <- FALSE
          if ((typeof(mzdiffs) != 'logical') & (TRUE %in% index)) {
            for (mzdiff in mzdiffs) {
              # calculate the difference between the mz and the mz entries with exactmass calculated,
              # if the difference matches the context-based mzdiffs, then the two mz are from the same
              # species
              if (TRUE %in% (abs(abs(mzWithEM - pl[i, 'mz']) - mzdiff) < (DeltaMZ /
                                                                          2))) {
                pl[i, 'exactmass'] <-
                  paste(EMList[(abs(abs(mzWithEM - pl[i, 'mz']) - mzdiff) < (DeltaMZ / 2))], collapse = ",")
                pl[i, 'adduct'] <- mzdiff
                # if the more possible context-based mzdiff can be matched, then stop searching for others
                flag <- TRUE
                break
              }
            }
          }
          if (!flag) {
            pl[i, 'exactmass'] <-
              as.character(round((pl[i, 'mz'] - 1.0072) * charge, 4))
          }
        }
      }
    }
  }
  return(pl)
}

mergePeaklist <- function(pl1, pl2, colname = FALSE, ppm = 20, Delta = TRUE, DeltaMZ = 0.01, rt = 20,
           exactmass = FALSE, adduct = FALSE) {
  # sort pl1 and pl2 based on mz first then based on rt, ascending
  # browser()
  pl1 <- pl1[order(pl1$mz, pl1$rt),]
  pl2 <- pl2[order(pl2$mz, pl2$rt),]
  nullflag = 0
  if (nrow(pl1) < 1) {
    nullflag = nullflag + 1
  }
  if (nrow(pl2) < 1) {
    nullflag = nullflag + 2
  }
  if (nullflag == 0) {

    if (colname != FALSE) {
      # add a new column to pl1, initiate with -1
      pl1[, colname] <- -1
    }
    i2 <- 1
    i1 <- 1
    nr1 <- nrow(pl1)
    while (i2 <= nrow(pl2)) {
      match <- 0

      while (i1 <= nr1) {
        mz1 <- pl1[i1, 'mz']
        mz2 <- pl2[i2, 'mz']
        rt1 <- pl1[i1, 'rt']
        rt2 <- pl2[i2, 'rt']
        if (Delta == TRUE) {
          check <- ((abs(mz1 - mz2) < (ppm * 1e-06 * mz1)) &
                        (abs(mz1 - mz2) < DeltaMZ) &
                        (abs(rt1 - rt2) < rt))
        }
        else {
          check <- ((abs(mz1 - mz2) < (ppm * 1e-06 * mz1)) &
                        (abs(rt1 - rt2) < rt))
        }
        if (check) {
          if (colname != FALSE) {
            pl1[i1, colname] <-
                pl2[i2, colname] #add the intensity value from pl2
          }
          else {
            i1maxo <- pl1[i1, 'maxo']
            i2maxo <- pl2[i2, 'maxo']
            if (length(i1maxo) < 1 | length(i2maxo) < 1) {
              print('maxo_not_found')
            }
            else {
              pl1[i1, 'maxo'] <- max(i1maxo, i2maxo)
            }
          }
          # combine exact mass info from different peaklist
          if (exactmass != FALSE) {
            exactmass1 <- unlist(strsplit(pl1[i1, 'exactmass'], split = ','))
            exactmass2 <-
                unlist(strsplit(pl2[i2, 'exactmass'], split = ','))
            pl1[i1, 'exactmass'] <-
                paste(unique(c(exactmass1, exactmass2)), collapse = ',')
          }
          # combine adduct info from different peaklist
          if (adduct != FALSE) {
            adduct1 <- pl1[i1, 'adduct']
            adduct2 <- pl2[i2, 'adduct']
            pl1[i1, 'adduct'] <-
                paste(adduct1, adduct2)
          }
          match <- 1
          break
        }
        if (Delta == FALSE) {
          if (pl1[i1, 'mz'] < (pl2[i2, 'mz'] * (1 + (1e-06) * ppm))) {
            i1 <- i1 + 1
          }
          else {
            break
          }
        }
        else {
          if (pl1[i1, 'mz'] < (pl2[i2, 'mz'] * (1 + (1e-06) * ppm)) &
                (pl1[i1, 'mz'] < (pl2[i2, 'mz'] + DeltaMZ))) {
            i1 <- i1 + 1
          }
          else {
            break
          }
        }
      }
      if (match < 1) {
        pl1 <- dplyr::bind_rows(pl1, pl2[i2,]) #pl1 must have colname
        #pl1 must have colname. bind_rows can handle lists with different columnnames
      }
      i2 <- i2 + 1
      i1 <-
          max(1, i1 - 10) #to avoid very similar m/z close together
    }
    return(pl1)
  }
  else if (nullflag == 1) {
    print("pl1_notfound")
    return(pl2)
  }
  else if (nullflag == 2) {
    print("pl2_not_found")
    return(pl1)
  }
  else {
    print("both_not_found")
    return(pl1)
  }
}

DBsearchMS <- function(pl, DB, types, ppm = 20, Delta = FALSE, DeltaMZ = 0.0) {
  pl[, 'match'] <- -1
  sqlFun <- function(t, DB, ppm = 20) {

    if (is.numeric(t)) {
      find_mass = t
    }
    else {
      stop(simpleError(paste(i, ' is not a number!')))
    }
    UPexactmass <- as.character(find_mass + find_mass * ppm * 1e-06)
    DOWNexactmass <- as.character(find_mass - find_mass * ppm * 1e-06)
    DeltaLow <- as.character(find_mass - DeltaMZ)
    DeltaHigh <- as.character(find_mass + DeltaMZ)
    if (Delta == TRUE) {
      SQLma <- paste("SELECT `ID` FROM ", as.character(types), " WHERE `pubchem_mass`<'", UPexactmass, "'AND `pubchem_mass`>'", DOWNexactmass, "'AND `pubchem_mass`>'", DeltaLow, "'AND `pubchem_mass`<'", DeltaHigh, "';", sep = "")
    }
    else {
      SQLma <- paste("SELECT `ID` FROM '", as.character(types), "' WHERE `pubchem_mass`<'", UPexactmass, "'AND `pubchem_mass`>'", DOWNexactmass, "';", sep = "")
    }

    res <- dbSendQuery(DB, SQLma)
    SQLresults = dbFetch(res)
    comID = SQLresults$ID
    dbClearResult(res)
    return(comID)
  }

  myFun <- function(x, y, ppm = 20) {
    hit <- unlist(lapply(x, sqlFun, y, ppm))
    return(hit)
  }

  for (i in 1:nrow(pl)) {
    exactmass <- as.numeric(unlist(strsplit(pl[i, 'exactmass'], split = ',')))
    hits <- lapply(exactmass, myFun, DB, ppm = ppm)
    if (length(unlist(hits)) > 0) {
      hits <- as.vector(unlist(hits))
      pl[i, 'match'] <- paste(hits, collapse = ',')
    }
  }
  return(pl)
}

rm(speak)
rm(speaklist)
rm(mg_speaklist)

load(MergeControlPeaklistFile)
if (looseCtrRemove) {
  load(MergeControlPeaklist2File)
}

rules <- data.frame(read.csv(rulesfile, sep = '\t'))

sfilelist <- list.files(path = file.path(getwd(), sfolder),
                        pattern = paste0('\\.', sfileExtension, '$'),
                        ignore.case = TRUE,
                        full.names = TRUE)
# sfielist and sampleName both list file with full.names, so the order is same
sampleName <-
  tools::file_path_sans_ext(list.files(path = file.path(getwd(), sfolder),
                                       pattern = paste0('\\.', sfileExtension, '$'),
                                       ignore.case = TRUE,
                                       full.names = TRUE))
sampleName <- unlist(lapply(sampleName, FUN = function(x) {
  words <- strsplit(x, '/')
  return(words[[1]][length(words[[1]])])
}))

if (length(sampleName) < 1) {
  print('No .mzXML Files Found!')
}


speak <- lapply(sfilelist, extractPeaklist, method = method, ppm = ppm, Delta = Delta, DeltaMZ = DeltaMZ,
                perfwhm = perfwhm, peakwidth = peakwidth, noise = snoise, snthresh = snthresh, rtCheck = rtCheck,
                rules = rules, calAdduct = calAdductS, remopeak = remopeak, Cutoff = Cutoff)
speaklist <- lapply(speak, getPeaklist, intval = 'maxo')

# if mzdiffCal is TRUE, will calculate adduct and in-source fragmentation based on context
if (mzdiffCal) {
  mzdiffs <- lapply(speaklist, calMZdiff, DeltaMZ)
} else {
  mzdiffs <- rep(FALSE, length(speaklist))
}

save(speak, file = 'speak.rData')
# filter out signals below intensity Cutoff
for (i in 1:length(speaklist)) {
  speaklist[[i]] <- speaklist[[i]][speaklist[[i]]$maxo > sThreshold,]
}

save(speaklist, file = 'before_remove_control.rdata')
print('speaklist saved before remove control')

# change the column name of 'maxo' to the sample name
for (i in seq_along(sampleName)) {
  colnames(speaklist[[i]])[colnames(speaklist[[i]]) == 'maxo'] <- sampleName[[i]]
}

# export all peaks, including those without isotopes detected, merge through 
# samples to generate one table first
mg_speaklist <- speaklist[[1]]
# move the intensity column to the last column, then add other sample intensity 
# to the end. easier to analyze
col <- which(colnames(mg_speaklist) == sampleName[1])
colnum <- length(colnames(mg_speaklist))
mg_speaklist <- mg_speaklist[, c(1:(col - 1), (col + 1):colnum, col)]
if (length(speaklist) > 1) {
  for (i in 2:length(speaklist)) {
    pl2 <- speaklist[[i]]
    mg_speaklist <- mergePeaklist(mg_speaklist, pl2, colname = sampleName[i],
                                  exactmass = FALSE, adduct = TRUE)
  }
}
mg_speaklist[is.na(mg_speaklist)] <- -1
# export table to csv file
folderNameString <- strsplit(sfolder, '/')[[1]] #strsplit gives a list containing vector
exportname <- paste0(folderNameString[length(folderNameString) - 1], '_',
                     folderNameString[length(folderNameString)])
write.csv(mg_speaklist, file.path(getwd(), sfolder, paste0(exportname, '_all_ions.csv')))

for (i in seq_along(speaklist)) {
  # recommend to use mzdiffs parameter only for mass features identified as [M],
  # not [M+1] or [M+2]. However the assigned exactmass of [M+1] or [M+2] can be
  # used for exactmass calculation for [M] with no exactmass assigned.
  index <- grep('.\\[M\\].', speaklist[[i]][, 'isotopes'])
  tryCatch({
    if (typeof(mzdiffs[[i]]) != 'logical') {
      speaklist[[i]] <-
        addExactMass(speaklist[[i]][index,],
                     speak[[i]],
                     mzdiffs = mzdiffs[[i]]$meandiff,
                     DeltaMZ = DeltaMZ)
    } else {
      speaklist[[i]] <-
        addExactMass(speaklist[[i]][index,],
                     speak[[i]],
                     mzdiffs = FALSE,
                     DeltaMZ = DeltaMZ)
    }
  },
  error = function(e) {
    message('Error, maybe sample files be input as controls by mistake')
    print(e)
    speaklist[[i]] <- c()
  })
}

#remove empty speaklist
index <- c()
for (i in seq_along(sampleName)) {
  if (nrow(speaklist[[i]]) > 0) {
    index <- c(index, i)
  }
}
speaklist <- speaklist[index]
sampleName <- sampleName[index]


#if no signals left, give error message and terminate
if (length(sampleName) < 1) {
  message('error, all signals exist in controls')
  quit()
}

# change the column name of of sample name to 'maxo', for removec
for (i in seq_along(sampleName)) {
  colnames(speaklist[[i]])[which(colnames(speaklist[[i]]) == sampleName[i])] <- 'maxo'
}

# remove signals in control samples
print('first round: removing signals from controls')
speaklist <-
  lapply(speaklist, function(x) {
    removeControl(mg_cpeaklist, x, ppm = ppm, Delta = Delta, DeltaMZ = DeltaMZ, rt = rtThrFilter, foldchange = foldchangeFilter,
                  isotope = isotopet)
  })
if (looseCtrRemove) {
  print('second round: removing signals from controls')
  speaklist <-
    lapply(speaklist, function(x) {
      removeControl(mg_cpeaklist2, x, ppm = ppm, Delta = Delta, DeltaMZ = DeltaMZ, rt = rtThrFilter, foldchange = foldchangeFilter,
                    isotope = isotopet)
    })
}

# mg_speaklist <- speaklist[[1]]
# filteredSpeaklist <- DBsearchMS(mg_speaklist, DB, types = DBtype, ppm = 20, Delta = Delta, DeltaMZ = DeltaMZ)

#remove empty speaklist
index <- c()
for (i in seq_along(sampleName)) {
  if (nrow(speaklist[[i]]) > 0) {
    index <- c(index, i)
  }
}
speaklist <- speaklist[index]
sampleName <- sampleName[index]

#if no signals left, give error message and terminate
if (length(sampleName) < 1) {
  message('error, all signals exist in controls')
  quit()
}
# change the column name of 'maxo' to the sample name
for (i in seq_along(sampleName)) {
  colnames(speaklist[[i]])[colnames(speaklist[[i]]) == 'maxo'] <- sampleName[[i]]
}

# merge all features into one table
mg_speaklist <- speaklist[[1]]
# move the intensity column to the last column, then add other sample intensity 
# to the end. easier to analyze
col <- which(colnames(mg_speaklist) == sampleName[1])


colnum <- length(colnames(mg_speaklist))
mg_speaklist <- mg_speaklist[, c(1:(col - 1), (col + 1):colnum, col)]
if (length(speaklist) > 1) {
  for (i in 2:length(speaklist)) {
    pl2 <- speaklist[[i]]
    mg_speaklist <- mergePeaklist(mg_speaklist, pl2, colname = sampleName[i],
                                  exactmass = TRUE, adduct = TRUE)
  }
}
mg_speaklist[is.na(mg_speaklist)] <- -1

write.csv(mg_speaklist, file.path(getwd(), sfolder, paste0(exportname, '_afterControlRemove.csv')))

all_hits <- c()

for (i in 1:length(tabel_list))
{
  filteredSpeaklist <- DBsearchMS(mg_speaklist, DB, types = tabel_list[i], ppm = 20, Delta = Delta, DeltaMZ = DeltaMZ)
  filteredSpeaklist[is.na(filteredSpeaklist)] <- -1
  all_hits <- c(all_hits, filteredSpeaklist)
}

write.csv(all_hits, file.path(getwd(), sfolder, paste0(exportname, "_", as.character(isotopet), "_", as.character(Delta), "_", as.character(DBtype), "_", as.character(remopeak), '_MS1sql.csv')))
