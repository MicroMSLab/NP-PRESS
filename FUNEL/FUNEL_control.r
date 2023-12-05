# CODE UTF-8
# LICENSE Apache-2.0
# AUTHOR 

library(xcms)
library(CluMSID)
library(CAMERA)
library(dplyr)
library(stringr)
library(rjson)

configs <- as.data.frame(fromJSON(file = ""))

setwd(configs$wd)
MergeControlPeaklistFile <- configs$MergeControlPeaklistFile
MergeControlPeaklist2File <- configs$MergeControlPeaklist2File
cfolder <- configs$ControlFolder
cfileExtension <- configs$ControlFileExtension
rulesfile <- configs$rulesfile
ctrThreshold <- as.double(configs$ControlThreshold)
method <- configs$method
ppm = as.double(configs$ppm)
DeltaMZ <- as.double(configs$DeltaMZ)
peakwidth = c(as.double(configs$lowPeakwidth), as.double(configs$highPeakwidth))
rtCheck = as.double(configs$rtCheck)
cnoise = as.double(configs$ControlNoise)
snthresh = as.double(configs$snthresh)
perfwhm = as.double(configs$perfwhm)
consecMissedLimit = as.double(configs$consecMissedLimit)
integrate = as.double(configs$integrate)
looseCtrRemove = as.logical(configs$looseCtrRemove)
polarity = as.character(configs$polarity)

if (DeltaMZ == 0) {
  Delta <- FALSE
} else {
  Delta <- TRUE
}

extractPeaklist <- function(filelist, ppm = 20, Delta = FALSE, DeltaMZ = 0.0, perfwhm = 1, peakwidth = c(5, 20),
                            noise = 10000, snthresh = 10, rtCheck = 30, rules = NULL,
                            calAdduct = TRUE, method = 'centWave',
                            consecMissedLimit = 1, integrate = 2, remopeak = FALSE, cutoff = 1000000) {
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
  if (remopeak == TRUE) {
    rindex = which(xs.s@groupInfo[, 'maxo'] < cutoff)
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
  # calculate rules-allowed adduct
  if (calAdduct) {
    xs.s <- groupCorr(xs.s, cor_eic_th = 0.75)
    xs.s <- findAdducts(xs.s, polarity = polarity, rules = rules)

  }

  return(xs.s)
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
      match <- 0 #

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
        #pl1 must have colname??? bind_rows can handle lists with different columnnames
      }
      i2 <- i2 + 1
      i1 <-
        max(1, i1 - 10) # to avoid very similar m/z close together
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

rules <- data.frame(read.csv(rulesfile, sep = '\t'))
# start process of control samples
# extract locations of all .mzXML files from the folder
cfilelist <- list.files(path = file.path(getwd(), cfolder),
                        pattern = paste0('\\.', cfileExtension, '$'), ignore.case = TRUE,
                        full.names = TRUE)
cpeak <- lapply(cfilelist, extractPeaklist, method = method, ppm = ppm, Delta = Delta, DeltaMZ = DeltaMZ,
                perfwhm = perfwhm, peakwidth = peakwidth, noise = cnoise, snthresh = snthresh, rtCheck = rtCheck,
                rules = rules, calAdduct = FALSE)

# when only getPeaklist applies to one sample xcmsSet instead of grouped samples,
# intval paremeter not working, though "maxo" selected, "into","intb" will be kept
cpeaklist <- lapply(cpeak, getPeaklist, intval = "maxo")
# merge peaks from control samples
mg_cpeaklist <- cpeaklist[[1]][(cpeaklist[[1]][, 'maxo'] > ctrThreshold),]
mg_cpeaklist <- mg_cpeaklist[, c('mz', 'rt', 'isotopes', 'maxo')]

if (length(cpeaklist) > 1) {
  for (i in 2:length(cpeaklist)) {
    pl2 <- cpeaklist[[i]][(cpeaklist[[i]][, 'maxo'] > ctrThreshold),]
    pl2 <- pl2[, c('mz', 'rt', 'isotopes', 'maxo')]
    mg_cpeaklist <- mergePeaklist(mg_cpeaklist, pl2)
  }
}
save(mg_cpeaklist, file = MergeControlPeaklistFile)

# to make sure all control peaks are removed, loose the extraction parameters,
# extract peaks again
if (looseCtrRemove) {
  cpeak2 <- lapply(cfilelist, extractPeaklist, method = method, ppm = ppm, Delta = Delta, DeltaMZ = DeltaMZ,
                   perfwhm = perfwhm, peakwidth = c(round(peakwidth[1] * 0.5), round(peakwidth[2] * 1.5)),
                   noise = cnoise / 10, snthresh = snthresh, rtCheck = rtCheck, rules = rules,
                   calAdduct = FALSE)
  cpeaklist2 <- lapply(cpeak2, getPeaklist, intval = "maxo")
  mg_cpeaklist2 <- cpeaklist2[[1]][(cpeaklist2[[1]][, 'maxo'] > ctrThreshold),]
  mg_cpeaklist2 <- mg_cpeaklist2[, c('mz', 'rt', 'isotopes', 'maxo')]
  if (length(cpeaklist2) > 1) {
    for (i in 2:length(cpeaklist2)) {
      pl2 <- cpeaklist2[[i]][(cpeaklist2[[i]][, 'maxo'] > ctrThreshold),]
      pl2 <- pl2[, c('mz', 'rt', 'isotopes', 'maxo')]
      mg_cpeaklist2 <- mergePeaklist(mg_cpeaklist2, pl2)
    }
  }
  save(mg_cpeaklist2, file = MergeControlPeaklist2File)
}

print('control_all_done,thanks')

