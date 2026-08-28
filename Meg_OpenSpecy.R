# https://cran.r-project.org/web/packages/OpenSpecy/vignettes/sop.html
#install.packages("OpenSpecy")

############################################
## Load OpenSpecy and reference libraries ##
############################################
## OpenSpecy provides libraries with both FTIR and Raman spectra, so
## we subset for only FTIR matches. We also load the BLoP library and 
## structure the object in an OpenSpecy-compatible way. We similarly
## construct libraries from the macro and cryo data.
library(OpenSpecy)
get_lib()

# Unit conversion function, transmittance to asorbance
trans.to.abs <- function(X){
  return(-log(X/100, base = 10))
}

# Load in the OpenSpecy library with only FTIR spectra
lib <- load_lib(type = "nobaseline")
ftir_ids <- which(lib$metadata$spectrum_type == "ftir")
lib_OS <- as_OpenSpecy(list(wavenumber = lib$wavenumber,
                            spectra = lib$spectra[,..ftir_ids],
                            metadata = lib$metadata[ftir_ids,]))

# Load the BLoP library of Milne and Rochman (2026)
fldr <- "BLoP SI Library Spectra Files copy/ATR-FTIR Files/BLoP .csv files - with FLoPP duplicates"
pths <- list.files(path = fldr, full.names = TRUE)
nms <- sapply(strsplit(list.files(path = fldr, full.names = FALSE), split = "[.]"), "[[", 1)
nms <- gsub(" ", "_", nms)

materials <- sapply(pths, function(x) read.csv(file = x, header = F)[,2])
materials <- apply(materials, 2, trans.to.abs)
materials <- as.data.frame(materials[-c(1:143, 1869),])
colnames(materials) <- nms

lib_blop <- as_OpenSpecy(x = round(read.csv(paste0(fldr, "/Bamboo Bioplastic 1. Brown Straw Fragment.csv"), 
                                            header = FALSE)[-c(1:143, 1869),1]), 
                         spectra = materials, 
                         metadata = NULL)

# Build the macro library
fldr <- "test_files_July20/macro"
pths <- list.files(path = fldr, full.names = TRUE)
nms <- sapply(strsplit(list.files(path = fldr, full.names = FALSE), split = "[.]"), "[[", 1)
nms <- gsub(" ", "_", nms)

macros <- as.data.frame(sapply(pths, 
                               function(x) read.csv(file = x, header = F)[,2]))
colnames(macros) <- nms

lib_macro <- as_OpenSpecy(x = read.csv(paste0(fldr, "/AMF white_macro.csv"), 
                                             header = FALSE)[,1], 
                         spectra = macros, 
                         metadata = NULL)

# Build the cryo library
fldr <- "test_files_July20/cryo"
pths <- list.files(path = fldr, full.names = TRUE)
nms <- sapply(strsplit(list.files(path = fldr, full.names = FALSE), split = "[.]"), "[[", 1)
nms <- gsub(" ", "_", nms)

cryos <- as.data.frame(sapply(pths, 
                               function(x) read.csv(file = x, header = F)[,2]))
colnames(cryos) <- nms

lib_cryo <- as_OpenSpecy(x = read.csv(paste0(fldr, "/AMF_cryo.csv"), 
                                             header = FALSE)[,1], 
                          spectra = cryos, 
                          metadata = NULL)


###############################
## Read in and organize data ##
###############################
dat <- read.csv(file = "test_files_July20/cryo/SB_cryo.CSV", header = F)
colnames(dat) <- c("wavenumber", "spectra")

## Convert transmittance to absorption (if needed)
#dat$spectra <- trans.to.abs(dat$spectra)


##############
## Raw data ##
##############
## Convert data into OpenSpecy object. This just restructures the data 
## into a more specific format the functions are expecting. The 
## check_OpenSpecy() function should return TRUE.
dat_raw <- as_OpenSpecy(x = dat)
check_OpenSpecy(dat_raw)

## Plot raw spectra
#plot(dat_raw)


##################
## Process data ##
##################
## Data are processed with the process_spec() function which has a 
## number of arguments to appropriately smooth each spectra.
## After adjusting the processing arguments we can compare the processed
## plot to the original. Make sure that the baseline is at zero, all
## peaks are captured, and extraneous noise is smoothed away. This may
## take a bit of trial and error and will require playing with the
## following arguments:
##
## (1) conform_spec_args - A list which changes the current spectral and
## wavelength values into a structure which will match the library.
## The range argument makes sure the wavelength ranges agree. The res
## argument specifies the step-size between wavelengths and is set to 6
## to match OpenSpecy or 1 to match BLoP.
##
## (2) restrict_range_args - A list indicating values to truncate the spectra.
##
## (3) smooth_intens_args - A list specifying how the spectral plot is smoothed.
## Larger polynomial values will capture more peaks and noise, smaller
## values with smooth more. Larger derivative values will pick out more
## individual peaks, but too large will find false peaks.
##
## (4) subtr_baseline_args - The degree argument specifies the polynomial
## structure used in determining the baseline for the curve. Larger
## values seem to increase how strongly values are pulled toward 0.

dat_OS <- process_spec(x = dat_raw,
                       conform_spec_args = list(range = lib_OS$wavenumber, 
                                                res = 6),
                       conform_spec = T,
                       adj_intens = T,
                       restrict_range = T,
                       flatten_range = T,
                       subtr_baseline = T,
                       smooth_intens = T,
                       smooth_intens_args = list(polynomial = 3, 
                                                 window = 11, 
                                                 derivative = 0),
                       make_rel = T)
plotly_spec(dat_OS, dat_raw)

dat_blop <- process_spec(x = dat_raw,
                         conform_spec_args = list(range = lib_blop$wavenumber,
                                                  res = 1),
                         conform_spec = T,
                         adj_intens = F,
                         restrict_range = F,
                         flatten_range = F,
                         subtr_baseline = F,
                         smooth_intens = F,
                         make_rel = F)
plotly_spec(dat_blop, dat_raw)


####################################################
## Compare processed spectra to reference library ##
####################################################
## Compare the data to each reference library. Use the processed data for 
## OpenSpecy as this library is built with processed data. Use raw data
## for other comparisions since the other libraries are built from raw data.
## The top_n argument specifies how many top matches to return. The 
## add_library_metadata argument indicates a file with information about 
## the matched spectra. Notice there may be some warnings. These reflect 
## the restricted range we created when processing the data. 
## The warning happens because there are no values outside of the
## restricted range where the library is looking to certain spectra.

# OpenSpecy library
matches_OS_raw <- match_spec(x = dat_OS, 
                             library = lib_OS, 
                             top_n = 5,
                             na.rm = T,
                             add_library_metadata = "sample_name")
matches_OS <- sort_by(matches_OS_raw[, c("object_id", "match_val", "spectrum_identity")],
                      ~ match_val, decreasing = T)

# BLoP library
matches_blop <- match_spec(x = dat_blop, 
                           library = lib_blop, 
                           top_n = 5,
                           na.rm = T)

# Macro library
matches_macro <- match_spec(x = dat_raw, 
                            library = lib_macro, 
                            top_n = 5,
                            na.rm = T)

# Cryo library
matches_cryo <- match_spec(x = dat_raw, 
                           library = lib_cryo, 
                           top_n = 5,
                           na.rm = T)

# Display match results
matches_OS
matches_blop
matches_macro
matches_cryo


## Plot the matches to assess the accuracy. We can inspect all the plots 
## and use visual inspection to help with a positive ID. 

# OpenSpecy plots
plotly_spec(dat_OS, filter_spec(lib_OS, logic = matches_OS_raw[[1,"library_id"]]))
plotly_spec(dat_OS, filter_spec(lib_OS, logic = matches_OS_raw[[2,"library_id"]]))
plotly_spec(dat_OS, filter_spec(lib_OS, logic = matches_OS_raw[[3,"library_id"]]))
plotly_spec(dat_OS, filter_spec(lib_OS, logic = matches_OS_raw[[4,"library_id"]]))
plotly_spec(dat_OS, filter_spec(lib_OS, logic = matches_OS_raw[[5,"library_id"]]))

# BLoP plots
plotly_spec(dat_blop, filter_spec(lib_blop, logic = "Mater-Bi_Bioplastic_1"))
plotly_spec(dat_blop, filter_spec(lib_blop, logic = "Mater-Bi_Bioplastic_2"))
plotly_spec(dat_blop, filter_spec(lib_blop, logic = "Flaxstic_Bioplastic_1"))
plotly_spec(dat_blop, filter_spec(lib_blop, logic = "PLA_Bioplastic_1"))

# Macro plots
plotly_spec(dat_raw, filter_spec(lib_macro, logic = "AMF_white_macro"))
plotly_spec(dat_raw, filter_spec(lib_macro, logic = "AMF_macro"))
plotly_spec(dat_raw, filter_spec(lib_macro, logic = "DS_macro"))
plotly_spec(dat_raw, filter_spec(lib_macro, logic = "P_macro"))
plotly_spec(dat_raw, filter_spec(lib_macro, logic = "PB_macro"))
plotly_spec(dat_raw, filter_spec(lib_macro, logic = "SB_macro"))

# Cryo plots
plotly_spec(dat_raw, filter_spec(lib_cryo, logic = ""))
plotly_spec(dat_raw, filter_spec(lib_cryo, logic = ""))
plotly_spec(dat_raw, filter_spec(lib_cryo, logic = ""))
plotly_spec(dat_raw, filter_spec(lib_cryo, logic = ""))



#######################
## Calculate indices ##
#######################
library(DescTools)
## The function below calculates the ratio of areas over two wavelength
## ranges from processed spectral data. Arguments lower1 and upper1 specify 
## the lower and upper ranges for the numerator, lower2 and upper2 do the
## same for the denominator. Argument dat is the spectral data.
calc.ind <- function(dat, lower1, upper1, lower2, upper2){
  ratios <- data.frame(lower1 = lower1,
                       upper1 = upper1,
                       lower2 = lower2,
                       upper2 = upper2,
                       ratio = NA)
  for (rw in 1:nrow(ratios)) {
    ## Calculate area under processed curve for each range.
    ## Uses default trapezoid interpolation method.
    area1 <- AUC(x = dat$wavenumber, y = dat$spectra$spectra,
                 from = ratios$lower1[rw], to = ratios$upper1[rw])
    area2 <- AUC(x = dat$wavenumber, y = dat$spectra$spectra,
                 from = ratios$lower2[rw], to = ratios$upper2[rw])
    
    ratios$ratio[rw] <- area1/area2
  }
  
  return(ratios)
}

## If you want to compute over multiple ranges, the above function is 
## easily vectorized. For each upper/lower argument, instead of an 
## integer enter a vector with the arguments for the first ratio 
## in the first spot, the second ratio in the second spot, etc.
calc.ind(dat = dat_OS, 
         lower1 = c(3100, 1600, 1000), 
         upper1 = c(3700, 1800, 1200),
         lower2 = c(2780, 2780, 2780), 
         upper2 = c(2970, 2970, 2970))



#################################################
## Example comparing raw and processed matches ##
#################################################
## As an example, we use PVC spectral data to demonstrate why processing
## is needed for OpenSpecy library comparison but not for BLoP, macro, nor 
## cryo library comparison


## Load the data as usual. P for OpenSpecy and AMF for BLoP
dat <- read.csv(file = "test_files_July20/cryo/P_cryo.CSV", header = F)
colnames(dat) <- c("wavenumber", "spectra")
dat_P <- as_OpenSpecy(x = dat)

dat <- read.csv(file = "test_files_July20/cryo/AMF_cryo.CSV", header = F)
colnames(dat) <- c("wavenumber", "spectra")
dat_AMF <- as_OpenSpecy(x = dat)

## Create raw versions of data
dat_P_raw <- process_spec(x = dat_P,
                        conform_spec_args = list(range = lib_OS$wavenumber, 
                                                 res = 6),
                        conform_spec = T,
                        adj_intens = F,
                        restrict_range = F,
                        flatten_range = F,
                        subtr_baseline = F,
                        smooth_intens = F,
                        make_rel = F)
dat_AMF_raw <- process_spec(x = dat_AMF,
                          conform_spec_args = list(range = lib_blop$wavenumber,
                                                   res = 1),
                          conform_spec = T,
                          adj_intens = F,
                          restrict_range = F,
                          flatten_range = F,
                          subtr_baseline = F,
                          smooth_intens = F,
                          make_rel = F)

## Create processed versions of data
dat_P_proc <- process_spec(x = dat_P,
                        conform_spec_args = list(range = lib_OS$wavenumber, 
                                                 res = 6),
                        conform_spec = T,
                        adj_intens = T,
                        restrict_range = T,
                        flatten_range = T,
                        subtr_baseline = T,
                        smooth_intens = T,
                        smooth_intens_args = list(polynomial = 3, 
                                                  window = 11, 
                                                  derivative = 0),
                        make_rel = F)
dat_AMF_proc <- process_spec(x = dat_AMF,
                          conform_spec_args = list(range = lib_blop$wavenumber, 
                                                   res = 1),
                          conform_spec = T,
                          adj_intens = T,
                          restrict_range = T,
                          flatten_range = T,
                          subtr_baseline = T,
                          smooth_intens = T,
                          smooth_intens_args = list(polynomial = 3, 
                                                    window = 11, 
                                                    derivative = 0),
                          make_rel = F)


## OpenSpecy matches
matches_P_raw <- match_spec(x = dat_P_raw, 
                              library = lib_OS, 
                              top_n = 5,
                              na.rm = T,
                              add_library_metadata = "sample_name")

matches_P_proc <- match_spec(x = dat_P_proc, 
                              library = lib_OS, 
                              top_n = 5,
                              na.rm = T,
                              add_library_metadata = "sample_name")

## BLoP matches
matches_AMF_raw <- match_spec(x = dat_AMF_raw, 
                            library = lib_blop, 
                            top_n = 5,
                            na.rm = T)

matches_AMF_proc <- match_spec(x = dat_AMF_proc, 
                            library = lib_blop, 
                            top_n = 5,
                            na.rm = T)

## Inspect match results
sort_by(matches_P_raw[, c("object_id", "match_val", "spectrum_identity")],
        ~ match_val, decreasing = T)
sort_by(matches_P_proc[, c("object_id", "match_val", "spectrum_identity")],
        ~ match_val, decreasing = T)

matches_AMF_raw
matches_AMF_proc

## Visually compare best matches
plotly_spec(dat_P_raw, filter_spec(lib_OS, logic = matches_P_proc[[1,"library_id"]]),
            line = list(color = "blue", width = 5), 
            line2 = list(dash = "dot", color = "red", width = 5),
            paper_bgcolor = "white", 
            plot_bgcolor = "white",
            font = list(color = "black", size = 20))
plotly_spec(dat_P_proc, filter_spec(lib_OS, logic = matches_P_proc[[1,"library_id"]]),
            line = list(color = "blue", width = 5), 
            line2 = list(dash = "dot", color = "red", width = 5),
            paper_bgcolor = "white", 
            plot_bgcolor = "white",
            font = list(color = "black", size = 20))

plotly_spec(dat_AMF_raw, filter_spec(lib_blop, logic = matches_AMF_raw[[1,"library_id"]]),
            line = list(color = "blue", width = 5), 
            line2 = list(dash = "dot", color = "red", width = 5),
            paper_bgcolor = "white", 
            plot_bgcolor = "white",
            font = list(color = "black", size = 20))
plotly_spec(dat_AMF_proc, filter_spec(lib_blop, logic = matches_AMF_raw[[1,"library_id"]]),
            line = list(color = "blue", width = 5), 
            line2 = list(dash = "dot", color = "red", width = 5),
            paper_bgcolor = "white", 
            plot_bgcolor = "white",
            font = list(color = "black", size = 20))

## Export P data
hold <- cbind.data.frame(wavenumber = filter_spec(lib_OS, logic = matches_P_proc[[1,"library_id"]])$wavenumber,
                         spectra = as.numeric(unlist(filter_spec(lib_OS, logic = matches_P_proc[[1,"library_id"]])$spectra)))
hold <- subset(hold, wavenumber %in% dat_P_raw$wavenumber)
                         
P_comp <- cbind.data.frame(wavenumber = dat_P_raw$wavenumber,
                           spectra_raw = unlist(dat_P_raw$spectra)/max(unlist(dat_P_raw$spectra), na.rm=T),
                           spectra_proc = unlist(dat_P_proc$spectra)/max(unlist(dat_P_proc$spectra), na.rm=T),
                           spectra_match = hold$spectra/max(hold$spectra, na.rm=T))
write.csv(P_comp, "P_matching_data.csv")


## Export AMF data
hold <- cbind.data.frame(wavenumber = filter_spec(lib_blop, logic = matches_AMF_raw[[1,"library_id"]])$wavenumber,
                         spectra = as.numeric(unlist(filter_spec(lib_blop, logic = matches_AMF_raw[[1,"library_id"]])$spectra)))

AMF_comp <- cbind.data.frame(wavenumber = dat_AMF_raw$wavenumber,
                             spectra_raw = unlist(dat_AMF_raw$spectra)/max(unlist(dat_AMF_raw$spectra), na.rm=T),
                             spectra_proc = unlist(dat_AMF_proc$spectra)/max(unlist(dat_AMF_proc$spectra), na.rm=T))
AMF_comp <- subset(AMF_comp, wavenumber %in% hold$wavenumber)
AMF_comp <- cbind.data.frame(AMF_comp, spectra_match = hold$spectra/max(hold$spectra, na.rm=T))
write.csv(AMF_comp, "AMF_matching_data.csv")



#################################
## Cryo vs degradation figures ##
#################################
## Load data
dat <- read.csv(file = "test_files_July20/cryo/DS_cryo.CSV", header = F)
colnames(dat) <- c("wavenumber", "spectra")
dat_raw_cryo <- as_OpenSpecy(x = dat)

dat <- read.csv(file = "test_files_July20/digestion/PLA 1.CSV", header = F)
colnames(dat) <- c("wavenumber", "spectra")
dat_raw_deg <- as_OpenSpecy(x = dat)

## Process data
dat_cryo <- process_spec(x = dat_raw_cryo,
                         conform_spec_args = list(range = lib_OS$wavenumber, 
                                                  res = 6),
                         conform_spec = T,
                         adj_intens = T,
                         restrict_range = T,
                         flatten_range = T,
                         subtr_baseline = T,
                         smooth_intens = T,
                         smooth_intens_args = list(polynomial = 3, 
                                                   window = 11, 
                                                   derivative = 0),
                         make_rel = F)
dat_deg <- process_spec(x = dat_raw_deg,
                        conform_spec_args = list(range = lib_OS$wavenumber, 
                                                 res = 6),
                        conform_spec = T,
                        adj_intens = T,
                        restrict_range = T,
                        flatten_range = T,
                        subtr_baseline = T,
                        smooth_intens = T,
                        smooth_intens_args = list(polynomial = 3, 
                                                  window = 11, 
                                                  derivative = 0),
                        make_rel = F)

## Export processed data
write.csv(cbind.data.frame(wavenumber = dat_cryo$wavenumber, spectra = dat_cryo$spectra), 
          "Processed_spectra/DS_cryo.csv", row.names = FALSE)
write.csv(cbind.data.frame(wavenumber = dat_deg$wavenumber, spectra = dat_deg$spectra), 
          "Processed_spectra/DS_degradation.csv", row.names = FALSE)

## Plot
plotly_spec(dat_cryo, dat_deg,
            line = list(color = "blue", width = 5), 
            line2 = list(dash = "dot", color = "red", width = 5),
            paper_bgcolor = "white", 
            plot_bgcolor = "white",
            font = list(color = "black", size = 20),
            make_rel = F)

