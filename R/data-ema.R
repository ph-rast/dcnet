##' Ecological Momentary Assment Data
##'
##' This data is similar to the Fitbit data described in O’Laughlin et al. (2020) and Williams et al. (2020).
##' It contains the jittered data (added random noise to the original data) for 38 participants with each 98 observations.
##' Each participant self-reported their daily score, ranging from 0-100, on the items: "Active", "Excited", and "Interested". Further, "Total Distance" records kilometers walked on each of the 98 days.
##'
##' @format A list containing 38 data frames, each representing a participant with 4 variables measured across 98 time points:
##' \describe{
##'   \item{active}{Score on the "Active" item, ranging from 0-100.}
##'   \item{excited}{Score on the "Excited" item, ranging from 0-100.}
##'   \item{totalDistance}{Total distance walked on a given day in kilometers (KM).}
##'   \item{interested}{Score on the "Interested" item, ranging from 0-100.}
##' }
"ema_fitbit"
