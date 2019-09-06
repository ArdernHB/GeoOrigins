
#' \emph{Rattus praetor} shape dataset
#'
#' A geometric morphometric dataset containing the XY coordinates for 24 2D landmarks
#' defining the shape of the first mandibular molar of 48 \emph{R. praetor}
#' specimens with corresponding latitude and longitude information. These data cover
#' approximately the entire range of this species. This speices has some geographic
#' structure. These data come from Hulme-Beaman et al. 2019.
#'
#' @format A list of 2 objects:
#' \describe{
#'   \item{Lat.Long}{Latitude longitude data for \emph{R. praetor} specimens in the form of a data.frame with 2 columns (Lat, Long) and 48 rows. Rownames are specimen IDs.}
#'   \item{LMs}{Landmark data of \emph{R. praetor} dental shape in the form of an array. The array has 24 rows representing landmarks, with 2 columns representing XY and 48 individuals.}
#' }
#' @source \url{https://doi.org/10.1016/j.mambio.2018.04.002}
"Rpraetor"


#' \emph{Microtus arvalis} shape dataset
#'
#' A geometric morphometric derived dataset from the data of Cucchi et al. 2014.
#' The dataset contains a data.frame for all 553 specimens, which includes relevant
#' information for each specimen. This dataset also contains a Procrustes distance
#' table calculated pairwise between all specimens. For convenience and as the trait
#' boundary correlation calculations in the \code{BoundaryFinder} function is time
#' consuming, two arrays are also provided containing the trait boundary correlation
#' scores for the total dataset and also the Orkney subset. The original landmark data
#' can be found in the supporting information of Cucchi et al. 2014.
#'
#' @format A list of 3 objects:
#' \describe{
#'   \item{Info}{Specimen information for \emph{M. arvalis} dataset in the form of a data.frame with 4 columns (Population, Geographic Division, Lat, Long) and 553 rows. Rownames are specimen IDs.}
#'   \item{Proc.Dist}{Square distance matrix with pairwise Procrustes Distances between all 553 \emph{M. arvalis} specimens.}
#'   \item{Orkney.Boundary}{A list including a dataframe with provenancing area (in m2) and an array of trait boundary finding correlation values for all 131 Orkney specimens for inputting into the \code{PlotBoundaries} function. This is provided as an example because the \code{BoundaryFinder} function can take a long time to compute.}
#'   \item{Total.Boundary}{A list including a dataframe with provenancing area (in m2) and an array of trait boundary finding correlation values for all 553 vole specimens for inputting into the \code{PlotBoundaries} function. This is provided as an example because the \code{BoundaryFinder} function can take a long time to compute.}
#' }
#' @source \url{https://doi.org/10.1111/evo.12476}
"Marvalis"


#' \emph{Fringilla teydea} song dataset
#'
#' A song based dataset constructed from dynamic time-warping to create dissimilarity
#' metrics. The datasets includes 4 objects: a specimens information data.frame, a
#' square dissimilarity matrix a kin to a distance matrix, a boundary correlation array
#' for the total dataset and a boundary correlation array for the NE subregion. The
#' boundary correlation arrays are provided for convenience and as the trait
#' boundary correlation calculations in the \code{BoundaryFinder} function is time consuming.
#'
#' @format A list of 4 objects:
#' \describe{
#'   \item{Info}{Specimen information for \emph{F. teydea} dataset in the form of a data.frame with 3 columns (ID, Lat, Long) and 116 rows.}
#'   \item{SongDisMat}{Square distance matrix with pairwise dynamic time-warping distance between all 116 \emph{F. teydea} specimens.}
#'   \item{Total.Boundary}{A list including a dataframe with provenancing area (in m2) and an array of trait boundary finding correlation values for all 116 Tenerifian specimens for inputting into the \code{PlotBoundaries} function. This is provided as an example because the \code{BoundaryFinder} function can take a long time to compute.}
#'   \item{NESubRegion.Boundary}{A list including a dataframe with provenancing area (in m2) and an array of trait boundary finding correlation values for all 29 NE subregion specimens for inputting into the \code{PlotBoundaries} function. This is provided as an example because the \code{BoundaryFinder} function can take a long time to compute.}
#' }
"Fteydea"

