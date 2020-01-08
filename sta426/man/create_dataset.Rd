\name{create_dataset}
\alias{create_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a simulated SingleCellExperiment object.
}
\description{
Function to create a simulated SingleCellExperiment object from real data.
}
\usage{
sce_sim = create_dataset(sce_prep2, n_comp, n_cells, kNN, kNN_subsample, n_sample, logFC, probs, library_size_muscat,verbose)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sce_prep}{SingleCellExperiment object prepared with muscat::prepSim()}
  \item{n_comp}{Numeric. Number of PC to average for creation of artifical cell}
  \item{n_cells}{Numeric or numeric vector. Indicates the number of cells per sample. If numeric each sample contains n_cells cells}
  \item{kNN}{k-nearest neighbours to compute}
  \item{kNN_subsample}{Size of subset of nearest neighbours to use for averaging}
  \item{logFC}{List of length 1 or 2. If length(logFC) == 1 a list of logFC changes for each gene has to be provided. If length(logFC) == 2 the first entry of the list is numeric and corresponds to the logFC of each gene (the same for all, see Details) and the second is the fraction of genes to be DE.}
  \item{library_size_muscat}{Either numeric or a vectore of length n_cell. If numeric == 0 the library size of the sampled cell is used. If numeric != 0 this library size is used for ALL cells. If a vector of library sizes is provided these are used.}
  \item{probs}{List of length 3 and corresponds to probs of muscat::simData. Containing probabilities of a cell belonging to each cluster, sample, and group, respectively.}
  \item{verbose}{Control amount of progress information printed. Either 1 or 2}
}
\details{
If logFC is numeric the all logFC are either +logFC or -logFC with probability 0.5.
}
\value{
SingleCellExperiment object.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Adriano Martinelli
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
muscat
}
\examples{
library(muscat)
library(sta426)
library(irlba)

data("sta426Data")
sce = logNormCounts(sce)

#Extract WT samples
sce_wt = sce[, colData(sce)$group_id == "ctrl"]

sce_prep <- prepSCE(sce_wt,
  cluster_id = "cluster_id",
  sample_id = "sample_id",
  group_id = "group_id",
  drop = FALSE)
sce_prep <- prepSim(sce_prep)
logcounts(sce_prep) = as.matrix(logcounts(sce_prep))
counts(sce_prep) = as.matrix(counts(sce_prep))

#Set params our simulation
n_comp = 10
n_cells = 100
kNN = 10
kNN_subsample = 5
n_sample = 2
logFC = list(magnitude = 4, proportion = 0.1)
n_group = 2
freq_cluster = table(colData(sce_prep)$cluster_id)
freq_cluster = freq_cluster / sum(freq_cluster)
probs = list(cluster = freq_cluster,
             sample = rep(1/n_sample, n_sample),
             group = rep(1/n_group, n_group))
verbose = 2
library_size_muscat = 0

sce_sim = create_dataset(sce_prep, n_comp, n_cells, kNN, kNN_subsample, n_sample, logFC, probs, library_size_muscat,verbose)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
