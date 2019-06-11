\name{ComBuilder}
\alias{ComBuilder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComBuilder(area,community,rank,loci,add_list = F, wiki_search = T,group = T, acc_only = F, summ_only = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{area}{
%%     an imported .sbh shape file for the area to be querid in gbif for occurance data
  \item{community}{
%%     A taxonomic discription for the occurance data wanted from NCBI (examples include: Animalia, Chordata, )
  \item{rank}{
%%     The taxonomic rank of the communitry discription
  \item{loci}{
%%     Specify the loci that will be searched on NCBI
  \item{add_list}{
%%     A specific list of species that needs to be added to the gbif occurance data before querying NCBI (e.g. an invasive species watchlist). If no list is required, keep as default F
  \item{wiki_search}{
%%     logical. If true, wikipedia is used as a check after gbif occurance data is obtained to ensure that all species listed match the community term
  \item{group}{
%%     logical. If true this outputs one large .fasta file containing all sequences found in the NCBI search. If false, a folder with small .fasta files for each species is the NCBI output instead.
  \item{acc_only}{
%%     logical. If true the NCBI search will output a series of accession numbers rather than sequences.
  \item{summ_only}{
%%     logical. If true the NCBI search will not look for accession numbers or sequences, and only outputs a dataframe with the number of sequences per species.


}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
