\name{Italy_Updated}
\docType{data}
\alias{Italy_Updated}
\title{Real GDP per Capita in Italy's 20 Administrative Regions}
\description{
  Annual real GDP per capita for the 20 Italian administrative regions,
  from 1951 through 2024. Each region contributes one observation per year.
  The ARDECO projections for 2025 and 2026 are excluded.
}
\usage{data("Italy_Updated")}
\format{
  A data frame with 1,480 rows and three columns.
  \describe{
    \item{year}{An ordered factor with 74 levels, from 1951 to 2024.}
    \item{gdp}{Numeric real GDP per capita in thousands of euros on the
      ARDECO 2020 reference scale.}
    \item{region}{The administrative region name, stored as a character
      vector linked by geographic code to the underlying source panel.}
  }
  Rows are ordered by year and then by the source geographic code.
}
\details{
  For 1951--1999, the historical reconstruction uses SVIMEZ chapter VI,
  Tables 3 and 4. The two SVIMEZ volume series, with reference years
  1970 and 2000, are linked at their common year, 1970. Real GDP per
  capita is calculated by multiplying nominal GDP per capita by the
  ratio of linked volume GDP to nominal GDP.

  One common multiplier (approximately 1.478612) places the historical
  values on the ARDECO reference scale. It is the median of the
  ARDECO-to-SVIMEZ real GDP per capita ratios for all 20 regions over
  1995--2009. A common multiplier preserves within-region growth rates
  and the cross-sectional shape of the historical data. For 2000--2024,
  the panel uses ARDECO series SOVGDP (EUR2020, NUTS 2024).
  ARDECO's Bolzano and Trento observations are combined into
  Trentino--Alto Adige using annual population weights from SNPTD.
  Values are divided by 1,000. Across regions at the 2000 handoff, the
  median absolute gap between the two linked sources is about 4.1 percent
  and the maximum about 10.8 percent.

  Giovanni Baiocchi compiled the antecedent \code{Italy} panel using
  CRENoS (1951--1963), Istat (1964--1995), and SVIMEZ (1996--1998).
  The original package object has 21 unlabelled observations per year.
  Its last annual observation is within 0.66 percent of a
  population-weighted reconstruction from the other 20 observations
  throughout 1951--1998, supporting its identification as the Italy-wide
  aggregate. The updated panel follows the 20-administrative-region
  convention and omits the national aggregate. The updated observations
  were independently reconstructed from SVIMEZ and ARDECO.
}
\source{
  Historical reconstruction: SVIMEZ, \emph{150 anni di statistiche
  italiane: Nord e Sud, 1861--2011}, chapter VI, Tables 3 and 4.
  \url{https://datimezzogiorno.svimez.it/en/the-data/svimez-statistics-1861-2011/}

  Modern data: European Commission Joint Research Centre, Annual
  Regional Database of the European Commission (ARDECO), SOVGDP
  (real GDP per capita) and SNPTD (annual average population).
  \url{https://territorial.ec.europa.eu/ardeco/viewer/SOVGDP}
}
\references{
  Baiocchi, G. (2006). \emph{Economic Applications of Nonparametric
  Methods}. Ph.D. thesis, University of York, chapter 2, p. 14,
  footnote 2. \url{https://etheses.whiterose.ac.uk/id/eprint/14117/}

  Giannola, A., Lepore, A., Padovani, R., Bianchi, L., and Miotti, D.
  (eds.) (2011). \emph{150 anni di statistiche italiane: Nord e Sud,
  1861--2011}. Bologna: il Mulino.
}
\examples{
data("Italy_Updated")
boxplot(gdp ~ year, data = Italy_Updated, xlab = "Year",
        ylab = "Real GDP per capita (thousands of 2020 euros)")
}
\keyword{datasets}
