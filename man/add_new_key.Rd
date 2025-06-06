% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_new_key.R
\name{add_new_key}
\alias{add_new_key}
\title{Add a new key to Project, Sample_note, or Material}
\usage{
add_new_key(
  key,
  variable,
  synology_path = "//Synology_m1/Synology_folder/AVENIO/"
)
}
\arguments{
\item{key}{\code{Character} of the key name you want to add as possible entry}

\item{variable}{\code{Character} of the variable name in the AVENIO_runs.xlsx
file you want to add the key to.}

\item{synology_path}{\code{Character} string with the full path to the directory
containing AVENIO_runs.xlsx & AVENIO_keys.rds. Default is
"//Synology_m1/Synology_folder/AVENIO/"}
}
\value{
Nothing is returned. The key is just added to the list of possible
entries
}
\description{
This function can add new possible keys for Project, Sample_note, or Material
in the AVENIO_runs.xlsx file.
}
\examples{
add_new_key(key = "test",
            variable = "Project")
}
