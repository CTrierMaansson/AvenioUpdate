# AvenioUpdate 
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FCTrierMaansson%2FAvenioUpdate&count_bg=%23DAFF3E&title_bg=%23031432&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

<img src="AvenioUpdate.png" width="200" align="right">

An R package to be used for updating and exploring NGS data generated with 
the AVENIO pipeline at AUH. 

**This package is only relevant for people working in the Boe Sandahl Sørensen**
**group at Aarhus University Hospital.**
The package requires the ability to access data only members of the BSS group
have the privileges to. 


## Installation

Use devtools to install the most recent version of AvenioUpdate from the GitHub repository.

```R
if (!require(devtools)) install.packages('devtools')

library(devtools)

devtools::install_github("CTrierMaansson/AvenioUpdate")

```

To check if the package has been installed correctly run the following 
command:

```R
library(AvenioUpdate)

result_stats_Info()
```

Which will output a `data.frame` containing information about some of the
information available using the `result_stats()` function.

## Quick guide

As of early 2025 the package contains two **main** functions 
(`add_run_to_list()` & `create_simple_output()`) which are used to update the
Avenio results and get an overview of the results for a single patient.
In addition, the package contains five **smaller** functions which are 
used to get some statistics on the data we have collected and to make sure the
results file is updated correctly. 

Example of how to use the package:

```R
# Reading the results file
results <- readRDS("/PATH/TO/AVENIO_results_patients.rds") 

# Example path to the specific run
test_path <- "/PATH/TO/AVENIO_results/" 

# Updating the results file
results <- add_run_to_list(master_list = results, 
                           Directory = test_path)
```

## Manual

A [html manual](https://github.com/CTrierMaansson/AvenioUpdate/blob/main/AvenioUpdate_manual.html) 
has been created to describe how to use the package and explain how the
individual functions work. 

The manual is rendered using 
[AvenioUpdate_manual.Rmd](https://github.com/CTrierMaansson/AvenioUpdate/blob/main/AvenioUpdate_manual.Rmd) 

## Issues

If you experience any issues when using the package I refer to the 
[issues page](https://github.com/CTrierMaansson/AvenioUpdate/issues) or
write Christoffer an e-mail at [ctm@clin.au.dk](ctm@clin.au.dk)

## Authors

- Christoffer Trier Maansson, 
Department of Clinical Biochemistry, Aarhus University Hospital, Denmark
    
- Lærke Rosenlund Nielsen, 
Department of Clinical Biochemistry, Aarhus University Hospital, Denmark