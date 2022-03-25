#!/usr/bin/env bash

R -e 'library(devtools); devtools::document("zaminfluence"); devtools::install_local("zaminfluence", force=TRUE, upgrade="never")'
