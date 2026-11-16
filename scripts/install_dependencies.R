packages<-c('ggplot2','dplyr','tidyr','readxl','ape','sf','patchwork','jsonlite','rnaturalearth','rnaturalearthdata','stringi','rmarkdown','flexdashboard','knitr','renv','digest')
missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
cat('Analysis and reporting dependencies available. sf also needs platform-specific GDAL/GEOS/PROJ libraries when compiled from source.\n')
