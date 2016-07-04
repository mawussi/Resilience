# load sources
file(
	GLOB
	source_files
	*.f
)

# library
add_library(fgmres
  ${source_files}
)

