# load sources
file(
	GLOB
	source_files
	*.f
)

# library
add_library(gmres
  ${source_files}
)

