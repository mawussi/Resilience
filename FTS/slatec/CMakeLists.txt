# load sources
file(
	GLOB_RECURSE
	source_files
	*.f
)

# library
add_library(
  slatec
  ${source_files}
)

