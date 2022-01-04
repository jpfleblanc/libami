find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx)

add_custom_target(Sphinx 
                  COMMAND
                  ${SPHINX_EXECUTABLE} -b html
				   # Tell Breathe where to find the Doxygen output
                  -Dbreathe_projects.Libami=${PROJECT_BINARY_DIR}/DOCS/xml
                  ${SPHINX_SOURCE} ${SPHINX_BUILD}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMENT "Generating documentation with Sphinx")
