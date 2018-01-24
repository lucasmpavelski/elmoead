include_directories(.)

file(GLOB header_files ./*.h)
aux_source_directory(. src_files)

add_library(tests STATIC
${src_files}
${header_files}
)

target_link_libraries(tests surrogate hv svm m)
