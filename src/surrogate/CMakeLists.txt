include_directories(.)

file(GLOB header_files           ./*.h        )

aux_source_directory(.         src_files          )
aux_source_directory(benchmark benchmark_src_files)
aux_source_directory(ml        ml_src_files       )
aux_source_directory(mo        mo_src_files       )
aux_source_directory(obj_dec   obj_dec_src_files  )
aux_source_directory(psp       psp_src_files      )
aux_source_directory(var       var_src_files      )
aux_source_directory(aos       aos_src_files      )
aux_source_directory(aux       aux_src_files      )
aux_source_directory(moead     moead_src_files    )


add_library(surrogate STATIC
${src_files}
${benchmark_src_files}
${ml_src_files}
${mo_src_files}
${obj_dec_src_files}
${psp_src_files}
${var_src_files}
${aos_src_files}
${aux_src_files}
${moead_src_files}
${header_files}
)

add_subdirectory(benchmark/Toolkit)
target_link_libraries(surrogate hv format m wfg_problems inireader)
