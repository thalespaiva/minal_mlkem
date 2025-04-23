add_compile_options(-Wall -Wextra -Wpedantic)
add_compile_options(-Wredundant-decls -Wshadow -Wpointer-arith)
add_compile_options(-fomit-frame-pointer)
add_compile_options(-march=native -mtune=native -O3 -fwrapv)

add_library(cpucycles cpucycles.c)

add_executable(minal_codes minal_codes.c)
target_link_libraries(minal_codes cpucycles)
