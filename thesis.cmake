set(CMAKE_C_COMPILER "C:/Program Files/LLVM/bin/clang.exe")
set(CMAKE_CXX_COMPILER "C:/Program Files/LLVM/bin/clang++.exe")
set(OPENMP_LIBRARIES "C:/Program Files/LLVM/lib")
set(OPENMP_INCLUDES "C:/Program Files/LLVM/include")

set(OPENMP_LIBRARIES "C:/Program Files/LLVM/lib")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Xclang -fopenmp")

link_directories(${OPENMP_LIBRARIES})
