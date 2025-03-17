cmake_minimum_required(VERSION 3.21)

project(md_tian VERSION 1.0
        DESCRIPTION "md_tian an MD code in Fortran for surface scattering"
        LANGUAGES Fortran)
enable_language(Fortran)

#set(FVERSION "-std=f95") use only for gfortran

if(INTEL)
    set(CMAKE_Fortran_COMPILER ifx)
else()
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

if(INTEL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -flto")
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -qmkl -xHost -fp-model strict -check bounds -g -traceback") # -Ofast -ipo
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -flto -xHost -qmkl -Ofast")
    endif()
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -parallel -parallel-source-info=1 ")
    endif()
else()
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -std=legacy")
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -Wall -Wconversion -g -Og -fbacktrace -fbounds-check") # -ffpe-trap=zero,invalid") #-flto
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -ffast-math -march=native") #-flto
    endif()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        set(libraries ${libraries} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
endif()

add_executable(md_tian.x
    # src/md_tian.f90
	src/atom_class.f90
	src/dtrnlspbc.f90
	src/evalulation.f90
	src/fit4_tian.f90
	src/force.f90
	src/mdalgo.f90
	src/md_init.f90
	src/model.f90
	src/nllsq_dbl.f
	src/open_file.f90
	src/output.f90
	src/useful_things.f90
)

target_link_libraries(md_tian.x ${libraries})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
