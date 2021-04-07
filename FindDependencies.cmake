# MKL - On Unix the compilervars.sh should be run to find the MKL libraries.
if(DEFINED ENV{MKLROOT})
    set(MKLROOT $ENV{MKLROOT} CACHE DOC "MKL root directory")
else()
    if(WIN32)
        set(MKLPATHS $ENV{ProgramFiles\(x86\)}/Intel* $ENV{PROGRAMFILES}/Intel* $ENV{SystemDrive} $ENV{SystemDrive}/Intel*)
        set(MKLSUFFIXES "compilers_and_libraries/windows")
    elseif(APPLE)
        set(MKLPATHS /opt/intel /intel /usr/local/intel /usr/local/opt/intel)
        set(MKLSUFFIXES "compilers_and_libraries/mac")
    else()
        set(MKLPATHS /opt/intel /intel /usr/local/intel /usr/local/opt/intel $ENV{HOME}/intel $ENV{HOME}/*/intel)
        set(MKLSUFFIXES "compilers_and_libraries/linux")
    endif()
    
    find_file(MKLROOT mkl
		PATHS ${MKLPATHS}
		PATH_SUFFIXES ${MKLSUFFIXES}
		DOC "MKL root directory")
endif()

if(DEFINED MKLROOT)
    set(MKL_INC ${MKLROOT}/include CACHE DOC "MKL include directory")
    set(MKL_LIB_DIR ${MKLROOT}/lib/intel64 CACHE DOC "MKL Library directory")
    
    find_library(MKL_OMP_LIB 
        NAMES iomp5 iomp5md
        PATHS ${MKLROOT}/../lib ${MKLROOT}/../compiler/lib
        PATH_SUFFIXES "intel64" "intel32"
        NO_DEFAULT_PATH
        DOC "MKL OMP Library")		
endif()

if(MKL_INC AND MKL_LIB_DIR AND MKL_OMP_LIB)
	option(USE_MKL "Required for pardiso and iterative solvers" ON)
else()
	option(USE_MKL "Required for pardiso and iterative solvers" OFF)
endif()

if(NOT USE_MKL)
    # OpenMP
    find_package(OpenMP)
endif()

# HYPRE
if(WIN32)
	find_path(HYPRE_INC HYPRE.h 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		PATH_SUFFIXES "include" "include/hypre" "src" "src/include" "src/hypre/include"
        DOC "HYPRE include directory")
	find_library(HYPRE_LIB HYPRE 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
        PATH_SUFFIXES "src" "src/build" "src/mbuild"
		DOC "HYPRE library path")
else()
	find_path(HYPRE_INC HYPRE_IJ_mv.h
        PATHS /opt/hypre* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "include" "include/hypre" "src" "src/include" "src/hypre/include"
		DOC "HYPRE include directory")
	find_library(HYPRE_LIB HYPRE 
        PATHS /opt/hypre* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "src" "src/build" "src/cbuild"
		DOC "HYPRE library path")
endif()	

if(HYPRE_INC AND HYPRE_LIB)		
	option(USE_HYPRE "Required for HYPRE solver" ON)
    mark_as_advanced(HYPRE_INC HYPRE_LIB)
else()
	option(USE_HYPRE "Required for HYPRE solver" OFF)
    mark_as_advanced(CLEAR HYPRE_INC HYPRE_LIB)
endif()

# MMG
if(WIN32)
	find_path(MMG_INC mmg/mmg3d/libmmg3d.h
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		PATH_SUFFIXES "include" "include/mmg*" "src" "build" "build/include"
        DOC "MMG include directory")
	find_library(MMG_LIB mmg3d 
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
        PATH_SUFFIXES "build/lib" "cmbuild/lib" "src/build/lib" "src/cmbuild/lib"
		DOC "MMG library path")
else()
	find_path(MMG_INC mmg/mmg3d/libmmg3d.h
        PATHS /opt/hypre* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "include" "include/mmg" "build" "build/include" "cbuild" "cbuild/include" "src" 
		DOC "MMG include directory")
	find_library(MMG_LIB mmg3d 
        PATHS /opt/mmg* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "build/lib" "cbuild/lib" "src/build/lib" "src/cbuild/lib"
		DOC "MMG library path")
endif()	

if(MMG_INC AND MMG_LIB)		
	option(USE_MMG "Required for MMG use" ON)
    mark_as_advanced(MMG_INC MMG_LIB)
else()
	option(USE_MMG "Required for MMG use" OFF)
    mark_as_advanced(CLEAR MMG_INC MMG_LIB)
endif()

# TETGEN
if(WIN32)
	find_path(TET_INC tetgen.h
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		PATH_SUFFIXES "include" "include/tetgen*" "src" "build" "build/include"
        DOC "TetGen include directory")
	find_library(TET_LIB 
        NAMES tet tetgen
        PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
        PATH_SUFFIXES "build/lib" "cmbuild/lib" "src/build/lib" "src/cmbuild/lib"
		DOC "TetGen library path")
else()
	find_path(TET_INC tetgen.h
        PATHS /opt/tetgen* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "include" "include/tetgen*" "src" "build" "build/include"
		DOC "TetGen include directory")
	find_library(TET_LIB 
        NAMES tet tetgen 
        PATHS /opt/tetgen* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "build/lib" "cmbuild/lib" "src/build/lib" "src/cmbuild/lib"
		DOC "TetGen library path")
endif()	

if(TET_INC AND TET_LIB)		
	option(USE_TETGEN "Required for adaptive remeshing" ON)
    mark_as_advanced(TET_INC TET_LIB)
else()
	option(USE_TETGEN "Required for adaptive remeshing" OFF)
    mark_as_advanced(CLEAR TET_INC TET_LIB)
endif()

# LEVMAR
if(WIN32)
	find_path(LEVMAR_INC levmar.h PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		DOC "Levmar include directory")
	find_library(LEVMAR_LIB levmar PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		DOC "Levmar library path")
else()
	find_path(LEVMAR_INC levmar.h PATHS /usr/local/ /opt/levmar* $ENV{HOME}/* $ENV{HOME}/*/*
		DOC "Levmar include directory")
	find_library(LEVMAR_LIB levmar PATHS /usr/local/ /opt/levmar* $ENV{HOME}/* $ENV{HOME}/*/*
        PATH_SUFFIXES "build" "cbuild" "cmbuild"
		DOC "Levmar library path")
endif()	

if(LEVMAR_INC AND LEVMAR_LIB)		
	option(USE_LEVMAR "Required for optimization in FEBio" ON)
    mark_as_advanced(LEVMAR_INC LEVMAR_LIB)
else()
	option(USE_LEVMAR "Required for optimization in FEBio" OFF)
    mark_as_advanced(CLEAR LEVMAR_INC LEVMAR_LIB)
endif()

# ZLIB
include(FindZLIB)

if(ZLIB_INCLUDE_DIR AND ZLIB_LIBRARY_RELEASE)		
	option(USE_ZLIB "Required for compressing xplt files" ON)
    mark_as_advanced(ZLIB_INCLUDE_DIR ZLIB_LIBRARY_RELEASE)
else()
	option(USE_ZLIB "Required for compressing xplt files" OFF)
    mark_as_advanced(CLEAR ZLIB_INCLUDE_DIR ZLIB_LIBRARY_RELEASE)
endif()

# GSL
if(WIN32)
	find_path(GSL_INC gsl/gsl_spline.h PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		DOC "GSL include directory")
	find_library(GSL_LIB gsl PATHS C::/Program\ Files/* $ENV{HOMEPATH}/* $ENV{HOMEPATH}/*/*
		DOC "GSL library path")
else()
	find_path(GSL_INC gsl/gsl_spline.h PATHS /usr/local/include /opt/gsl* $ENV{HOME}/* $ENV{HOME}/*/*
		DOC "GSL include directory")
	find_library(GSL_LIB NAMES libgsl.a gsl PATHS /usr/local/lib64 /usr/local/lib /opt/gsl* $ENV{HOME}/* $ENV{HOME}/*/*
		DOC "GSL library path")
endif()	

if(GSL_INC AND GSL_LIB)		
	option(USE_GSL "Required for strain energy calculation for certain materials." ON)
    mark_as_advanced(GSL_INC GSL_LIB)
else()
	option(USE_GSL "Required for strain energy calculation for certain materials." OFF)
    mark_as_advanced(CLEAR GSL_INC GSL_LIB)
endif()

if(WIN32)
    set(VCPKG_GSL false CACHE BOOL DOC "Use Vcpkg to link GSL")
endif()
