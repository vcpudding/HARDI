# This project is designed to be built outside the Insight source tree.
PROJECT(NAMIC_RESEARCH)

#SET(CMAKE_CXX_FLAGS "-g")

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


SET(NAMIC_INCLUDE_DIR
${CMAKE_CURRENT_SOURCE_DIR}/namicModules
)

#Include Directories
INCLUDE_DIRECTORIES(${NAMIC_INCLUDE_DIR}/DiffusionTensorStatistics .)
ADD_EXECUTABLE(tensorFilter  tensorFilter.cxx)
ADD_EXECUTABLE(dwiFilter  dwiFilter.cxx)

#Libraries
TARGET_LINK_LIBRARIES(tensorFilter ITKCommon ITKIO)
TARGET_LINK_LIBRARIES(dwiFilter ITKCommon ITKIO)

