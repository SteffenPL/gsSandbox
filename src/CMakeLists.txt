cmake_minimum_required( VERSION 2.6 )

project( gsSandbox )

# boost is not essetial for this project, it is just used to create folders
find_package( Boost COMPONENTS filesystem system REQUIRED )
include_directories ( ${Boost_INCLUDE_DIR})


set( SOURCES
gsOptParametrization.hpp
algorithms.cpp
)

set( HEADERS

gsOptParametrization.h
algorithms.h
)

add_executable( genParam genParam.cpp ${SOURCES} ${HEADERS} )
target_link_libraries( genParam gismo ${IPOPT_LIBRARIES}  ${Boost_LIBRARIES} )

add_executable( testParam testParam.cpp ${SOURCES} ${HEADERS} )
target_link_libraries( testParam gismo ${Boost_LIBRARIES} )

add_executable( optParam optParam.cpp  ${SOURCES} ${HEADERS} )
target_link_libraries( optParam gismo ${IPOPT_LIBRARIES} ${Boost_LIBRARIES} )

add_executable( sandbox sandbox.cpp ${SOURCES} ${HEADERS})
target_link_libraries( sandbox gismo ${IPOPT_LIBRARIES}  )

