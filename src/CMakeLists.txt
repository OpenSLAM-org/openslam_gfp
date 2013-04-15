#==========================================================================================================
 
SET(gefp_cl_SRC
    gefp_cl.cpp)

SET(gefp_cl_onequery_SRC
    gefp_cl_onequery.cpp)

SET(libgefp_SRC
    gefp_engine.cpp)

#~ target
ADD_LIBRARY(gefp ${libgefp_SRC})
ADD_EXECUTABLE(gefp_cl ${gefp_cl_SRC})
TARGET_LINK_LIBRARIES(gefp_cl gefp)
ADD_EXECUTABLE(gefp_cl_onequery ${gefp_cl_onequery_SRC})
TARGET_LINK_LIBRARIES(gefp_cl_onequery gefp)


#~ install
install (TARGETS gefp DESTINATION /usr/lib)
install (FILES gefp_engine.hpp DESTINATION /usr/include/gefp)
