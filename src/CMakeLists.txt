
add_library(libsignalserver STATIC)

target_sources(libsignalserver
    PRIVATE
    inputs.cc
    outputs.cc
    image.cc
    image-ppm.cc
    tiles.cc
    models/cost.cc
    models/ecc33.cc
    models/egli.cc
    models/ericsson.cc
    models/hata.cc
    models/itwom3.0.cc
    models/los.cc
    models/sui.cc
    PUBLIC
        FILE_SET HEADERS
    FILES
    common.h
    image.hh
    image-ppm.hh
    inputs.hh
    main.hh
    outputs.hh
    tiles.hh
    models/cost.hh
    models/ecc33.hh
    models/egli.hh
    models/ericsson.hh
    models/fspl.hh
    models/hata.hh
    models/itwom3.0.hh
    models/los.hh
    models/pel.hh
    models/soil.hh
    models/sui.hh
    models/thread_pool.hh
)

target_link_libraries(libsignalserver PUBLIC
    BZip2::BZip2
    ZLIB::ZLIB
)

target_compile_options(libsignalserver PUBLIC
    "-ffast-math"
)
# target_compile_definitions(libsignalserver CROPPING)

# CFLAGS    = -ffast-math -DCROPPING
# CXXFLAGS	= -Wall -O3 -s -ffast-math
# LIBS		= -lm -lpthread -ldl

add_executable(signalserver
    main.cc
)

target_link_libraries(signalserver PUBLIC
  libsignalserver
)

add_executable(testmodels
    models/testmodels.cc
)

target_link_libraries(testmodels PUBLIC
  libsignalserver
)

if(ENABLE_CLANG_TIDY)
    AddClangTidy(libsignalserver)
    AddClangTidy(signalserver)
    AddClangTidy(testmodels)
endif()


## TODO: Create links to signalserverHD and signalserverLIDAR

## TODO: Add Tests
add_executable(gtestmodels)

target_sources(gtestmodels
    PRIVATE
        models/gtest_itwom.cc
)

target_link_libraries(gtestmodels PUBLIC
    libsignalserver
    GTest::gtest
    GTest::gtest_main
)

add_test(NAME models_test COMMAND gtestmodels)

include(GoogleTest)
gtest_discover_tests(gtestmodels)
