# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/common/home/twm70/Repos/matter/_deps/googletest-src"
  "/common/home/twm70/Repos/matter/_deps/googletest-build"
  "/common/home/twm70/Repos/matter/_deps/googletest-subbuild/googletest-populate-prefix"
  "/common/home/twm70/Repos/matter/_deps/googletest-subbuild/googletest-populate-prefix/tmp"
  "/common/home/twm70/Repos/matter/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
  "/common/home/twm70/Repos/matter/_deps/googletest-subbuild/googletest-populate-prefix/src"
  "/common/home/twm70/Repos/matter/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/common/home/twm70/Repos/matter/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/common/home/twm70/Repos/matter/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
