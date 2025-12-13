# CPMAddPackage(
#   NAME oneTBB
#   GIT_REPOSITORY "https://github.com/oneapi-src/oneTBB.git"
#   GIT_TAG "master"
#   OPTIONS
#   "TBB_TETS OFF"
# )
find_package(TBB REQUIRED)
