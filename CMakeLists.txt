cmake_minimum_required(VERSION 3.10)
project(DeCasteljauAlgorithm)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Trouve Eigen3
find_package(Eigen3 REQUIRED)

# Trouve Boost (pas besoin de composants, car header-only)
find_package(Boost REQUIRED)

# Ajoute les répertoires d'include
include_directories(
    include
    ${Boost_INCLUDE_DIRS}  # ← Pour que Boost::math soit trouvé
)

# Déclare l'exécutable
add_executable(DeCasteljauAlgorithm main.cpp src/Pyramide.cpp)

# Lie Eigen3 (pas besoin pour Boost)
target_link_libraries(DeCasteljauAlgorithm Eigen3::Eigen)
