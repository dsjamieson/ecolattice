
sed '/#include <mpi.h>/d' src/parallel/simulation.h > src/series/simulation.h
sed '/MPI_Finalize()/d' src/parallel/simulation.cpp > src/series/simulation.cpp
sed '/MPI_Finalize()/d' src/parallel/getparameters.cpp > src/series/getparameters.cpp
sed '/MPI_Finalize()/d' src/parallel/competition.cpp > src/series/competition.cpp
sed '/MPI_Finalize()/d' src/parallel/restart.cpp > src/series/restart.cpp

