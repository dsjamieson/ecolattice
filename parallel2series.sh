
sed '/#include <mpi.h>/d' src/mpi/simulation.h > src/series/simulation.h
sed '/MPI_Finalize()/d' src/mpi/simulation.cpp > src/series/simulation.cpp
sed '/MPI_Finalize()/d' src/mpi/getparameters.cpp > src/series/getparameters.cpp
sed '/MPI_Finalize()/d' src/mpi/competition.cpp > src/series/competition.cpp
sed '/MPI_Finalize()/d' src/mpi/restart.cpp > src/series/restart.cpp

