#include <mpi.h>

#ifndef MPIX_H
#define	MPIX_H

/**
 * This is a stub function of the original provided by MPIX for the Blue Gene.
 * As a stub it will set pset_Comm as MPI_COMM_WORLD.\n
 * \n
 * Original function has the following description:\n
 * This function is a collective operation that creates a set of communicators (each node
 * seeing only one), where no two nodes in a given communicator are part of the same pset
 * (all have different I/O Nodes), see Figure 7-3 on page 78. The most common use for this
 * function is to coordinate access to the outside world to maximize the number of I/O Nodes.
 * For example, an application that has an extremely high bandwidth per node requirement
 * can run both MPIX_Pset_same_comm_create() and MPIX_Pset_diff_comm_create().
 * Nodes without rank 0 in MPIX_Pset_same_comm_create() can sleep, leaving those with
 * rank 0 independent and parallel access to the functional Ethernet. Those nodes all belong
 * to the same communicator from MPIX_Pset_diff_comm_create(), allowing them to use
 * that communicator instead of MPI_COMM_WORLD for group communication or
 * coordination.
 */
int MPIX_Pset_same_comm_create(MPI_Comm *pset_comm)
{
    int mpiRank=-1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if (mpiRank==0)
    	printf("Using the stub function MPIX_Pset_same_comm_create to simulate a BG implementation...\n");

    *pset_comm = MPI_COMM_WORLD;
    return MPI_SUCCESS;
}
#endif	/* MPIX_H */
