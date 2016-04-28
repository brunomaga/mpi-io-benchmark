#include <stdio.h>
#include "mpi.h"

#define OUTPUT_FOLDER "./"
#define MB(x) (x*1024*1024)

#ifdef DONT_USE_MPIX
#include <mpix.h>
#else
#include "mpix.stub.h"
#endif


int getDetailsForIOgroup(int *fileNumber_ptr, MPI_Comm *sameIOcomm_ptr);

int main(int argc, char* argv[])
{
    int mpiRank=-1, mpiSize=-1; //Rank and Size in MPI_COMM_WORLD
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    int IoNodeNumber=-1; //id of this process' IO node
    MPI_Comm mpiComm_IO = MPI_COMM_NULL; //MPI_Comm of this process' IO node
    int numberOfIoNodes = getDetailsForIOgroup(&IoNodeNumber, &mpiComm_IO);
    if (mpiRank==0) printf("Number of ranks:%d\n", mpiSize);
    if (mpiRank==0) printf("Number of IO nodes:%d\n", numberOfIoNodes);
    if (mpiRank==0) printf("MPI_WTIME_IS_GLOBAL=%d\n",MPI_WTIME_IS_GLOBAL);

    int mpiRank_IO=-1, mpiSize_IO=-1;
    MPI_Comm_rank(mpiComm_IO, &mpiRank_IO);
    MPI_Comm_size(mpiComm_IO, &mpiSize_IO);

    int err=-1; //MPI error return variable
    double t1=-1,t=-1, t_sum=-1, t_min=-1, t_max=-1; //used for timing

    //write data
    for (long long size = MB(0.1); size < MB(1024); size*=10)
    {
      if (mpiRank==0) printf("\nSIZE = %.3f MB\n", (double) size/MB(1));

      //open file
      MPI_File file; //output file
      char filename[2048]; //filename to write to (one per IO node)
      sprintf(filename, "%s/output.%d", OUTPUT_FOLDER, IoNodeNumber);
      err = MPI_File_open(mpiComm_IO, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
      if (err!=MPI_SUCCESS) printf("Could not open %s for writing. Error %d\n", filename, err);

      unsigned char * buffer = new unsigned char[size];

      //#### MPI_File_write_at
      t1 = MPI_Wtime(); 
      long long * sizes = new long long[mpiSize_IO];
      long long offset = 0;
      MPI_Allgather(&size, 1, MPI_LONG_LONG, sizes, 1, MPI_LONG_LONG, mpiComm_IO);
      for (int i=0; i<mpiRank_IO; i++) offset+=sizes[i];
      MPI_File_write_at(file, offset, buffer, size, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
      MPI_Barrier(MPI_COMM_WORLD);
      t = MPI_Wtime()-t1; 
      if (err!=MPI_SUCCESS) printf("Could not perform MPI_File_write_at. Error %d\n", err);
      //else printf("rank %d: MPI_File_write_at : %.2f\n", mpiRank, t );
      MPI_Allreduce(&t,&t_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&t,&t_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&t,&t_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (mpiRank==0) printf("# MPI_File_write_at: avg %.2f, min %.2f, max %.2f, speed %.2f MB/s\n", t_sum/mpiSize, t_min, t_max, size*mpiSize/MB(1)/t_max);

      //#### MPI_File_write_at
      t1 = MPI_Wtime(); 
      MPI_Allgather(&size, 1, MPI_LONG_LONG, sizes, 1, MPI_LONG_LONG, mpiComm_IO);
      for (int i=0; i<mpiRank_IO; i++) offset+=sizes[i];
      MPI_File_write_at_all(file, offset, buffer, size, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
      MPI_Barrier(MPI_COMM_WORLD);
      t = MPI_Wtime()-t1; 
      if (err!=MPI_SUCCESS) printf("Could not perform MPI_File_write_at_all. Error %d\n", err);
      //else printf("rank %d: MPI_File_write_at : %.2f\n", mpiRank, t );
      MPI_Allreduce(&t,&t_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&t,&t_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&t,&t_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (mpiRank==0) printf("# MPI_File_write_at_all: avg %.2f, min %.2f, max %.2f, speed %.2f MB/s\n", t_sum/mpiSize, t_min, t_max, size*mpiSize/MB(1)/t_max);

     
      //#### MPI_File_write_all
      t1 = MPI_Wtime(); 
      MPI_File_write_all(file, buffer, size, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
      MPI_Barrier(MPI_COMM_WORLD);
      t = MPI_Wtime()-t1; 
      if (err!=MPI_SUCCESS) printf("Could not perform MPI_File_write_all. Error %d\n", err);
      //else printf("rank %d: MPI_File_write_all: %.2f\n", mpiRank, t);
      MPI_Allreduce(&t,&t_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&t,&t_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&t,&t_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (mpiRank==0) printf("# MPI_File_write_all: avg %.2f, min %.2f, max %.2f, speed %.2f MB/s\n", t_sum/mpiSize, t_min, t_max, size*mpiSize/MB(1)/t_max);

      //clean up
      delete [] sizes;
      delete [] buffer;
      err =  MPI_File_close(&file);
      if (err != MPI_SUCCESS) printf("Could not close file. Error %d\n", err);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpiRank==0) printf("Done!\n");
    MPI_Finalize();

    return 0;
}

/// Returns (in the arguments) the number of IO node and the MPI_Comm of the IO node network; as return value, the total number of IO nodes;
int getDetailsForIOgroup(int *fileNumber_ptr, MPI_Comm *sameIOcomm_ptr)
{
    int globalMpiRank=-1;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalMpiRank);

    int fileNumber=-1;//number of this CPU's mpi file
    int numberIOnodes = -1; //number of IOnodes/files

    //1. Split the Comm World into sub groups of nodes facing each IO node
    MPI_Comm sameIOcomm;
    int mpiResult = MPIX_Pset_same_comm_create(&sameIOcomm);
    if (globalMpiRank == 0 && (mpiResult != MPI_SUCCESS || sameIOcomm == MPI_COMM_NULL))
        printf("Could not run MPIX_Pset_same_comm_create in order to get the IO node behind each CPU. FLAG53.\n");

    //2. get a "global fileNumber" for the files to be written [0... Nfiles-1]
    //2.1 gets the lowest index of its new IO node group
    int sameIOminRank, sameIOrank, sameIOsize;
    MPI_Comm_rank(sameIOcomm, &sameIOrank);
    MPI_Comm_size(sameIOcomm, &sameIOsize);
    MPI_Allreduce(&globalMpiRank, &sameIOminRank, 1, MPI_INT, MPI_MIN, sameIOcomm);

    //2.2 On each (IO node) CPUs group they will determine a order for files:
    //each min rank of each IO group will get its order

    MPI_Comm minsRanksComm;
    int groupCode = globalMpiRank == sameIOminRank ? 1 : 0;
    MPI_Comm_split(MPI_COMM_WORLD, groupCode, globalMpiRank, &minsRanksComm);

    if (groupCode == 1) //if he's part of the minimal ranks...
    {
        MPI_Comm_size(minsRanksComm, &numberIOnodes);

        //All minimals get the other minimals rank:
        int * otherMinRanks = new int [numberIOnodes];
        MPI_Allgather(&globalMpiRank, 1, MPI_INT, otherMinRanks, 1, MPI_INT, minsRanksComm);

        //the calculate which fileNumbet they write to [0,1,2, .. Nfiles-1]
        for (int i = 0; i < numberIOnodes; i++)
            if (globalMpiRank == otherMinRanks[i])
                fileNumber = i;

        //They send that fileNumber to all other elements in same IO node
        for (int i = 0; i < sameIOsize; i++)
            if (i != sameIOrank)
            {
                MPI_Send(&fileNumber, 1, MPI_INT, i, 1234, sameIOcomm);
                MPI_Send(&numberIOnodes, 1, MPI_INT, i, 1235, sameIOcomm);
            }
        delete [] otherMinRanks; otherMinRanks=NULL;
    } 
    else //non minimal nodes will just receive file number to write to
    {
        MPI_Recv(&fileNumber, 1, MPI_INT, MPI_ANY_SOURCE, 1234, sameIOcomm, MPI_STATUS_IGNORE);
        MPI_Recv(&numberIOnodes, 1, MPI_INT, MPI_ANY_SOURCE, 1235, sameIOcomm, MPI_STATUS_IGNORE);
    }
    MPI_Comm_free(&minsRanksComm);

    *sameIOcomm_ptr = sameIOcomm;
    *fileNumber_ptr = fileNumber;

    return numberIOnodes;
}

