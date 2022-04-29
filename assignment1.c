#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#define B 2
#define EPS 1e-6

void allocateMatrixMemory();

void allocateCpairArrayMemory();

void insertValueToMatrix();

void insertValueToCpairArray();

void resetAllCpairResult();

double calculateTime(struct timeval start_time, struct timeval end_time);

void sequenceAlgorithm();

void baselineAlgorithm();

void *baselineCalculator(void *threadarg);

void efficientAlgorithm();

void *efficientCalculator(void *threadarg);

void storeToEfficientResult();

int checkCorrectness(double *CpairResult);


struct Cpair {
    int start;
    int end;
};

struct baselineThreadData {
    int id;
    int workLoad;
    int startIndex;
};

struct efficientThreadData {
    int id;
    int workLoadOfCpair;
    int numberOfDiagonalBlock;
    int numberOfNormalBlock;
    int numberOfRemainCpair;
    int *diagonalArray;
    int *normalArray;
    int *remainArray;
};

int N, M, T;
int state = 0;
int numberOfCpair;
double **matrix;
struct Cpair *cpairArray;
int *diagonalBlock;
int *normal;
int *remain;
double **efficientResultMatrix;
double *sequenceCpairResult;
double *baselineCpairResult;
double *efficientCpairResult;


int main(int argc, char *argv[]) {
    if (argc == 4 || argc == 5) {
        N = atoi(argv[1]);
        M = atoi(argv[2]);
        T = atoi(argv[3]);
        numberOfCpair = (N * (N + 1) / 2);

        insertValueToMatrix();
        insertValueToCpairArray();

        sequenceCpairResult = (double *) malloc(numberOfCpair * sizeof(double));
        baselineCpairResult = (double *) malloc(numberOfCpair * sizeof(double));
        efficientCpairResult = (double *) malloc(numberOfCpair * sizeof(double));
        resetAllCpairResult();
        if (argc == 4) {
            printf("N= %d; M= %d; threads= %d; state= default\n\n", N, M, T);
        } else {
            state = atoi(argv[4]);
//            printf("N= %d; M= %d; threads= %d; state= chart\n\n", N, M, T);
        }
    } else {
        printf("Usage: %s N M T\n\n"
               " N: matrix N length\n"
               " M: matrix M length\n"
               " T: number of thread\n", argv[0]);
        return -2;
    }
    if (state == 0) {
        sequenceAlgorithm();
    }
    baselineAlgorithm();

    efficientAlgorithm();

    pthread_exit(NULL);

}

void allocateMatrixMemory() {
    double *matrixSpace;
    double *testMatrixSpace;
    matrix = (double **) malloc(N * sizeof(double *));
    matrixSpace = (double *) malloc(N * M * sizeof(double));
    for (int i = 0; i < N; ++i) {
        matrix[i] = &(matrixSpace[i * M]);
    }

    efficientResultMatrix = (double **) malloc(N * sizeof(double *));
    testMatrixSpace = (double *) malloc(N * N * sizeof(double));
    for (int i = 0; i < N; ++i) {
        efficientResultMatrix[i] = &(testMatrixSpace[i * N]);
    }

}

void allocateCpairArrayMemory() {
    int blockSide = (N / B);
    int remainSide = (N % B);
    int numberOfNormalBlock = (blockSide - 1) * (blockSide) / 2;
    int numberOfRemainCpair = 0;

    if (remainSide != 0) {
        for (int i = 0; i < B - 1; ++i) {
            numberOfRemainCpair += N - (B - 2);
        }
    }
//    workloadRecord = malloc(T * sizeof(int));
    cpairArray = malloc(numberOfCpair * sizeof(struct Cpair));
    diagonalBlock = (int *) malloc(blockSide * sizeof(int));
//    diagonalR = (int *) malloc(N * sizeof(int));
    normal = (int *) malloc(numberOfNormalBlock * sizeof(int));
    if (numberOfRemainCpair > 0) {
        remain = (int *) malloc(numberOfRemainCpair * sizeof(int));
    } else {
        remain = NULL;
    }

}

void insertValueToMatrix() {
    allocateMatrixMemory();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            matrix[i][j] = (double) rand() / RAND_MAX;
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            efficientResultMatrix[i][j] = 0.0;
        }
    }
    ////check matrix
//    for (int i = 0; i < N; ++i) {
//        for (int j = 0; j < M; ++j) {
//            printf("%f \t",matrix[i][j]);
//        }
//        printf(" \n");
//    }

//    for (int i = 0; i < N; ++i) {
//        for (int j = 0; j < M; ++j) {
//            printf("%f \t",efficientResultMatrix[i][j]);
//        }
//        printf(" \n");
//    }
}

void insertValueToCpairArray() {
    allocateCpairArrayMemory();
//    int diagonalRIndex = 0;
    int diagonalBlockIndex = 0;
    int normalIndex = 0;
    int remainIndex = 0;
    for (int index = 0; index < numberOfCpair; ++index) {
        for (int i = 0; i < N; ++i) {
            for (int j = i; j < N; ++j) {
                cpairArray[index].start = i;
                cpairArray[index].end = j;
                if (j < N - B + 1) {
                    if (i == j && i % B == 0) {
                        diagonalBlock[diagonalBlockIndex] = index;
                        diagonalBlockIndex++;
                    } else if (i % B == 0 && j % B == 0) {
                        normal[normalIndex] = index;
                        normalIndex++;
                    }
                } else {
                    if (remain != NULL) {
                        remain[remainIndex] = index;
                        remainIndex++;
                    }
                }
                index++;
            };
        }
//// check CpairArray diagonalBlock normal and remain
//        printf("Cpair:\n");
//        for (int i = 0; i < numberOfCpair; ++i) {
//            printf("C%d %d\t", cpairArray[i].start, cpairArray[i].end);
//        }
//        printf(" \n");
//        printf("diagonalBlock:\n");
//        for (int i = 0; i < diagonalBlockIndex; ++i) {
//            printf("C%d %d\t", cpairArray[diagonalBlock[i]].start, cpairArray[diagonalBlock[i]].end);
//        }
//        printf(" \n");
//        printf("normal:\n");
//        for (int i = 0; i < normalIndex; ++i) {
//            printf("C%d %d\t", cpairArray[normal[i]].start, cpairArray[normal[i]].end);
//        }
//        printf(" \n");
//        printf("remain:\n");
//        for (int i = 0; i < remainIndex; ++i) {
//            printf("C%d %d\t", cpairArray[remain[i]].start, cpairArray[remain[i]].end);
//        }
    }

    ////check  cpairArray; format: start end
//    for (int i = 0; i < numberOfCpair; ++i) {
//        printf("C%d %d\t", cpairArray[i].start, cpairArray[i].end);
//    }
}

void resetAllCpairResult() {
    for (int i = 0; i < numberOfCpair; ++i) {
        sequenceCpairResult[i] = 0.0;
    }

    for (int i = 0; i < numberOfCpair; ++i) {
        baselineCpairResult[i] = 0.0;
    }

    for (int i = 0; i < numberOfCpair; ++i) {
        efficientCpairResult[i] = 0.0;
    }

    ////check sequenceCpairResult
//    for (int i = 0; i < numberOfCpair; ++i) {
//        printf("%d \t",sequenceCpairResult[i]);
//    }

    ////check baselineCpairResult
//    for (int i = 0; i < numberOfCpair; ++i) {
//        printf("%d \t",baselineCpairResult[i]);
//    }

    ////check pairEfficientResult
//    for (int i = 0; i < numberOfCpair; ++i) {
//        printf("%d \t",efficientCpairResult[i]);
//    }
}

double calculateTime(struct timeval start_time, struct timeval end_time) {
    long second = end_time.tv_sec - start_time.tv_sec;
    long microseconds = end_time.tv_usec - start_time.tv_usec;
    return second + 1e-6 * microseconds;
}

//// codes of Sequence Algorithm
void sequenceAlgorithm() {
    if (state == 0) {
        printf("\n-------------------------------------------\n");
        printf("%30s\n", "Sequence Algorithm");
        printf("-------------------------------------------\n");
    }
    struct timeval start_time, end_time;
    gettimeofday(&start_time, 0);

    for (int i = 0; i < numberOfCpair; ++i) {
        int start = cpairArray[i].start;
        int end = cpairArray[i].end;
        double sum = 0.0;
        for (int j = 0; j < M; ++j) {
            sum += matrix[start][j] * matrix[end][j];
        }
        sequenceCpairResult[i] = sum;
    }
    gettimeofday(&end_time, 0);
    double elapsed = calculateTime(start_time, end_time);
    if (state == 0) {
        printf("Sequence algorithm took %f second to complete.\n\n", elapsed);
    }
    ////check sequenceCpairResult
//    for (int i = 0; i < numberOfCpair; ++i) {
//
//        if (i % 7 == 0) {
//            printf("\n");
//        }
//        printf("%f \t", sequenceCpairResult[i]);
//
//    }
//    printf("\n");

}

//// codes of Baseline Algorithm
void baselineAlgorithm() {
    if (state == 0) {
        printf("\n-------------------------------------------\n");
        printf("%30s\n", "Baseline Algorithm");
        printf("-------------------------------------------\n");
    }
    int remainWork = numberOfCpair % T;
    int res;
    struct timeval start_time, end_time;
//start time
    gettimeofday(&start_time, 0);

    //create thread
    pthread_t pthread[T];
    struct baselineThreadData threadDataArray[T];

    int cur = 0;
    for (int i = 0; i < T; ++i) {
        threadDataArray[i].id = i;

        if (remainWork > 0) {
            threadDataArray[i].workLoad = numberOfCpair / T + 1;
            remainWork--;
        } else {
            threadDataArray[i].workLoad = numberOfCpair / T;
        }

        threadDataArray[i].startIndex = cur;
//        threadDataArray[i].baselineCpairResult = baselineCpairResult;

        cur += threadDataArray[i].workLoad;

        if (state == 0) {
            printf("Creating thread %d\n", i);
        }
        res = pthread_create(&pthread[i], NULL, baselineCalculator, &threadDataArray[i]);
        if (res) {
            printf("Pthread_create() error code: %d\n", res);
            exit(-1);
        }
    }

    for (int i = 0; i < T; ++i) {
        res = pthread_join(pthread[i], NULL);
        if (res) {
            printf("Pthread_join() error code: %d\n", res);
            exit(-1);
        }
    }

    gettimeofday(&end_time, 0);
    double elapsed = calculateTime(start_time, end_time);
    if (state == 0) {
        printf("\nCheck correctness: ");
        if (checkCorrectness(baselineCpairResult) == 1) {
            printf("True\n");
        } else {
            printf("False\n");
        }

        printf("All threads have been completed.\n");
        printf("Baseline algorithm took %f second to complete.\n\n", elapsed);
    } else if(state==1){
        printf("%-13f\t",elapsed);
    }


////check baselineCpairResult
//    for (int i = 0; i < numberOfCpair; ++i) {
//        if (i % 7 == 0) {
//            printf("\n");
//        }
//
//        printf("%f\t", baselineCpairResult[i]);
//    }


}

void *baselineCalculator(void *threadarg) {
    int id;
    int workLoad;
    int startIndex;
    struct baselineThreadData *cur_thread;


    if (state == 0) {
        sleep(1);
    }

    cur_thread = (struct baselineThreadData *) threadarg;
    id = cur_thread->id;
    workLoad = cur_thread->workLoad;
    startIndex = cur_thread->startIndex;

    for (int i = startIndex; i < startIndex + workLoad; ++i) {
        double sum = 0.0;
        for (int j = 0; j < M; ++j) {
            sum += matrix[cpairArray[i].start][j] * matrix[cpairArray[i].end][j];
        }
        baselineCpairResult[i] = sum;
    }
    if (state == 0) {
        printf("Thread %d completed %d workloads\n", id, workLoad);
    }
    pthread_exit(NULL);
}

//// codes of Efficient Algorithm
void efficientAlgorithm() {
    if (state == 0) {
        printf("\n-------------------------------------------\n");
        printf("%30s\n", "Efficient Algorithm");
        printf("-------------------------------------------\n");
    }
    int res;
    struct timeval start_time, end_time;
    int blockSide = (N / B);
    int remainSide = (N % B);
    int numberOfDiagonalBlock = blockSide;
    int numberOfNormalBlock = (blockSide - 1) * (blockSide) / 2;
    int numberOfRemainCpair = 0;
    gettimeofday(&start_time, 0);

    //create thread
    pthread_t efficientPthread[T];
    struct efficientThreadData efficientThreadDataArray[T];

    int remainDiagonalBlock = numberOfDiagonalBlock % T;
    int remainNormalBlock = numberOfNormalBlock % T;

    if (remainSide != 0) {
        for (int i = 0; i < B - 1; ++i) {
            numberOfRemainCpair += N - (B - 2);
        }
    }

    int remainCpair = numberOfRemainCpair % T;

    int cursorId = 0;
    //load balance
    for (int i = 0; i < T; ++i) {
        efficientThreadDataArray[i].id = i;
        if (remainDiagonalBlock > 0) {
            efficientThreadDataArray[i].numberOfDiagonalBlock = numberOfDiagonalBlock / T + 1;
            remainDiagonalBlock--;
            cursorId = i;
        } else {
            efficientThreadDataArray[i].numberOfDiagonalBlock = numberOfDiagonalBlock / T;
        }
    }

    int cursorId2 = 0;
    for (int i = cursorId + 1; i < T + cursorId + 1; ++i) {
        if (remainNormalBlock > 0) {
            efficientThreadDataArray[i % T].numberOfNormalBlock = numberOfNormalBlock / T + 1;
            remainNormalBlock--;
            cursorId2 = i % T;
        } else {
            efficientThreadDataArray[i % T].numberOfNormalBlock = numberOfNormalBlock / T;
        }

    }


    for (int i = cursorId2 + 1; i < T + cursorId2 + 1; ++i) {
        if (remainCpair > 0) {
            efficientThreadDataArray[i % T].numberOfRemainCpair = numberOfRemainCpair / T + 1;
            remainCpair--;
            cursorId = i;
        } else {
            efficientThreadDataArray[i % T].numberOfRemainCpair = numberOfRemainCpair / T;
        }
    }

    int workLoadCpairCursor = 0;
    for (int i = 0; i < T; ++i) {
        efficientThreadDataArray[i].workLoadOfCpair =
                efficientThreadDataArray[i].numberOfDiagonalBlock * ((1 + B) * B / 2)
                + efficientThreadDataArray[i].numberOfNormalBlock * B * B +
                efficientThreadDataArray[i].numberOfRemainCpair;
        workLoadCpairCursor += efficientThreadDataArray[i].workLoadOfCpair;
    }

    ////check the load
//    for (int i = 0; i < T; ++i) {
//        struct efficientThreadData e=efficientThreadDataArray[i];
//
//        printf("\nThread %d:\nworkLoad: %d \n diagonalBlock: %d \n normal: %d "
//               "\n remain: %d \n",e.id,e.workLoadOfCpair,e.numberOfDiagonalBlock,e.numberOfNormalBlock,e.numberOfRemainCpair);
//    }

    for (int i = 0; i < T; ++i) {
        if (efficientThreadDataArray[i].numberOfDiagonalBlock > 0) {
            efficientThreadDataArray[i].diagonalArray = (int *) malloc(
                    efficientThreadDataArray[i].numberOfDiagonalBlock * sizeof(int));
        } else {
            efficientThreadDataArray[i].diagonalArray = NULL;
        }
        if (efficientThreadDataArray[i].numberOfNormalBlock > 0) {
            efficientThreadDataArray[i].normalArray = (int *) malloc(
                    efficientThreadDataArray[i].numberOfNormalBlock * sizeof(int));
        } else {
            efficientThreadDataArray[i].normalArray = NULL;
        }
        if (efficientThreadDataArray[i].numberOfRemainCpair > 0) {
            efficientThreadDataArray[i].remainArray = (int *) malloc(
                    efficientThreadDataArray[i].numberOfRemainCpair * sizeof(int));
        } else {
            efficientThreadDataArray[i].remainArray = NULL;
        }
    }

    ////check the distribute of remain
//    printf("remain:\n");
//    for (int i = 0; i < numberOfRemainCpair; ++i) {
//        printf("C%d%d\t",cpairArray[remain[i]].start,cpairArray[remain[i]].end );
//    }

//allocate work
    int diagonalCursor = 0;
    int normalCursor = 0;
    int remainCursor = 0;
    for (int i = 0; i < T; ++i) {
        int lengthOfDiagonal = efficientThreadDataArray[i].numberOfDiagonalBlock;
        int lengthOfnormal = efficientThreadDataArray[i].numberOfNormalBlock;
        int lengthOfRemain = efficientThreadDataArray[i].numberOfRemainCpair;
//        printf("%d\n",efficientThreadDataArray[i].numberOfDiagonalBlock);
        for (int j = 0; j < lengthOfDiagonal; ++j) {
            efficientThreadDataArray[i].diagonalArray[j] = diagonalBlock[diagonalCursor];
            diagonalCursor++;
        }

        for (int j = 0; j < lengthOfnormal; ++j) {
            efficientThreadDataArray[i].normalArray[j] = normal[normalCursor];
            normalCursor++;
        }

        if (lengthOfRemain > 0) {
            for (int j = 0; j < efficientThreadDataArray[i].numberOfRemainCpair; ++j) {
                efficientThreadDataArray[i].remainArray[j] = remain[remainCursor];
                remainCursor++;
            }
        } else {
            efficientThreadDataArray[i].remainArray = NULL;
        }
    }

    ////check the state of threads
//    printf("%d \n", numberOfCpair);
//    for (int i = 0; i < T; ++i) {
//        printf("\n\nThread %d: %d\n", efficientThreadDataArray[i].id, efficientThreadDataArray[i].workLoadOfCpair);
//        printf("%d %d %d\n", efficientThreadDataArray[i].numberOfDiagonalBlock,
//               efficientThreadDataArray[i].numberOfNormalBlock, efficientThreadDataArray[i].numberOfRemainCpair);
//        printf("diagonalBlock: %d\n", efficientThreadDataArray[i].numberOfDiagonalBlock);
//        for (int j = 0; j < efficientThreadDataArray[i].numberOfDiagonalBlock; ++j) {
//            int num = efficientThreadDataArray[i].diagonalArray[j];
//            printf("C%d %d\t ", cpairArray[num].start, cpairArray[num].end);
//        }
//        printf(" \n");
//        printf("normal: %d\n", efficientThreadDataArray[i].numberOfNormalBlock);
//        for (int j = 0; j < efficientThreadDataArray[i].numberOfNormalBlock; ++j) {
//            int num = efficientThreadDataArray[i].normalArray[j];
//            printf("C%d %d\t ", cpairArray[num].start, cpairArray[num].end);
//        }
//        printf(" \n");
//        printf("remain: %d\n", efficientThreadDataArray[i].numberOfRemainCpair);
//        if (remain != NULL && efficientThreadDataArray[i].remainArray != NULL) {
//            for (int j = 0; j < efficientThreadDataArray[i].numberOfRemainCpair; ++j) {
//                int num = efficientThreadDataArray[i].remainArray[j];
//                printf("C%d %d\t ", cpairArray[num].start, cpairArray[num].end);
//            }
//        }
//    }
//    printf(" \n");

    for (int i = 0; i < T; ++i) {
        if (state == 0) {
            printf("Creating thread %d\n", i);
        }
        res = pthread_create(&efficientPthread[i], NULL, efficientCalculator, (void *) &efficientThreadDataArray[i]);
        if (res) {
            printf("Pthread_create() error code: %d\n", res);
            exit(-1);
        }
    }
    for (int i = 0; i < T; i++) {
        res = pthread_join(efficientPthread[i], NULL);
        if (res) {
            printf("ERROR; return code from pthread_join() is %d\n", res);
            exit(-1);
        }
    }

    storeToEfficientResult();

    gettimeofday(&end_time, 0);
    double elapsed = calculateTime(start_time, end_time);
    if (state == 0) {
        printf("\nCheck correctness: ");
        if (checkCorrectness(baselineCpairResult) == 1) {
            printf("True\n");
        } else {
            printf("False\n");
        }
        printf("All threads have been completed.\n");
        printf("Efficient algorithm took %f second to complete.\n\n", elapsed);
    }else if(state==1){
        printf("%-13f\n",elapsed);
    }

    ////check efficientCpairResult
//    for (int i = 0; i < numberOfCpair; ++i) {
//        if (i % 7 == 0) {
//            printf("\n");
//        }
//
//        printf("%f\t", efficientCpairResult[i]);
//    }
    pthread_exit(NULL);
}


void *efficientCalculator(void *threadarg) {
    int id;
    int workLoadOfCpair;
    int numberOfDiagonalBlock;
    int numberOfNormalBlock;
    int numberOfRemainCpair;
    int *diagonalArray;
    int *normalArray;
    int *remainArray;

    struct efficientThreadData *cur_thread;
    if (state == 0) {
        sleep(1);
    }

    cur_thread = (struct efficientThreadData *) threadarg;
    id = cur_thread->id;
    workLoadOfCpair = cur_thread->workLoadOfCpair;
    numberOfDiagonalBlock = cur_thread->numberOfDiagonalBlock;
    numberOfNormalBlock = cur_thread->numberOfNormalBlock;
    numberOfRemainCpair = cur_thread->numberOfRemainCpair;
    diagonalArray = cur_thread->diagonalArray;
    normalArray = cur_thread->normalArray;
    remainArray = cur_thread->remainArray;

    ////calculate diagonalBlock
    for (int i = 0; i < numberOfDiagonalBlock; ++i) {
        struct Cpair temp = cpairArray[diagonalArray[i]];
        int startRow = temp.start;
        int endRow = temp.end;

        double sum[(B + 1) * B / 2] = {0.0, 0.0, 0.0};
        int matrixRowCursor = 0;
        for (int matrixRowIndex = 0; matrixRowIndex < M / 5; matrixRowIndex += 5) {
            sum[0] +=
                    matrix[startRow][matrixRowIndex] * matrix[endRow][matrixRowIndex] +
                    matrix[startRow][matrixRowIndex + 1] * matrix[endRow][matrixRowIndex + 1] +
                    matrix[startRow][matrixRowIndex + 2] * matrix[endRow][matrixRowIndex + 2] +
                    matrix[startRow][matrixRowIndex + 3] * matrix[endRow][matrixRowIndex + 3] +
                    matrix[startRow][matrixRowIndex + 4] * matrix[endRow][matrixRowIndex + 4];
            sum[1] +=
                    matrix[startRow][matrixRowIndex] * matrix[endRow + 1][matrixRowIndex] +
                    matrix[startRow][matrixRowIndex + 1] * matrix[endRow + 1][matrixRowIndex + 1] +
                    matrix[startRow][matrixRowIndex + 2] * matrix[endRow + 1][matrixRowIndex + 2] +
                    matrix[startRow][matrixRowIndex + 3] * matrix[endRow + 1][matrixRowIndex + 3] +
                    matrix[startRow][matrixRowIndex + 4] * matrix[endRow + 1][matrixRowIndex + 4];
            sum[2] +=
                    matrix[startRow + 1][matrixRowIndex] * matrix[endRow + 1][matrixRowIndex] +
                    matrix[startRow + 1][matrixRowIndex + 1] * matrix[endRow + 1][matrixRowIndex + 1] +
                    matrix[startRow + 1][matrixRowIndex + 2] * matrix[endRow + 1][matrixRowIndex + 2] +
                    matrix[startRow + 1][matrixRowIndex + 3] * matrix[endRow + 1][matrixRowIndex + 3] +
                    matrix[startRow + 1][matrixRowIndex + 4] * matrix[endRow + 1][matrixRowIndex + 4];
            matrixRowCursor += 5;
        }
        for (int matrixRowIndex = matrixRowCursor;
             matrixRowIndex < matrixRowCursor + M % 5; matrixRowIndex++) {
            sum[0] += matrix[startRow][matrixRowIndex] * matrix[endRow][matrixRowIndex];
            sum[1] += matrix[startRow][matrixRowIndex] * matrix[endRow + 1][matrixRowIndex];
            sum[2] += matrix[startRow + 1][matrixRowIndex] * matrix[endRow + 1][matrixRowIndex];
        }
        efficientResultMatrix[startRow][endRow] = sum[0];
        efficientResultMatrix[startRow][endRow + 1] = sum[1];
        efficientResultMatrix[startRow + 1][endRow + 1] = sum[2];
    }


    ////calculate normal
    for (int i = 0; i < numberOfNormalBlock; ++i) {
        struct Cpair temp = cpairArray[normalArray[i]];
        int startRow = temp.start;
        int endRow = temp.end;

        double sum[2][2] = {{0.0, 0.0},
                            {0.0, 0.0}};
        int matrixRowCursor = 0;
        for (int matrixRowIndex = 0; matrixRowIndex < M / 5; matrixRowIndex += 5) {
            sum[0][0] += matrix[startRow][matrixRowIndex] * matrix[endRow][matrixRowIndex] +
                         matrix[startRow][matrixRowIndex + 1] * matrix[endRow][matrixRowIndex + 1] +
                         matrix[startRow][matrixRowIndex + 2] * matrix[endRow][matrixRowIndex + 2] +
                         matrix[startRow][matrixRowIndex + 3] * matrix[endRow][matrixRowIndex + 3] +
                         matrix[startRow][matrixRowIndex + 4] * matrix[endRow][matrixRowIndex + 4];

            sum[0][1] += matrix[startRow][matrixRowIndex] * matrix[endRow + 1][matrixRowIndex] +
                         matrix[startRow][matrixRowIndex + 1] * matrix[endRow + 1][matrixRowIndex + 1] +
                         matrix[startRow][matrixRowIndex + 2] * matrix[endRow + 1][matrixRowIndex + 2] +
                         matrix[startRow][matrixRowIndex + 3] * matrix[endRow + 1][matrixRowIndex + 3] +
                         matrix[startRow][matrixRowIndex + 4] * matrix[endRow + 1][matrixRowIndex + 4];

            sum[1][0] += matrix[startRow + 1][matrixRowIndex] * matrix[endRow][matrixRowIndex] +
                         matrix[startRow + 1][matrixRowIndex + 1] * matrix[endRow][matrixRowIndex + 1] +
                         matrix[startRow + 1][matrixRowIndex + 2] * matrix[endRow][matrixRowIndex + 2] +
                         matrix[startRow + 1][matrixRowIndex + 3] * matrix[endRow][matrixRowIndex + 3] +
                         matrix[startRow + 1][matrixRowIndex + 4] * matrix[endRow][matrixRowIndex + 4];

            sum[1][1] += matrix[startRow + 1][matrixRowIndex] * matrix[endRow + 1][matrixRowIndex] +
                         matrix[startRow + 1][matrixRowIndex + 1] * matrix[endRow + 1][matrixRowIndex + 1] +
                         matrix[startRow + 1][matrixRowIndex + 2] * matrix[endRow + 1][matrixRowIndex + 2] +
                         matrix[startRow + 1][matrixRowIndex + 3] * matrix[endRow + 1][matrixRowIndex + 3] +
                         matrix[startRow + 1][matrixRowIndex + 4] * matrix[endRow + 1][matrixRowIndex + 4];

            matrixRowCursor += 5;
        }

        for (int matrixRowIndex = matrixRowCursor; matrixRowIndex < matrixRowCursor + M % 5; matrixRowIndex++) {
            sum[0][0] += matrix[startRow][matrixRowIndex] * matrix[endRow][matrixRowIndex];
            sum[0][1] += matrix[startRow][matrixRowIndex] * matrix[endRow + 1][matrixRowIndex];
            sum[1][0] += matrix[startRow + 1][matrixRowIndex] * matrix[endRow][matrixRowIndex];
            sum[1][1] += matrix[startRow + 1][matrixRowIndex] * matrix[endRow + 1][matrixRowIndex];
        }
        efficientResultMatrix[startRow][endRow] = sum[0][0];
        efficientResultMatrix[startRow][endRow + 1] = sum[0][1];
        efficientResultMatrix[startRow + 1][endRow] = sum[1][0];
        efficientResultMatrix[startRow + 1][endRow + 1] = sum[0][0];
    }

    ////calculate Cpair
    for (int i = 0; i < numberOfRemainCpair; ++i) {
        int startRow = cpairArray[remainArray[i]].start;
        int endRow = cpairArray[remainArray[i]].end;
        double sum = 0;
        int matrixRowCursor = 0;

        for (int matrixRowIndex = 0; matrixRowIndex < M / 5; matrixRowIndex += 5) {
            sum += matrix[startRow][matrixRowIndex] * matrix[endRow][matrixRowIndex] +
                   matrix[startRow][matrixRowIndex + 1] * matrix[endRow][matrixRowIndex + 1] +
                   matrix[startRow][matrixRowIndex + 2] * matrix[endRow][matrixRowIndex + 2] +
                   matrix[startRow][matrixRowIndex + 3] * matrix[endRow][matrixRowIndex + 3] +
                   matrix[startRow][matrixRowIndex + 4] * matrix[endRow][matrixRowIndex + 4];
            matrixRowCursor += 5;
        }
        for (int matrixRowIndex = matrixRowCursor; matrixRowIndex < matrixRowCursor + M % 5; matrixRowIndex++) {
            sum += matrix[startRow][matrixRowIndex] * matrix[endRow][matrixRowIndex];
        }
        efficientResultMatrix[startRow][endRow] = sum;
    }
    if (state == 0) {
        printf("Thread %d completed %d workloads\n", id, workLoadOfCpair);
    }
    pthread_exit(NULL);

}

void storeToEfficientResult() {
    ////check efficientResultMatrix
//    for (int i = 0; i < N; ++i) {
//        for (int j = 0; j < N; ++j) {
//            printf("%f\t", efficientResultMatrix[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
    ////check initial efficientCpairResult
    //    for (int i = 0; i < numberOfCpair; ++i) {
//        printf("%f\t",efficientCpairResult[i]);
//    }

    int cursor = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i <= j) {
                efficientCpairResult[cursor] = efficientResultMatrix[i][j];
                cursor++;
            }
        }
    }
}

int checkCorrectness(double *CpairResult) {
    for (int i = 0; i < numberOfCpair; ++i) {
        if (fabs(CpairResult[i] - sequenceCpairResult[i]) > EPS) {
            return 0;
        }
    }
    return 1;
}

//
// Created by lo2 on 4/3/2022.
//