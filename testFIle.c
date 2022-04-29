//#include <stdlib.h>
//#include <stdio.h>
//#include <pthread.h>
//#include <unistd.h>
//#include <time.h>
//
//
//void setB();
//
//void allocateMatrixMemory();
//
//void allocateCpairArrayMemory();
//
//void insertValueToMatrix();
//
//void insertValueToCpairArray();
//
//void resetAllCpairResult(double *sequenceCpairResult, double *baselineCpairResult, double *cpairEfficientResult);
//
//double calculateTime(struct timeval start_time, struct timeval end_time);
//
//void sequenceAlgorithm(double *sequenceCpairResult);
//
//void baselineAlgorithm(double *baselineCpairResult);
//
//void *baselineCalculator(void *threadarg);
//
//void efficientAlgorithm(double *efficientCpairResult);
//
//void *efficientCalculator(void *threadarg);
//
//
//struct Cpair {
//    int start;
//    int end;
//};
//
//struct baselineThreadData {
//    int id;
//    int workLoad;
//    int startIndex;
//    double *baselineCpairResult;
//};
//
//struct efficientThreadData {
//    int id;
//    int workLoadOfCpair;
//    int numberOfDiagonalBlock;
//    int numberOfNormalBlock;
//    int numberOfRemainCpair;
//    double *efficientCpairResult;
//};
//
//int row, column, T;
//int B;
//int numberOfCpair;
//double **matrix;
//struct Cpair *cpairArray;
//int *diagonal;
//
//int main(int argc, char *argv[]) {
//    double *sequenceCpairResult;
//    double *baselineCpairResult;
//    double *efficientCpairResult;
//
//    if (argc == 4) {
//        row = atoi(argv[1]);
//        column = atoi(argv[2]);
//        T = atoi(argv[3]);
//        numberOfCpair = (row * (row + 1) / 2);
//        setB();
//
//        diagonal = (int *) malloc(row * sizeof(int));
//        insertValueToMatrix();
//        insertValueToCpairArray();
//
//        sequenceCpairResult = (double *) malloc(numberOfCpair * sizeof(double));
//        baselineCpairResult = (double *) malloc(numberOfCpair * sizeof(double));
//        efficientCpairResult = (double *) malloc(numberOfCpair * sizeof(double));
//        resetAllCpairResult(sequenceCpairResult, baselineCpairResult, efficientCpairResult);
//
//        printf("row= %d; column= %d; threads= %d \n\n", row, column, T);
//    } else {
//        printf("Usage: %s N M T\n\n"
//               " N: matrix row length\n"
//               " M: matrix column length\n"
//               " T: number of thread\n", argv[0]);
//        return -2;
//    }
//
////    sequenceAlgorithm(sequenceCpairResult);
//
////    baselineAlgorithm(baselineCpairResult);
//
//
//
//    efficientAlgorithm(efficientCpairResult);
//
//    pthread_exit(NULL);
//
//}
//
//void setB() {
////    int b=B;
////    int minRemain=row%(b*T);
////
////    while (b<T&&row/(b*T)>0){
////        int curRemain=row%(b*T);
////        if(curRemain<=minRemain){
////            minRemain=curRemain;
////            B=b;
////        }
////        b++;
////    }
//    B = 2;
//
//    ////check B
////    printf("%d\n",B);
//}
//
//void allocateMatrixMemory() {
//    double *matrixSpace;
//    matrix = (double **) malloc(row * sizeof(double *));
//    matrixSpace = (double *) malloc(row * column * sizeof(double));
//    for (int i = 0; i < row; ++i) {
//        matrix[i] = &(matrixSpace[i * column]);
//    }
//}
//
//void allocateCpairArrayMemory() {
//    cpairArray = malloc(numberOfCpair * sizeof(struct Cpair));
//}
//
//void insertValueToMatrix() {
//    allocateMatrixMemory();
//    for (int i = 0; i < row; ++i) {
//        for (int j = 0; j < column; ++j) {
//            matrix[i][j] = (double) rand() / RAND_MAX;
//        }
//    }
//
//    ////check matrix
////    for (int i = 0; i < row; ++i) {
////        for (int j = 0; j < column; ++j) {
////            printf("%f \t",matrix[i][j]);
////        }
////        printf(" \n");
////    }
//}
//
//void insertValueToCpairArray() {
//    allocateCpairArrayMemory();
//    int diagonalIndex = 0;
//
//    for (int index = 0; index < numberOfCpair; ++index) {
//        for (int i = 0; i < row; ++i) {
//            for (int j = i; j < row; ++j) {
//                cpairArray[index].start = i;
//                cpairArray[index].end = j;
//
//                if (i == j) {
//                    diagonal[diagonalIndex] = index;
//                    diagonalIndex++;
//                }
//                index++;
//            };
//        }
//    }
//
//    ////check  cpairArray; format: start end
////    for (int i = 0; i < numberOfCpair; ++i) {
////        printf("%d %d\t", cpairArray[i].start, cpairArray[i].end);
////    }
//}
//
//void resetAllCpairResult(double *sequenceCpairResult, double *baselineCpairResult, double *cpairEfficientResult) {
//    for (int i = 0; i < (row * (row + 1) / 2); ++i) {
//        sequenceCpairResult[i] = 0.0;
//    }
//
//    for (int i = 0; i < (row * (row + 1) / 2); ++i) {
//        baselineCpairResult[i] = 0.0;
//    }
//
//    for (int i = 0; i < (row * (row + 1) / 2); ++i) {
//        cpairEfficientResult[i] = 0.0;
//    }
//
//    ////check sequenceCpairResult
////    for (int i = 0; i < numberOfCpair; ++i) {
////        printf("%d \t",sequenceCpairResult[i]);
////    }
//
//    ////check baselineCpairResult
////    for (int i = 0; i < numberOfCpair; ++i) {
////        printf("%d \t",baselineCpairResult[i]);
////    }
//
//    ////check pairEfficientResult
////    for (int i = 0; i < numberOfCpair; ++i) {
////        printf("%d \t",cpairEfficientResult[i]);
////    }
//}
//
//double calculateTime(struct timeval start_time, struct timeval end_time) {
//    long second = end_time.tv_sec - start_time.tv_sec;
//    long microseconds = end_time.tv_usec - start_time.tv_usec;
//    return second + 1e-6 * microseconds;
//}
//
////// codes of Sequence Algorithm
//void sequenceAlgorithm(double *sequenceCpairResult) {
//    printf("\n-------------------------------------------\n");
//    printf("% 30s\n", "Sequence Algorithm");
//    printf("-------------------------------------------\n");
//
//    struct timeval start_time, end_time;
//    mingw_gettimeofday(&start_time, 0);
//
//    for (int i = 0; i < numberOfCpair; ++i) {
//        int start = cpairArray[i].start;
//        int end = cpairArray[i].end;
//        double sum = 0.0;
//        for (int j = 0; j < column; ++j) {
//            sum += matrix[start][j] * matrix[end][j];
//        }
//        sequenceCpairResult[i] = sum;
//    }
//    mingw_gettimeofday(&end_time, 0);
//    double elapsed = calculateTime(start_time, end_time);
//
//    printf("Sequence algorithm took %f second to complete.\n\n", elapsed);
//
//    ////check sequenceCpairResult
//
////    for (int i = 0; i < numberOfCpair; ++i) {
////        printf("%f \t", sequenceCpairResult[i]);
////    }
//
//}
//
////// codes of Baseline Algorithm
//void baselineAlgorithm(double *baselineCpairResult) {
//    printf("\n-------------------------------------------\n");
//    printf("% 30s\n", "Baseline Algorithm");
//    printf("-------------------------------------------\n");
//    int remainWork = numberOfCpair % T;
//    int res;
//    struct timeval start_time, end_time;
////start time
//    mingw_gettimeofday(&start_time, 0);
//
//    //create thread
//    pthread_t pthread[T];
//    struct baselineThreadData threadDataArray[T];
//
//    int cur = 0;
//    for (int i = 0; i < T; ++i) {
//        threadDataArray[i].id = i;
//
//        if (remainWork > 0) {
//            threadDataArray[i].workLoad = numberOfCpair / T + 1;
//            remainWork--;
//        } else {
//            threadDataArray[i].workLoad = numberOfCpair / T;
//        }
//
//        threadDataArray[i].startIndex = cur;
//        threadDataArray[i].baselineCpairResult = baselineCpairResult;
//
//        cur += threadDataArray[i].workLoad;
//
//
//        printf("Creating thread %d\n", i);
//        res = pthread_create(&pthread[i], NULL, baselineCalculator, &threadDataArray[i]);
//        if (res) {
//            printf("Pthread_create() error code: %d\n", res);
//            exit(-1);
//        }
//    }
//
//    for (int i = 0; i < T; ++i) {
//        res = pthread_join(pthread[i], NULL);
//        if (res) {
//            printf("Pthread_join() error code: %d\n", res);
//            exit(-1);
//        }
//    }
//
//    mingw_gettimeofday(&end_time, 0);
//    double elapsed = calculateTime(start_time, end_time);
//
//    printf("All threads have been completed.\n");
//    printf("Baseline algorithm took %f second to complete.\n\n", elapsed);
//
//////check baselineCpairResult
////    for (int i = 0; i < numberOfCpair; ++i) {
////        printf("%f\t",baselineCpairResult[i]);
////    }
//
//
//}
//
//void *baselineCalculator(void *threadarg) {
//    int id;
//    int workLoad;
//    int startIndex;
//    double *baselineCpairResult;
//    struct baselineThreadData *cur_thread;
//
//
//    sleep(1);
//    cur_thread = (struct baselineThreadData *) threadarg;
//    id = cur_thread->id;
//    workLoad = cur_thread->workLoad;
//    startIndex = cur_thread->startIndex;
//    baselineCpairResult = cur_thread->baselineCpairResult;
//
//    for (int i = startIndex; i < startIndex + workLoad; ++i) {
//        double sum = 0.0;
//        for (int j = 0; j < column; ++j) {
//            sum += matrix[cpairArray[i].start][j] * matrix[cpairArray[i].end][j];
//        }
//        baselineCpairResult[i] = sum;
//    }
//    printf("Thread %d completed %d workloads\n", id, workLoad);
//}
//
////// codes of Efficient Algorithm
//void efficientAlgorithm(double *efficientCpairResult) {
//    printf("\n-------------------------------------------\n");
//    printf("% 30s\n", "Efficient Algorithm");
//    printf("-------------------------------------------\n");
//    int res;
//    struct timeval start_time, end_time;
//
//    int blockSide = (row / B);
//    int remainSide = (row % B);
//    int numberOfDiagonalBlock = blockSide;
//    int numberOfNormalBlock = (blockSide - 1) * (blockSide) / 2;
//    int numberOfRemainCpair = 0;
//    int numberOfDivideBlock = blockSide * (blockSide + 1) / 2;
//    mingw_gettimeofday(&start_time, 0);
//
//    //create thread
//    pthread_t efficientPthread[T];
//    struct efficientThreadData efficientThreadDataArray[T];
//
//    int remainDiagonalBlock = numberOfDiagonalBlock % T;
//    int remainNormalBlock = numberOfNormalBlock % T;
//    int remainCpair = numberOfRemainCpair % T;
//
//    int maxWorkload = 0;
//
//    int cursor=INT_MAX;
//
//    if(remainSide!=0) {
//        for (int i = row; i > row - B + 1; --i) {
//            numberOfRemainCpair+=i;
//        }
//    }
//
//    for (int i = 0; i < T; ++i) {
//        efficientThreadDataArray[i].id = i;
//        efficientThreadDataArray[i].efficientCpairResult = efficientCpairResult;
//    }
//
//    for (int i = 0; i < T; ++i) {
//        if (remainDiagonalBlock > 0) {
//            efficientThreadDataArray[i].numberOfDiagonalBlock = numberOfDiagonalBlock / T + 1;
//            remainDiagonalBlock--;
//        } else {
//            if(cursor>i){
//                cursor=i;
//            }
//            efficientThreadDataArray[i].numberOfDiagonalBlock = numberOfDiagonalBlock / T;
//        }
//    }
//
//    for (int i = cursor; i < T+cursor; ++i) {
//        if (remainNormalBlock > 0) {
//            efficientThreadDataArray[i%T].numberOfNormalBlock = numberOfNormalBlock / T + 1;
//            remainNormalBlock--;
//        } else {
//            if(cursor>i){
//                cursor=i;
//            }
//            efficientThreadDataArray[i%T].numberOfNormalBlock = numberOfNormalBlock / T;
//        }
//
//    }
//
//    for (int i = cursor; i < T+cursor; ++i) {
//        if (remainCpair > 0) {
//            efficientThreadDataArray[i%T].numberOfRemainCpair = numberOfRemainCpair / T + 1;
//            remainNormalBlock--;
//        } else {
//            efficientThreadDataArray[i%T].numberOfRemainCpair = numberOfRemainCpair / T;
//        }
//    }
//
//    for (int i = 0; i < T; ++i) {
//        efficientThreadDataArray[i].workLoadOfCpair =
//                efficientThreadDataArray[i].numberOfDiagonalBlock * ((1 + B) * B / 2)
//                + efficientThreadDataArray[i].numberOfNormalBlock * B * B+efficientThreadDataArray[i].numberOfRemainCpair;
//    }
//
//
////    for (int i = 0; i < T; ++i) {
////        printf("%d\t%d\t%d\t%d\n" ,efficientThreadDataArray[i].workLoadOfCpair,efficientThreadDataArray[i].numberOfDiagonalBlock * ((1 + B) * B / 2),
////               efficientThreadDataArray[i].numberOfNormalBlock*B * B,efficientThreadDataArray[i].numberOfRemainCpair);
////    }
////    printf("%d",row*(row+1)/2);
//
//
//    for (int i = 0; i < T; ++i) {
//        printf("%d \n%d %d %d %d\n",i,efficientThreadDataArray[i].workLoadOfCpair,efficientThreadDataArray[i].numberOfNormalBlock,efficientThreadDataArray[i].numberOfDiagonalBlock);
//    }
//
//    for (int i = 0; i < T; ++i) {
//        printf("Creating thread %d\n", i);
//        res = pthread_create(&efficientPthread[i], NULL, efficientCalculator, &efficientThreadDataArray[i]);
//        if (res) {
//            printf("Pthread_create() error code: %d\n", res);
//            exit(-1);
//        }
//    }
//
////    for (int i = 0; i < T; ++i) {
////        if (remainCpair > 0) {
////            efficientThreadDataArray[i].numberOfRemainCpair = numberOfCpair / T + 1;
////            remainCpair--;
////        } else {
////            efficientThreadDataArray[i].numberOfRemainCpair = numberOfCpair / T;
////        }
////    }
//
//
//
////        for (int i = 0; i < row; ++i) {
////            printf("%d\t", diagonal[i]);
////        }
////        printf(" \n");
////        for (int i = 0; i < numberOfCpair; ++i) {
////            printf("%d %d\t", cpairArray[i].start, cpairArray[i].end);
////        }
//}
//
//
//void *efficientCalculator(void *threadarg) {
//    int id;
//    int workLoadOfCpair;
//    int numberOfDiagonalBlock;
//    int numberOfNormalBlock;
//    int numberOfRemainBlock;
//    double *efficientCpairResult;
//
//    struct efficientThreadData *cur_thread;
//
//    sleep(1);
//    cur_thread=(struct efficientThreadData *)threadarg;
//    printf("a");
//    id=cur_thread->id;
//    printf("b");
//    workLoadOfCpair=cur_thread->workLoadOfCpair;
//
//    printf("%d %d",id ,workLoadOfCpair);
//
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
