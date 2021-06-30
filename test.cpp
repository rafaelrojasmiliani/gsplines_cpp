#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void displaySubset(int subSet[], int size) {
  printf("One possible subset is:\n");
  for (int i = 0; i < size; i++) {
    printf("%d ", subSet[i]);
  }
  printf("\n");
}

void subsetSum(int set[], int subSet[], int n, int subSize, int total,
               int nodeCount, int sum) {
  if (total == sum) {
    displaySubset(subSet, subSize); // print the subset
    subsetSum(set, subSet, n, subSize - 1, total - set[nodeCount],
              nodeCount + 1, sum); // for other subsets
    return;
  } else {
    for (int i = nodeCount; i < n; i++) { // find node along breadth
      subSet[subSize] = set[i];
      subsetSum(set, subSet, n, subSize + 1, total + set[i], i + 1,
                sum); // do for next node in depth
    }
  }
}

void findSubset(int set[], int size, int sum) {
  int *subSet = (int *)malloc(
      sizeof(int) * size); // create subset array to pass parameter of subsetSum
  subsetSum(set, subSet, size, 0, 0, 0, sum);
  free(subSet);
}

int main() {
  clock_t t;
  t = clock();
  int size = 10;

  int weights[size];
  for (int i = 0; i < size; i++) {
    weights[i] = i + 1;
  }

  findSubset(weights, size, 10);
  t = clock() - t;
  double time_taken = ((double)t) / CLOCKS_PER_SEC; // in seconds

  printf("\nThis probelm took %f seconds to execute \n\n", time_taken);
}
