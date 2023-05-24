//
// Created by Dennis-Florian Herr on 13/06/2022.
//

#include <condition_variable>
#include <deque>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "Utility.h"

#define NUM_THREADS 32

std::thread threads[NUM_THREADS];

struct Problem {
  Sha1Hash sha1_hash;
  int problemNum;
};

// TO-DO: implement a thread-safe queue
// tip: use a condition variable to make threads wait when the queue is empty
class ProblemQueue {
public:
  void push(Problem problem) {
    {
      std::lock_guard<std::mutex> lock(mutex);
      problemQueue.push_back(problem);
    }
    cv.notify_one();
  }

  Problem pop() {
    std::unique_lock<std::mutex> lock(mutex);
    while (problemQueue.empty()) {
      cv.wait(lock);
    }
    Problem p = problemQueue.front();
    problemQueue.pop_front();
    return p;
  }

  bool empty() { return problemQueue.empty(); }

private:
  std::deque<Problem> problemQueue;
  std::mutex mutex;
  std::condition_variable cv;
};

ProblemQueue problemQueue;

// generate numProblems sha1 hashes with leadingZerosProblem leading zero bits
// This method is intentionally compute intense so you can already start working
// on solving problems while more problems are generated
void generateProblem(int seed, int numProblems, int leadingZerosProblem) {
  srand(seed);

  for (int i = 0; i < numProblems; i++) {
    std::string base = std::to_string(rand()) + std::to_string(rand());
    Sha1Hash hash = Utility::sha1(base);
    do {
      // we keep hashing ourself until we find the desired amount of leading
      // zeros
      hash = Utility::sha1(hash);
    } while (Utility::count_leading_zero_bits(hash) < leadingZerosProblem);
    problemQueue.push(Problem{hash, i});
  }

  for (int i = 1; i < NUM_THREADS; i++) {
    problemQueue.push(Problem{Sha1Hash(), -1});
  }
}

// This method repeatedly hashes itself until the required amount of leading
// zero bits is found
Sha1Hash findSolutionHash(Sha1Hash hash, int leadingZerosSolution) {
  do {
    // we keep hashing ourself until we find the desired amount of leading zeros
    hash = Utility::sha1(hash);
  } while (Utility::count_leading_zero_bits(hash) < leadingZerosSolution);

  return hash;
}

void parallel_work(std::vector<Sha1Hash> &solutionHashes,
                   int leadingZerosSolution) {
  while (true) {
    Problem p = problemQueue.pop();
    if (p.problemNum == -1) {
      break;
    }
    solutionHashes[p.problemNum] =
        findSolutionHash(p.sha1_hash, leadingZerosSolution);
  }
}

int main(int argc, char *argv[]) {
  int leadingZerosProblem = 8;
  int leadingZerosSolution = 11;
  int numProblems = 10000;

  // Not interesting for parallelization
  Utility::parse_input(numProblems, leadingZerosProblem, leadingZerosSolution,
                       argc, argv);
  std::vector<Sha1Hash> solutionHashes(numProblems);

  unsigned int seed = Utility::readInput();

  // TO-DO: generate problems in another thread and work on solving them while
  // generation continues
  threads[0] =
      std::thread(generateProblem, seed, numProblems, leadingZerosProblem);

  //   while (!problemQueue.empty()) {
  //     Problem p = problemQueue.pop();
  //     solutionHashes[p.problemNum] =
  //         findSolutionHash(p.sha1_hash, leadingZerosSolution);
  //   }

  for (int i = 1; i < NUM_THREADS; i++) {
    threads[i] = std::thread(parallel_work, std::ref(solutionHashes),
                             leadingZerosSolution);
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    threads[i].join();
  }

  Sha1Hash solution;
  // guarantee initial solution hash data is zero
  memset(solution.data, 0, SHA1_BYTES);
  // this doesn't need parallelization. it's neglectibly fast
  for (int i = 0; i < numProblems; i++) {
    solution = Utility::sha1(solution, solutionHashes[i]);
  }

  Utility::printHash(solution);
  printf("DONE\n");
  return 0;
}
