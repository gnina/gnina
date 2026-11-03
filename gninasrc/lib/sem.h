// Counting semaphore — portable implementation using std::mutex +
// std::condition_variable. POSIX sem_init is deprecated and broken on macOS
// (returns ENOSYS), so we avoid it entirely.
#pragma once
#include <mutex>
#include <condition_variable>

struct sem {
  sem() : count(0) {}

  void wait() {
    std::unique_lock<std::mutex> lk(mtx);
    cv.wait(lk, [this]{ return count > 0; });
    --count;
  }

  void signal() {
    std::unique_lock<std::mutex> lk(mtx);
    ++count;
    cv.notify_one();
  }

  int value() {
    std::unique_lock<std::mutex> lk(mtx);
    return count;
  }

private:
  std::mutex mtx;
  std::condition_variable cv;
  int count;
};
