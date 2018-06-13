#include <semaphore.h>

struct sem {
    sem();
    ~sem();

    void wait();
    void signal();

  private:
    sem_t pthread_sem;
};

sem::sem() {
  sem_init(&pthread_sem, 0, 0);
}

sem::~sem() {
  sem_destroy(&pthread_sem);
}

void sem::wait() {
  int ret = sem_wait(&pthread_sem);
  // We don't expect to be interrupted by signals.
  assert(ret == 0);
}

void sem::signal() {
  int ret = sem_post(&pthread_sem);
  // We don't expect to overflow the signal count.
  assert(ret == 0);
}
