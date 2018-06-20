#ifndef __LOOP_TIMER_H
#define __LOOP_TIMER_H

#include <boost/timer/timer.hpp>

struct loop_timer : boost::timer::cpu_timer {
    /* TODO: switch to gnu++11 and use i = 0*/
    int i;

    loop_timer()
        : i(0) {
    }

    void resume(void) {
      if (!i++)
        boost::timer::cpu_timer::start();
      else
        boost::timer::cpu_timer::resume();
    }

    ~loop_timer() {
      boost::timer::cpu_times t = elapsed();
      std::cout << "Loops:" << i << " wall:" << t.wall << " avg:" << t.wall / i
          << "\n";
    }
};

#endif
