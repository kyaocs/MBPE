//
//  Timer.h
//  maximal_balanced_kplex_enum
//
//  Created by kai on 2023/3/21.
//

#ifndef TIMER_H_
#define TIMER_H_

#include <cstdlib>
#include <sys/time.h>

class Timer {
public:
    Timer() { m_start = timestamp(); }
    void restart() { m_start = timestamp(); }
    long long elapsed() { return timestamp() - m_start; }
    
private:
    long long m_start;
    
    // Returns a timestamp ('now') in microseconds
    long long timestamp() {
        struct timeval tp;
        gettimeofday(&tp, nullptr);
        return ((long long)(tp.tv_sec))*1000000 + tp.tv_usec;
    }
};

#endif /* TIMER_H_ */


