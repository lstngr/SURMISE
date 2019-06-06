/**
 * @file timer_hpp_ I time you
 */
#ifndef SURMISE_TIMER_HPP_
#define SURMISE_TIMER_HPP_

#include "serrors.hpp"
#include "iomanager.hpp"

enum STimerNames {
    T_ITER=0,
    T_DISTR,
    T_FORCES,
    T_EVOL,
    T_LEAFSYNC,
    T_TREEUPDATE,
    T_OUTPUT,
    T_COUNT
};

/** @brief A class measuring the almighty time.
 */
class Timer {
    public:
        Timer( unsigned count=T_COUNT );
        ~Timer();
        SError StartTimer( unsigned idx );
        SError StopTimer(  unsigned idx );
        bool IsActive( unsigned idx ) const;
        unsigned GetNumber() const;
        double GetTime( unsigned idx ) const;
    private:
        friend class IOManager;
        unsigned counters_;
        double* timers_;
        double* starts_;
        // double* ends_; ?
};

#endif // SURMISE_TIMER_HPP_
