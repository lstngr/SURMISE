#include <iostream>
#include <mpi.h>
#include "timer.hpp"

Timer::Timer( unsigned count )
    :counters_(count), timers_(NULL), starts_(NULL)
{
    if( not count ){
        std::cout << "Timer is fucked up." << std::endl;
        return;
    }
    timers_ = new double[counters_];
    starts_ = new double[counters_];
}

Timer::~Timer() {
    if( timers_!=NULL ) {
        delete[] timers_;
        timers_=NULL;
    }
    if( starts_ != NULL ) {
        delete[] starts_;
        starts_ = NULL;
    }
}

SError Timer::StartTimer( unsigned idx ) {
    timers_[idx] = 0.0;
    starts_[idx] = MPI_Wtime();
    return E_SUCCESS;
}

SError Timer::StopTimer( unsigned idx ) {
    timers_[idx] = MPI_Wtime() - starts_[idx];
    starts_[idx] = 0.0;
    return E_SUCCESS;
}

bool Timer::IsActive( unsigned idx ) const {
    return starts_[idx]==0.0;
}

unsigned Timer::GetNumber() const {
    return counters_;
}

double Timer::GetTime( unsigned idx ) const {
    if( not IsActive(idx) )
        return timers_[idx];
    return MPI_Wtime()-starts_[idx];
}
