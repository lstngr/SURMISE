#include <iostream>
#include <iomanip>
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

SError Timer::SetTimer( unsigned idx, const double& val ) {
    StopTimer( idx );
    timers_[idx] = val;
    return E_SUCCESS;
}

SError Timer::StopTimer( unsigned idx ) {
    timers_[idx] = MPI_Wtime() - starts_[idx];
    starts_[idx] = 0.0;
    return E_SUCCESS;
}

bool Timer::IsActive( unsigned idx ) const {
    return timers_[idx]==0.0;
}

unsigned Timer::GetNumber() const {
    return counters_;
}

double Timer::GetTime( unsigned idx ) const {
    if( not IsActive(idx) )
        return timers_[idx];
    return MPI_Wtime()-starts_[idx];
}

std::ostream& operator<<( std::ostream& os, const Timer& t ) {
    int os_prec( os.precision() );
    os << "Timers (in seconds)" << std::endl << std::setprecision(8)
        << std::setw(25) << "Iteration: "
        << std::setw(10) << t.GetTime( T_ITER ) << std::endl
        << std::setw(25) << "Tree Distribution: "
        << std::setw(10) << t.GetTime( T_DISTR ) << std::endl
        << std::setw(25) << "Force Computation: "
        << std::setw(10) << t.GetTime( T_FORCES ) << std::endl
        << std::setw(25) << "Time Evolution: "
        << std::setw(10) << t.GetTime( T_EVOL ) << std::endl
        << std::setw(25) << "Leaf Synchronisation: "
        << std::setw(10) << t.GetTime( T_LEAFSYNC ) << std::endl
        << std::setw(25) << "Tree Update: "
        << std::setw(10) << t.GetTime( T_TREEUPDATE ) << std::endl
        << std::setw(25) << "File Outputs: "
        << std::setw(10) << t.GetTime( T_OUTPUT ) << std::endl;
    os << std::setprecision( os_prec );
    return os;
}
