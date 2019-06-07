#include <iostream>
#include <iomanip>
#include <mpi.h>
#include "timer.hpp"

/** Timer constructor expecting a count for dynamic allocation calls.
 * @param[in] count Number of timers to set up. Usually, an enumeration will be
 * queried to fetch this number (default argument).
 */
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

/** Frees the memory allocated for the timers.
 */
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

/** Starts the timer with the requested index. The timer storage array is
 * cleared at the index location.
 * @param[in] idx Index of the timer to start
 */
SError Timer::StartTimer( unsigned idx ) {
    timers_[idx] = 0.0;
    starts_[idx] = MPI_Wtime();
    return E_SUCCESS;
}

/** Overrides calls to start and stop and assigns a given value to the requested
 * timer. The timer is left stopped.
 * @param[in] idx Index of the timer.
 * @param[in] val Value to store in the timer's slot.
 */
SError Timer::SetTimer( unsigned idx, const double& val ) {
    StopTimer( idx );
    timers_[idx] = val;
    return E_SUCCESS;
}

/** Stops a running timer. The difference between the actual time and the
 * starting time is computed, and the start time is cleared.
 * @param[in] idx Index of the timer.
 */
SError Timer::StopTimer( unsigned idx ) {
    timers_[idx] = MPI_Wtime() - starts_[idx];
    starts_[idx] = 0.0;
    return E_SUCCESS;
}

/** Checks if a timer is active (if it's start time is set).
 * @param[in] idx Timer index.
 * @return True if the timer is running.
 */
bool Timer::IsActive( unsigned idx ) const {
    return timers_[idx]==0.0;
}

/** Returns the number of timers are handled.
 */
unsigned Timer::GetNumber() const {
    return counters_;
}

/** Returns the requested timer's time.
 * If the timer is running, this time is computed on the fly.
 * @param[in] idx Index of the timer.
 * @return Time elapsed since the timer started.
 */
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
