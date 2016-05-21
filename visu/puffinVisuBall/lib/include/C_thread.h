#ifndef C_THREAD_H
#define C_THREAD_H
//http://tommed.tumblr.com/post/90393657/multithreading-in-c-using-posix-threads-tutorial
//http://en.wikipedia.org/wiki/POSIX_Threads

#include <pthread.h>

#ifdef ERROR_HANDLING_THREAD
    #include <errno.h> /**errno*/
    #include <stdlib.h> /**EXIT_FAILURE*/
    #include <stdio.h> /**perror*/
#endif

#ifdef VERBOSE_THREAD
    #include <iostream>
    using namespace std;
#endif

#ifdef ERROR_HANDLING_THREAD
    #define handle_error_en(en, msg) \
            do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)
#endif

class C_thread
{
    public:
        /** Default constructor */
        C_thread();
        /** Default destructor */
        virtual ~C_thread();

        //methods that implement theads.
        void start(); //start a thread that run the method execute
        virtual void execute(); //method to be overload: the working method
        bool wait_for_exit(); //method that waits till the end of the thread launch by start
        //bool waitForAllThreadToEndUp(threadStruct* THREAD);
        //threadStruct* launchThreads(unsigned long nbThread, voidStruct* t);
    private:
        pthread_t handle; //thread handler: created by start and used by wait_for_exit (no use to be in public scope)
        pthread_attr_t *attrp; //thread attribute
};










#endif // C_THREAD_H
