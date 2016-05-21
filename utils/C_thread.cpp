#include <C_thread.h>

C_thread::C_thread()
{
    //ctor
    attrp = NULL;
    //handle
}

C_thread::~C_thread()
{
    //dtor
}


extern "C"
{
    // this C function will be used to receive the thread and pass it back to the Thread instance
    void* thread_catch_fun(void* arg)
    {
        C_thread* t = static_cast<C_thread*>(arg);
        t->execute();
        return 0;
    }
}

// method which starts the new thread
void C_thread::start()
{
    ///create a thread that run the c-function thread_catch_fun which run the method execute
    attrp = NULL;
    #ifndef ERROR_HANDLING_THREAD
    pthread_create(&handle, attrp/*0*/, thread_catch_fun, this);
    #else
    int s = pthread_create(&handle, attrp/*0*/, thread_catch_fun, this);
    if (s != 0)
        handle_error_en(s, "pthread_create");
    #endif
}

// code which will be run on the new thread
void C_thread::execute(/*you should not add any argument here. You will not be able to pass it properly. You can use attibute of you inherited class*/)
{
    #ifdef VERBOSE_THREAD
    cout << "Thread:Hello From a new Thread!" << endl;
    cout << "You should overload me!!! unless you're big fat lazy guy who... sorry, don't pay attetion to it, it's another thread."  << endl;
    #endif
}


// wait until this thread has finished executing
bool C_thread::wait_for_exit()
{
    //free thread attribute (no effect on running thread: http://linux.die.net/man/3/pthread_attr_destroy)
    if(attrp!=NULL)
    {
        pthread_attr_destroy(attrp);
    }

    //wait for thread to send an exit signal
    int rc = pthread_join(handle, NULL); //
    if(rc)
    {
        #ifdef VERBOSE_THREAD
        cout << "Error:unable to join," << rc << endl;
        #endif
        return false;
    }
    return true;
}

