#ifndef C_PROCESS_H
#define C_PROCESS_H

#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <sys/types.h> /* for pid_t */
#include <sys/wait.h> /* for wait */



class C_process
{
    public:
        /** Default constructor */
        C_process(unsigned long NbCore);
        /** Default destructor */
        virtual ~C_process();
        void launch(unsigned long arg, void** myArgs);
    protected:
        unsigned long nbCore;
        virtual void run(void* myArg) = 0;
    private:

};

#endif // C_PROCESS_H
