#include <C_process.h>

C_process::C_process(unsigned long NbCore)
{
    //ctor
    nbCore = NbCore;
}

C_process::~C_process()
{
    //dtor
}

void C_process::launch(unsigned long arg, void** myArgs)
{
    if(arg>nbCore)
    {
        //set a pid storage
        pid_t pids[nbCore];

        //init
        for (unsigned long i = 0; i < nbCore; ++i) {
            if ((pids[i] = fork()) < 0) {
                perror("fork");
                abort();
            } else if (pids[i] == 0)
            {
                run(myArgs[i]);
                exit(0);
            }
        }

        //loop
        int status;
        pid_t pid;
        unsigned long idx = 0;
        bool cond;
        for(unsigned long k=nbCore ; k<arg ; k++)
        {
            pid = wait(&status);
            printf("Child with PID %ld exited with status 0x%x.\n", (long)pid, status);

            //find idx of pid in pids
            cond = false;
            for(unsigned long g=0 ; g<nbCore ; g++)
            {
                if(pids[g]==pid)
                {
                    cond = true;
                    idx = g;
                    g=nbCore;
                }
            }
            if(!cond)
            {
                return;
            }

            //create new process and store pid at the right place
            if((pids[idx] = fork()) < 0)
            {
                perror("fork");
                abort();
            } else if (pids[idx] == 0)
            {
                run(myArgs[k]);
                exit(0);
            }
        }

        /* Wait for children to exit. */
        unsigned long n = nbCore;
        while (n > 0) {
          pid = wait(&status);
          printf("Child with PID %ld exited with status 0x%x.\n", (long)pid, status);
          --n;
        }
    }
    else
    {
        //set a pid storage
        pid_t pids[arg];

        //init
        for (unsigned long i = 0; i < arg; ++i) {
            if ((pids[i] = fork()) < 0) {
                perror("fork");
                abort();
            } else if (pids[i] == 0)
            {
                run(myArgs[i]);
                exit(0);
            }
        }

        int status;
        pid_t pid;

        /* Wait for children to exit. */
        unsigned long n = arg;
        while (n > 0) {
          pid = wait(&status);
          printf("Child with PID %ld exited with status 0x%x.\n", (long)pid, status);
          --n;
        }
    }

    return;
}
