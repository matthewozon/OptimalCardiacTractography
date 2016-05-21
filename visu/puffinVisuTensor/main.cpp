#include <iostream>
#include <C_visu_IO.h>
#include <C_visu_draw.h>


using namespace std;


int main(int argc, char **argv)
{

    if(argc==5)
    {
        C_visu_IO* ioPuffin = new C_visu_IO();
        ioPuffin->m_FAth = atof(argv[1]);
        ioPuffin->Zmin  = atof(argv[2]);
        ioPuffin->Zmax  = atof(argv[3]);
        bool canDraw = false;
        if(ioPuffin->load(argv[4], M_TENSOR_ANA))
        {
            canDraw = true;
        }
        else
        {
            canDraw = false;
        }



        if(canDraw)
        {
            //draw it
            C_visu_draw* drawPuffin = new C_visu_draw();
            drawPuffin->m_FAth = atof(argv[1]);
            drawPuffin->draw(ioPuffin->data, 1, 0.6, 1.0, 1.0, 1.0);
            //cout << "I am here!" << endl;
            delete drawPuffin;
        }
        else
        {
            cout << "nothing can be drawn... together :-)" << endl;
        }
        //cout << "I am there!" << endl;
        delete ioPuffin;
    }
    else
    {
        cout << "command " << argv[0] << " should be used as: " << argv[0] << " FA threshold <double> Zmin <double> Zmax <double>" << endl;
        return -1;
    }
    return 0;
}
