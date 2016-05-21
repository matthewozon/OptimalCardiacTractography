#include <iostream>
#include <C_visu_IO.h>
#include <C_visu_draw.h>


using namespace std;


int main(int argc, char **argv)
{

    if(argc==6)
    {
        C_visu_IO* ioPuffin = new C_visu_IO();
        bool canDraw = false;
        if(ioPuffin->load(argv[5], M_FIBERS /*_RAW_FIBERS*/))
        {
            //select the fiber in the given range of detail
            canDraw = ioPuffin->selectFiberDetail(atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]));
        }
        else
        {
            canDraw = false;
        }



        if(canDraw)
        {
            //draw it
            C_visu_draw* drawPuffin = new C_visu_draw();
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
        cout << "function: " << argv[0] << " must be used as: " << argv[0] << " x <double> y <double> z <double> R <double> fiber_file <string>" << endl;
        return -1;
    }
    return 0;
}
