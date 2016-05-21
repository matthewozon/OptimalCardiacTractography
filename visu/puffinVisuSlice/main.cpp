#include <iostream>
#include <C_visu_IO.h>
#include <C_visu_draw.h>


using namespace std;


int main(int argc, char **argv)
{

    if(argc==5)
    {
        C_visu_IO* ioPuffin = new C_visu_IO();
        if(ioPuffin->load(argv[4], M_FIBERS /*_RAW_FIBERS*/,atof(argv[1]), atof(argv[2])))
        {
            //draw it
            C_visu_draw* drawPuffin = new C_visu_draw();
            drawPuffin->draw(ioPuffin->data, 1, atof(argv[3]), 1.0, 1.0, 1.0);
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
        cout << "function: " << argv[0] << " must be used as: " << argv[0] << " zmin <double> zmax <double> fiber radius <double> fiber_file <string>" << endl;
        return -1;
    }
    return 0;
}
