#include <iostream>
#include <C_visu_IO.h>
#include <C_visu_draw.h>
#include <sstream>


using namespace std;


int main(int argc, char **argv)
{

    if(argc==5)
    {
        C_visu_IO* ioPuffin = new C_visu_IO();
        bool canDraw = false;
        if(ioPuffin->load(argv[4], M_FIBERS /*M_RAW_FIBERS*/))
        {
            //select the fiber in the given range of detail
            if(false)
            {
                canDraw = ioPuffin->selectFiberDetail(atof(argv[1]), atof(argv[2]));
            }
            else
            {
                canDraw = ioPuffin->selectFiberDetailConvexHullTotalPerimeter(atof(argv[1]), atof(argv[2]));
                std::ostringstream buffer;
                buffer << "savedFiber_min_"<< argv[1] << "_max_" << argv[2];
                ioPuffin->saveFiber(buffer.str());
            }
        }
        else
        {
            canDraw = false;
        }



        if(canDraw)
        {
            //draw it
            C_visu_draw* drawPuffin = new C_visu_draw();
            drawPuffin->draw(ioPuffin->data, 1, atof(argv[3])/*0.35*0.6*/, 1.0, 1.0, 1.0);
            //ioPuffin->saveFiber("fiber_threshold_bulk");
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
        cout << "function: " << argv[0] << " must be used as: " << argv[0] << " detail min <double> detail max <double> fiber radius<double> fiber_file <string>" << endl;
        return -1;
    }
    return 0;
}
