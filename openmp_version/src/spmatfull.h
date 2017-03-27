#ifndef SPMATFULL_H
#define SPMATFULL_H

#include <string>

#include "spmat.h"

using namespace std;

class SpMatFull : public SpMat
{
protected:
    int *m_m;
public:
    enum {
        EMPTY = 0,
        BLUE,
        RED
    };
    //Ctor
    SpMatFull(int R = 0, int C = 0) : SpMat(R,C){
        if (R != 0 && C != 0) {
            m_m = new int[R*C];
        } else {
            m_m = nullptr;
        }
    }
    //Dtor
    ~SpMatFull() {
        delete[] m_m;
    }

    //Utilities
    int index(int r, int c) {return r*m_nC+c;}

    //Sets & Gets
    virtual void setN0el(int i) {} //No function. This function is neeeded on in the sparse case

    //Initialization
    virtual void readInit(const string& filename);

    //Recognize:
    virtual string whoRU() {
        return "Full";
    }

    //Moves
    virtual bool moveBlue();
    virtual bool moveRed ();
    virtual bool move    (){
        bool mB = moveBlue();
        bool mR = moveRed();
        return mB || mR;
    }

    //Print
    virtual void oprint(ostream &o) {
        for (int r = 0; r < m_nR; r++) {
            for (int c = 0; c < m_nC-1; c++) { //Notice nC-1
                o << m_m[r*m_nC+c] << ",";
            }//end_for_c
            o << m_m[(r+1)*m_nC-1] << "\n";
        }//end_for_r
    }//end_OPRINT

};

#endif // SPMATFULL_H
