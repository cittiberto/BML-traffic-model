#ifndef SPMATVWBBS_H
#define SPMATVWBBS_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "spmat.h"
#include "vectbeg.h"

class SpMatVwBBS : public SpMat
{
private:
    vector<VectBeg> m_bl;
    vector<VectBeg> m_rd;
    int             n0el;   //Estimation of the #Not-0 element on a row/col
public:
    //Ctor
    SpMatVwBBS(int R = 0, int C = 0, int n0el = -1) : SpMat(R,C), m_bl(), m_rd(), n0el(n0el){
        if (C!=0)     m_bl.reserve(C);
        if (R!=0)     m_rd.reserve(R);
        else if(C!=0) m_rd.reserve(C);
    }

    //Initialization
    virtual void readInit(const string& filename);

    //Sets
    void setN0el(int n) {n0el = n;}

    //Recognize
    virtual string whoRU() {
        return "VwBBS";
    }
    //Overloads
    bool isThere(bool color,int r, int c) {
        if (color)  //blue
            return m_bl[c].isThere(r);
        else
            return m_rd[r].isThere(c);
    }
    bool isThereBS(bool color,int r, int c) {
        if (color)  //blue
            return m_bl[c].isThereBS(r);
        else
            return m_rd[r].isThereBS(c);
    }

    //Moves
            bool basicMove(bool color);
    virtual bool moveBlue()             {return basicMove(1);}
    virtual bool moveRed ()             {return basicMove(0);}
    virtual bool move    (){
        bool mB = basicMove(1);     //moveBlue
        bool mR = basicMove(0);     //moveRed
        return (mB || mR);
    }

            void oprint_norm(ostream &o);
            void oprint_fast(ostream &o);
    virtual void oprint     (ostream &o) {

                if (m_empty) {
                    for (int r = 0; r < m_nR; r++) {
                        for (int c = 0; c < m_nC-1; c++)
                            o << "0,";
                        o << "0\n";
                    }//end_for_R
                } else {

                    ///////
                    bool fastprint = 1;
                    ///////

                    if (fastprint)  oprint_fast(o);
                    else            oprint_norm(o);
                }


    }//end_OPRINT
};

#endif // SPMATVWBBS_H
