#include "spmatfull.h"

#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <omp.h>

void SpMatFull::readInit(const string& filename) {
    int r(0), c(0);
    ifstream    file(filename);
    string      tmpline;
//    setm(m_nR,m_nC);
    if (file) {
        string       anelem;
        stringstream aline;
        getline(file, tmpline);  //First row is already been read
        tmpline.clear();
        while ( getline(file, tmpline) ) {
            //Parse the document and store the line
            c = 0;
            aline << tmpline;
                while ( getline(aline, anelem, ',') ) {
                    //Parse the line
                    if (anelem.compare("0") == 0) {             //Empty location
                        m_full = false;
                        m_m[r*m_nC+c] =     EMPTY;
                    } else {
                        m_empty = false;
                        if (anelem.compare("1") == 0){          //Blue block
                            m_allR = false;
                            m_m[r*m_nC+c] = BLUE;
                        } else if (anelem.compare("2") == 0) {  //Red block
                            m_allB = false;
                            m_m[r*m_nC+c] = RED;
                        } else {                                //Not an allowed format
                            //Throw an exception
                            throw("Not allowed char");
                        }
                    }//end_if_NOT_empty
                    c++;
                }//end_while_parseALine
                if (c!=0 && c != m_nC) {                                //Different number of col in this row
                    //Throw an exception
                    throw("Not same number of col");
                }
            aline.clear();
            r++;
            tmpline.clear();
        }//end_while_getline
        file.close();
    } else throw("File not found");    //end_if_file
}//end_readInit

bool SpMatFull::moveBlue() {
    bool moved(0),  //At least one movement has been done
            mv1;    //Take track of the movement of the first of the col
#pragma omp parallel for schedule(guided,5)
    for (int c = 0; c < m_nC; c++) {
        bool mv1(false);
        for (int r = 0; r < m_nR; r++) {
            int nxr = (r+1)%m_nR;
            if (   m_m[r*m_nC+c]   != BLUE  //Curr not blue
                || m_m[nxr*m_nC+c] != EMPTY //Next not empty
                || (r == m_nR-1 && mv1)     //Last and first moved
                ) {
                //Here, we do NOT move
                continue;
            } else {
                //Here, we MOVE
                //Updates
                moved = true;
                if (r == 0)
                    mv1 = true;
                //Free & Move
                m_m[r  *m_nC+c] = EMPTY;
                m_m[nxr*m_nC+c] = BLUE;
                r++;    //Avoid finding again block just moved
            }//end_if_move
        }//end_for_r
    }//end_for_c
    return moved;
}//end_MOVEBLUE
bool SpMatFull::moveRed() {
    bool moved(0),  //At least one movement has been done
            mv1;    //Take track of the movement of the first of the row
    #pragma omp parallel for schedule(guided,5)
    for (int r = 0; r < m_nR; r++) {
        bool mv1(false);
        for (int c = 0; c < m_nC; c++) {
            int nxc = (c+1)%m_nC;
            if (   m_m[r*m_nC+c]   != RED   //Curr not red
                || m_m[r*m_nC+nxc] != EMPTY //Next not empty
                || (c == m_nC-1 && mv1)     //Last and first moved
               ) {
                //Here, we do NOT move
                continue;
            } else {
                //Here, we MOVE
                //Updates
                moved = true;
                if (c == 0)
                    mv1 = true;
                //Free & Move
                m_m[r*m_nC+c]   = EMPTY;
                m_m[r*m_nC+nxc] = RED;
                c++;    //Avoid finding again block just moved
            }//end_if_move
        }//end_for_c
    }//end_for_r
    return moved;
}//end_MOVERED

