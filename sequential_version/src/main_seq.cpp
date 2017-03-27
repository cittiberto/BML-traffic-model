#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <array>
#include <cmath>
#include <stdlib.h>
//#include <algorithm>
//#include <ctime>
//#include <chrono>
//#include <thread>


#include "spmatfull.h"
#include "spmatvwbbs.h"

using namespace std;

void firstRead  (const string & filename, array<int,4>& ret, vector<int> &steps, int treshold = -1, int trR = -1);
void eval       (SpMat *m_m, int N);
void perform    (SpMat *m_m, int N);

int main()
{

    string filename("problem.csv");

    array<int,4> ret = {0,0,0,0};
    vector<int> vect;

    //Creating a pointer to an abstract object Matrix
    SpMat* m;

    try {
        //A first read of the matrix to extract some useful infos
        firstRead(filename,ret,vect,-1,10);

        if (ret[1] == 0) /*Empty matrix*/ throw("No matrix");
        double est;
        /* According to some rough BigO estimations, iterations are faster with sparse method if
         * 1<d*ln(d*M) where d is the density and M is the dimension of the matrix
         * We'll use this to choose the right method to use. However, we'll set treshold
         * lower than 1 in order to take into account the higher time needed when printing
         * in sparse method
         */
        est = (double(ret[3])/double(ret[1]))*log(double(ret[3]));


        //Test which model should be better
        if (   (ret[1]*ret[2] < 40000)     //Small matrix, smaller than 200x200
            || (est > 0.85)                 //High density
              )
        {
            //Full
            m  = new SpMatFull(ret[1],ret[2]);
            //Second read and actual initialization
            m->readInit(filename);
        } else {
            //Sparse
            m  = new SpMatVwBBS(ret[1],ret[2]);
            m->setN0el(ret[3]);
            //Second read and actual initialization
            m->readInit(filename);
        }
        m->copyS(vect);

        //Perform everything
//        perform(m, ret[0]);
        eval(m,ret[0]);

        delete m;

    } catch (...) {
        //Just finish (that's what it has been said on the forum)
        delete m;
        exit(EXIT_FAILURE);
    }

    return 0;
}

void firstRead(const string & filename, array<int,4>& ret, vector<int>& steps, int treshold, int trR) {
    /* Inputs:
     * [string] filename: .csv to be read
     * [array]  ret     : cf here below
     * [vector] steps   : steps at which the matrix should be printed
     * [int]    treshold: if the max iteration to be printed is less than this,
     *                      break the process (if max_iter is very small, reading
     *                      two times would mean spend an amount of time comparable to the whole computational time)
     *                      -1 means that the matrix is read
     * [int]    trR     : max num of rows to read (if trR great, we'll have an accurate approx of the density but this costs time)
     *                      -1 means that all the rows are read
     */

    /* ret:
     *  - [0] = max N_iterations
     *  - [1] = N_Row
     *  - [2] = N_Col
     *  - [3] = N_nn, non-null elements per row, used to compute an extimation of the density
     */

    //Be sure to read at least one line
    if (trR == 0) trR = 1;

    double valD;
    int r(0), c(0), nn(0), rr(-1), stdC(c);
    ifstream    file(filename);
    string      tmpline;
    if (file) {
        string       anelem;
        stringstream aline;

        //First line
        if (getline(file, tmpline)) {
            aline << tmpline;
            while ( getline(aline, anelem, ',') ) {
                //Parse the line
                stringstream tmp(anelem);
                tmp >> valD;
                if (floor(valD) == valD) /*Integer*/ steps.push_back(floor(valD));
                tmp.clear();
            }//end_WHILE

            //Sort the steps...
            sort(steps.begin(),steps.end(),[](int i, int j) -> bool {return i<j;});
            //...erase duplicates
            steps.erase( unique( steps.begin(), steps.end() ), steps.end() );

            ret[0] = steps.back();
        }//end_if_GETLINE
        aline.clear();

        if ( treshold == -1 || ret[0] >= treshold ) {
            while ( getline(file, tmpline) ) {
                r++;
                if (trR == -1 || r <= trR) { //Read each element of the line
                    //Parse the document and store the line
                    aline << tmpline;
                    c = 0;
                    while ( getline(aline, anelem, ',') ) {
                        c++;
                        if(   anelem.compare("1") == 0
                           || anelem.compare("2") == 0)
                            nn++;
                        else if (anelem.compare("0") != 0)
                            throw("Not allowed char");
                    }//end_while_ANELEM

                    if(r==1) stdC=c;
                    else if(c!=0 && c!=stdC) throw("Not same number of col");

                    aline.clear();
                }//end_if_trR
                if (c==0) /*Problem with the last line*/ r--;
                rr = r;     //#rows actually entirely read (can be <trR)
            }//end_while_GETLINE
            //Build ret
            ret[1] = r;
            ret[2] = stdC;
            ret[3] = nn/rr;

            file.close();
        }//end_if_TRESHOLD
    } else throw("File not found");    //end_if_file
}//END_FIRST_READ

void eval(SpMat *m_m, int N) {
    //N : max #iteration to print
    bool moved(true), moved_blue(true), moved_red(true);
    int n(1);
    int s(0), toPrint(m_m->getAStep(s));
    int tp = m_m->getType();

    if (toPrint == 0) {    //Deal with exception: print before any moves
        ofstream outfile("0.csv");
        outfile << *m_m;
        outfile.close();
        toPrint = m_m->getAStep(++s);
    }

    if ( tp > 0) {                  //Not Empty nor Full
        while (moved && n<=N) {
            if (    n%2             //Time for blue
                 && tp != 2) {      //There are blues
                moved_red  = true;  //A new cycle starts, reset
                moved_blue = m_m->moveBlue();
            } else if ( tp != 1){   //There are reds
                moved_red  = m_m->moveRed();
            }//end_if

            moved = moved_blue || moved_red;

            if ( n == toPrint) {
                ofstream outfile(to_string(n)+".csv");
                outfile << *m_m;
                outfile.close();
                toPrint = m_m->getAStep(++s);
            }
            n++;
        }//end_while

        if (!moved) {
            //Here if we haven't reached N, but no more moves are possible
            for (; s < static_cast<int>(m_m->getSize()); s++) {
                ofstream outfile(to_string(m_m->getAStep(s))+".csv");
                outfile << *m_m;
                outfile.close();
            }//end_for
        }//end_if_moved
    }//end_if_tp
}//END_EVAL


void perform(SpMat *m_m, int N) {
    //N : max #iteration to print
    bool moved(true);
    int n(1);
    int s(0), toPrint(m_m->getAStep(s));
    int tp = m_m->getType();

    if (toPrint == 0) {    //Deal with exception: print before any moves
        ofstream outfile("0.csv");
        outfile << *m_m;
        outfile.close();
        toPrint = m_m->getAStep(++s);
    }

    //Select a specific way of moving accordingly to the type
    switch (tp)
    {
    case -1://Empty ; no moves needed, just print
        break;

    case 0: //Full; no moves needed, just print
        break;

    case 1: //Only Blue
        while (moved && n <= N) {
            moved = m_m->moveBlue();
            if ( n == toPrint) {
                ofstream outfile(to_string(n)+".csv");
                outfile << *m_m;
                outfile.close();
                toPrint = m_m->getAStep(++s);
            }
            n++;
        }//end_while
        break;

    case 2: //Only Red
        while (moved && n <= N) {
            moved = m_m->moveRed();
            if (n == toPrint) {
                ofstream outfile(to_string(n)+".csv");
                outfile << *m_m;
                outfile.close();
                toPrint = m_m->getAStep(++s);
            }
            n++;
        }//end_while
        break;

    default: //Generic
        while (moved && n <= N) {
            moved = m_m->move();
            if (n == toPrint) {
                ofstream outfile(to_string(n)+".csv");
                outfile << *m_m;
                outfile.close();
                toPrint = m_m->getAStep(++s);
            }
            n++;
        }//end_while
        break;
    }//end_switch

    if (!moved) {
        //Here if we haven't reached N, but no more moves are possible
        for (; s < static_cast<int>(m_m->getSize()); s++) {
            ofstream outfile(to_string(m_m->getAStep(s))+".csv");
            outfile << *m_m;
            outfile.close();
        }//end_for
    }//end_if_nN
}//END_PERFORM


