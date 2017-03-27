#ifndef SPMAT_H
#define SPMAT_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;


class SpMat
{
protected:
    int m_nR, m_nC;
    bool m_full, m_empty;
    bool m_allR, m_allB;
    vector<int> m_steps;

public:
    //Ctor
    SpMat(int R = 0, int C = 0) : m_nR(R), m_nC(C),
                                  m_full(true), m_empty(true),
                                  m_allR(true), m_allB(true),
                                  m_steps() {}
    //Destructor
    virtual ~SpMat() {}
    //Gets & Sets
    int         getnR()     {return m_nR;}
    int         getnC()     {return m_nC;}

    void    setnR(int nR)           {m_nR = nR;   return;}
    void    setnC(int nC)           {m_nC = nC;   return;}
    void    copyS(vector<int>& s)   {m_steps = s; return;}
    virtual void setN0el(int i) = 0;

    //Initialization
    virtual void readInit (const string& filename) = 0;
    void sortUniqueSteps() {
        //Sort the steps...
        sort(m_steps.begin(),m_steps.end(),[](int i, int j) -> bool {return i<j;});
        //...erase duplicates
        m_steps.erase( unique( m_steps.begin(), m_steps.end() ), m_steps.end() );
    }

    //Moves
    virtual bool moveBlue() = 0;
    virtual bool moveRed () = 0;
    virtual bool move()     = 0;

    //Recognize
    virtual string whoRU() = 0;

    //I/O functions
    virtual void oprint(ostream &o) = 0;
    friend ostream& operator<< (ostream &o, SpMat &m) {
        m.oprint(o);
        return o;
    }


    //Utilities
    bool getFull   ()        {return m_full;}
    bool getEmpty  ()        {return m_empty;}
    bool getAllR   ()        {return m_allR;}
    bool getAllB   ()        {return m_allB;}
    int& getMaxStep()        {return m_steps.back();}
    int  getSize   ()        {return m_steps.size();}
    int  getAStep  (int i)   {return m_steps[i];}

    int getType    () {
        if      (m_empty)       return -1;
        else if (m_full)        return 0;
        else if (m_allB)        return 1;
        else if (m_allR)        return 2;
        else /*generic matrix*/ return 3;
    }


};

#endif // SPMAT_H
