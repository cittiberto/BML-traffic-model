#ifndef VECTBEG_H
#define VECTBEG_H

#include <vector>
#include <algorithm>
#include <Eigen/Sparse>

using namespace std;

class VectBeg
{
protected:
    vector<int> m_v;    //Positions of the cars in acendent order
    int m_b;            //Position (in the vector here above) of the frist car of the row/col
public:
    //Ctor
    VectBeg() : m_v(), m_b(0) {}

    //"Overloads"
    int                     size()              {return m_v.size();}
    vector<int>::iterator   beginT()            {return m_v.begin();}
    vector<int>::iterator   endT  ()            {return m_v.end();}
    bool                    empty()             {return m_v.empty();}
    int                     operator[](int i)   {return m_v[i];}
    void                    push_back(int val)  {m_v.push_back(val);}
    void                    reserve(int i)      {m_v.reserve(i);}

    //Gets & Sets
    void setithvalue (int i, int nwval) {
        m_v[i] = nwval;
    }

    //Find
    bool isThere(int i) {
        return (find(m_v.begin(),m_v.end(),i) != m_v.end());
    }
    bool isThereBS(int i) {
        //Binary search
        int sz(static_cast<int>(m_v.size())),
                min(m_b),
                max(m_b+sz-1);
        while (min < max) {
            int mid((min+max)/2);
            if ( mid < max) {
                if (m_v[mid%sz] < i) {
                    min = mid + 1;
                } else {
                    max = mid;
                }
            } else {
                break;
            }
        }//end_while
        if ( (min == max) && (m_v[min%sz] == i) )
            return true;
        else
            return false;
    }

    void addTriplets(vector<Eigen::Triplet<int>>& elem, int color, int pr_index) {
        //Color: 1 == Blue => pr_index=col
        //       2 == Red  => pr_index=row
        if (color == 1) {
            for(auto it = m_v.begin(); it != m_v.end(); it++)
                elem.push_back(Eigen::Triplet<int>(*it,pr_index,1));
        } else {
            for(auto it = m_v.begin(); it != m_v.end(); it++)
                elem.push_back(Eigen::Triplet<int>(pr_index,*it,2));
        }
    }

    friend class SpMatVwBBS;
};

#endif // VECTBEG_H
