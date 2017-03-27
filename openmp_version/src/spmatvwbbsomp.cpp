#include "spmatvwbbs.h"

#include <omp.h>
#include <ostream>
#include <unordered_map>
#include <Eigen/Sparse>

void SpMatVwBBS::readInit(const string &filename) {
    /* nest = number of estimated elements in a row/col
     * It'll be used to call a reserve on the VectBeg (avoid the time-consuming reallocation)
     * In a VectBeg there will be blocks of only one type, so their number should be half of nest
     * This means that calling a reserve w/ nest should allocate enough memory
     */
    int r(0), c(0);

    ifstream    file(filename);
    string      tmpline;
    if (file) {
        string       anelem;
        stringstream aline;
        getline(file, tmpline); //First row is already been read
        while ( getline(file, tmpline) ) {
            //Parse the document and store the line
            c = 0;
            aline << tmpline;

                    m_rd.push_back(VectBeg());
                    if (n0el > 0)
                        m_rd[r].reserve(n0el);
                while ( getline(aline, anelem, ',') ) {
                    if (c >= static_cast<int>(m_bl.size())){ //We should insert an element for the new col
                        m_bl.push_back(VectBeg());
                        if (n0el > 0)
                            m_bl[c].reserve(n0el);
                    }
                    //Parse the line
                    if (anelem.compare("0") == 0) {             //Empty location
                        m_full = false;
                    } else {
                        m_empty = false;
                        if (anelem.compare("1") == 0){          //Blue block
                            m_allR = false;
                            //Principal index = Col
                            m_bl[c].push_back(r);
                        } else if (anelem.compare("2") == 0){   //Red block
                            m_allB = false;
                            //Principal index = Row
                            m_rd[r].push_back(c);
                        } else {
                            //Throw an exception
                            throw("Not allowed char");
                        }
                    }//end_if_NOT_empty
                    c++;
                }//end_while_parseALine
                if (c!=0 && c != m_nC) {
                    //Throw an exception
                    throw("Not the same number of col");
                }
            aline.clear();
            r++;
        }//end_while_getline
        file.close();
    } else throw("File not found");
}//end_readInit

bool SpMatVwBBS::basicMove(bool color) {
    vector<VectBeg>* vv;
    vector<VectBeg>* notv;
    bool moved(0);
    int nV, nT;
    // color == 1 iff blue iff col-based Remember: t is the principal index
    if (color) {vv = &m_bl; notv = &m_rd; nV = m_nR; nT = m_nC;}
    else       {vv = &m_rd; notv = &m_bl; nV = m_nC; nT = m_nR;}

#pragma omp parallel for schedule(dynamic)
    for (int t = 0; t < nT; t++) {
        //Alternative version
//        vector<VectBeg>::iterator itT = vv->begin();
//        advance(itT,t);
        auto itT = &((*vv)[t]); //Pointer to VectBeg
        if (!itT->empty()) {
            bool mv1(false);
            int s(0), sz = itT->size();
            while (s < sz) {
                int idx = (itT->m_b+s++)%sz,    //Indexes and s update
                    nxi = (idx+1)%sz;
                int v   = (*itT)[idx],          //Values
                    nxv = (v+1)%nV;
                if (    (*itT)[nxi] == nxv                      //Next engaged, same type
                     || (v == nV-1 && mv1)                      //Last & first moved
                     || (nxv < static_cast<int>(notv->size())
                         && (*notv)[nxv].isThereBS(t))          //Next engaged, other type
                    ) {
                    //Here we do NOT move
                    continue;
                } else {
                    //Here we move
                    moved = true;
                    if (v == 0)
                        mv1 = true;
                    if (nxv == 0)
                        itT->m_b = idx;
                    itT->setithvalue(idx,nxv);
                }
            }//end_while
        }//end_if_NOT_empty
    }//end_for_itT
    return moved;
}//end_BASICMOVE

void SpMatVwBBS::oprint_norm(ostream &o) {
int tp = getType();
if (m_nR*m_nC > 2500000) {

    //Parallel version

    //"ORDERED" VERSION
    switch(tp)
    {
    case -1://Empty; only zeros
        for (int r = 0; r < m_nR; r++) {
            for (int c = 0; c < m_nC-1; c++)
                o << "0,";
            o << "0\n";
        }//end_for_R

        break;

    case 1: //Only Blues
        #pragma omp parallel for ordered schedule(dynamic)
        for (int r = 0; r < m_nR; r++) {
            stringstream ss;
            for (int c = 0; c < m_nC-1; c++)
                ss << m_bl[c].isThereBS(r) << ",";
            ss << m_bl[m_nC-1].isThereBS(r) << "\n";
            #pragma omp ordered
            o << ss.str();
        }//end_for_R

        break;

    case 2: //Only Reds
        #pragma omp parallel for ordered schedule(dynamic)
        for (int r = 0; r < m_nR; r++) {
            stringstream ss;
            for (int c = 0; c < m_nC-1; c++)
                ss << ( m_rd[r].isThereBS(c) ? 2 : 0 ) << ",";
            ss << ( m_rd[r].isThereBS(m_nC-1) ? 2 : 0 ) << "\n";
            #pragma omp ordered
            o << ss.str();
        }//end_for_R

        break;

    default://Full or Generic
        if (m_allB) {       //It can be full of Blues...
            for (int r = 0; r < m_nR; r++) {
                for (int c = 0; c < m_nC-1; c++)
                    o << "1,";
                o << "1\n";
            }//end_for_R

            break;
        } else if (m_allR) {//...or of Reds
            for (int r = 0; r < m_nR; r++) {
                for (int c = 0; c < m_nC-1; c++)
                    o << "2,";
                o << "2\n";
             }//end_for_R

            break;
        } else {            //Generic
            #pragma omp parallel for ordered schedule(dynamic, 5)
            for (int r = 0; r < m_nR; r++) {
                stringstream ss;
                for (int c = 0; c < m_nC; c++) {
                    if (m_rd[r].isThereBS(c)) {
                        //It's a RED
                        ss << 2;
                    } else if (m_bl[c].isThereBS(r)) {
                        //It's a BLUE
                        ss << 1;
                    } else {
                        //It's a ZERO
                        ss << 0;
                    }//end_if_COLOR
                    if (c != m_nC-1)
                        ss << ",";
                }//end_for_C
                ss << "\n";
                #pragma omp ordered
                {
                o << ss.str();
                }
            }//end_for_R
            break;
        }
        break;
    }//end_switch_TP
} else {
    //Sequential version

    switch(tp)
    {
    case -1://Empty; only zeros
        for (int r = 0; r < m_nR; r++) {
            for (int c = 0; c < m_nC-1; c++)
                o << "0,";
            o << "0\n";
        }//end_for_R
        break;

    case 1: //Only Blues
        for (int r = 0; r < m_nR; r++) {
            for (int c = 0; c < m_nC-1; c++)
                o << m_bl[c].isThereBS(r) << ",";
            o << m_bl[m_nC-1].isThereBS(r) << "\n";
        }//end_for_R
        break;

    case 2: //Only Reds
        for (int r = 0; r < m_nR; r++) {
            for (int c = 0; c < m_nC-1; c++)
                o << ( m_rd[r].isThereBS(c) ? 2 : 0 ) << ",";
            o << ( m_rd[r].isThereBS(m_nC-1) ? 2 : 0 ) << "\n";
        }//end_for_R
        break;

    default://Full or Generic
        if (m_allB) {       //It can be full of Blues...
            for (int r = 0; r < m_nR; r++) {
                for (int c = 0; c < m_nC-1; c++)
                    o << "1,";
                o << "1\n";
            }//end_for_R
            break;
        } else if (m_allR) {//...or of Reds
            for (int r = 0; r < m_nR; r++) {
                for (int c = 0; c < m_nC-1; c++)
                    o << "2,";
                o << "2\n";
            }//end_for_R
            break;
        } else {            //Generic
            for (int r = 0; r < m_nR; r++) {
                for (int c = 0; c < m_nC; c++) {
                    if (m_rd[r].isThereBS(c)) {
                        //It's a RED
                        o << 2;
                    } else if (m_bl[c].isThereBS(r)) {
                        //It's a BLUE
                        o << 1;
                    } else {
                        //It's a ZERO
                        o << 0;
                    }//end_if_COLOR
                    if (c != m_nC-1)
                        o << ",";
                }//end_for_C
                o << "\n";
            }//end_for_R
            break;
        }
        break;
    }//end_switch_TP
}
}//end_OPRINT
void SpMatVwBBS::oprint_fast(ostream &o) {

    Eigen::SparseMatrix<int> cop(m_nR,m_nC);
    cop.reserve(m_nR*n0el);
    Eigen::MatrixXi dcop(m_nR,m_nC);
    vector<vector<Eigen::Triplet<int>>> elem;

    int nt;
#pragma omp parallel
    {
    int id = omp_get_thread_num();
    nt = omp_get_num_threads();
#pragma omp master //blocking
        {
            elem.reserve(nt);
            vector<Eigen::Triplet<int>> tmp;
            tmp.reserve(m_nR*n0el);
            elem.push_back(tmp);
            for (int i=1; i<nt; i++) {
                vector<Eigen::Triplet<int>> tmp1;
                tmp1.reserve(int(m_nR/nt*n0el));
                elem.push_back(tmp1);
            }
        }
#pragma omp barrier
#pragma omp for nowait schedule(dynamic) //blocking not necessary
    for(int i = 0; i < m_nC; i++)
        m_bl[i].addTriplets(elem[id],1,i); //BLUE => 1
#pragma omp for        schedule(dynamic)
    for(int i = 0; i < m_nR; i++)
        m_rd[i].addTriplets(elem[id],2,i); //RED  => 2
#pragma omp master
    {
    for(int i=1; i<nt; i++)
        elem[0].insert(elem[0].end(),elem[i].begin(),elem[i].end());

    cop.setFromTriplets(elem[0].begin(), elem[0].end());

    dcop = Eigen::MatrixXi(cop);    //Conversion to dense type in order to eploit the operator<<

    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");

    o << dcop.format(CSVFormat) << "\n";
    }
    }//END_OMP_PARALLEL

}



//////////////////////////////////////////////////////////////////
///////////           ALTERNATIVE   VERSION            ///////////
//////////////////////////////////////////////////////////////////

 //Parallel version
//int numT;
//int tp = getType();
//unordered_map<int,string> os;
////    vector<int> ids;
//bool para = false;
//stringstream op;

//#pragma omp parallel private(op)
//{
//    int tID = omp_get_thread_num();
//    numT = omp_get_num_threads();
//    int chunk = int(m_nR/numT);
//    if (chunk != 0) {
//        //This means m_nR > numT
//        para = true;
//        int r_begin = chunk*tID,
//            r_end   = ( (tID==(numT-1)) ? m_nR:(chunk*(tID+1)) ); //Last thread take all the remaining rows

//        switch(tp)
//        {
//        case -1://Empty; only zeros
//            for ( int r=r_begin; r<r_end; r++) {
//                for (int c = 0; c < m_nC-1; c++)
//                    op << "0,";
//                op << 0 << endl;
//            }//end_for_R
//            break;

//        case 1: //Only Blues
//            for ( int r=r_begin; r<r_end; r++) {
//                for (int c = 0; c < m_nC-1; c++)
//                    op << m_bl[c].isThereBS(r) << ",";
//                op << m_bl[m_nC-1].isThereBS(r) << endl;
//            }//end_for_R
//            break;

//        case 2: //Only Reds
//            for ( int r=r_begin; r<r_end; r++) {
//                for (int c = 0; c < m_nC-1; c++)
//                    op << ( m_rd[r].isThereBS(c) ? 2 : 0 ) << ",";
//                op << ( m_rd[r].isThereBS(m_nC-1) ? 2 : 0 ) << endl;
//            }//end_for_R
//            break;

//        default://Full or Generic
//            if (m_allB) {       //It can be full of Blues...
//                for ( int r=r_begin; r<r_end; r++) {
//                    for (int c = 0; c < m_nC-1; c++)
//                        op << "1,";
//                    op << 1 << endl;
//                }//end_for_R
//                break;
//            } else if (m_allR) {//...or of Reds
//                for ( int r=r_begin; r<r_end; r++) {
//                    for (int c = 0; c < m_nC-1; c++)
//                        op << "2,";
//                    op << 2 << endl;
//                }//end_for_R
//                break;
//            } else {            //Generic
//                for ( int r=r_begin; r<r_end; r++) {
//                    for (int c = 0; c < m_nC; c++) {
//                        if (m_rd[r].isThereBS(c)) {
//                            //It's a RED
//                            op << 2;
//                        } else if (m_bl[c].isThereBS(r)) {
//                            //It's a BLUE
//                            op << 1;
//                        } else {
//                            //It's a ZERO
//                            op << 0;
//                        }//end_if_COLOR
//                        if (c != m_nC-1)
//                            op << ",";
//                    }//end_for_C
//                    op << endl;
//                }//end_for_R
//                break;
//            }
//            break;
//        }//end_switch_TP

//        os.insert(make_pair(tID,op.str()));

//    #pragma omp single
//        {
//            for (int i = 0; i<numT; i++) {
//                o << os[i];
//            }
//        }

//    }//end_if_CHUNK
//    else
//    {
//    #pragma omp master
//        {
//            switch(tp)
//            {
//            case -1://Empty; only zeros
//                for (int r = 0; r < m_nR; r++) {
//                    for (int c = 0; c < m_nC-1; c++)
//                        o << "0,";
//                    o << 0 << endl;
//                }//end_for_R
//                break;

//            case 1: //Only Blues
//                for (int r = 0; r < m_nR; r++) {
//                    for (int c = 0; c < m_nC-1; c++)
//                        o << m_bl[c].isThereBS(r) << ",";
//                    o << m_bl[m_nC-1].isThereBS(r) << endl;
//                }//end_for_R
//                break;

//            case 2: //Only Reds
//                for (int r = 0; r < m_nR; r++) {
//                    for (int c = 0; c < m_nC-1; c++)
//                        o << ( m_rd[r].isThereBS(c) ? 2 : 0 ) << ",";
//                    o << ( m_rd[r].isThereBS(m_nC-1) ? 2 : 0 ) << endl;
//                }//end_for_R
//                break;

//            default://Full or Generic
//                if (m_allB) {       //It can be full of Blues...
//                    for (int r = 0; r < m_nR; r++) {
//                        for (int c = 0; c < m_nC-1; c++)
//                            o << "1,";
//                        o << 1 << endl;
//                    }//end_for_R
//                    break;
//                } else if (m_allR) {//...or of Reds
//                    for (int r = 0; r < m_nR; r++) {
//                        for (int c = 0; c < m_nC-1; c++)
//                            o << "2,";
//                        o << 2 << endl;
//                    }//end_for_R
//                    break;
//                } else {            //Generic
//                    for (int r = 0; r < m_nR; r++) {
//                        for (int c = 0; c < m_nC; c++) {
//                            if (m_rd[r].isThereBS(c)) {
//                                //It's a RED
//                                o << 2;
//                            } else if (m_bl[c].isThereBS(r)) {
//                                //It's a BLUE
//                                o << 1;
//                            } else {
//                                //It's a ZERO
//                                o << 0;
//                            }//end_if_COLOR
//                            if (c != m_nC-1)
//                                o << ",";
//                        }//end_for_C
//                        o << endl;
//                    }//end_for_R
//                    break;
//                }
//                break;
//            }//end_switch_TP
//        }//END_OF_MASTER
//    }//END_OF_IFELSE_CHUNK
//}//end_PRAGMA_PARALLEL

