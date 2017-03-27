#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <list>
#include <algorithm>


using namespace std;

# include "mpi.h"



// Questa versione MPI è strutturata come segue:
// 1) Nel main vengono lette dimensioni e densità della matrice
// 2) In base alla densità viene scelto l'algoritmo per la matrice SPARSA o quella PIENA
// 3) Nelle due funzioni, si legge la matrice e la si divide tra i vari processori.
// 4) Vengono fatti gli spostamenti per ogni processore separatamente, scambiando le frontiere dei blocchi ad ogni iterazione (lunghezza del messaggio = lunghezza della linea)
// 5) La matrice viene infine stampata su file da tutti i processori in parallelo con la funzione MPI_File_write

// NB Per la matrice sparsa abbiamo usato delle liste. Per poter passare le righe tra un processore all'altro, le trasformiamo in un vettore di interi nella forma (10 * posizione + colore), in cui colore è o 1 (per i blu) o 2 (per i rossi).






// Classe per le liste di oggetti della verione sparsa
class car{
public:
    car(int p=0, int c=0): pos(p), col(c), moved(false) {};
    int pos;
    int col;
    bool moved;
};

// Funzione per leggere le dimensioni
void read_dim(string s_filename,  int *numb , vector<int> &steps);

// Funzione per leggere la matrice se piena
void read_mat(string s_filename, int *vec , int *numb, vector<int> &steps);

// Funzione per leggere la matrice se sparsa
void read_mat_list(string s_filename, list<car>& l,int *numb, vector<int> &steps);

// Funzione per fare tutti gli spostamenti e gli output per la matrice SPARSA
void mpi_sparse(int argc, char **argv, int lun, int lar, string s_filename);

// Funzione per fare tutti gli spostamenti e gli output per la matrice PIENA
void mpi_full(int argc, char **argv, int lun, int lar,string s_filename);

// Funzione per fare tutti gli spostamenti e gli output per la matrice PIENA se lavora un solo processore
void mpi_full_oneproc(int *numbers, vector<int> &steps, int *vec);



int main(int argc, char **argv) {
    
    try{
        ///// 1.  initialize MPI.
        MPI::Init ( argc, argv );
        int id = MPI::COMM_WORLD.Get_rank ( );
        int p = MPI::COMM_WORLD.Get_size ( );
        
        //  MPI::COMM_WORLD.Set_errhandler ( MPI::ERRORS_THROW_EXCEPTIONS );
        
       // if(id==0) cout << "Number of Proc: " << p << endl;
        
        // TIME
//        clock_t start, end;
//        if (id==0) start = clock();
//        
        int *numbers_dim;
        vector<int>  steps;
        int lar, lun, dens, FulSpar;
        
        string s_filename = "problem.csv";
        
        
        ///// 2.  READ dimensions from the root
        if ( id == 0 ){
            numbers_dim = new int[3];
            read_dim(s_filename, numbers_dim, steps);
            
            lun=numbers_dim[0];
            lar=numbers_dim[1];
            dens=numbers_dim[2];
            
            // Chose between SPARSE (1) or FULL (0)
            // FulSpar= ( (dens<8 && lar*lun>10000) ? 1 : 0);
            FulSpar= ( (dens<8 ) ? 1 : 0);
            
            delete[] numbers_dim;
        }
        
        ///// 3. broadcast dimension and sparse indicator
        MPI::COMM_WORLD.Bcast(&FulSpar, 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&lun, 1, MPI::INT, 0);
        MPI::COMM_WORLD.Bcast(&lar, 1, MPI::INT, 0);
        


        
        ///// 4.  SPARSE or FULL
        
        
        if(FulSpar) mpi_sparse(argc, argv, lun, lar, s_filename);
        
        else mpi_full(argc, argv, lun, lar, s_filename);
        
        
        
        
        // final TIME
//        if (id==0) {
//            end = clock();
//            std::cout << "Process took " << (double(end - start) / CLOCKS_PER_SEC) << "seconds" << '\n';
//        }
//        
        MPI::Finalize ( );
        
    } catch (...) {
        
        MPI::COMM_WORLD.Abort (-1) ;
        
        MPI::Finalize ( );
        exit(EXIT_FAILURE);
    }
    
    
    return 0;
} //end


void mpi_sparse(int argc, char **argv, int lun, int lar,string s_filename){
    
    
    ////// INITIALIZATION //////
    
    // reinizializzazione MPI se serve
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized) MPI::Init ( argc, argv );
    
    int id = MPI::COMM_WORLD.Get_rank ( );
    int p = MPI::COMM_WORLD.Get_size ( );
    
    if(id==0) cout << "SPARSE" << endl;
    
    MPI::Status status;
    MPI_Status status1;
    MPI_Offset offset;
    MPI_File   file;
    
    // variabili COMUNICAZIONE
    int chunk_size, dest, tag, num_workers, ini_root, lun_chunk, size_steps, ini_root_chunk, root_chunk_size;
    int *numbers, *vec, *vec_proc, *row_bef, *row_second, *row_last, *row_aft, *prova, *vec_rec,*vec_send;
    int size_vector_trasf1, size_vector_trasf2;
    
    // variabili EXPORT
    ofstream myfile;
    char *filename;
    string ss_in_string;
    int msgsize;
    
    // variabili matrice
    int dens, full, empty, allbl, allrd, iterazioni_max, s2pr(0), nxst;
    vector<int>  steps;
    numbers = new int[7];
    list<car> l, l_tot;
    
    int position, color;
    
    list<car>::iterator it_send1;
    list<car>::iterator it_send2;
    list<car>::iterator  it_next;
    
    ////// READ MATRIX //////
    
    if(id==0){
        //inizializzo
        // vec = new int[lun*lar];
        //leggo
        
        //read_mat(s_filename, vec, numbers, steps);
        read_mat_list(s_filename, l_tot, numbers, steps);
        
        cout << "legge matrice " << endl;
        
        it_send1=l_tot.begin();

        
        //steps
        size_steps=steps.size();
        iterazioni_max=steps.back();
        nxst=steps[0];
        cout << "passa i parametri" << endl;
        
    }//end_if_ID0 (initialization)
    
    
    // passo i parametri a tutti i processori
    MPI::COMM_WORLD.Bcast(numbers, 7, MPI::INT, 0);
    
    MPI::COMM_WORLD.Bcast(&iterazioni_max, 1, MPI::INT, 0);
    
    lun=numbers[0];         lar=numbers[1];
    full  = numbers[3];     empty = numbers[4];
    allbl = numbers[5];     allrd = numbers[6];
    
    
    // numero di PROCESSI
    num_workers=p;
    if (lun<=4) num_workers=1;
    
    // chunk SIZE
    chunk_size= lun / num_workers;
    chunk_size*= lar;
    
    
    
    //if full or empty, print
    
    
    if ( (full!=0) || (empty!=0)) {
        if(id==0){
            for (;s2pr<size_steps;s2pr++) {
                nxst  = steps[s2pr];
                string file_exp = to_string(nxst) + ".csv";
                ofstream myfile(file_exp);
                
                string ss_in_string;
                it_next=l_tot.begin();
                
                for (int i=0; i<lun; i++) {
                    for(int j=0; j< lar ; j++){
                        
                        if(it_next->pos == (j + lar*i) && it_next!=l_tot.end() ){
                            ss_in_string.append( to_string( it_next->col ) );
                            it_next++;
                        }
                        else ss_in_string.append("0");
                        
                        if(j<lar-1) ss_in_string.append(",");
                    } ss_in_string.append("\n");
                }
                
                myfile << ss_in_string;
                myfile.close();
            }
        }
    }
    else{
        ////// dividing the matrix and broadcast of the CHUNKS
        
        if ( id == 0 ){
            cout << "non è piena" << endl;
            it_send1=l_tot.begin();
            
            while (it_send1!=l_tot.end() ){
                it_send1++;
            }
            
            
            if (num_workers>1) {
                vector<int> vector_trasf1, vector_trasf2;
                vector_trasf1.reserve(chunk_size/10);
                vector_trasf2.reserve(chunk_size/10);
                
                //SEND chunks
                // send to 1
                
                dest = 1;
                tag = 0;
                
                //trasformo la lista in vettore
                
                //ultima riga
                it_send1=l_tot.end();
                
                //sposto l'iteratore alla penultima linea
                it_send1--;
                while (it_send1->pos >=  (lun-1)*lar && it_send1!=l_tot.begin() ) it_send1--;
                it_send1++;
                
                while (it_send1->pos <  lun*lar && it_send1!=l_tot.end() ){
                    
                    position= it_send1->pos - (lun-1)*lar;
                    
                    vector_trasf1.push_back( 10*(position) + it_send1->col );
                    it_send1++;
                }
                
                //resto
                it_send1=l_tot.begin();
                
                while (it_send1->pos <  (chunk_size+lar)  && it_send1!=l_tot.end() ){
                    
                    position= it_send1->pos + lar;
                    
                    vector_trasf1.push_back( 10*(position) + it_send1->col );
                    it_send1++;
                    
                }
                
                size_vector_trasf1= vector_trasf1.size();
                
                
                //mando
                MPI::COMM_WORLD.Send ( &size_vector_trasf1, 1 , MPI::INT, dest, 670 );
                MPI::COMM_WORLD.Send ( &vector_trasf1[0], size_vector_trasf1, MPI::INT, dest, tag );
                

                
                
                //send to others
                for (int i = 2; i < num_workers; i++ )
                {
                    dest = i;
                    tag = chunk_size*(i-1);
                    
                    //sposto l'iteratore alla linea -1 del chunk corrispondente
                    it_send1--;
                    while (it_send1->pos >=  tag-lar && it_send1!=l_tot.begin() ) it_send1--;
                    it_send1++;
                    
                    while (it_send1->pos <  tag+chunk_size+lar  && it_send1!=l_tot.end() ){
                        
                        position= it_send1->pos - (chunk_size*(i-1)-lar);
                        
                        vector_trasf2.push_back( 10*(position) + it_send1->col );
                        it_send1++;
                        
                    }
                    size_vector_trasf2= vector_trasf2.size();
                    
                    
                    MPI::COMM_WORLD.Send ( &size_vector_trasf2, 1 , MPI::INT, dest, 670 );
                    MPI::COMM_WORLD.Send ( &vector_trasf2[0], size_vector_trasf2, MPI::INT, dest, tag );
                    
                    
                    
                }
                
                int tag_root = chunk_size*(num_workers-1);
                tag= chunk_size*(num_workers-1);
                
                // chunk for the root
                
                ini_root_chunk= chunk_size*(num_workers-1)-lar;
                
                root_chunk_size= lar*lun - chunk_size*(num_workers-1);
                
                //  vec_proc = new int[root_chunk_size+(2*lar)];
                
                l.clear();
                
                //metto l'ultimo chunk in l e poi la prima linea
                it_send1--;
                while (it_send1->pos >=  tag-lar && it_send1!=l_tot.begin() ) it_send1--;
                it_send1++;
                
                while (it_send1->pos <  lar*lun  && it_send1!=l_tot.end() ){
                    
                    position= it_send1->pos - (tag-lar); // -(tag-lar)
                    
                    car t(position, it_send1->col);
                    l.push_back(t);
                    it_send1++;
                }
                
                it_send1=l_tot.begin();
                
                
                while (it_send1->pos <  lar && it_send1!=l_tot.end() ){
                    
                    position= it_send1->pos + root_chunk_size + lar;
                    
                    car t(position, it_send1->col);
                    l.push_back(t);
                    it_send1++;
                }
                
                 it_send2=l.begin();
                

                
            }
            else{ //se num_workers è 1, lo faccio in seriale cambiando le frontiere con sè stesso. Quindi vec_proc è uguale a vec con ultima linea all'inizio e prima linea alla fine
                
                root_chunk_size=lar*lun;
                
                l=l_tot;
                
            }
            
            
        } //end if id==0
        else if (id<num_workers)    /////  Ricevo i vari chunk i tutti i processori divversi da root
        {
            
            
            //ricevo
            MPI::COMM_WORLD.Recv ( &size_vector_trasf1, 1 , MPI::INT, 0 , 670, status );
            
            row_bef = new int[size_vector_trasf1];
            
            MPI::COMM_WORLD.Recv ( row_bef, size_vector_trasf1 , MPI::INT, 0 , MPI::ANY_TAG, status );
            tag = status.Get_tag();
            
            
            //trasformo e inserisco in l
            for (int i = 0; i<size_vector_trasf1 ; i++ ) {
                
                position= (int) row_bef[i] / 10;
                color= row_bef[i] - 10*position;
                
                car t(position, color);
                l.push_back(t);
            }
            delete[] row_bef;
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        
        
        /////   MOUVEMENT and EXPORT /////
        
        //variabili per passaggi liste
        list<car>::iterator p_contr;
        
        
        vector<int> vector_send1, vector_send2;
        vector_send1.reserve(lar);
        vector_send2.reserve(lar);
        
        int size_vec_send1, size_vec_rec1, size_vec_send2, size_vec_rec2;
        
        // variabili varie
        if(id==0) lun_chunk=(root_chunk_size/lar)+2;
        else if (id<num_workers) lun_chunk=(chunk_size/lar)+2;
        
        if (num_workers>1) {
            row_bef = new int[lar];
            row_aft = new int[lar];
            row_second = new int[lar];
            row_last = new int[lar];
        }
        
        int proc_prec, proc_succ;
        
        // scambi tra processori
        if(id==0){
            proc_succ=1;
            proc_prec=(num_workers-1);
        }
        else if (id<num_workers-1){
            proc_succ=id+1;
            proc_prec=id-1;
        }
        else if (id==num_workers-1){
            proc_succ=0;
            proc_prec=id-1;
        }
        
        
        
        // PRINT at 0 if necessary
        if(id==0){
            if(nxst==0){
                if (id==0)   cout << "EXPORT all'iterazione 0" << endl;

                // EXPORT 0.csv
                string file_exp = to_string(nxst) + ".csv";
                ofstream samefile(file_exp);
                string ss_in_string;
                

                it_next=l_tot.begin();
                while(it_next->pos < lar && it_next!=l_tot.end()) it_next++;
                
                for (int i=0; i<lun; i++) {
                    for(int j=0; j< lar ; j++){
                        position = j + lar*i;
                        if(it_next->pos == position && it_next!=l_tot.end() ){
                            ss_in_string.append( to_string( it_next->col ) );
                            it_next++;
                        }
                        else ss_in_string.append("0");
                        
                        if(j<lar-1) ss_in_string.append(",");
                    } ss_in_string.append("\n");
                }
                
                
                samefile << ss_in_string;
                samefile.close();
                s2pr++;
                nxst  = steps[s2pr];
            }
        }
        // mando a tutti il prossimo step da stampare
        MPI::COMM_WORLD.Bcast(&nxst, 1, MPI::INT, 0);
        
        car temp;
        list<car>::iterator s_row = l.begin();
        vector<bool>  spost_ini(lar, false);
        
        
        int cont(1);
        int nexti;
        
        
        ///// BEGIN MOUVEMENTS ////
        
        if (id<num_workers){
            
            while (cont<=iterazioni_max) {
                
                if (cont%2==1){
                    
                    // per ogni iterazione:
                    
                    // 1) sposto blu in ogni chunk
                    // 2) passo l'ultima riga (row_aft) al proc successivo (ultimo manda a primo)
                    // 3) ricevo vettore e lo sovrascrivo su prima riga (row_bef)
                    // 4) passo seconda riga (row_second) al proc precedente (ultimo manda a primo)
                    // 5) ricevo vettore e lo sovrascrivo su ultima riga (row_last)
                    // 6) sposto rossi in ogni chunk
                    
                    // if there are some blues
                    if(allrd==0){
                        
                        // I) BLUE
                        
                        // 1) move
                        
                        if (num_workers>1) {
                            p_contr = l.begin();
                            for (list<car>::iterator it=l.begin(); it != l.end(); it++) {
                                
                                if(it->col == 1){ // if BLU
                                    
                                    if (it->pos /lar < lun_chunk-1) { //non ultima riga   // DA OTTIMIZZARE CAMBIANDO IL FOR IN WHILE
                                        
                                        if (it->moved==true) {
                                            it->moved=false;}
                                        
                                        else{
                                            // sposto contr se è dietro it
                                            if (p_contr->pos < it->pos) {
                                                p_contr= ++it;
                                                it--;
                                            }
                                            
                                            // sposto p_contr
                                            while (p_contr->pos < (it->pos +lar) && p_contr!=l.end() ) { p_contr++; }
                                            if (p_contr->pos !=(it->pos +lar) ) {
                                                
                                                temp = *it;
                                                temp.pos+=lar; temp.moved=true;
                                                l.insert(p_contr ,temp);
                                                it = l.erase(it);
                                                it--;
                                            }
                                            
                                        }//end moved==false
                                    }//end ultima riga
                                }// end  ifblu
                            }// end for BLU
                            
                            vector_send1.clear();
                            vector_send2.clear();
                            
                            // 2) SEND penultima riga
                            
                            // trasformo in vettore
                            it_send1=l.begin();
                            while (it_send1->pos <  (lun_chunk-2)*lar && it_send1!=l.end() ) it_send1++;
                            
                            while (it_send1->pos <  (lun_chunk-1)*lar && it_send1!=l.end() ){
                                
                                position= it_send1->pos - (lun_chunk-2)*lar;
                                
                                vector_send1.push_back( 10*(position) + it_send1->col );
                                it_send1++;
                            }
                            
                            size_vec_send1= vector_send1.size();
                            
                            //mando
                            MPI::COMM_WORLD.Send ( &size_vec_send1, 1 , MPI::INT, proc_succ, 666 );
                            MPI::COMM_WORLD.Send ( &vector_send1[0], size_vec_send1, MPI::INT, proc_succ, 667 );
                            
                            
                            // 3) SEND seconda riga
                            
                            // trasformo in vettore
                            it_send2=l.begin();
                            while (it_send2->pos <  lar && it_send2!=l.end() ) it_send2++;
                            
                            while (it_send2->pos <  2*lar && it_send2!=l.end() ){
                                
                                vector_send2.push_back( 10*(it_send2->pos) + it_send2->col );
                                it_send2++;
                            }
                            
                            size_vec_send2= vector_send2.size();
                            
                            
                            //mando
                            MPI::COMM_WORLD.Send ( &size_vec_send2, 1 , MPI::INT, proc_prec, 668 );
                            MPI::COMM_WORLD.Send ( &vector_send2[0], size_vec_send2, MPI::INT, proc_prec, 669 );
                            
                            // 4) RECV in prima
                            //erase prima riga
                            while (l.begin()->pos <  lar && l.begin()!=l.end()) {
                                l.pop_front();
                            }
                            
                            //ricevo
                            MPI::COMM_WORLD.Recv ( &size_vec_rec1, 1 , MPI::INT, proc_prec , 666, status );
                            MPI::COMM_WORLD.Recv ( row_bef, size_vec_rec1 , MPI::INT, proc_prec , 667, status );
                            
                            //trasformo e inserisco all'inizio
                            for (int i = size_vec_rec1-1; i>=0 ; i-- ) {
                                position= (int) row_bef[i] / 10;
                                color= row_bef[i] - 10*position;
                                
                                car t(position, color);
                                l.push_front(t);
                            }
                            
                            // 5) RECV in ultima
                            //erase ultima riga
                            while (it_send1->pos <  lun_chunk*lar && it_send1!=l.end() ) {
                                it_send1 = l.erase(it_send1);
                            }
                            
                            //ricevo
                            MPI::COMM_WORLD.Recv ( &size_vec_rec2, 1 , MPI::INT, proc_succ , 668, status );
                            MPI::COMM_WORLD.Recv ( row_last, size_vec_rec2 , MPI::INT, proc_succ , 669, status );
                            
                            //trasformo e inserisco alla fine
                            for (int i = 0; i<size_vec_rec2 ; i++ ) {
                                position= (int) row_last[i] / 10;
                                color= row_last[i] - 10*position;
                                
                                position = position + (lun_chunk-2)*lar;
                                
                                car t(position, color);
                                l.push_back(t);
                            }
                        } //end if num_workers>1
                        
                        else{ //se il processore 0 è da solo, faccio seriale
                            
                            // riazzero spot_ini
                            for (int i=0; i<spost_ini.size(); i++)  spost_ini[i]=false;
                            
                            p_contr = l.begin();
                            for (list<car>::iterator it=l.begin(); it != l.end(); it++) {
                                
                                if(it->col == 1){ // if BLU
                                    if (it->moved==true) {
                                        it->moved=false;}
                                    else{
                                        
                                        if (it->pos /lar != lun-1) { //non ultima riga
                                            // sposto contr se è dietro it
                                            if (p_contr->pos < it->pos) {
                                                p_contr= ++it;
                                                it--;
                                            }
                                            // sposto p_contr
                                            while (p_contr->pos < (it->pos +lar) && p_contr!=l.end() ) { p_contr++; }
                                            if (p_contr->pos !=(it->pos +lar) ) {
                                                
                                                if (it->pos < lar) spost_ini[it->pos]=true; //nel vettore spostini metto true nella posizione da cui ho spostato un elemento
                                                
                                                temp = *it;
                                                temp.pos+=lar; temp.moved=true;
                                                l.insert(p_contr ,temp);
                                                it = l.erase(it);
                                                it--;
                                            }
                                        } else{ //ultima riga
                                            //aggiorno p_contr
                                            if (p_contr->pos/lar!=0) { p_contr = l.begin(); }
                                            
                                            while (p_contr->pos < (it->pos %lar) && p_contr->pos <lar  && p_contr!=l.end() ) { p_contr++; }
                                            if (p_contr->pos !=(it->pos %lar) && !spost_ini[it->pos %lar]) {
                                                
                                                temp = *it;
                                                temp.pos=temp.pos%lar;
                                                l.insert(p_contr ,temp);
                                                it = l.erase(it);
                                                it--;
                                            }
                                        }
                                        
                                    }//end moved==false
                                }// end  ifblu
                            }// end for BLU
                        } //end if there is only one proc
                        
                        
                    }//end if some blues
                    
                }
                else{ //cont%1==0
                    
                    // if there are some reds
                    if (allbl == 0) {
                        // II) RED
                        
                        int spost;
                        
                            s_row = l.begin();
                            for (list<car>::iterator it=l.begin(); it != l.end(); it++) {
                                
                                if (it->pos/lar > s_row->pos/lar) {
                                    if (s_row->moved == true) s_row->moved = false;
                                    s_row=it;
                                }
                                if(it->col == 2){
                                    // aggiornamento inizio riga
                                    if (it->pos %lar != lar-1) {  // se non fine riga
                                        it_next= ++it;
                                        it--;
                                        if (it_next->pos != it->pos +1) {
                                            if (it->pos % lar ==0) {
                                                it->moved=true;
                                            }
                                            it->pos ++;
                                        }
                                    } else{ //se fine riga
                                        if (s_row->moved == true) {
                                            s_row->moved = false;
                                        }else if (s_row->pos %lar !=0) {
                                            //sposto
                                            temp = *it;
                                            temp.pos-=temp.pos%lar;
                                            l.insert(s_row ,temp);
                                            s_row--;
                                            it = l.erase(it);
                                            it--;
                                        } //end if
                                    }//end else fine riga
                                }//end if red
                            } //end for RED
                            if (s_row->moved == true) s_row->moved = false;
                    } //end if some reds
                    
                } //end red or blue mouvement
                
                MPI_Barrier(MPI_COMM_WORLD); // serve?

                //// PRINT NOW
                if (cont == nxst) {
                    
                    if (id==0)   cout << "EXPORT all'iterazione " << cont << endl;
                    
                    // nome file in string
                    string file_exp = to_string(nxst) + ".csv";
                    
                    // nome file in char
                    filename= new char[file_exp.size()];
                    strcpy(filename, file_exp.c_str());
                    
                    
                    string ss_in_string;
                    it_next=l.begin();
                    while(it_next->pos < lar && it_next!=l.end()) it_next++;
                    
                    for (int i=1; i<lun_chunk-1; i++) {
                        for(int j=0; j< lar ; j++){
                            position = j + lar*i;
                            if(it_next->pos == position && it_next!=l.end() ){
                                ss_in_string.append( to_string( it_next->col ) );
                                it_next++;
                            }
                            else ss_in_string.append("0");
                            
                            if(j<lar-1) ss_in_string.append(",");
                        } ss_in_string.append("\n");
                    }
                    
                    
                    msgsize=ss_in_string.size();
                    
                    char ss_in_char[msgsize];
                    strcpy(ss_in_char, ss_in_string.c_str() );
                    
                    if(id==0) offset = (chunk_size*2*(num_workers-1));
                    else offset = (chunk_size*2*(id-1));
                    
                    MPI_File_open(MPI_COMM_WORLD, filename , MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
                    
                    MPI_File_seek(file, offset, MPI_SEEK_SET);
                    MPI_File_write(file, ss_in_char, msgsize, MPI_CHAR, &status1);
                    MPI_File_close(&file);
                    
                    delete[] filename;
                    
                    if(id==0){
                        s2pr++;
                        nxst  = steps[s2pr];
                    }
                    
                    MPI::COMM_WORLD.Bcast(&nxst, 1, MPI::INT, 0);
                    
                } // end print
                cont++;

            } // end while iterazioni
        } // end if id < num_workers
        
        
        
    }//END_IF_NOT-FULL/NOT-EMPTY_ln83
    
    //Free the memory
    delete[] numbers;
    
    if ( (full==0) && (empty==0)) {
        
        if (num_workers>1) {
            delete[] row_bef;
            delete[] row_aft;
            delete[] row_second;
            delete[] row_last;
        }
    }
}


void mpi_full(int argc, char **argv, int lun, int lar,string s_filename){
    
    // reinizializzazione MPi se serve
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized) MPI::Init ( argc, argv );
    
    int id = MPI::COMM_WORLD.Get_rank ( );
    int p = MPI::COMM_WORLD.Get_size ( );
    
    if(id==0) cout << "FULL" << endl;
    
    
    MPI::Status status;
    MPI_Status status1;
    MPI_Offset offset;
    MPI_File   file;
    
    // variabili COMUNICAZIONE
    int chunk_size, dest, tag, num_workers, ini_root, lun_chunk, size_steps, ini_root_chunk, root_chunk_size;
    int *numbers, *vec, *vec_proc, *row_bef, *row_second, *row_last, *row_aft, *prova, *vec_rec,*vec_send;
    
    // variabili EXPORT
    ofstream myfile;
    char *filename;
    //char *ss_in_char;
    string ss_in_string;
    int msgsize;
    
    // variabili matrice
    int dens, full, empty, allbl, allrd, iterazioni_max, s2pr(0), nxst;
    vector<int>  steps;
    numbers = new int[7];
    
    
    ////// READ MATRIX //////
    
    if(id==0){
        //inizializzo
        vec = new int[lun*lar];
        //leggo
        read_mat(s_filename, vec, numbers, steps);
        //steps
        size_steps=steps.size();
        iterazioni_max=steps.back();
        nxst=steps[0];
    }//end_if_ID0 (initialization)
    
    MPI::COMM_WORLD.Bcast(numbers, 7, MPI::INT, 0);
    
    MPI::COMM_WORLD.Bcast(&iterazioni_max, 1, MPI::INT, 0);
    
    lun=numbers[0];         lar=numbers[1];
    full  = numbers[3];     empty = numbers[4];
    allbl = numbers[5];     allrd = numbers[6];
    
    
    num_workers=p;
    if (lun<=4) num_workers=1;
    
    //se non ci son processori, finisci
    if(num_workers==0){ return;
    }
    //se 1 processore, seriale
    else if(num_workers==1){
        
        if(id==0) mpi_full_oneproc(numbers, steps, vec);
    
    }
    else{ //num_workers >1
        
        steps.push_back(iterazioni_max+1) ;
        
        ///// 3. passo i parametri a tutti i processori
        
        // chunk SIZE
        chunk_size= lun / num_workers;
        chunk_size*= lar;
        
        
        
        //if full or empty, print
        ////// dividing the matrix and broadcast of the CHUNKS
        
        if ( id == 0 ){

            //SEND chunks
            // send to 1
            dest = 1;
            tag = 0;
            
            vec_send = new int[chunk_size+(2*lar)];
            for (int i = 0; i < lar; i++ ) vec_send[i] = vec[lar*(lun-1) +i];
            for (int i = lar; i < chunk_size+(2*lar); i++ ) vec_send[i] = vec[i-lar];
            
            MPI::COMM_WORLD.Send ( vec_send, chunk_size+(2*lar), MPI::INT, dest, tag );
            
            
            //send to others
            for (int i = 2; i < num_workers; i++ )
            {
                dest = i;
                tag = chunk_size*(i-1);
                MPI::COMM_WORLD.Send ( vec+(tag-lar), chunk_size+(2*lar), MPI::INT, dest, tag );
            }
            
            int tag_root = chunk_size*(num_workers-1);
            tag= chunk_size*(num_workers-1);
            
            
            // chunk del root
            
            ini_root_chunk= chunk_size*(num_workers-1)-lar;
            root_chunk_size= lar*lun - chunk_size*(num_workers-1);
            vec_proc = new int[root_chunk_size+(2*lar)];
            int ind=0;
            for (int i = ini_root_chunk; i < lar*lun; i++ ) {
                vec_proc[ind] = vec[i];
                ind++;
            }
            for (int i = 0; i < lar; i++ ){
                vec_proc[ind] = vec[i];
                ind++;
            }
            delete[] vec_send;
            
        } //end if id==0
        ///// 5. Ricevo i vari chunk i tutti i processori divversi da root
        else if (id<num_workers)
        {
            
            vec_proc = new int[chunk_size+(2*lar)];
            
            MPI::COMM_WORLD.Recv ( vec_proc,  chunk_size+(2*lar) , MPI::INT, 0 , MPI::ANY_TAG, status );
            tag = status.Get_tag();
        }
        
        MPI::COMM_WORLD.Barrier();
        
        
        /////   MOUVEMENT and EXPORT /////
        
        //variabili per passaggi liste
        if(id==0) lun_chunk=(root_chunk_size/lar)+2;
        else if (id<num_workers) lun_chunk=(chunk_size/lar)+2;
        
        
        int proc_prec, proc_succ;
        
        row_bef = new int[lar];
        row_aft = new int[lar];
        row_second = new int[lar];
        row_last = new int[lar];
        
        // scambi tra processori
        if(id==0){
            proc_succ=1;
            proc_prec=(num_workers-1);
        }
        else if (id<num_workers-1){
            proc_succ=id+1;
            proc_prec=id-1;
        }
        else if (id==num_workers-1){
            proc_succ=0;
            proc_prec=id-1;
        }
        
        
        
        // PRINT at 0 if necessary
        if(id==0){
            if(nxst==0){
                // EXPORT 0.csv
                if (id==0)   cout << "EXPORT all'iterazione 0" << endl;
                
                string file_exp = to_string(nxst) + ".csv";
                ofstream samefile(file_exp);
                string ss_in_string;
                for (int i=0; i<lun; i++) {
                    for(int j=0; j< lar ; j++){
                        ss_in_string.append( to_string(vec[j + lar*i]) );
                        
                        if(j<lar-1) ss_in_string.append(",");
                    } ss_in_string.append("\n");
                }
                samefile << ss_in_string;
                samefile.close();
                s2pr++;
                nxst  = steps[s2pr];
            }
        }
        // mando a tutti il prossimo step da stampare
        MPI::COMM_WORLD.Bcast(&nxst, 1, MPI::INT, 0);
        
        bool moved_local(true), moved_local_blu(true), moved_local_red(true), moved_global(true);
        int cont(1);
        
        
        ///// BEGIN MOUVEMENTS ////
        
        if (id<num_workers){
            
            while (moved_global && cont<=iterazioni_max) {
                
                if (cont%2==1){
                    // per ogni iterazione:
                    
                    // 1) sposto blu in ogni chunk
                    // 2) passo l'ultima riga (row_aft) al proc successivo (ultimo manda a primo)
                    // 3) ricevo vettore e lo sovrascrivo su prima riga (row_bef)
                    // 4) passo seconda riga (row_second) al proc precedente (ultimo manda a primo)
                    // 5) ricevo vettore e lo sovrascrivo su ultima riga (row_last)
                    // 6) sposto rossi in ogni chunk
                    
                    
                    // if there are some blues
                    if(allrd==0){
                        // I) BLUE
                        
                        // 1) move
                        int nexti;
                        
                        moved_local_blu=false;
                        for (int col=0; col<lar; col++) {
                            for (int i=col; i<(lun_chunk-1)*lar; i=i+lar) { //vado da col fino alla penultima riga del vettore del chunk
                                if (vec_proc[i]==1) {
                                    nexti=i+lar;
                                    if(vec_proc[nexti]==0){
                                        vec_proc[nexti]=1;
                                        vec_proc[i]=0;
                                        i=i+lar;
                                        if(!moved_local_blu) moved_local_blu=true;
                                    }
                                }//end if blue
                            } // end row blue
                        } //end col blue
                        
                        
                        
                        // 2) SEND penultima e RECV in prima
                        for (int i = 0; i < lar; i++ ) row_aft[i]=vec_proc[i + (lun_chunk-2)*lar];
                        
                        MPI::COMM_WORLD.Sendrecv( row_aft, lar, MPI::INT, proc_succ, 666, row_bef, lar, MPI::INT, proc_prec, 666, status);
                        
                        for (int i = 0; i < lar; i++ ) vec_proc[i]=row_bef[i];
                        
                        // 3) SEND seconda e RECV in ultima
                        for (int i = 0; i < lar; i++ ) row_second[i]=vec_proc[i + lar];
                        MPI::COMM_WORLD.Sendrecv( row_second, lar, MPI::INT, proc_prec, 667, row_last, lar, MPI::INT, proc_succ, 667, status);
                        for (int i = 0; i < lar; i++ ) vec_proc[i + (lun_chunk-1)*lar]=row_last[i];
                        
                        
                        
                        
                        
                    }//end_if_NOT_ALL_RED
                    
                }
                else{   //cont%2==0
                    
                    // if there are some reds
                    if (allbl == 0) {
                        // II) RED
                        
                        int spost;
                        int nexti;
                        
                        moved_local_red=false;
                        
                        for (int row=0; row<(lun_chunk); row++) {
                            spost=false;
                            for (int i=(row*lar); i<(row*lar)+lar; i++) {
                                
                                if (vec_proc[i]==2) {
                                    
                                    nexti=i+1;
                                    
                                    if ( (i%lar)== (lar-1) ) {
                                        if(!spost) nexti=(row*lar); else {nexti=i; spost = false;}
                                    }
                                    if(vec_proc[nexti]==0){
                                        vec_proc[nexti]=2;
                                        vec_proc[i]=0;
                                        if(!moved_local_red) moved_local_red=true;
                                        if(i==(row*lar)) spost=true;
                                        if(nexti%lar!=0) i++;
                                    }
                                } //end if 2
                            } //end for i
                        } //end red
                        
                    } //end if some reds
                    
                }
                
                //// PRINT NOW
                if (cont == nxst) {
                    
                    if (id==0)   cout << "EXPORT all'iterazione " << cont << endl;
                    
                    
                    // nome file in string
                    string file_exp = to_string(nxst) + ".csv";
                    
                    
                    // nome file in char
                    filename= new char[file_exp.size()];
                    strcpy(filename, file_exp.c_str());
                    
                    string ss_in_string;
                    for (int i=1; i<lun_chunk-1; i++) {
                        for(int j=0; j< lar ; j++){
                            ss_in_string.append( to_string(vec_proc[j + lar*i]) );
                            
                            if(j<lar-1) ss_in_string.append(",");
                        } ss_in_string.append("\n");
                        
                    }
                    msgsize=ss_in_string.size();
                    
                    char ss_in_char[msgsize];
                    
                    strcpy(ss_in_char, ss_in_string.c_str() );
                    
                    
                    if(id==0) offset = (chunk_size*2*(num_workers-1));
                    else offset = (chunk_size*2*(id-1));
                    
                    
                    MPI_File_open(MPI_COMM_WORLD, filename , MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
                    
                    MPI_File_seek(file, offset, MPI_SEEK_SET);
                    MPI_File_write(file, ss_in_char, msgsize, MPI_CHAR, &status1);
                    MPI_File_close(&file);
                    
                    delete[] filename;
                    
                    if(id==0){
                        s2pr++;
                        nxst  = steps[s2pr];
                    }
                    
                    MPI::COMM_WORLD.Bcast(&nxst, 1, MPI::INT, 0);
                    
                } // end print
                
                
                
                if(cont%2==1){
                    
                    moved_local= (moved_local_blu || moved_local_red);
                    
                    if(id==0) moved_global=moved_local;
                    
                    //// gestione moved_local e moved_global
                    if(num_workers>1 ){
                        if(id>0){
                            //mando a root tutti i moved_local
                            MPI::COMM_WORLD.Send ( &moved_local, 1, MPI::BOOL, 0, 777 );
                        }
                        else if(id==0){
                            //ricevo tutti i moved_local
                            bool buffer(true);
                            for (int i=1; i<num_workers; i++) {
                                MPI::COMM_WORLD.Recv ( &buffer,  1 , MPI::BOOL, MPI::ANY_SOURCE, 777, status );
                                //aggiorno il moved_global
                                if(!buffer) moved_global= false;
                            }
                        }
                        
                        //broadcasto il moved_global
                        MPI::COMM_WORLD.Bcast(&moved_global, 1, MPI::BOOL, 0);
                    }
                    
                    if(!moved_global) break;
                }
                
                
                cont++;
            } // end while iterazioni
            
        } //end if id<num_workers
        
        MPI::COMM_WORLD.Barrier(); // serve?
        
        
        ////// other PRINTS if moved_smth==false
        
        
        
        if (!moved_global) {
            
            while(nxst<=iterazioni_max) {
                
                if (id==0)  cout << "EXPORT all'iterazione finale " << nxst << endl;
                
                // nome file in string
                string file_exp = to_string(nxst) + ".csv";
                
                
                // nome file in char
                filename= new char[file_exp.size()];
                strcpy(filename, file_exp.c_str());
                
                string ss_in_string;
                for (int i=1; i<lun_chunk-1; i++) {
                    for(int j=0; j< lar ; j++){
                        ss_in_string.append( to_string(vec_proc[j + lar*i]) );
                        
                        if(j<lar-1) ss_in_string.append(",");
                    } ss_in_string.append("\n");
                    
                }
                msgsize=ss_in_string.size();
                
                char ss_in_char[msgsize];
                strcpy(ss_in_char, ss_in_string.c_str() );
                
                if(id==0) offset = (chunk_size*2*(num_workers-1));
                else offset = (chunk_size*2*(id-1));
                
                MPI_File_open(MPI_COMM_WORLD, filename , MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
                
                MPI_File_seek(file, offset, MPI_SEEK_SET);
                MPI_File_write(file, ss_in_char, msgsize, MPI_CHAR, &status1);
                MPI_File_close(&file);
                
                delete[] filename;
                
                if(id==0){
                    s2pr++;
                    nxst  = steps[s2pr];
                }
                
                MPI::COMM_WORLD.Bcast(&nxst, 1, MPI::INT, 0);
                
                
            }//end while
        } // end print
        
        
        
        /////// Free the memory
        delete[] numbers;
        if ( id == 0 ) delete[] vec;
        delete[] vec_proc;
        delete[] row_bef;
        delete[] row_aft;
        delete[] row_second;
        delete[] row_last;
        
    }//end  num_workers>1
}


void mpi_full_oneproc(int *numbers, vector<int> &steps, int *vec){
    
    cout << "ONLY PROC 0" << endl;
    
    int chunk_size, num_workers,lun_chunk, size_steps, root_chunk_size;
    int *vec_proc;
    
    
    // variabili EXPORT
    ofstream myfile;
    char *filename;
    //char *ss_in_char;
    string ss_in_string;
    int msgsize;
    
    // variabili matrice
    int dens, full, empty, allbl, allrd, iterazioni_max, s2pr(0), nxst, lar, lun;
    
    
    size_steps=steps.size();
    iterazioni_max=steps.back();
    nxst=steps[0];
    
    lun=numbers[0];         lar=numbers[1];
    full  = numbers[3];     empty = numbers[4];
    allbl = numbers[5];     allrd = numbers[6];
    
    num_workers=1;
    chunk_size= lun / num_workers;
    chunk_size*= lar;
    
    //inizializzazione mat
    
    root_chunk_size=lar*lun;
    vec_proc = new int[root_chunk_size+(2*lar)];
    
    for (int i = 0; i < lar; i++ ) vec_proc[i]=vec[i + lar*(lun-1)];
    for (int i = lar; i < lar*(lun+1); i++ ) vec_proc[i]=vec[i-lar];
    for (int i = lar*(lun+1); i < lar*(lun+2); i++ ) vec_proc[i]=vec[i-lar*(lun+1)];
    

    lun_chunk=(root_chunk_size/lar)+2;
    
    bool moved_global(true), moved_local_blu(true), moved_local_red(true);
    
    // PRINT at 0 if necessary
    
    if(nxst==0){
        // EXPORT 0.csv
          cout << "EXPORT all'iterazione 0" << endl;
        
        string file_exp = to_string(nxst) + ".csv";
        ofstream samefile(file_exp);
        string ss_in_string;
        for (int i=0; i<lun; i++) {
            for(int j=0; j< lar ; j++){
                ss_in_string.append( to_string(vec[j + lar*i]) );
                
                if(j<lar-1) ss_in_string.append(",");
            } ss_in_string.append("\n");
        }
        samefile << ss_in_string;
        samefile.close();
        s2pr++;
        nxst  = steps[s2pr];
    }
    
    
    //MOUVEMENTS
    int cont(1);
    while (moved_global && cont<=iterazioni_max) {
        
        if (cont%2==1){
            
            if(allrd==0){
                // I) BLUE
                
                // 1) move
                int nexti;
                
                moved_local_blu=false;
                for (int col=0; col<lar; col++) {
                    for (int i=col; i<(lun_chunk-1)*lar; i=i+lar) { //vado da col fino alla penultima riga del vettore del chunk
                        if (vec_proc[i]==1) {

                            nexti=i+lar;
                            if(vec_proc[nexti]==0){
                                vec_proc[nexti]=1;
                                vec_proc[i]=0;
                                i=i+lar;
                                if(!moved_local_blu) moved_local_blu=true;
                            }
                        }//end if blue
                    } // end row blue
                } //end col blue
                

                //aggiorno ultima e prima
                for (int i = 0; i < lar; i++ ) vec_proc[i]=vec_proc[i + (lun_chunk-2)*lar];
                for (int i = 0; i < lar; i++ ) vec_proc[i + (lun_chunk-1)*lar]=vec_proc[i + lar];
                
            }//end_if_NOT_ALL_RED
            
        }
        //cont%2==0
        else{

            // if there are some reds
            if (allbl == 0) {
                // II) RED
                
                int spost;
                int nexti;
                
                moved_local_red=false;
                
                for (int row=0; row<(lun_chunk); row++) {
                    spost=false;
                    for (int i=(row*lar); i<(row*lar)+lar; i++) {
                        
                        if (vec_proc[i]==2) {
                            
                            nexti=i+1;
                            
                            if ( (i%lar)== (lar-1) ) {
                                if(!spost) nexti=(row*lar); else {nexti=i; spost = false;}
                            }
                            if(vec_proc[nexti]==0){
                                vec_proc[nexti]=2;
                                vec_proc[i]=0;
                                if(!moved_local_red) moved_local_red=true;
                                if(i==(row*lar)) spost=true;
                                if(nexti%lar!=0) i++;
                            }
                        } //end if 2
                    } //end for i
                } //end red
                
            } //end if some reds
            
        }
        
        //// PRINT NOW
        if (cont == nxst) {
            
            cout << "EXPORT all'iterazione " << cont << endl;
            
            string file_exp = to_string(nxst) + ".csv";
            ofstream samefile(file_exp);
            string ss_in_string;
            for (int i=1; i<lun+1; i++) {
                for(int j=0; j< lar ; j++){
                    ss_in_string.append( to_string(vec_proc[j + lar*i]) );
                    
                    if(j<lar-1) ss_in_string.append(",");
                } ss_in_string.append("\n");
            }
            samefile << ss_in_string;
            samefile.close();
            s2pr++;
            nxst  = steps[s2pr];
            
        } // end print
        
        if(cont%2==1){
            moved_global= (moved_local_blu || moved_local_red);
            if(!moved_global) break;
        }
        
        cont++;
    } // end while iterazioni
    
    
    if (!moved_global) {
        
        while(s2pr<size_steps) {
            
            cout << "EXPORT all'iterazione finale " << nxst << endl;
            
            string file_exp = to_string(nxst) + ".csv";
            ofstream samefile(file_exp);
            string ss_in_string;
            for (int i=1; i<lun+1; i++) {
                for(int j=0; j< lar ; j++){
                    ss_in_string.append( to_string(vec_proc[j + lar*i]) );
                    
                    if(j<lar-1) ss_in_string.append(",");
                } ss_in_string.append("\n");
            }
            samefile << ss_in_string;
            samefile.close();
            s2pr++;
            nxst  = steps[s2pr];
            
        }
    } //end final prints
    delete[] numbers;
    delete[] vec;
    
    delete[] vec_proc;
    
}




void read_mat(string s_filename, int *vec,int *numb, vector<int> &steps){
    
    
    double val(-1);
    int  e(0), r(0), c(0), stdc;
    int Full(0), Empty(0);
    int cont(0);
    int allbl(1), allrd(1);
    
    const char *filenamec = s_filename.c_str();
    
    ifstream file;
    file.open(filenamec);
    
    string      tmpline;
    if (file) {
        string       anelem;
        getline(file, tmpline);
        
        
        stringstream steps_line(tmpline);
        
        while ( getline(steps_line, anelem, ',') ) {
            //converte in int
            stringstream convertor(anelem);
            convertor >> val;
            convertor.clear();
            if (int(val) == val) /*Integer*/ steps.push_back(int(val));
            //salva
        }//end while
        
        //Sort the steps...
        sort(steps.begin(),steps.end(),[](int i, int j) -> bool {return i<j;});
        //...erase duplicates
        steps.erase( unique( steps.begin(), steps.end() ), steps.end() );
        
        
        
        
        while ( getline(file, tmpline) ) {
            //Parse all the rows
            r++;
            
            stringstream aline(tmpline);
            
            c=0;
            while ( getline(aline, anelem, ',') ) {
                c++;
                if (anelem.compare("0") == 0) {             //Empty location
                    vec[cont++] = 0;
                } else if (anelem.compare("1") == 0) {
                    e++;
                    allrd = 0;
                    vec[cont++] = 1;
                } else if (anelem.compare("2") == 0) {
                    e++;
                    allbl = 0;
                    vec[cont++] = 2;
                } else //Throw an exception
                    throw("Not allowed character");
            }//end_while_parseALine
            if (r==1) stdc = c;
            else if(stdc != c) throw("Not the same number of col");
        }//end_while_getline
        
    } else throw("no file");    //end_if_file
    if(e==0)   Empty = 1;
    if(e==r*c) Full = 1;
    
    
    
    numb[0]=r;
    numb[1]=c;
    numb[2]=e;
    numb[3]=Full;
    numb[4]=Empty;
    numb[5]=allbl;
    numb[6]=allrd;
    
    
    file.close();
    
}



void read_dim(string s_filename, int *numb , vector<int> &steps){
    
    int val(-1);
    int  e(0), r(0), c(0), stdc(0);
    int max_cont=0;
    const char *filenamec = s_filename.c_str();
    ifstream  file;
    file.open(filenamec);
    string      tmpline;
    
    if (file) {
        string       anelem;
        getline(file, tmpline);
        
        stringstream steps_line(tmpline);
        
        while ( getline(steps_line, anelem, ',') ) {
            //converte in int
            stringstream convertor(anelem);
            convertor >> val;
            convertor.clear();
            if (int(val) == val) /*Integer*/ steps.push_back(int(val));
            //salva
        }//end while
        
        //Sort the steps...
        sort(steps.begin(),steps.end(),[](int i, int j) -> bool {return i<j;});
        //...erase duplicates
        steps.erase( unique( steps.begin(), steps.end() ), steps.end() );
        
        // 1 first row
        if(getline(file, tmpline)){
            r++;
            
            stringstream first_line(tmpline);
            
            while ( getline(first_line, anelem, ',') ) {
                c++;
                if (   anelem.compare("1") == 0
                    || (anelem.compare("2") == 0))          //Blue block
                    e++;
                else if (anelem.compare("0") != 0)
                    throw("Not allowed character");
            }//end_while_parseALine
            stdc = c; //Set the reference for columns
            
        }
        
        // 2 after fisrt row to max_cont
        max_cont=min(c/5, 10);
        
        for (int i=0; i<max_cont; i++){
            if(getline(file, tmpline)){
                r++; c=0;
                stringstream aline(tmpline);
                
                while ( getline(aline, anelem, ',') ) {
                    c++;
                    if (   anelem.compare("1") == 0
                        || (anelem.compare("2") == 0))
                        e++;
                    else if (anelem.compare("0") != 0)
                        throw("Not an allowed character");
                }//end_while_parseALine
                if (stdc != c){
                    cout << "Not the same number of col" << endl;
                    throw("Not the same number of col");
                }
            }
        }
        
        // 3 after max_cont
        while ( getline(file, tmpline) ) {
            r++;
        }
        
        if(r==0) throw("no matrix");
        
    } else throw("no file");    //end_if_file
    
    numb[0]=r;
    numb[1]=c;
    numb[2]= (100*e) / ((1+max_cont)*c);
    
    
    file.close();
    
    
}

void read_mat_list(string s_filename, list<car>& l,int *numb, vector<int> &steps){
    
    
    double val(-1);
    int  e(0), r(0), c(0), stdc;
    int Full(0), Empty(0);
    int cont(0);
    int allbl(1), allrd(1);
    
    const char *filenamec = s_filename.c_str();
    
    
    ifstream file;
    file.open(filenamec);
    
    string      tmpline;
    if (file) {
        string       anelem;
        getline(file, tmpline);
        
        
        stringstream steps_line(tmpline);
        
        while ( getline(steps_line, anelem, ',') ) {
            //converte in int
            stringstream convertor(anelem);
            convertor >> val;
            convertor.clear();
            if (int(val) == val) /*Integer*/ steps.push_back(int(val));
            //salva
        }//end while
        
        //Sort the steps...
        sort(steps.begin(),steps.end(),[](int i, int j) -> bool {return i<j;});
        //...erase duplicates
        steps.erase( unique( steps.begin(), steps.end() ), steps.end() );
        
        
        
        
        while ( getline(file, tmpline) ) {
            //Parse all the rows
            r++;
            
            stringstream aline(tmpline);
            
            c=0;
            while ( getline(aline, anelem, ',') ) {
                c++;
                if (anelem.compare("0") == 0) {             //Empty location
                    cont++;
                } else if (anelem.compare("1") == 0) {
                    e++;
                    allrd = 0;
                    car t(cont++, 1);
                    l.push_back(t);
                } else if (anelem.compare("2") == 0) {
                    e++;
                    allbl = 0;
                    car t(cont++, 2);
                    l.push_back(t);
                } else //Throw an exception
                    throw("Not allowed character");
            }//end_while_parseALine
            if (r==1) stdc = c;
            else if(stdc != c) throw("Not the same number of col");
        }//end_while_getline
        
    } else throw("no file");    //end_if_file
    if(e==0)   Empty = 1;
    if(e==r*c) Full = 1;
    
    
    
    numb[0]=r;
    numb[1]=c;
    numb[2]=e;
    numb[3]=Full;
    numb[4]=Empty;
    numb[5]=allbl;
    numb[6]=allrd;
    
    
    file.close();
    
}


