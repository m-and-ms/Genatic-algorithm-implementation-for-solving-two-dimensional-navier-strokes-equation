#ifndef GENATIC3_H
#define GENATIC3_H
const int  NP= 100;
const int n= 100;
const int sz= 50;
#include <cstdlib>
#include <ctime>
#include<math.h>
#include<bits/stdc++.h>
#include<vector>
#include<iostream>
#include<random>
#define pop_size 100
#define gen_sz 100


using namespace std ;




class cross{

public:
    cross();

    float getrand();
    float s1[50][50];

    float s2[sz][sz];

    float s3[sz][sz];

    float child1[sz][sz];

    float child2[sz][sz];

    float child3[sz][sz];

    float child4[sz][sz];

    float child5[sz][sz];

    float child6[sz][sz];




    float cross_children[sz][sz][100];
    float cross_children1[sz][sz][100];
    float cross_children2[sz][sz][100];
    void cross_child(float p1[sz][sz],float p2[sz][sz]);
    void cross_child1(float p1[sz][sz],float p2[sz][sz]);

    void cross_child2(float p1[sz][sz],float p2[sz][sz]);



    float passing_children[sz][sz][sz];


    void gets3();
    void gets2();
    void gets1();
    void getchild(float p1[sz][sz],float p2[sz][sz]);






 };




class imagration:public cross{

public:
    imagration();
void getG1();
void getG2();
void getG4();
float getrandom();
float G1[sz][sz];
float G2[sz][sz];
float G4[sz][sz];
void getO(float OI[sz][sz]);
void getO1(float OI[sz][sz]);
void getO2(float OI[sz][sz]);
float residue;

float Oimegrats2[sz][sz];
float Oimegrants3[sz][sz];
float Oiemgrants[sz][sz];
float best_fit[sz][sz];






 float getrand3();












 float imgration_children[sz][sz][sz];
 float imgration_children2[sz][sz][sz];
 float imgration_children3[sz][sz][sz];






void terminaate(float bestfit[sz][sz],float fitthreshold,float);






void bypass(float OI[sz][sz]);
void bypass1(float OI[sz][sz]);

void bypass2(float OI[sz][sz]);


void getn1();
   void getn2();
   void getn3();
   void getn4();

   float n1[sz][sz];
   float n2[sz][sz];
   float n3[sz][sz];
   float n4[sz][sz];














};



class menah:public imagration{

public:

    menah();

   float eval_fitness(float);
   float calc_sum_res(int ,float [sz][sz][gen_sz]);
   float RandomFloat(float,float);
   float fx(float,float);
   float rand_between_neg_pos(int num);


   float generate_soll_1(float x,float y);
   float generate_soll_2(float x,float y);
   float generate_soll_3(float x,float y);
   float generate_soll_4(float x,float y);
// u soll is for all the generation similar to the children arrays at mine
   float u_soll1[sz][sz][gen_sz];
   float v_soll1[sz][sz][gen_sz];
   float res_soll1[sz][sz][gen_sz];

   float u_soll2[sz][sz][gen_sz];
   float v_soll2[sz][sz][gen_sz];
   float res_soll2[sz][sz][gen_sz];

   float u_soll3[sz][sz][gen_sz];
   float v_soll3[sz][sz][gen_sz];
   float res_soll3[sz][sz][gen_sz];

   float u_soll4[sz][sz][gen_sz];
   float v_soll4[sz][sz][gen_sz];
   float res_soll4[sz][sz][gen_sz];
   float fitness_soll1 ;
   float fitness_soll2 ;
   float fitness_soll3 ;
   float fitness_soll4 ;
   float total_fitness ;
   float raverage();

   float J();
   void formgenrate();
   float generations[sz][sz][gen_sz];





   float c_x(void);

   float sigma_x(void);
   float c_y (void );
   float sigma_y(void) ;
   float q_x (void);
   float q_y (void);
   float res;




};


class mutation:public menah{
public :
    mutation();



    float li ();

    float cx();
    float cy();
    float sigmax();
    float sigmay();
    float w(float);























   //us are the n1,n2,n3,n4 at mine

       float rave;
        float u_soll1[sz][sz];
          float u_soll2[sz][sz];
            float u_soll3[sz][sz];
              float u_soll4[sz][sz];

              float mutations_ij[sz][sz];// is the value of the muttation mixing curve

              float mutations_ij1[sz][sz];

              float mutations_ij2[sz][sz];





              float m1(float x,float y);




              float m2(float x,float y);
              float m3();
              void eval_m();
              void eval_m1();

              void eval_m2();

              void mutate_getchild(float arr[sz][sz],float);




};




















#endif // GENATIC3_H

