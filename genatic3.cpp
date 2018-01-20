

#include<iostream>
#include<bits/stdc++.h>
#include"genatic3.h"
#include<algorithm>
#define NP 100
#define n 100












cross::cross(){

}


float cross ::getrand(){
    float minx=0.0;
    float maxx=1.0;
    float miny=0.0;
    float maxy=1.0;
    float cx= (0.9*(maxx-minx))+ static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(((minx+0.1*(maxx-minx))-(0.9*(maxx-minx)))+1)));
    float cy=(0.9*(maxy-miny)) + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((miny+0.1*(maxy-miny))-(0.9*(maxy-miny)))));
    float sigmax=-1;
    float sigmay=1;
    float qxy=-0.2 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.2-(-0.2))));
    int qy=-1;
    int qx=1;


    return cx,cy,qx,qy,qxy,sigmay,sigmax;



}


void cross ::gets1(){

        float cx,cy,qx,qy,qxy,sigmay,sigmax =getrand();


        float i,j=0;
        int k,h;

        while(i<1){
k=i*10;
            while(j<1){
h=j*10;
       this->s1[k][h]=0.5*(1+tanh((qx*(i-cx)/sigmax)+(qy*(j-cy)/sigmay)));

        }

            j+=0.1;

        }
    i+=0.1;


}




void cross::gets2(){
    float cx,cy,qx,qy,qxy,sigmay,sigmax =getrand();
    float i,j=0;
    int k,h;
    while(i<1){
k=i*10;

        while(j<1){

h=j*10;




    this->s2[k][h]=0.5*(1+tanh((qx*(i-cx)/sigmax)+((qxy*(i-cx)*(j/10-cy))/(sigmax*sigmay))+(qy*(j-cy)/sigmay)));



    j+=0.1;

}
   i+=0.1;

    }
}

   void cross::gets3(){

       float cx,cy,qx,qy,qxy,sigmay,sigmax =getrand();

   float i=0;
   float j=0;
   int k,h;
   while(i<1){
       k=i*10;
       while(j<1){
           h=j*10;


       this->s3[k][h]=(0.5*((1+tanh(qx*(i-cx)/sigmax))+tanh(qy*(j-cy)/sigmay)));


           j+=0.1;

       }
 i+=0.1;
}


   }


   void cross::getchild(float p1[sz][sz],float p2[sz][sz]){
       for(int i=0;i<9;i++){
           for(int j=0;j<9;j++){


           this->child1[i][j]=((s1[i][j])*p1[i][j])+((1-(s1[i][j])*p2[i][j]));
           this->child2[i][j]=(1-(s1[i][j])*p1[i][j])+ (s1[i][j])*p2[i][j];

   }
   }


       for(int i=0;i<9;i++){
          for(int j=0;j<9;j++){
           this->child3[i][j]=(s2[i][j])*p1[i][j]+((1-(s2[i][j])*p2[i][j]));
           this->child4[i][j]=((1-(s2[i][j])*p1[i][j])+ (s2[i][j])*p2[i][j]);
   }

      }

       for(int i=0;i<9;i++){
          for(int j=0;j<9;j++){

              this->child5[i][j]=((s3[i][j])*p1[i][j])+((1-(s3[i][j])*p2[i][j]));
              this->child6[i][j]=((1-(s3[i][j])*p1[i][j])+ (s3[i][j])*p2[i][j]);

          }
}




    }





   void cross::cross_child(float p1[sz][sz],float p2[sz][sz]){

       for(int i=0;i<NP/4;i++){

           this->getchild(p1,p2);



     for(int j=0;j<9;j++){

         for(int k=0;k<9;k++){

             (this->cross_children)[j][k][i]=this->child1[j][k];

             (this->cross_children1)[j][k][i]=this->child3[j][k];






             (this->cross_children2)[j][k][i]=this->child5[j][k];



         }


     }


   }
}









   imagration::imagration():cross(){

   }
   float imagration::getrandom(){
       int maxx=1.0;
       int mn=-1.0;
       float qxy2=-1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1-(-1))));
       //float qxy2=rand()%(maxx-mn+1)+min;

       float mcd=-0.25 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.75-(-0.25))));


       //float mcd=rand() %(max2-mn2+1.0)+mn2;
       float mdc4=rand() %1;
       float cx,cy,qx,qy,qxy,sigmay,sigmax =getrand();
       return cx,cy,qx,qy,qxy2,sigmay,sigmax,mdc4;

   }



   void imagration::getG1(){
       float cx,cy,qx,qy,qxy2,sigmay,sigmax,mdc4=getrandom();


       float i=0;
       float j=0;
       int k,h;
       while(i<1){
           k=i*10;
           while(j<1){
               h=j*10;


    float mcd=-0.25 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.75-(-0.25))));
    this->G1[k][h]=mcd+0.5*(1+tanh((qx*(i-cx)/sigmax)+((qxy2*(i-cx)*(j-cy))/sigmax*sigmay)+qy*(j-cy)/sigmay));

               j+=0.1;

           }
     i+=0.1;
    }


   }





   void imagration::getG2(){
       float cx,cy,qx,qy,qxy2,sigmay,sigmax,mdc4=getrandom();
       float i=0;
       float j=0;
       int k,h;
       while(i<=1){
           k=i*10;
           while(j<=1){
               h=j*10;
               float mdc=-0.25 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.75-(-0.25))));
               this->G2[k][h]=mdc+0.5*(1+tanh(((qx*(i-cx))/sigmax)))*tanh(qy*(j-cy)/sigmay);









               j+=0.1;

           }
     i+=0.1;
    }


   }



       void imagration::getG4(){

           float cx,cy,qx,qy,qxy2,sigmay,sigmax,mdc4=getrandom();
           float i=0;
           float j=0;
           int k,h;
           while(i<=1){
               k=i*10;
               while(j<=1){
                   h=j*10;
                   float mdc4= static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1)));


                   this->G4[k][h]=mdc4;







                   j+=0.1;

               }
         i+=0.1;
        }


       }

void imagration::getO(float OI[50][50]){
//oi is the best fit individual so far



     float f= (-(this->residue)) + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((this->residue)-(-(this->residue)))));
     float i=0;
     float j=0;
     int k,h;
     while(i<1){
         k=i*10;
         while(j<1){
             h=j*10;



             this->Oiemgrants[k][h]=OI[k][h]+(f*(this->G1)[k][h]);


 j+=0.1;

         }
   i+=0.1;
  }


 }



void imagration::getO1(float OI[50][50]){
//oi is the best fit individual so far



     float f= (-(this->residue)) + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((this->residue)-(-(this->residue)))));
     float i=0;
     float j=0;
     int k,h;
     while(i<1){
         k=i*10;
         while(j<1){
             h=j*10;



 this->Oimegrats2[k][h]=OI[k][h]+(f*(this->G2)[k][h]);

 j+=0.1;

         }
   i+=0.1;
  }


 }

void imagration::getO2(float OI[50][50]){
//oi is the best fit individual so far



     float f= (-(this->residue)) + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((this->residue)-(-(this->residue)))));
     float i=0;
     float j=0;
     int k,h;
     while(i<1){
         k=i*10;
         while(j<1){
             h=j*10;



             this->Oimegrants3[k][h]=OI[k][h]+(f*(this->G4)[k][h]);


 j+=0.1;

         }
   i+=0.1;
  }


 }





 float imagration::getrand3(){
     float xmin=0;
     float xmax=1;
     float ymin=0;
     float ymax=1;


     float k = 1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((NP/5)-(1))));

     float B1 = -2 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((2)-(-2))));
     float B2= -2 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((2)-(-2))));
     float B3 = -0.25 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.25)-(-0.25))));
     float cxx = xmin + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((xmax)-(xmin))));
     float cyy=ymin + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((ymax)-(ymin))));

     float sigmaxx =  static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(xmax-xmin)));

     float sigmayy = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(ymax-ymin)));
     float l = -0.5 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.5)-(-0.5))));

     float qxx =  static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));

     float qyy =  static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
     return k,B1,B2,B3,sigmaxx,sigmayy,l,qxx,qyy,cxx,cyy;



 }

 void imagration::getn1(){



     float i=0;
     float j=0;
     int k,h;
     while(i<1){
         k=i*10;
         while(j<1){
             h=j*10;
             float nn=0;
              float kk,B1,B2,B3,sigmaxx,sigmayy,l,qxx,qyy,cxx,cyy=getrand3();
                 for(int ll=0;ll<kk;ll++){
                      kk,B1,B2,B3,sigmaxx,sigmayy,l,qxx,qyy,cxx,cyy=getrand3();

                     nn+=B1*exp(-(1/2*(1-pow(l,2)))*((pow(((i)-cxx),2)/pow(sigmaxx,2))-(2*l*((i)-cxx)*((j)-cyy))/(sigmaxx*sigmayy)+(pow((j)-cyy,2)/pow(sigmayy,2))));
                       }





         (this->n1)[k][h]=nn;

         }

                 j+=0.1;

         }
   i+=0.1;
  }




     void imagration::getn2(){
         int kk=NP/5;
         float i=0;
         float j=0;
         int k,h;
         while(i<1){
             k=i*10;
             while(j<1){
                 h=j*10;
                  float n2c=0;
                  for(int LL=0;LL<kk;LL++){

                   float kk,B1,B2,B3,sigmaxx,sigmayy,l,qxx,qyy,cxx,cyy=getrand3();
                   n2c+=(B2*0.5*(1+tanh((qxx*(i-cxx)/sigmaxx)+((2*l*(i-cxx)*(j-cyy))/sigmaxx*sigmayy)+(qyy*(j-cyy)/sigmayy))));



}

(this->n2)[k][h]=n2c;
                 j+=0.1;

             }
       i+=0.1;
      }





}

void imagration::getn3(){
    float kk,B1,B2,B3,sigmaxx,sigmayy,l,qxx,qyy,cxx,cyy=getrand3();
    float i=0;
    float j=0;
    int k,h;
    while(i<1){
        k=i*10;
        while(j<1){
            h=j*10;

n3[k][h]=B3+ 0.5*(1+ tanh(qxx*(i-cxx)/sigmaxx)*tanh(qyy*(j-cyy)/sigmayy));





            j+=0.1;
}
  i+=0.1;
 }

}




void imagration::getn4(){
    for(int i=0;i<9;i++){


        for(int j=0;j<9;j++){





            (this->n4)[i][j]=(this->n2)[i][j]+(this->n3)[i][j]+(this->n1)[i][j];


}


 }

}



void imagration::bypass(float OI[sz][sz]){
    for(int i=0;i<NP/5;i++){
         this->getO(OI);



        for(int j=0;j<9;j++){
            for(int k=0;k<9;k++){



                (this->imgration_children)[j][k][i]=this->Oiemgrants[j][k];


            }


        }
}
}
        //this->getO1(OI);
        //this->getO2(OI);





void imagration::bypass1(float OI[sz][sz]){
    for(int i=0;i<NP/5;i++){
         this->getO1(OI);
        for(int j=0;j<9;j++){
            for(int k=0;k<9;k++){



                (this->imgration_children2)[j][k][i]=this->Oimegrats2[j][k];


            }


        }
}
}

void imagration::bypass2(float OI[sz][sz]){
    for(int i=0;i<NP/5;i++){
         this->getO2(OI);
        for(int j=0;j<9;j++){
            for(int k=0;k<9;k++){



                (this->imgration_children3)[j][k][i]=this->Oimegrants3[j][k];


            }


        }
}
}

//  imagration::imagration():cross(){


menah::menah():imagration(){


}


float menah::eval_fitness(float residue){







    float fitness=1/(1+residue);


    return fitness;

}



float menah::calc_sum_res(int k, float arr[sz][sz][gen_sz]){
    float sum =0.0 ;
    for(int j=0; j<=9; j++)
    {
        for(int i=0; i<=9; i++)
        {
            sum +=pow(arr[j][i][k],2);
        }
    }
    sum =sqrt(sum);
    return sum;



}





float menah:: rand_between_neg_pos(int num){

int randomNum = rand() %((2*num)+1) + (-num);
return randomNum;


}

float menah::RandomFloat(float min, float max){

    assert(max > min);
    float random = ((float) rand()) / (float) RAND_MAX;


    float range = max - min;
    return (random*range) + min;

 }


float menah::c_x(void){
    float cx =rand()/(double)RAND_MAX;
    return cx;
}
float menah::sigma_x(void){
    float sigma_x =rand()/(double)RAND_MAX;
    return sigma_x ;}

float menah::c_y (void ){
   float  c_y =rand()/(double)RAND_MAX;
    return c_y ;
}
float menah::sigma_y(void) {
    float sigma_y =rand()/(double)RAND_MAX;
    return sigma_y;
}
float menah::q_x (void){

    return static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
}
float menah::q_y (void){ return static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));}
float res;

float menah::J(){ return 0.5 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.5)-(-0.5))));}






//float exact_u[sz][sz] ;
//float exact_v[sz][sz];
float menah::fx(float x, float y){
    float val=((8*pow(x,2)*pow((x-1),2)*(y-1))
               +(4*pow(x,2)*((2*y)-1)*pow((x-1),2))+(8*pow(x,2)*y*pow((x-1),2))
               +(4*y*((2*y)-1)*pow((x-1),2)*(y-1))+(4*pow(x,2)*y*((2*y)-1)*( y - 1))
               +( 8* x * y *(( 2 * x ) - 2 )*(( 2 * y )-1)*( y - 1 )));

    return val;


 }
float menah::generate_soll_1(float x, float y){
    float B1 =RandomFloat(-2,2) ;

    float c1 = ( (B1) * exp( ((-1)/( (2) * ( 1 - pow( J(), 2 ) ) )) * ( ( pow( ( x - c_x() ), 2 ) / ( pow( sigma_x(), 2 ) ) ) - ( ( 2 * J() * ( x - c_x() ) * ( y - c_y() ) ) / ( sigma_x() * sigma_y() ) ) + ( ( y - c_y() ) / pow( sigma_y(), 2 ) ) ) ) );

    return c1;


 }

float menah::generate_soll_2(float x, float y)
{

    float B2 =RandomFloat(-2,2) ;

    float c2 = (B2) * (0.5) * (1 + tanh( ( ( (q_x()) * ( x - c_x() ) ) / ( sigma_x() ) ) +( ( 2 * J() * ( x - c_x() ) * ( y - c_y() ) ) / ( sigma_x() * sigma_y() ) ) + ( ( q_y() * ( y - c_y() ))/(sigma_y()))));

    return c2 ;
}



float menah::generate_soll_3(float x, float y)
{
    float B3= RandomFloat(-0.25,0.25);
    float c3 = (B3) + ( (0.5)*((1) + ( tanh( ((q_x()) * (x - c_x())) / ( sigma_x()) ) * ( tanh(((q_y())*(y - c_y())) / (sigma_y())) ) )) );
    return c3;
}


float menah::generate_soll_4(float x, float y)
{
    float c4 =generate_soll_1(x,y) +generate_soll_2(x,y)+generate_soll_3(x,y);
    return c4;
}



mutation::mutation():menah(){


}
float mdc2 =-0.25 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.75)-(-0.25))));
float mutation::li(){return rand() / double(RAND_MAX);}

float qx =1;
float qy =-1;
float qxy =-0.2 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.2)-(-0.2))));
float mdc1 = -0.25 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.75)-(-0.25))));
//float mdc2 =-0.25 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.75)-(-0.25))));
float mdc3  =-0.5 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.5)-(-0.5))));
int xmin=0;
int xmax=1;
int delta_x =xmax-xmin ;



float mutation::cx(){return rand() / double(RAND_MAX);}
float mutation::cy(){return rand() / double(RAND_MAX);}
float mutation::sigmax(){return rand() / double(RAND_MAX);}
float mutation::sigmay(){return rand() / double(RAND_MAX);}
//float mutation::w =RandomFloat(-int rave,int rave);




float mutation::m1(float x, float y){


    float val_1 =mdc1 + exp((-1/1-li())*((pow(2,(x-cx()))/sigmax()*sigmax())
                                    -((-2l*(x-cx())*(y-cy()))/sigmax()*sigmay())+((pow(2,(y-cy())))/sigmay()*sigmay())));
    return val_1;
    }


//this is different form my mixing functions because it takes x ,y discrete values not a loop and returns a sigle value at a time so we will call it in a loop to fill the 2 dimensional grid of the mixing function
 float mutation::m2(float x, float y){



    float val_2=mdc2
             + 0.5*(1+tanh((qx*(x -cx())/sigmax())+(qxy*(x-cx())*((y-cy())/sigmax()*sigmay()))+(qy*(x -cy())/sigmay()))) ;

     return val_2;

 }

float mutation::m3(){
    float mdc3=-0.5 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((0.5)-(-0.5))));




    float val_3 = mdc3;
    return val_3 ;
}


float mutation::w(float ravg){

    float w=-ravg + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/((ravg)-(-ravg))));
    return w;

}




void mutation::eval_m(){
    for(float j=0;j<=1;j+=0.1){
     for(float i =0;i<= 1;i+=1){
         int x=j*10;
         int y =i*10;
         float mutation_val = m1(i,j) ;
         this->mutations_ij[x][y] =mutation_val ;
     }
    }
}



void mutation::eval_m1(){
    for(float j=0;j<=1;j+=0.1){
     for(float i =0;i<= 1;i+=1){
         int x=j*10;
         int y =i*10;
         float mutation_val = m2(i,j) ;
         this->mutations_ij1[x][y] =mutation_val ;
     }
    }
}




void mutation::eval_m2(){
    for(float j=0;j<=1;j+=0.1){
     for(float i =0;i<= 1;i+=1){
         int x=j*10;
         int y =i*10;
         float mutation_val = m3() ;
         this->mutations_ij2[x][y] =mutation_val ;
     }
    }
}
//gets muataion child which is similar as my omigrents geto
void mutation:: mutate_getchild(float mutation_curve[sz][sz],float ravg){

    //mutation curve which is the n1,n2,n3 we take it as input



    //this is not similar as mine its a generic one in which we will call sevreal time to get output of mutation children one time with mut_c as n1
    //the other time as n2

    for (int i =0;i<=9;i++){
        for(int j =0;j<=9 ;j++){
           mutation_curve[j][j]=mutation_curve[j][i]+(w(ravg)*(this->mutations_ij[i][j]));
        }
    }
  }









