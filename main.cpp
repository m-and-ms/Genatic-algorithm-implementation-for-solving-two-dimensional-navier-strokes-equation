
#include"genatic3.h"
using namespace std;



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
float mennah[sz][sz][20];
float fitness_soll1 ;
float fitness_soll2 ;
float fitness_soll3 ;
float fitness_soll4 ;
float total_fitness ;
float exact_u[sz][sz] ;
float exact_v[sz][sz];

/**
float RandomFloat(float min, float max)
{

    assert(max > min);
    float random = ((float) rand()) / (float) RAND_MAX;


    float range = max - min;
    return (random*range) + min;
}

float generate_soll_1(float x,float y)
{
    float B1 =RandomFloat(-2,2) ;

    float c1 = ( (B1) * exp( ((-1)/( (2) * ( 1 - pow( J, 2 ) ) )) * ( ( pow( ( x - c_x ), 2 ) / ( pow( sigma_x, 2 ) ) ) - ( ( 2 * J * ( x - c_x ) * ( y - c_y ) ) / ( sigma_x * sigma_y ) ) + ( ( y - c_y ) / pow( sigma_y, 2 ) ) ) ) );

    return c1;
}
float generate_soll_2(float x,float y)
{

    float B2 =RandomFloat(-2,2) ;

    float c2 = (B2) * (0.5) * (1 + tanh( ( ( (q_x) * ( x - c_x ) ) / ( sigma_x ) ) +( ( 2 * J * ( x - c_x ) * ( y - c_y ) ) / ( sigma_x * sigma_y ) ) + ( ( q_y * ( y - c_y ))/(sigma_y))));

    return c2 ;
}

float generate_soll_3(float x, float y )
{
    float B3= RandomFloat(-0.25,0.25);
    float c3 = (B3) + ( (0.5)*((1) + ( tanh( ((q_x) * (x - c_x)) / ( sigma_x) ) * ( tanh(((q_y)*(y - c_y)) / (sigma_y)) ) )) );
    return c3;
}

float generate_soll_4(float x,float y)
{
    float c4 =generate_soll_1(x,y) +generate_soll_2(x,y)+generate_soll_3(x,y);
    return c4;
}



**/








int main()
{

   menah m1;
  mutation  mutate1;
    cross c1solution;
    imagration img1;

    cout<<"menaa ya 3beeta"<<endl;
    cout<<"may gamda gdn :)"<<endl;

    for(int i=0; i<100; i++)
    {
        for(float y=0.1; y<1; y +=0.1)
        {
            for(float x=0.1; x<1; x +=0.1)
            {

                int k =x*10 ;
                int j =y*10;
                float u =-2*x*y*(x-1)*(y-1)*x*(x-1)*(2*y -1);
                float v =2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x -1);
                exact_u[j][k] =u ;
                exact_v[j][k] =v ;
            }
        }
    }

    for(int i=0; i<100; i++)
    {
        for(float y=0.1; y<1; y +=0.1)
        {
            for(float x=0.1; x<1; x +=0.1)
            {



                int k =x*10 ;
                int j =y*10;




                u_soll1[j][k][i] =m1.generate_soll_1(x,y);


                u_soll2[j][k][i] =m1.generate_soll_2(x,y);

                u_soll3[j][k][i] =m1.generate_soll_3(x,y);

                u_soll4[j][k][i] =m1.generate_soll_4(x,y);
            }
        }
    }


    for(int i=0; i<100; i++)
    {
        for(float y=0.1; y<1; y +=0.1)
        {
            for(float x=0.1; x<1; x +=0.1)
            {



               int k =x*10 ;
                int j =y*10;
                //cout<<exact_u[j][k]-(u_soll1[j][k-1]+u_soll1[j][k+1]+u_soll1[j+1][k]+u_soll1[j-1][k]+(fx(x,y)/4)*pow(0.1,2))<<endl;
                res_soll1[j][k][i] =exact_u[j][k]-(u_soll1[j][k-1][i]+u_soll1[j][k+1][i]+u_soll1[j+1][k][i]+u_soll1[j-1][k][i]+(m1.fx(x,y)/4)*pow(0.1,2))                                        ;
                //cout<<res_soll1[j][k][i]<<endl;
                res_soll2[j][k][i] =exact_u[j][k]-(u_soll2[j][k-1][i]+u_soll2[j][k+1][i]+u_soll2[j+1][k][i]+u_soll2[j-1][k][i]+(m1.fx(x,y)/4)*pow(0.1,2))
                                ;
                res_soll3[j][k][i] =exact_u[j][k]-(u_soll3[j][k-1][i]+u_soll3[j][k+1][i]+u_soll3[j+1][k][i]+u_soll3[j-1][k][i]+(m1.fx(x,y)/4)*pow(0.1,2))
                                 ;
                res_soll4[j][k][i] =exact_u[j][k]-(u_soll4[j][k-1][i]+u_soll4[j][k+1][i]+u_soll4[j+1][k][i]+u_soll4[j-1][k][i]+(m1.fx(x,y)/4)*pow(0.1,2))
                                 ;

            }
        }

    }
    cout<<u_soll1[9][9][1]<<endl;
    cout<<u_soll2[9][9][1]<<endl;
    cout<<res_soll3[9][9][1]<<endl;
    cout<<res_soll4[9][9][1]<<endl;
    cout<<res_soll4[9][9][1]<<endl;

for(int i=0 ;i<gen_sz/4;i++){


        float res_1 =2*(m1.calc_sum_res(i,res_soll1));
        float res_2 =2*(m1.calc_sum_res(i,res_soll2));
        float res_3 =2*(m1.calc_sum_res(i,res_soll3));
        float res_4 =2*(m1.calc_sum_res(i,res_soll4));
        float res=res_1+res_2+res_3+res_4;

 fitness_soll1 =m1.eval_fitness(res_1);
 fitness_soll2 =m1.eval_fitness(res_2);
 fitness_soll3 =m1.eval_fitness(res_3);

 fitness_soll4 =m1.eval_fitness(res_4);
 total_fitness  =fitness_soll1 +fitness_soll2 +fitness_soll3+fitness_soll4;
 cout<<"hi"<< " "<<total_fitness<<"d" ;
 float most_fit=0;
 for(int k=0;k<=4;k++){
     most_fit =fitness_soll1 ;
     if( fitness_soll2>most_fit){
       most_fit =fitness_soll2;

     }
       if( fitness_soll3>fitness_soll2){
          most_fit =fitness_soll3;
       }
       if( fitness_soll4>fitness_soll3){
          most_fit =fitness_soll4;
 }
 }





float p1[sz][sz];
 float p2[sz][sz];



for(int j=0;j<9;j++){
    for(int k=0;k<9;k++){


 float p1[sz][sz];
  float p2[sz][sz];


p1[j][k]=u_soll1[j][k][i] ;

 p2[j][k]=u_soll2[j][k][i] ;



    c1solution.cross_child(p1,p2);



u_soll1[j][k][i] =(c1solution.cross_children)[j][k][i];

u_soll2[j][k][i]=(c1solution.cross_children1)[j][k][i];
u_soll3[j][k][i]=(c1solution.cross_children2)[j][k][i];

}
}
}          //kda el members ra7o el 3 -d array el esmo cross children1
            //cross children 2 and cross children





//segementation fault happens due to too much memory usage  :)


}



