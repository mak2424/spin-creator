#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <random>
#include <ctime>

//#define DBL_EPSILON 2.2204460492503131e-16
#define DBL_EPSILON 1e-12
//typedef unsigned long long ulong;

using namespace std;

bool is_equal(double a, double b) { return fabs(a-b) <= DBL_EPSILON * fmax(fabs(a),fabs(b)); }

#define NumMC 30000 //* количество спинов, перенести!!!!!!!!!!!!!!!!!!!!

#define Tmin 0.0001
#define Tmax 10.01
float Tstep = 0.01;

double C1=0;


int main()
{
   srand(5);
   default_random_engine generator; //тип генератора случайных чисел
   
   ofstream foutc("coord.dat"); //координаты точек
   ofstream fouts("spins.dat"); //магнитные моменты
   ofstream fout_centers("centers.dat"); //центры гексагонов
   
   vector <double> coorm; // array for scaled lattice
   vector <double> coor(24);//array of coordinates of unit cell число координат
   vector <double> mm(24); //array of magnetic moments of unit cel  число магнитных моментов
   vector <double> mxmy; //array of magnetic moments of scaled lattice
   vector <double> hex_centers(8);
   vector <double> temp_vector;
   
   double sqrt_3 = sqrt(3);
   double a = sqrt_3/2; //параметр решетки
   
   unsigned int total_num_of_spins_in_core = 6; //количество спинов в ядре
   unsigned int total_num_of_spins_in_boundaries = 24; //количество спинов на границе = 24
   unsigned int total_num_of_neighbors = 14; //количество соседей для каждого спина (для проверки) = 14
   
   double interaction_radius_squared_for_core = 3/4. + 0.0001; //радиус взаимодействия в квадрате (от центра гексагона), 
   //используется для определения ядра и границы
   double interaction_radius_squared_for_border = 27/4. + 0.0001; //(2*a)^2, радиус взаимодействия в квадрате (от центра гексагона), 
   //используется для определения границы (+ 0.0001)
   ///double interaction_radius_squared_for_border = 3. - 0.0001; 
   double interaction_radius_squared_for_neighbors = 3. + 0.0001; //(2*a)^2, радиус взаимодействия в квадрате (от спина до спина), 
   //используется для определения соседей
   ///double interaction_radius_squared_for_neighbors = 3/4. + 0.0001;
   cout << "interaction_radius_squared = " << interaction_radius_squared_for_core << endl;
   
   double period_X=4*a;
   double period_Y=3;
   
   
   hex_centers[0]=sqrt_3/2;        hex_centers[1]=0.75;
   hex_centers[2]=(3*sqrt_3)/2;    hex_centers[3]=0.75;
   hex_centers[4]=sqrt_3;          hex_centers[5]=2.25;
   hex_centers[6]=2*sqrt_3;        hex_centers[7]=2.25;
   
   
   //hc
   //      X            Y
   //1.73205 = sqrt_3
   //1.29904 = (3*sqrt_3)/4
   //0.433013 = sqrt_3/4
   //0.866025 = sqrt_3/2
   //2.59808 = (3*sqrt_3)/2
   //3.03109 = (7*sqrt_3)/4
   //2.16506 = (5*sqrt_3)/4
   
   
   coor[0]=sqrt_3;             coor[1]=0.75;
   coor[2]=(3*sqrt_3)/4;       coor[3]=1.5;
   coor[4]=sqrt_3/4;           coor[5]=1.5;
   coor[6]=0;                  coor[7]=0.75;
   coor[8]=sqrt_3/4;           coor[9]=0;
   coor[10]=(3*sqrt_3)/4;      coor[11]=0;
   
   coor[12]=(7*sqrt_3)/4;      coor[13]=1.5;
   coor[14]=(5*sqrt_3)/4;      coor[15]=1.5;
   coor[16]=(5*sqrt_3)/4;      coor[17]=0;
   coor[18]=(7*sqrt_3)/4;      coor[19]=0;
   coor[20]=sqrt_3/2;          coor[21]=2.25;
   coor[22]=(3*sqrt_3)/2;      coor[23]=2.25;
   
   
   //unit cell magnetic moments
   
   //hc
   mm[0]=0;            mm[1]=1; //1.75
   mm[2]=-sqrt_3/2;    mm[3]=0.5; //1.25
   mm[4]=sqrt_3/2;     mm[5]=0.5; //0.25
   mm[6]=0;            mm[7]=1;  //-0.25
   mm[8]=-sqrt_3/2;    mm[9]=0.5; //0.25
   mm[10]=sqrt_3/2;    mm[11]=0.5; //1.25
   
   mm[12]=-sqrt_3/2;   mm[13]=0.5; //1.25
   mm[14]=sqrt_3/2;    mm[15]=0.5; //0.25
   mm[16]=-sqrt_3/2;   mm[17]=0.5; //0.25
   mm[18]=sqrt_3/2;    mm[19]=0.5; //1.25
   mm[20]=0;           mm[21]=1; //1.75
   mm[22]=0;           mm[23]=1; //1.75
   
   cout<<"Number of spins in the unit cell = "<<coor.size()/2<<endl;
   
   cout<<"Period X of lattice = "<<period_X<<endl;
   cout<<"Period Y of lattice = "<<period_Y<<endl;
   
   /*
    //check of coordinates
    cout<<"{";
    for(unsigned int i=0; i<coor.size(); i+=2)
    {
        cout<<"{"<<coor[i]<<","<<coor[i+1]<<"},";
    }
    cout<<"}\n";
    
    //check of magnetic moments
    cout<<"{";
    for(unsigned int i=0; i<mm.size(); i+=2)
    {
        cout<<"{"<<mm[i]<<","<<mm[i+1]<<"},";
    }
    cout<<"}\n";
    //*/
   
   
   
   int numuc = 5;
   
   cout<<"Input Number of Unit Cels in linear dimension "<<endl;
   cin>>numuc;
   cout<<"Lenear size of system in unit cell's spins = "<<numuc*coor.size()/2<<" spins"<<endl;
   cout<<"Number of spins in system in linear = "<<numuc*coor.size()/2<<endl;
   
   numuc=numuc*numuc;
   cout<<"Total number of unit cells = "<<numuc<<endl;
   
   if(sqrt(numuc)==1)
      coorm=coor;
   vector<double> repit;
   vector<double> repitm;
   
   //cout<<"hex_centers.size()_0 = "<<hex_centers.size()<<endl;
   
   //масштабирование по Y центров гексагонов
   for(unsigned int y=0;y<hex_centers.size(); y+=2)
   {
      for(int j=0; j<sqrt(numuc); j++ )
      {
         if(hex_centers[y]!=hex_centers[y+2])
         {
            temp_vector.push_back(hex_centers[y]);
            temp_vector.push_back(hex_centers[y+1]+period_Y*j);
         }
         else
         {
            repit.push_back(hex_centers[y]);
            repit.push_back(hex_centers[y+1]);
            
            while(hex_centers[y]==hex_centers[y+2])
            {
               y+=2;	
               repit.push_back(hex_centers[y]);
               repit.push_back(hex_centers[y+1]);
            }
            
            for(unsigned int i=0; i<repit.size(); i++)
            {
               temp_vector.push_back(repit[i]);
            }
            
            for(int k=j+1; k<sqrt(numuc); k++)
            {
               for(unsigned int i=0; i<repit.size(); i+=2)
               { 
                  temp_vector.push_back(repit[i]);  temp_vector.push_back(repit[i+1]+period_Y*k);
               }
            }
            j=sqrt(numuc);
            repit.erase(repit.begin(),repit.end());
         }
      }
   }
   hex_centers=temp_vector;
   //cout<<"hex_centers.size()_1 = "<<hex_centers.size()<<endl;
   
   //scaling over y and ordering of array of coorm and magnetic moments
   for(unsigned int y=0;y<coor.size(); y+=2)
   {
      for(int j=0; j<sqrt(numuc); j++ )
      {
         if(coor[y]!=coor[y+2])
         {
            coorm.push_back (coor[y]);
            coorm.push_back (coor[y+1]+period_Y*j);
            
            mxmy.push_back(mm[y]);
            mxmy.push_back(mm[y+1]); 
            
         }
         else
         {  repit.push_back (coor[y]);
            repit.push_back (coor[y+1]);
            repitm.push_back(mm[y]);
            repitm.push_back(mm[y+1]);
            
            while(coor[y]==coor[y+2])
            {
               y+=2;	
               repit.push_back (coor[y]);
               repit.push_back (coor[y+1]);
               repitm.push_back(mm[y]);
               repitm.push_back(mm[y+1]);
            }
            
            for(unsigned int i=0; i<repit.size(); i++)
            {
               coorm.push_back(repit[i]);
               mxmy.push_back(repitm[i]);
            }
            for(int k=j+1; k<sqrt(numuc); k++)
               for(unsigned int i=0; i<repit.size(); i+=2)
               { 
                  coorm.push_back(repit[i]);  coorm.push_back(repit[i+1]+period_Y*k);
                  mxmy.push_back(repitm[i]);  mxmy.push_back(repitm[i+1]);
               }
            j=sqrt(numuc);
            repit.erase(repit.begin(),repit.end());
            repitm.erase(repitm.begin(),repitm.end());
         }
      }
   }
   
   //scaling over x
   coor=coorm;
   mm=mxmy;
   
   for(int j=1; j<sqrt(numuc); j++ )
   {
      for(unsigned int x=0;x<coor.size(); x+=2)
      {
         coorm.push_back(coor[x]+period_X*j);
         coorm.push_back(coor[x+1]);
         
         mxmy.push_back(mm[x]);
         mxmy.push_back(mm[x+1]); 
      }
      for(unsigned int x=0;x<hex_centers.size(); x+=2)
      {
         temp_vector.push_back(hex_centers[x]+period_X*j);
         temp_vector.push_back(hex_centers[x+1]);
      }
   }
   hex_centers=temp_vector;
   //cout<<"hex_centers.size()_2 = "<<hex_centers.size()<<endl;
   
   //запись координат в файл
   foutc<<"{";
   for(unsigned int i=0; i<coorm.size(); i+=2)
   {
      if(i!=coorm.size()-2)
         foutc<<"{"<<coorm[i]<<","<<coorm[i+1]<<"},";
      else
         foutc<<"{"<<coorm[i]<<","<<coorm[i+1]<<"}";
   }
   foutc<<"}\n\n";
   foutc.close();
   
   cout<<"  Number of spins in system = "<<(coorm.size()+1)/2<<endl;
   
   //запись магнитных моментов в файл
   fouts<<"{";
   for(unsigned int i=0; i<coorm.size(); i+=2)
   {
      if(i!=coorm.size()-2)
         fouts<<"{{"<<coorm[i]-mxmy[i]/2<<","<<coorm[i+1]-mxmy[i+1]/2<<"},{"<<coorm[i]+mxmy[i]/2<<","<<coorm[i+1]+mxmy[i+1]/2<<"}},";
      else
         fouts<<"{{"<<coorm[i]-mxmy[i]/2<<","<<coorm[i+1]-mxmy[i+1]/2<<"},{"<<coorm[i]+mxmy[i]/2<<","<<coorm[i+1]+mxmy[i+1]/2<<"}}";
      
   }
   fouts<<"}\n";
   fouts.close();
   
   
   fout_centers<<"{";
   for(unsigned int i=0; i<hex_centers.size(); i+=2)
   {
      if(i!=hex_centers.size()-2)
         fout_centers<<"{"<<hex_centers[i]<<","<<hex_centers[i+1]<<"},";
      else
         fout_centers<<"{"<<hex_centers[i]<<","<<hex_centers[i+1]<<"}";
      
   }
   fout_centers<<"}\n\n";
   fout_centers.close();        
   
   
   double r, max_X_of_centers, max_Y_of_centers, min_Y_of_centers; 
   
   max_X_of_centers=hex_centers[0];
   max_Y_of_centers=hex_centers[1];
   min_Y_of_centers=hex_centers[1];
   for(unsigned int i=2; i<hex_centers.size(); i+=2)
   {
      if(max_X_of_centers<hex_centers[i]) max_X_of_centers = hex_centers[i];
      if(min_Y_of_centers>hex_centers[i+1]) min_Y_of_centers = hex_centers[i+1];
      if(max_Y_of_centers<hex_centers[i+1]) max_Y_of_centers = hex_centers[i+1];
   }
   //cout<<"max_X = "<<max_X_of_centers<<endl<<endl;
   
   ///вывод всех спинов на экран
   //for(unsigned int i=0; i<coorm.size(); i+=2) cout<<"["<<i/2<<"] "<<coorm[i]<<", "<<coorm[i+1]<<endl;
   
   
   
   ///
   /// \определение ядер и границ для каждого гексагона
   ///
   vector <vector <unsigned int> > hex_array(hex_centers.size()/2);
   vector <vector <unsigned int> > bound_array(hex_centers.size()/2);
   
   map <double,map <double,unsigned int>> temp_core;
   map <double,map <double,unsigned int>> temp_bound;
   
   cout<<"Number of hexagons = " << hex_array.size() <<endl<<endl;
   
   unsigned int counter_spins_in_core=0;
   unsigned int counter_spins_on_bound=0;
   
   cout<<"determining the cores and borders...\n";
   
   for(unsigned int x=0; x<hex_centers.size(); x+=2)
   {
      for(unsigned int i=0; i<coorm.size(); i+=2)
      {
         //определяем спины, которые входят в гексагон и границу
         r=(hex_centers[x]-coorm[i])*(hex_centers[x]-coorm[i])+
               (hex_centers[x+1]-coorm[i+1])*(hex_centers[x+1]-coorm[i+1]);
         
         //для гексагона
         if(r<=interaction_radius_squared_for_core)
         {
            hex_array[x/2].push_back(i);
            temp_core[coorm[i]][coorm[i+1]] = i;
            counter_spins_in_core++;
            //cout<<"hex["<<x/2<<"] = ";
            //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
            if(counter_spins_in_core == total_num_of_spins_in_core && counter_spins_on_bound == total_num_of_spins_in_boundaries)
            {
               //cout<<"hex["<<x/2<<"] = ";
               //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
               ///break;
            }
         }
         else
         {
            //для границы
            if(r<=interaction_radius_squared_for_border)
            {
               bound_array[x/2].push_back(i);
               temp_bound[coorm[i]][coorm[i+1]] = i;
               counter_spins_on_bound++;
               if(counter_spins_in_core == total_num_of_spins_in_core && counter_spins_on_bound == total_num_of_spins_in_boundaries)
               {
                  //cout<<"hex["<<x/2<<"] = ";
                  //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                  ///break;
               }
            }
         }
         
         if(i==coorm.size()-2 && (counter_spins_in_core<total_num_of_spins_in_core || counter_spins_on_bound<total_num_of_spins_in_boundaries))
         {
            for(unsigned int ii=0; ii<coorm.size(); ii+=2)
            {
               //сдвигаем по X
               r=(hex_centers[x]-(max_X_of_centers+coorm[ii]))*(hex_centers[x]-(max_X_of_centers+coorm[ii]))+
                     (hex_centers[x+1]-coorm[ii+1])*(hex_centers[x+1]-coorm[ii+1]);
               
               if(r<=interaction_radius_squared_for_core)
               {
                  hex_array[x/2].push_back(ii);
                  temp_core[max_X_of_centers+coorm[ii]][coorm[ii+1]] = ii;
                  counter_spins_in_core++;
                  if(counter_spins_in_core == total_num_of_spins_in_core && counter_spins_on_bound == total_num_of_spins_in_boundaries)
                  {
                     //cout<<"hex["<<x/2<<"] = ";
                     //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                     ///break;
                  }
               }
               else
               {
                  if(r<=interaction_radius_squared_for_border)
                  {
                     bound_array[x/2].push_back(ii);
                     temp_bound[max_X_of_centers+coorm[ii]][coorm[ii+1]] = ii;
                     counter_spins_on_bound++;
                     if(counter_spins_in_core == total_num_of_spins_in_core && counter_spins_on_bound == total_num_of_spins_in_boundaries)
                     {
                        //cout<<"hex["<<x/2<<"] = ";
                        //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                        ///break;
                     }
                  }
               }
               
               if(ii==coorm.size()-2 && (counter_spins_in_core<total_num_of_spins_in_core || counter_spins_on_bound<total_num_of_spins_in_boundaries))
               {
                  for(unsigned int iii=0; iii<coorm.size(); iii+=2)
                  {
                     //сдвигаем по Y
                     r=(hex_centers[x]-coorm[iii])*(hex_centers[x]-coorm[iii])+
                           (hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iii+1]))*(hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iii+1]));
                     
                     if(r<=interaction_radius_squared_for_core)
                     {
                        hex_array[x/2].push_back(iii);
                        temp_core[coorm[iii]][max_Y_of_centers+min_Y_of_centers+coorm[iii+1]] = iii;
                        counter_spins_in_core++;
                        if(counter_spins_in_core == total_num_of_spins_in_core && counter_spins_on_bound == total_num_of_spins_in_boundaries)
                        {
                           //cout<<"hex["<<x/2<<"] = ";
                           //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                           ///break;
                        }
                     }
                     else
                     {
                        if(r<=interaction_radius_squared_for_border)
                        {
                           bound_array[x/2].push_back(iii);
                           temp_bound[coorm[iii]][max_Y_of_centers+min_Y_of_centers+coorm[iii+1]] = iii;
                           counter_spins_on_bound++;
                           if(counter_spins_in_core == total_num_of_spins_in_core && counter_spins_on_bound == total_num_of_spins_in_boundaries)
                           {
                              //cout<<"hex["<<x/2<<"] = ";
                              //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                              ///break;
                           }
                        }
                     }
                     
                     if(iii==coorm.size()-2 && (counter_spins_in_core<total_num_of_spins_in_core || counter_spins_on_bound<total_num_of_spins_in_boundaries))
                     {
                        for(unsigned int iiii=0; iiii<coorm.size(); iiii+=2)
                        {
                           //сдвигаем по диагонали (X и Y)
                           r=(hex_centers[x]-(max_X_of_centers+coorm[iiii]))*(hex_centers[x]-(max_X_of_centers+coorm[iiii]))+
                                 (hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iiii+1]))*(hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iiii+1]));
                           
                           if(r<=interaction_radius_squared_for_core)
                           {
                              hex_array[x/2].push_back(iiii);
                              temp_core[max_X_of_centers+coorm[iiii]][max_Y_of_centers+min_Y_of_centers+coorm[iiii+1]] = iiii;
                              counter_spins_in_core++;
                              if(counter_spins_in_core == total_num_of_spins_in_core && counter_spins_on_bound == total_num_of_spins_in_boundaries)
                              {
                                 //cout<<"hex["<<x/2<<"] = ";
                                 //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                 ///break;
                              }
                           }
                           else
                           {
                              if(r<=interaction_radius_squared_for_border)
                              {
                                 bound_array[x/2].push_back(iiii);
                                 temp_bound[max_X_of_centers+coorm[iiii]][max_Y_of_centers+min_Y_of_centers+coorm[iiii+1]] = iiii;
                                 counter_spins_on_bound++;
                                 if(counter_spins_in_core == total_num_of_spins_in_core && counter_spins_on_bound == total_num_of_spins_in_boundaries)
                                 {
                                    //cout<<"hex["<<x/2<<"] = ";
                                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                    ///break;
                                 }
                              }
                           }
                           
                           if(iiii==coorm.size()-2 && counter_spins_on_bound<total_num_of_spins_in_boundaries)
                           {
                              for(unsigned int iiiii=0; iiiii<coorm.size(); iiiii+=2)
                              {
                                 //сдвигаем по -X
                                 r=(hex_centers[x]-(-max_X_of_centers+coorm[iiiii]))*(hex_centers[x]-(-max_X_of_centers+coorm[iiiii]))+
                                       (hex_centers[x+1]-coorm[iiiii+1])*(hex_centers[x+1]-coorm[iiiii+1]);
                                 
                                 if(r>interaction_radius_squared_for_core && r<=interaction_radius_squared_for_border)
                                 {
                                    bound_array[x/2].push_back(iiiii);
                                    temp_bound[-max_X_of_centers+coorm[iiiii]][coorm[iiiii+1]] = iiiii;
                                    counter_spins_on_bound++;
                                    if(counter_spins_on_bound == total_num_of_spins_in_boundaries)
                                    {
                                       //cout<<"hex["<<x/2<<"] = ";
                                       //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                       ///break;
                                    }
                                 }
                                 
                                 if(iiiii==coorm.size()-2 && counter_spins_on_bound<total_num_of_spins_in_boundaries)
                                 {
                                    for(unsigned int iiiiii=0; iiiiii<coorm.size(); iiiiii+=2)
                                    {
                                       //сдвигаем по -Y
                                       r=(hex_centers[x]-coorm[iiiiii])*(hex_centers[x]-coorm[iiiiii])+
                                             (hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiii+1]))*(hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiii+1]));
                                       
                                       if(r>interaction_radius_squared_for_core && r<=interaction_radius_squared_for_border)
                                       {
                                          bound_array[x/2].push_back(iiiiii);
                                          temp_bound[coorm[iiiiii]][-max_Y_of_centers-min_Y_of_centers+coorm[iiiiii+1]] = iiiiii;
                                          counter_spins_on_bound++;
                                          if(counter_spins_on_bound == total_num_of_spins_in_boundaries)
                                          {
                                             //cout<<"hex["<<x/2<<"] = ";
                                             //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                             ///break;
                                          }
                                       }
                                       
                                       if(iiiiii==coorm.size()-2 && counter_spins_on_bound<total_num_of_spins_in_boundaries)
                                       {
                                          for(unsigned int iiiiiii=0; iiiiiii<coorm.size(); iiiiiii+=2)
                                          {
                                             //сдвигаем по диагонали (-X и -Y)
                                             r=(hex_centers[x]-(-max_X_of_centers+coorm[iiiiiii]))*(hex_centers[x]-(-max_X_of_centers+coorm[iiiiiii]))+
                                                   (hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiii+1]))*(hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiii+1]));
                                             
                                             if(r>interaction_radius_squared_for_core && r<=interaction_radius_squared_for_border)
                                             {
                                                bound_array[x/2].push_back(iiiiiii);
                                                temp_bound[-max_X_of_centers+coorm[iiiiiii]][-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiii+1]] = iiiiiii;
                                                counter_spins_on_bound++;
                                                if(counter_spins_on_bound == total_num_of_spins_in_boundaries)
                                                {
                                                   //cout<<"hex["<<x/2<<"] = ";
                                                   //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                   ///break;
                                                }
                                             }
                                             
                                             if(iiiiiii==coorm.size()-2 && counter_spins_on_bound<total_num_of_spins_in_boundaries)
                                             {
                                                for(unsigned int iiiiiiii=0; iiiiiiii<coorm.size(); iiiiiiii+=2)
                                                {
                                                   //сдвигаем по диагонали (-X и Y)
                                                   r=(hex_centers[x]-(-max_X_of_centers+coorm[iiiiiiii]))*(hex_centers[x]-(-max_X_of_centers+coorm[iiiiiiii]))+
                                                         (hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iiiiiiii+1]))*(hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iiiiiiii+1]));
                                                   
                                                   if(r>interaction_radius_squared_for_core && r<=interaction_radius_squared_for_border)
                                                   {
                                                      bound_array[x/2].push_back(iiiiiiii);
                                                      temp_bound[-max_X_of_centers+coorm[iiiiiiii]][max_Y_of_centers+min_Y_of_centers+coorm[iiiiiiii+1]] = iiiiiiii;
                                                      counter_spins_on_bound++;
                                                      if(counter_spins_on_bound == total_num_of_spins_in_boundaries)
                                                      {
                                                         //cout<<"hex["<<k<<"] = ";
                                                         //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                         ///break;
                                                      }
                                                   }
                                                   
                                                   if(iiiiiiii==coorm.size()-2 && counter_spins_on_bound<total_num_of_spins_in_boundaries)
                                                   {
                                                      for(unsigned int iiiiiiiii=0; iiiiiiiii<coorm.size(); iiiiiiiii+=2)
                                                      {
                                                         //сдвигаем по диагонали (X и -Y)
                                                         r=(hex_centers[x]-(max_X_of_centers+coorm[iiiiiiiii]))*(hex_centers[x]-(max_X_of_centers+coorm[iiiiiiiii]))+
                                                               (hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiiiii+1]))*(hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiiiii+1]));
                                                         
                                                         if(r>interaction_radius_squared_for_core && r<=interaction_radius_squared_for_border)
                                                         {
                                                            bound_array[x/2].push_back(iiiiiiiii);
                                                            temp_bound[max_X_of_centers+coorm[iiiiiiiii]][-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiiiii+1]] = iiiiiiiii;
                                                            counter_spins_on_bound++;
                                                            if(counter_spins_on_bound == total_num_of_spins_in_boundaries)
                                                            {
                                                               //cout<<"hex["<<k<<"] = ";
                                                               //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                               ///break;
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      counter_spins_in_core=0;
      counter_spins_on_bound=0;
      
      //проверка размера
      //cout<<"temp_core.size() = "<<temp_core.size()<<endl;
      //cout<<"temp_bound.size() = "<<temp_bound.size()<<endl;
      //cout<<"hex_array["<<x/2<<"].size() = "<<hex_array[x/2].size()<<endl;
      
      int k = 0;
      for (auto it : temp_core)
      {
         for (auto it2 : it.second)
         {
            //cout<<"spin "<<it2.second/2<<": ["; //номер спина
            //cout<<it.first<<", "; //координата X
            //cout<<it2.first<<"]\n"; //координата Y
            hex_array[x/2][k] = it2.second;
            k++;
         }
      }
      
      //cout<<"\nboundary:\n";
      k = 0;
      for (auto it : temp_bound)
      {
         for (auto it2 : it.second)
         {
            //cout<<"spin "<<it2.second/2<<": ["; //номер спина
            //cout<<it.first<<", "; //координата X
            //cout<<it2.first<<"]\n"; //координата Y
            bound_array[x/2][k] = it2.second;
            k++;
         }
      }
      
      //очистка
      temp_core.erase(temp_core.begin(),temp_core.end());
      temp_bound.erase(temp_bound.begin(),temp_bound.end());
      
      //cout<<"hex "<<x/2<<" ["<<hex_centers[x]<<","<<hex_centers[x+1]<<"] = ";
      //cout<<hex_array[x/2].size()/2<<endl;
      //cout<<bound_array[x/2].size()/2<<endl;
      
      //вывод номеров и координат спинов гексагона
      //if(x/2==6)
      /*
        {
            cout<<"hex "<<x/2<<" ["<<hex_centers[x]<<","<<hex_centers[x+1]<<"] = ";
            cout<<hex_array[x/2].size()<<endl;
            for(unsigned int i=0; i<hex_array[x/2].size(); i++)
            {
                cout<<"spin in core: "<<i<<" ("<<hex_array[x/2][i]/2<<") "<<" ["<<coorm[hex_array[x/2][i]]<<", "<<coorm[hex_array[x/2][i]+1]<<"]\n";
            }
        }
        system("pause");
        //*/
      
      //проверка на правильное количество спинов в каждом гексагоне
      if(hex_array[x/2].size()!=total_num_of_spins_in_core)
      {
         cout<<endl<<"Wrong hex size!!!!!!!!!!!!!!!!!!!!!"<<endl;
         cout<<"hex "<<x/2<<" ["<<hex_centers[x]<<", "<<hex_centers[x+1]<<"] = ";
         cout<<hex_array[x/2].size()<<endl<<endl;
      }
      
      //проверка на правильное количество спинов на границе каждого гексагона
      if(bound_array[x/2].size()!=total_num_of_spins_in_boundaries)
      {
         cout<<endl<<"Wrong border size!!!!!!!!!!!!!!!!!!!!!"<<endl;
         cout<<"hex "<<x/2<<" ["<<hex_centers[x]<<", "<<hex_centers[x+1]<<"] = ";
         cout<<bound_array[x/2].size()<<endl<<endl;
      }
      
      ///вывод на экран всех спинов, входящих в гексагон
      //for(unsigned int c=0; c<hex_array[x/2].size(); c+=2) cout<<hex_array[x/2][c]<<", "<<hex_array[x/2][c+1]<<endl;
   }
   
   
   
   ///
   /// определение соседей для каждого спина
   /// 
   cout << "searching neighbors...\n";
   
   struct spin_struct
   {
      unsigned int num;
      double x;
      double y;
      double mx;
      double my;
      
      //конструктор
      spin_struct(unsigned int num, double x, double y, double mx, double my){spin_struct::num = num; spin_struct::x = x; spin_struct::y = y; spin_struct::mx = mx; spin_struct::my = my; }
   };
   
   vector <vector <spin_struct> > neighbors(coorm.size()/2); //array of neighbors
   unsigned int counter_of_neighbors=0;
   
   cout<<"\nneighbors";
   //поиск соседей
   for(unsigned int i=0; i<coorm.size(); i+=2)
   {
      for(unsigned int j=0; j<coorm.size(); j+=2)
      {
         r=(coorm[i]-coorm[j])*(coorm[i]-coorm[j])+
               (coorm[i+1]-coorm[j+1])*(coorm[i+1]-coorm[j+1]);
         
         if(i!=j && r<=interaction_radius_squared_for_neighbors)
         {
            neighbors[i/2].push_back(spin_struct(j,coorm[j],coorm[j+1],mxmy[j],mxmy[j+1]));
            counter_of_neighbors++;
            if(counter_of_neighbors == total_num_of_neighbors)
            {
               //cout<<"\nneighbors["<<i/2<<"] = ";
               //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
               ///break;
            }
         }
         
         if(j==coorm.size()-2 && counter_of_neighbors<total_num_of_neighbors)
         {
            for(unsigned int jj=0; jj<coorm.size(); jj+=2)
            {
               //сдвигаем по X
               r=(coorm[i]-(max_X_of_centers+coorm[jj]))*(coorm[i]-(max_X_of_centers+coorm[jj]))+
                     (coorm[i+1]-coorm[jj+1])*(coorm[i+1]-coorm[jj+1]);
               
               if(i!=jj && r<=interaction_radius_squared_for_neighbors)
               {
                  neighbors[i/2].push_back(spin_struct(jj,max_X_of_centers+coorm[jj],coorm[jj+1],mxmy[jj],mxmy[jj+1]));
                  counter_of_neighbors++;
                  if(counter_of_neighbors == total_num_of_neighbors)
                  {
                     //cout<<"\nneighbors["<<i/2<<"] = ";
                     //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                     ///break;
                  }
               }
               
               if(jj==coorm.size()-2 && counter_of_neighbors<total_num_of_neighbors)
               {
                  for(unsigned int jjj=0; jjj<coorm.size(); jjj+=2)
                  {
                     //сдвигаем по -X
                     r=(coorm[i]-(-max_X_of_centers+coorm[jjj]))*(coorm[i]-(-max_X_of_centers+coorm[jjj]))+
                           (coorm[i+1]-coorm[jjj+1])*(coorm[i+1]-coorm[jjj+1]);
                     
                     if(i!=jjj && r<=interaction_radius_squared_for_neighbors)
                     {
                        neighbors[i/2].push_back(spin_struct(jjj, -max_X_of_centers+coorm[jjj], coorm[jjj+1], mxmy[jjj], mxmy[jjj+1]));
                        counter_of_neighbors++;
                        if(counter_of_neighbors == total_num_of_neighbors)
                        {
                           //cout<<"\nneighbors["<<i/2<<"] = ";
                           //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                           ///break;
                        }
                     }
                     
                     if(jjj==coorm.size()-2 && counter_of_neighbors<total_num_of_neighbors)
                     {
                        for(unsigned int jjjj=0; jjjj<coorm.size(); jjjj+=2)
                        {
                           //сдвигаем по Y
                           r=(coorm[i]-coorm[jjjj])*(coorm[i]-coorm[jjjj])+
                                 (coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjj+1]))*(coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjj+1]));
                           
                           if(i!=jjjj && r<=interaction_radius_squared_for_neighbors)
                           {
                              neighbors[i/2].push_back(spin_struct(jjjj, coorm[jjjj], max_Y_of_centers+min_Y_of_centers+coorm[jjjj+1], mxmy[jjjj], mxmy[jjjj+1]));
                              counter_of_neighbors++;
                              if(counter_of_neighbors == total_num_of_neighbors)
                              {
                                 //cout<<"\nneighbors["<<i/2<<"] = ";
                                 //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                 ///break;
                              }
                           }
                           
                           
                           if(jjjj==coorm.size()-2 && counter_of_neighbors<total_num_of_neighbors)
                           {
                              for(unsigned int jjjjj=0; jjjjj<coorm.size(); jjjjj+=2)
                              {
                                 //сдвигаем по -Y
                                 r=(coorm[i]-coorm[jjjjj])*(coorm[i]-coorm[jjjjj])+
                                       (coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjj+1]))*(coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjj+1]));
                                 
                                 if(i!=jjjjj && r<=interaction_radius_squared_for_neighbors)
                                 {
                                    neighbors[i/2].push_back(spin_struct(jjjjj, coorm[jjjjj], -max_Y_of_centers-min_Y_of_centers+coorm[jjjjj+1], mxmy[jjjjj], mxmy[jjjjj+1]));
                                    counter_of_neighbors++;
                                    if(counter_of_neighbors == total_num_of_neighbors)
                                    {
                                       //cout<<"\nneighbors["<<i/2<<"] = ";
                                       //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                       ///break;
                                    }
                                 }
                                 
                                 if(jjjjj==coorm.size()-2 && counter_of_neighbors<total_num_of_neighbors)
                                 {
                                    for(unsigned int jjjjjj=0; jjjjjj<coorm.size(); jjjjjj+=2)
                                    {
                                       //сдвигаем по диагонали (X и Y)
                                       r=(coorm[i]-(max_X_of_centers+coorm[jjjjjj]))*(coorm[i]-(max_X_of_centers+coorm[jjjjjj]))+
                                             (coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjjjj+1]))*(coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjjjj+1]));
                                       
                                       if(i!=jjjjjj && r<=interaction_radius_squared_for_neighbors)
                                       {
                                          neighbors[i/2].push_back(spin_struct(jjjjjj, max_X_of_centers+coorm[jjjjjj], max_Y_of_centers+min_Y_of_centers+coorm[jjjjjj+1], mxmy[jjjjjj], mxmy[jjjjjj+1]));
                                          counter_of_neighbors++;
                                          if(counter_of_neighbors == total_num_of_neighbors)
                                          {
                                             //cout<<"\nneighbors["<<i/2<<"] = ";
                                             //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                             ///break;
                                          }
                                       }
                                       
                                       if(jjjjjj==coorm.size()-2 && counter_of_neighbors<total_num_of_neighbors)
                                       {
                                          for(unsigned int jjjjjjj=0; jjjjjjj<coorm.size(); jjjjjjj+=2)
                                          {
                                             //сдвигаем по диагонали (-X и -Y)
                                             r=(coorm[i]-(-max_X_of_centers+coorm[jjjjjjj]))*(coorm[i]-(-max_X_of_centers+coorm[jjjjjjj]))+
                                                   (coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjj+1]))*(coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjj+1]));
                                             
                                             if(i!=jjjjjjj && r<=interaction_radius_squared_for_neighbors)
                                             {
                                                neighbors[i/2].push_back(spin_struct(jjjjjjj, -max_X_of_centers+coorm[jjjjjjj], -max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjj+1], mxmy[jjjjjjj], mxmy[jjjjjjj+1]));
                                                counter_of_neighbors++;
                                                if(counter_of_neighbors == total_num_of_neighbors)
                                                {
                                                   //cout<<"\nneighbors["<<i/2<<"] = ";
                                                   //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                                   ///break;
                                                }
                                             }
                                             
                                             if(jjjjjjj==coorm.size()-2 && counter_of_neighbors<total_num_of_neighbors)
                                             {
                                                for(unsigned int jjjjjjjj=0; jjjjjjjj<coorm.size(); jjjjjjjj+=2)
                                                {
                                                   //сдвигаем по диагонали (X и -Y)
                                                   r=(coorm[i]-(max_X_of_centers+coorm[jjjjjjjj]))*(coorm[i]-(max_X_of_centers+coorm[jjjjjjjj]))+
                                                         (coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjjj+1]))*(coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjjj+1]));
                                                   
                                                   if(i!=jjjjjjjj && r<=interaction_radius_squared_for_neighbors)
                                                   {
                                                      neighbors[i/2].push_back(spin_struct(jjjjjjjj, max_X_of_centers+coorm[jjjjjjjj], -max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjjj+1], mxmy[jjjjjjjj], mxmy[jjjjjjjj+1]));
                                                      counter_of_neighbors++;
                                                      if(counter_of_neighbors == total_num_of_neighbors)
                                                      {
                                                         //cout<<"\nneighbors["<<i/2<<"] = ";
                                                         //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                                         ///break;
                                                      }
                                                   }
                                                   
                                                   if(jjjjjjjj==coorm.size()-2 && counter_of_neighbors<total_num_of_neighbors)
                                                   {
                                                      for(unsigned int jjjjjjjjj=0; jjjjjjjjj<coorm.size(); jjjjjjjjj+=2)
                                                      {
                                                         //сдвигаем по диагонали (-X и Y)
                                                         r=(coorm[i]-(-max_X_of_centers+coorm[jjjjjjjjj]))*(coorm[i]-(-max_X_of_centers+coorm[jjjjjjjjj]))+
                                                               (coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjjjjjjj+1]))*(coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjjjjjjj+1]));
                                                         
                                                         if(i!=jjjjjjjjj && r<=interaction_radius_squared_for_neighbors)
                                                         {
                                                            neighbors[i/2].push_back(spin_struct(jjjjjjjjj, -max_X_of_centers+coorm[jjjjjjjjj], max_Y_of_centers+min_Y_of_centers+coorm[jjjjjjjjj+1], mxmy[jjjjjjjjj], mxmy[jjjjjjjjj+1]));
                                                            counter_of_neighbors++;
                                                            if(counter_of_neighbors == total_num_of_neighbors)
                                                            {
                                                               //cout<<"\nneighbors["<<i/2<<"] = ";
                                                               //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                                               ///break;
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      
      counter_of_neighbors=0;
      
      //cout << i/2 << ": " << neighbors[i/2].size() << endl;
      
      //проверка на правильное количество соседей
      if(neighbors[i/2].size()!=total_num_of_neighbors)
      {
         cout<<endl<<"Wrong number of neighbors!!!!!!!!!!!!!!!!!!!!!"<<endl;
         cout<<"spin "<<i/2<<" ["<<coorm[i]<<", "<<coorm[i+1]<<"] = ";
         cout<<neighbors[i/2].size()<<endl<<endl;
      }
   }
   
   cout<<" OK!\n";
   
   //проверяем спин на правильность соседей
   /*
    int s=38;
    cout<<"spin "<<s<<" ["<<coorm[s*2]<<", "<<coorm[s*2+1]<<"] = ";
    cout<<neighbors[s].size()<<" neighbors"<<endl<<endl;
    for(unsigned int i=0; i<neighbors[s].size(); ++i)
    {
        //выводим номера соседних спинов
        cout<<i<<": "<< neighbors[s][i].number/2;
        //выводим координаты соседних спинов
        cout<<" ["<<neighbors[s][i].x<<", "<<neighbors[s][i].y<<"]";
        cout<<" (["<<coorm[neighbors[s][i].number]<<", "<<coorm[neighbors[s][i].number+1]<<"])\n";
    }
    system("pause");
    //*/
   
   
   
   ///
   /// расчет свойств для блока спинов (ядро+граница)
   /// 
   cout<<"sample calculation...\n";
   
   //структура номер, координаты, магнитные моменты спина
   struct num_xy_MxMy_struct{
      unsigned int num;
      double x;
      double y;
      double mx;
      double my;
   };
   
   
   //
   //сохраняем отдельно блок для перебора, за образец берем блок с центром #6
   vector <num_xy_MxMy_struct> sample_core(total_num_of_spins_in_core);
   vector <num_xy_MxMy_struct> sample_bound(total_num_of_spins_in_boundaries);
   
   //cout<<"core #6: \n";
   for(unsigned int i=0; i<hex_array[6].size(); ++i)
   {
      sample_core[i].num = hex_array[6][i];
      sample_core[i].x = coorm[sample_core[i].num];
      sample_core[i].y = coorm[sample_core[i].num+1];
      sample_core[i].mx = mxmy[sample_core[i].num];
      sample_core[i].my = mxmy[sample_core[i].num+1];
      
      //cout<<sample_core[i].num/2<<", ";
      //cout<<"# "<<sample_core[i].num/2<<endl;
      //cout<<"coord: "<<sample_core[i].x<<", "<<sample_core[i].y<<endl;
      //cout<<"mm: "<<sample_core[i].mx<<", "<<sample_core[i].my<<endl;
   }
   
   //cout<<endl;
   
   //cout<<"bound #6: \n";
   for(unsigned int i=0; i<bound_array[6].size(); ++i)
   {
      sample_bound[i].num = bound_array[6][i];
      sample_bound[i].x = coorm[sample_bound[i].num];
      sample_bound[i].y = coorm[sample_bound[i].num+1];
      sample_bound[i].mx = mxmy[sample_bound[i].num];
      sample_bound[i].my = mxmy[sample_bound[i].num+1];
      
      //cout<<sample_bound[i].num/2<<", ";
      //cout<<"# "<<sample_bound[i].num/2<<endl;
      //cout<<"coord: "<<sample_bound[i].x<<", "<<sample_bound[i].y<<endl;
      //cout<<"mm: "<<sample_bound[i].mx<<", "<<sample_bound[i].my<<endl;
   }
   
   cout<<endl;
   
   //system("pause");
   //
   
   
   /** перебор всех конфигураций границ и блока, сохранение E и M в массив --->>>> */
   
   unsigned int states_amount_of_core = 1<<total_num_of_spins_in_core; //количество конфигураций ядра
   unsigned int boundaries_amount = 1<<total_num_of_spins_in_boundaries; //количество вариантов границ
   
   vector <double> Emin_array (boundaries_amount); //минимальные значания энергий для каждой границы
   
   //структура хранит энергию с учетом соседей, намагниченность по осям и конфигурации
   struct EM_state_struct{
      double E;
      double Mx;
      double My;
      vector <int> state_array;
      
      //конструктор
      EM_state_struct(double E, double Mx, double My){EM_state_struct::E = E; EM_state_struct::Mx = Mx; EM_state_struct::My = My;}
   };
   
   //структура хранит энергию с учетом соседей и намагниченность по осям
   struct EM_struct{
      double E;
      double Mx;
      double My;
   };
   
   //
   //массив границ, энергий, намагниченности и соответствующих состояний
   //boundaries_EM_core[boundaries_amount][num_of_unique_EMxMy].{E,Mx,My,state_array[]}
   vector < vector <EM_state_struct> > boundaries_EM_core (boundaries_amount); 
   
   //массив границ, состояний и их энергий и намагниченности
   //boundaries_core_EM[boundaries_amount][states_amount_of_core].{E,Mx,My}
   vector < vector <EM_struct> > boundaries_core_EM (boundaries_amount, vector <EM_struct> (states_amount_of_core)); 
   
   vector < vector <int> > core_array (states_amount_of_core, vector <int> (total_num_of_spins_in_core)); //хранит +1 и -1
   vector <int> boundaries_array (total_num_of_spins_in_boundaries); //хранит +1 и -1
   
   double Mx, Mx_1, Mx_2;
   double My, My_1, My_2;
   double E_block_tot = 0; //энергия всего блока (ядро + граница)
   
   unsigned int count=0;
   
   //перебор конфигураций ядра побитовым сдвигом
   for(unsigned int state_num=0; state_num<states_amount_of_core; ++state_num)
   {
      ///cout << "State: " << state_num << " ( ";
      count=0;
      
      for(int spin_num_in_core=(total_num_of_spins_in_core-1); spin_num_in_core>=0; --spin_num_in_core)
      {
         if(state_num&(1<<spin_num_in_core))
            core_array[state_num][count]=1;
         else
            core_array[state_num][count]=-1;
         
         ///cout << core_array[state_num][count] << " "; //вывод конфигурации на экран
         count++;
      }
      ///cout<<")"<<endl;
   }
   
   unsigned int num_of_unique_EM = 0;
   
   
   double X, Y;
   double mx_i, mx_j, my_i, my_j;
   
   ///int kk=0;
   
   bool flag=0;
   
   //перебор границ
   for(unsigned int boundary_num=0; boundary_num<boundaries_amount; ++boundary_num)
   {
      ///cout << "\nBoundary: " << boundary_num << endl;
      
      Mx_1=0; //проверить подсчет намагниченности!!!!!!!!!!!!!!!
      My_1=0; //проверить подсчет намагниченности!!!!!!!!!!!!!!!
      count=0;
      
      //перебор границ побитовым сдвигом
      ///cout << "\nBoundary: " << boundary_num << "( ";
      for(int spin_num_in_boundary=(total_num_of_spins_in_boundaries-1); spin_num_in_boundary>=0; --spin_num_in_boundary)
      {   
         if(boundary_num&(1<<spin_num_in_boundary)) 
            boundaries_array[count]=1;
         else 
            boundaries_array[count]=-1;
         
         ///cout << boundaries_array [count] << " "; //вывод границы на экран
         
         //Mx_1+=boundaries_array[count] * sample_bound[count].mx; //проверить!!!!!!!!!!!!
         //My_1+=boundaries_array[count] * sample_bound[count].my; //проверить!!!!!!!!!!!!
         
         //заменяем магнитные моменты на границе
         mxmy[sample_bound[count].num]=boundaries_array[count] * sample_bound[count].mx;
         Mx_1+=mxmy[sample_bound[count].num];
         mxmy[sample_bound[count].num+1]=boundaries_array[count] * sample_bound[count].my;
         My_1+=mxmy[sample_bound[count].num+1];
         
         count++;
      }
      ///cout<<")"<<endl;
      
      ///cout<<"My_1 = "<<My_1<<endl;
      
      //перебор состояний в рамках границ и подсчет энергии и намагниченности
      for(unsigned int state_num=0; state_num<states_amount_of_core; ++state_num)
      {
         Mx_2=0;
         My_2=0;
         E_block_tot = 0;
         
         //заменяем магнитные моменты в ядре
         for(unsigned int spin_num_in_core=0; spin_num_in_core<total_num_of_spins_in_core; ++spin_num_in_core)
         {
            mxmy[sample_core[spin_num_in_core].num]=sample_core[spin_num_in_core].mx * core_array[state_num][spin_num_in_core];
            mxmy[sample_core[spin_num_in_core].num+1]=sample_core[spin_num_in_core].my * core_array[state_num][spin_num_in_core];
         }
         
         //подсчет энергии и намагниченности в блоке
         ///cout << "State: " << state_num << " ( ";
         for(unsigned int spin_num_in_core=0; spin_num_in_core<total_num_of_spins_in_core; ++spin_num_in_core)
         {
            ///cout << core_array[state_num][spin_num_in_core] << " "; //вывод конфигурации на экран
            
            //cout<<"spin_num_in_core = "<<spin_num_in_core<<": "<<sample_core[spin_num_in_core].num/2<<endl;
            
            //mx_i = sample_core[spin_num_in_core].mx * core_array[state_num][spin_num_in_core];
            //my_i = sample_core[spin_num_in_core].my * core_array[state_num][spin_num_in_core];
            mx_i = mxmy[sample_core[spin_num_in_core].num];
            my_i = mxmy[sample_core[spin_num_in_core].num+1];
            
            ///cout<<"s"<<sample_core[spin_num_in_core].num/2<<"  mx_i = "<<mx_i<<", my_i = "<<my_i<<endl;
            ///system("pause");
            
            for(unsigned int neigh=0; neigh<neighbors[sample_core[spin_num_in_core].num/2].size(); ++neigh)
            {
               //cout<<neighbors[sample_core[spin_num_in_core].num/2][neigh].number/2<<", ";
               flag = 0;
               for(unsigned int i=0; i<spin_num_in_core; i++)
               {
                  if (neighbors[sample_core[spin_num_in_core].num/2][neigh].num == sample_core[i].num)
                     flag=1;
               }
               if(flag==0)
               {
                  mx_j = mxmy[neighbors[sample_core[spin_num_in_core].num/2][neigh].num];//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  my_j = mxmy[neighbors[sample_core[spin_num_in_core].num/2][neigh].num+1];
                  
                  ///cout<<"s"<<neighbors[sample_core[spin_num_in_core].num/2][neigh].num/2<<"  mx_j = "<<mx_j<<", my_j = "<<my_j<<endl;
                  
                  //энергия ядра с учетом границ 
                  X=(sample_core[spin_num_in_core].x-neighbors[sample_core[spin_num_in_core].num/2][neigh].x);
                  Y=(sample_core[spin_num_in_core].y-neighbors[sample_core[spin_num_in_core].num/2][neigh].y);
                  
                  E_block_tot += (mx_i*mx_j + my_i*my_j)/
                        sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y) -
                        3*(mx_i*X+my_i*Y)*(mx_j*X+my_j*Y)/
                        sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
               }
            }
            
            //намагниченность ядра
            Mx_2 += sample_core[spin_num_in_core].mx * core_array[state_num][spin_num_in_core]; //проверить!!!!!!!!
            My_2 += sample_core[spin_num_in_core].my * core_array[state_num][spin_num_in_core]; //проверить!!!!!!!!
         }
         
         //cout<<"E_block_tot = "<<E_block_tot<<endl;
         ///cout<<")"<<endl;
         ///system("pause");
         
         //намагниченность блока
         Mx = Mx_1+Mx_2; //проверить намагниченность!!!!!!!!!!!!
         My = My_1+My_2; //проверить намагниченность!!!!!!!!!!!!
         
         
         //!заменяем на 0 очень маленькие значения!!!!!
         if(fabs(Mx)<1e-13) Mx = 0; 
         if(fabs(My)<1e-13) My = 0; 
         if(fabs(E_block_tot)<1e-13) E_block_tot = 0; 
         
         //cout<<"My_2 = "<<My_2<<endl;
         
         //cout<<"E_block_tot = "<<E_block_tot<<", Mx= "<<Mx<<", My= "<<My<<endl;
         
         ///M и E блока посчитано, далее заполняем массив с учетом вырождений
         
         //system("pause");
         
         if(state_num == 0)
         {
            //ищем минимальную энергию для границы
            //запоминаем первую энергию
            Emin_array[boundary_num]=E_block_tot;
            
            boundaries_EM_core[boundary_num].push_back(EM_state_struct(E_block_tot,Mx,My));
            boundaries_EM_core[boundary_num][0].state_array.push_back(0);
            
            //заполняем массив границ и соответствующие блоки
            boundaries_core_EM[boundary_num][state_num].E = E_block_tot;
            boundaries_core_EM[boundary_num][state_num].Mx = Mx;
            boundaries_core_EM[boundary_num][state_num].My = My;
            
            
            //cout<<"E_block_tot = "<<E_block_tot<<", Mx= "<<Mx<<", My= "<<My<<endl;
            //cout<<"E= "<<boundaries_EM_core[boundary_num][boundaries_EM_core[boundary_num].size()-1].E<<endl;
            //system("pause");
         }
         else
         {
            //ищем минимальную энергию для границы
            if(Emin_array[boundary_num] > E_block_tot)
               Emin_array[boundary_num] = E_block_tot;
            
            //заполняем массив границ и соответствующие блоки
            boundaries_core_EM[boundary_num][state_num].E = E_block_tot;
            boundaries_core_EM[boundary_num][state_num].Mx = Mx;
            boundaries_core_EM[boundary_num][state_num].My = My;
            
            for(unsigned int ii_4=0; ii_4<boundaries_EM_core[boundary_num].size(); ++ii_4)
            {
               //cout<<"qq1"<<endl;
               //cout<<"E_block_tot = "<<E_block_tot<<", Mx= "<<Mx<<", My= "<<My<<endl;
               
               //!исправить!!!!!!! работает с ошибками!!!!!!!!
               //проверить!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               //if(boundaries_EM_core[boundary_num][ii_4].E == E_block_tot && boundaries_EM_core[boundary_num][ii_4].Mx == Mx && boundaries_EM_core[boundary_num][ii_4].My == My)
               if(is_equal(boundaries_EM_core[boundary_num][ii_4].E, E_block_tot) && is_equal(boundaries_EM_core[boundary_num][ii_4].Mx, Mx) && is_equal(boundaries_EM_core[boundary_num][ii_4].My, My))
               {
                  //cout<<"qq2"<<endl;
                  boundaries_EM_core[boundary_num][ii_4].state_array.push_back(state_num);
                  //cout<<"E_block_tot = "<<E_block_tot<<", Mx= "<<Mx<<", My= "<<My<<endl;
                  //cout << "E1= " <<boundaries_EM_core[boundary_num][boundaries_EM_core[boundary_num].size()-1].E<<endl;
                  //system("pause");
                  break;
               }
               else if(ii_4 == boundaries_EM_core[boundary_num].size()-1)
               {
                  //cout<<"qq3"<<endl;
                  boundaries_EM_core[boundary_num].push_back(EM_state_struct(E_block_tot,Mx,My));
                  boundaries_EM_core[boundary_num][ii_4+1].state_array.push_back(state_num);
                  ii_4++;
                  //cout<<"E_block_tot = "<<E_block_tot<<", Mx= "<<Mx<<", My= "<<My<<endl;
                  //cout<<"E2= "<<boundaries_EM_core[boundary_num][boundaries_EM_core[boundary_num].size()-1].E<<endl;
                  //system("pause");
               }
               //system("pause");
               ///cout << "Mx = " << states_array[ii][ii_4].Mx << endl;
               ///cout << "My = " << states_array[ii][ii_4].My << endl;
               ///cout << "E = " << states_array[ii][ii_4].E << endl;
            }
         }
         
         ///system("pause");
         //здесь заканчивается перебор состояний в рамках границ
      }
      
      ///cout<<"Emin_array[boundary_num] = "<<Emin_array[boundary_num]<<endl;
      
      num_of_unique_EM+=boundaries_EM_core[boundary_num].size();
      
      /*
       if(boundaries_EM_core[boundary_num].size() == states_amount_of_core)
       {
          kk++;
          
          cout << "\nBoundary: " << boundary_num << endl;
          cout << "Number of unique E&Mx&My = " << boundaries_EM_core[boundary_num].size() << endl;
          
          //вывод E и M с учетом кратности вырождения
          for(unsigned int ii=0; ii<boundaries_EM_core[boundary_num].size(); ++ii)
          {
             cout<<"________\n";
             cout<<"E = "<<boundaries_EM_core[boundary_num][ii].E;
             cout<<", Mx = "<<boundaries_EM_core[boundary_num][ii].Mx;
             cout<<", My = "<<boundaries_EM_core[boundary_num][ii].My<<endl;
             cout<<"States: { ";
             for(unsigned int jj=0; jj<boundaries_EM_core[boundary_num][ii].state_array.size(); ++jj){
                if(jj<boundaries_EM_core[boundary_num][ii].state_array.size()-1)
                   cout<<boundaries_EM_core[boundary_num][ii].state_array[jj]<<", ";
                else
                   cout<<boundaries_EM_core[boundary_num][ii].state_array[jj]<<" ";
             }
             cout<<"}";
             cout<<endl;
          }
          
          cout << "\nBoundary: " << boundary_num << endl;
          cout << "Number of unique E&Mx&My = " << boundaries_EM_core[boundary_num].size() << endl;
          
          cout<<"________\n";
          
       }
       system("pause");
       //*/
      
      //здесь заканчивается перебор границ
   }
   
   //cout<<"kk = "<<kk<<endl;
   
   cout << "Number of boundaries = " << boundaries_EM_core.size() << endl;
   //cout << "Number of states with unique E&M = " << boundaries_EM_core[0][0].state_array.size() << endl;
   //cout << "Mx = " << boundaries_EM_core[0][0].Mx << endl;
   //cout << "My = " << boundaries_EM_core[0][0].My << endl;
   //cout << "E = " << boundaries_EM_core[0][0].E << endl;
   ///cout << "The number of all pairs EM = "<<boundaries_amount*states_amount<<endl;
   ///cout << "The number of unique pairs EM = "<<num_of_unique_EM<<endl;
   
   ///system("pause");
   
   
   /** 
     * на выходе получаем массив boundaries_EM_core[][].{E,Mx,My,state_array[]}, 
     * который хранит все конфигурации границ, структуру уникальных значений E,Mx,My, и соответствующие им конфигурации блоков
     */
   /** <<<<----- здесь заканчивается перебор всех конфигураций границ и блока ----- */
   
   
   
   
   
   /** --------------- подставляем блок в систему --------------->>>> */
   
   
   ofstream C_data("C.txt"); //файл теплоемкости
   ofstream chi_x_data("chi_x.txt"); //файл магнитной восприимчивости по оси x
   ofstream chi_y_data("chi_y.txt"); //файл магнитной восприимчивости по оси y
   
   ofstream E_aver_data("E_aver.txt"); //файл средней энергии
   ofstream E2_aver_data("E2_aver.txt"); //файл средней энергии в квадрате
   ofstream E4_aver_data("E4_aver.txt"); //файл средней энергии в 4 степени
   
   ofstream Mx_aver_data("Mx_aver.txt"); //файл средней намагниченности по оси X
   ofstream Mx2_aver_data("Mx2_aver.txt"); //файл средней намагниченности в квадрате по оси X
   ofstream Mx4_aver_data("Mx4_aver.txt"); //файл средней намагниченности в 4 степени по оси X
   
   ofstream My_aver_data("My_aver.txt"); //файл средней намагниченности по оси Y
   ofstream My2_aver_data("My2_aver.txt"); //файл средней намагниченности в квадрате по оси Y
   ofstream My4_aver_data("My4_aver.txt"); //файл средней намагниченности в 4 степени по оси Y
   
   
   double E_aver = 0, E2_aver = 0, E4_aver = 0; //средняя энергия, ее квадрат и 4 степень
   double Mx_pr_aver = 0, Mx_pr2_aver = 0, Mx_pr4_aver = 0; //средняя проекция намагниченности по оси X и ее квадрат и 4 степень
   double My_pr_aver = 0, My_pr2_aver = 0, My_pr4_aver = 0; //средняя проекция намагниченности оси Y, ее квадрат и 4 степень
   
   
   unsigned int set_of_states, total_set_of_states=0;
   
   uniform_real_distribution<double> distribution_real(0,1); //вещественное равномерное распределение
   
   int rand_state=0; //случайная конфигурация с одинаковыми параметрами EM
   
   double Z = 0; //статсумма
   double r0_1; //число от 0 до 1
   
   double E_aver_i = 0;
   double E2_aver_i = 0;
   double E4_aver_i = 0;
   
   double Mx_aver_i = 0;
   double Mx2_aver_i = 0;
   double Mx4_aver_i = 0;
   
   double My_aver_i = 0;
   double My2_aver_i = 0;
   double My4_aver_i = 0;
   
   double interval=0;//интервалы вероятностей
   
   unsigned int rand_hex_num; //случайный гексагон
   
   unsigned int count_of_hex = hex_array.size();
   
   
   double E_sys = 0, Mx_sys = 0, My_sys = 0;
   
   //расчитываем начальную энергию системы E_sys и намагниченность Mx и My
   for(unsigned int i=0; i<coorm.size(); i+=2)
   {
      mx_i = mxmy[i];
      my_i = mxmy[i+1];
      
      for(unsigned int neigh=0; neigh<neighbors[i/2].size(); ++neigh)
      //for(unsigned int j=i+2; j<coorm.size(); j+=2)
      {
         //cout<<neighbors[i/2][neigh].number/2<<", ";
         if(neighbors[i/2][neigh].num > i)
         {
            mx_j = mxmy[neighbors[i/2][neigh].num];
            my_j = mxmy[neighbors[i/2][neigh].num+1];
            
            X = coorm[i]-neighbors[i/2][neigh].x;
            Y = coorm[i+1]-neighbors[i/2][neigh].y;
            
            E_sys+= (mx_i*mx_j + my_i*my_j)/
                  sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y) -
                  3*(mx_i*X+my_i*Y)*(mx_j*X+my_j*Y)/
                  sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
         }
      }
      
      Mx_sys += mx_i;
      My_sys += my_i;
   }
   
   cout<<"E_sys = "<<E_sys<<", Mx_sys = "<<Mx_sys<<", My_sys = "<<My_sys<<endl;
   ///system("pause");
   
   
   //проход блока по системе
   ///for(float T = Tmin; T<Tmax; T+=Tstep) //цикл по температуре
   for(float T = Tmax; T>Tmin; T-=Tstep)
   {
      if(T<0.2) Tstep=0.001;
      else 
         if(T<1.1) Tstep=0.01;
      else Tstep=0.1;
      ///if(T>2.1 && T<2.5) Tstep=0.01;
      ///else Tstep=0.1;
      cout << "T = " << T << endl << "----------"<<endl;
      
      E_aver = 0;
      E2_aver = 0;
      Mx_pr_aver = 0;
      Mx_pr2_aver = 0;
      My_pr_aver = 0;
      My_pr2_aver = 0;
      
      /**
      //heating (прогрев)
      for(unsigned int MCS = 0; MCS<1000*coorm.size()/2; ++MCS)
      {
         rand_hex_num = rand()%count_of_hex;
         
         unsigned int boundary_dec=0; //конфигурация границы в десятичной системе
         int bit;
         for (unsigned int ii=0; ii<total_num_of_spins_in_boundaries; ++ii)
         {
            if (mxmy[bound_array[rand_hex_num][ii]+1]<0) 
               bit = 0;
            else 
               bit = 1;
            
            boundary_dec = (boundary_dec | bit);
            
            if (ii<total_num_of_spins_in_boundaries-1) 
               boundary_dec = boundary_dec<<1;
         }
         
         unsigned int core_dec=0; //конфигурация ядра в десятичной системе
         for (unsigned int ii=0; ii<total_num_of_spins_in_core; ++ii)
         {
            if (mxmy[hex_array[rand_hex_num][ii]+1]<0) 
               bit = 0;
            else 
               bit = 1;
            
            core_dec = (core_dec | bit);
            
            if (ii<total_num_of_spins_in_core-1)
               core_dec = core_dec<<1;
         }
         
         unsigned int num_of_unique_EM_in_boundary = boundaries_EM_core[boundary_dec].size(); //g
         
         double P[num_of_unique_EM_in_boundary]; //массив вероятностей энергий
         
         E_sys  -= boundaries_core_EM[boundary_dec][core_dec].E;
         Mx_sys -= boundaries_core_EM[boundary_dec][core_dec].Mx;
         My_sys -= boundaries_core_EM[boundary_dec][core_dec].My;
         
         Z = 0; //статсумма
         
         //считаем статсумму
         for(set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; ++set_of_states)
         {
            //вырождение * exp
            P[set_of_states] = boundaries_EM_core[boundary_dec][set_of_states].state_array.size()*
                  exp((double)(-(boundaries_EM_core[boundary_dec][set_of_states].E-Emin_array[boundary_dec])/T));
            Z += P[set_of_states];
         }
         
         r0_1 = distribution_real(generator);
         
         //перебираем интервалы вероятностей
         interval=0;
         flag=0;
         
         //считаем вероятности
         for(set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; ++set_of_states)
         {
            P[set_of_states] /= Z;
            interval+=P[set_of_states];
            if(r0_1<=interval && flag==0)
            {
               total_set_of_states = set_of_states;
               flag = 1;
            }
         }
         
         Mx_sys += boundaries_EM_core[boundary_dec][total_set_of_states].Mx;
         My_sys += boundaries_EM_core[boundary_dec][total_set_of_states].My;
         E_sys += boundaries_EM_core[boundary_dec][total_set_of_states].E;
         
         //выбираем случайную конфигурацию с одинаковыми параметрами
         rand_state = rand() % boundaries_EM_core[boundary_dec][total_set_of_states].state_array.size();
         
         //запоминаем конфигурацию
         rand_state = boundaries_EM_core[boundary_dec][total_set_of_states].state_array[rand_state];
         
         //переводим в двоичную систему и в +-1
         int value;
         count=0;
         for(bit=total_num_of_spins_in_core-1; bit>=0; --bit)
         {
            value = (1 & rand_state >> bit);
            
            if(value==0)
            {
               mxmy[hex_array[rand_hex_num][count]] = -sample_core[count].mx;
               mxmy[hex_array[rand_hex_num][count]+1] = -sample_core[count].my;
            }
            else
            {
               mxmy[hex_array[rand_hex_num][count]] = sample_core[count].mx;
               mxmy[hex_array[rand_hex_num][count]+1] = sample_core[count].my;
            }
            
            count++;
         }
         //здесь заканчивается прогрев 
      }
      cout<<"E_sys after heating = "<<E_sys<<", Mx_sys = "<<Mx_sys<<", My_sys = "<<My_sys<<endl;
      
      //*/
      
      
      
      /**
      /////////////////////////////////////////////////////////////////////////////////
      //!пересчитываем энергию системы E_sys и намагниченность Mx и My для актуализации
      E_sys = 0, Mx_sys = 0, My_sys = 0;
      
      for(unsigned int i=0; i<coorm.size(); i+=2)
      {
         mx_i = mxmy[i];
         my_i = mxmy[i+1];
         
         for(unsigned int neigh=0; neigh<neighbors[i/2].size(); ++neigh)
         {
            if(neighbors[i/2][neigh].num > i)
            {
               mx_j = mxmy[neighbors[i/2][neigh].num];
               my_j = mxmy[neighbors[i/2][neigh].num+1];
               
               X = coorm[i]-neighbors[i/2][neigh].x;
               Y = coorm[i+1]-neighbors[i/2][neigh].y;
               
               E_sys+= (mx_i*mx_j + my_i*my_j)/
                     sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y) -
                     3*(mx_i*X+my_i*Y)*(mx_j*X+my_j*Y)/
                     sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
            }
         }
         
         Mx_sys += mx_i;
         My_sys += my_i;
      }
      
      cout<<"correct E_sys = "<<E_sys<<", Mx_sys = "<<Mx_sys<<", My_sys = "<<My_sys<<endl;
      /////////////////////////////////////////////////////////////////////////////
      //*/
      
      //Monte Carlo
      for(unsigned int MCS = 0; MCS<NumMC*coorm.size()/2; ++MCS)
      {
         rand_hex_num = rand()%count_of_hex;
         
         //bound_array[номер гексагона][порядковый номер спина на границе] = номер спина в системе;
         
         unsigned int boundary_dec=0; //конфигурация границы в десятичной системе
         int bit;
         for (unsigned int ii=0; ii<total_num_of_spins_in_boundaries; ++ii)
         {
            if (mxmy[bound_array[rand_hex_num][ii]+1]<0) 
               bit = 0;
            else 
               bit = 1;
            
            boundary_dec = (boundary_dec | bit);
            
            if (ii<total_num_of_spins_in_boundaries-1) 
               boundary_dec = boundary_dec<<1;
         }
         ///cout<<"boundary_dec = "<<boundary_dec<<endl;
         
         
         //hex_array[номер гексагона][порядковый номер спина в гексагоне] = номер спина в системе;
         
         unsigned int core_dec=0; //конфигурация ядра в десятичной системе
         for (unsigned int ii=0; ii<total_num_of_spins_in_core; ++ii)
         {
            if (mxmy[hex_array[rand_hex_num][ii]+1]<0) 
               bit = 0;
            else 
               bit = 1;
            
            core_dec = (core_dec | bit);
            
            if (ii<total_num_of_spins_in_core-1)
               core_dec = core_dec<<1;
         }
         
         unsigned int num_of_unique_EM_in_boundary = boundaries_EM_core[boundary_dec].size(); //g
         
         double P[num_of_unique_EM_in_boundary]; //массив вероятностей энергий
         
         ///cout<<"E0 = "<<boundaries_core_EM[boundary_dec][core_dec].E<<", Mx0 = "<<boundaries_core_EM[boundary_dec][core_dec].Mx<<", My0 = "<<boundaries_core_EM[boundary_dec][core_dec].My<<endl;
         ///cout << "Emin = " << Emin_array[boundary_dec] << endl;
         
         E_sys  -= boundaries_core_EM[boundary_dec][core_dec].E;
         Mx_sys -= boundaries_core_EM[boundary_dec][core_dec].Mx;
         My_sys -= boundaries_core_EM[boundary_dec][core_dec].My;
         ///cout << "boundaries_core_EM[boundary_dec][core_dec].E = " << boundaries_core_EM[boundary_dec][core_dec].E << endl;
         ///cout << "E_sys before = " << E_sys << endl;
         
         Z = 0; //статсумма
         
         E_aver_i = 0;
         E2_aver_i = 0;
         E4_aver_i = 0;
         
         Mx_aver_i = 0;
         Mx2_aver_i = 0;
         Mx4_aver_i = 0;
         
         My_aver_i = 0;
         My2_aver_i = 0;
         My4_aver_i = 0;
         
         //считаем статсумму
         for(set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; ++set_of_states)
         {
            //вырождение * exp
            P[set_of_states] = boundaries_EM_core[boundary_dec][set_of_states].state_array.size()*
                  exp((double)(-(boundaries_EM_core[boundary_dec][set_of_states].E-Emin_array[boundary_dec])/T));
            Z += P[set_of_states];
            ///cout<<exp((double)(-(boundaries_EM_core[boundary_dec][set_of_states].E-Emin)/T))<<endl;
            /// 
            /// ПЕРЕНЕСТИ СЮДА РАСЧЕТ СРЕДНИХ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         }
         ///system("pause");
         
         ///cout << "Z = " << Z << endl;
         
         r0_1 = distribution_real(generator);
         ///cout<<"r0_1 = "<<r0_1<<endl;
         
         //перебираем интервалы вероятностей
         interval=0;
         flag=0;
         ///double sum=0;
         
         //cout<<"E_aver_i = "<<E_aver_i<<endl;
         
         double E_temp, Mx_temp, My_temp;
         
         //считаем вероятности
         for(set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; ++set_of_states)
         {
            P[set_of_states] /= Z;
            ///sum+=P[set_of_states];
            ///cout << "P"<<set_of_states<<"= " << P[set_of_states] << endl;
            
            /// 
            E_temp = E_sys+boundaries_EM_core[boundary_dec][set_of_states].E;
            Mx_temp = Mx_sys+boundaries_EM_core[boundary_dec][set_of_states].Mx;
            My_temp = My_sys+boundaries_EM_core[boundary_dec][set_of_states].My;
            /// 
            
            E_aver_i += P[set_of_states]*E_temp;
            Mx_aver_i += P[set_of_states]*Mx_temp;
            My_aver_i += P[set_of_states]*My_temp;
            
            ///cout<<"E_aver_i = "<<E_aver_i<<endl;
            
            E2_aver_i += P[set_of_states]*E_temp*E_temp;
            Mx2_aver_i += P[set_of_states]*Mx_temp*Mx_temp;
            My2_aver_i += P[set_of_states]*My_temp*My_temp;
            
            E4_aver_i += P[set_of_states]*E_temp*E_temp*E_temp*E_temp;
            Mx4_aver_i += P[set_of_states]*Mx_temp*Mx_temp*Mx_temp*Mx_temp;
            My4_aver_i += P[set_of_states]*My_temp*My_temp*My_temp*My_temp;
            
            interval+=P[set_of_states];
            if(r0_1<=interval && flag==0)
            {
               total_set_of_states = set_of_states;
               flag = 1;
            }
         }
         ///cout<<"E_aver_i = "<<E_aver_i<<endl;
         
         ///cout<<"P_sum = " << sum << endl;
         
         
         //средние
         E_aver += (E_aver_i - E_aver) / (double)(MCS + 1);
         //E2_aver += (E_aver_i*E_aver_i - E2_aver) / (double)(MCS + 1);
         E2_aver += (E2_aver_i - E2_aver) / (double)(MCS + 1);
         E4_aver += (E4_aver_i - E4_aver) / (double)(MCS + 1);
         
         Mx_pr_aver += (abs(Mx_aver_i) - Mx_pr_aver) / (double)(MCS + 1);
         //Mx_pr2_aver += ((Mx_aver_i*Mx_aver_i) - Mx_pr2_aver) / (double)(MCS + 1);
         Mx_pr2_aver += (Mx2_aver_i - Mx_pr2_aver) / (double)(MCS + 1);
         Mx_pr4_aver += (Mx4_aver_i - Mx_pr4_aver) / (double)(MCS + 1);
         
         My_pr_aver += (abs(My_aver_i) - My_pr_aver) / (double)(MCS + 1);
         My_pr2_aver += (My2_aver_i - My_pr2_aver) / (double)(MCS + 1);
         My_pr4_aver += (My4_aver_i - My_pr4_aver) / (double)(MCS + 1);
         
         ///cout<<"E_aver = "<<E_aver<<endl;
         ///system("pause");
         
         ///cout<<"r0_1<="<<interval<<endl;
         
         Mx_sys += boundaries_EM_core[boundary_dec][total_set_of_states].Mx;
         My_sys += boundaries_EM_core[boundary_dec][total_set_of_states].My;
         E_sys += boundaries_EM_core[boundary_dec][total_set_of_states].E;
         
         ///cout << "boundaries_EM_core[boundary_dec][total_set_of_states].E = " << boundaries_EM_core[boundary_dec][total_set_of_states].E << endl;
         ///cout << "E_sys after = " << E_sys << endl<<endl;
         
         ///system("pause");
         
         //cout<<"E_sys = "<<E_sys<<endl;
         
         //выбираем случайную конфигурацию с одинаковыми параметрами
         rand_state = rand() % boundaries_EM_core[boundary_dec][total_set_of_states].state_array.size();
         
         //запоминаем конфигурацию
         rand_state = boundaries_EM_core[boundary_dec][total_set_of_states].state_array[rand_state];
         
         
         
         //переводим в двоичную систему и в +-1
         ///cout<<"rand_state = "<<rand_state<<endl;
         
         ///ПЕРЕДЕЛАЛ!!! ПРОВЕРИТЬ!!!
         int value;
         count=0;
         for(bit=total_num_of_spins_in_core-1; bit>=0; --bit)
         {
            value = (1 & rand_state >> bit);
            
            
            //hex_array[номер гексагона][порядковый номер спина в гексагоне] = номер спина в системе;
            //mxmy[hex_array[rand_hex_num][count]] = sample_core[spin_num_in_core].mx;//x
            //mxmy[hex_array[rand_hex_num][count]+1] = sample_core[spin_num_in_core].my;//y
            
            if(value==0)
            {
               mxmy[hex_array[rand_hex_num][count]] = -sample_core[count].mx;
               mxmy[hex_array[rand_hex_num][count]+1] = -sample_core[count].my;
            }
            else
            {
               mxmy[hex_array[rand_hex_num][count]] = sample_core[count].mx;
               mxmy[hex_array[rand_hex_num][count]+1] = sample_core[count].my;
            }
            
            ///cout << core_temp[count]<< " ";
            count++;
         }
         ///cout<<endl;
         //здесь заканчивается проход блока по системе 
      }
      
      
      //cout<<"Mx_pr_aver = "<<Mx_pr_aver<<endl<<"My_pr_aver = "<<My_pr_aver<<endl;
      
      cout << "E_aver = " << E_aver << endl;
      
      double C = ((E2_aver - E_aver*E_aver)/(T*T))/(double)(coorm.size()/2);
      cout << "C = " << C << endl;
      C_data << T << "\t" << C << endl;
      //if(C1>C) Tstep /= 1.1; 
      //else Tstep *= 1.1; 
      ///C1=C;
      
      double chi_x = ((Mx_pr2_aver - Mx_pr_aver*Mx_pr_aver)/T) / (double)(coorm.size()/2);
      cout << "chi_x = " << chi_x << endl;
      chi_x_data << T << "\t" << chi_x << endl;
      
      double chi_y = ((My_pr2_aver - My_pr_aver*My_pr_aver)/T) / (double)(coorm.size()/2);
      cout << "chi_y = " << chi_y << endl << endl;
      chi_y_data << T << "\t" << chi_y << endl;
      
      E_aver_data << T << "\t" << E_aver << endl;
      E2_aver_data << T << "\t" << E2_aver << endl;
      E4_aver_data << T << "\t" << E4_aver << endl;
      
      Mx_aver_data << T << "\t" << Mx_pr_aver << endl;
      Mx2_aver_data << T << "\t" << Mx_pr2_aver << endl;
      Mx4_aver_data << T << "\t" << Mx_pr4_aver << endl;
      
      My_aver_data << T << "\t" << My_pr_aver << endl;
      My2_aver_data << T << "\t" << My_pr2_aver << endl;
      My4_aver_data << T << "\t" << My_pr4_aver << endl;
      
      ///system("pause");
      
      //здесь заканчивается перебор температур 
   }
   
   //записываем последнюю конфигурацию в файл
   ofstream GS_conf("GS_conf.mfsys");
   GS_conf<<"[header]"<<endl;
   GS_conf<<"dimensions=2"<<endl;
   GS_conf<<"size="<<coorm.size()/2<<endl;
   GS_conf<<"state=";
   for(unsigned int i=0; i<coorm.size(); i+=2) GS_conf<<"0";
   GS_conf<<endl;
   GS_conf<<"[parts]"<<endl;
   
   for(unsigned int i=0; i<coorm.size(); i+=2)
   {
      GS_conf<<i/2<<"\t"; //порядковый номер
      GS_conf<<coorm[i]<<"\t"<<coorm[i+1]<<"\t"<<"0"<<"\t"; //координаты x,y,z
      GS_conf<<mxmy[i]<<"\t"<<mxmy[i+1]<<"\t"<<"0"<<"\t"; //магнитные моменты mx,my,mz
      GS_conf<<"0"<<endl;
   }
   GS_conf.close();
   
   
   C_data.close();
   chi_x_data.close();
   chi_y_data.close();
   
   E_aver_data.close();
   E2_aver_data.close();
   E4_aver_data.close();
   
   Mx_aver_data.close();
   Mx2_aver_data.close();
   Mx4_aver_data.close();
   
   My_aver_data.close();
   My2_aver_data.close();
   My4_aver_data.close();
   
   /** <<<<--------------------------------------------------- */
   
   cout << "\n" << (double)clock()/CLOCKS_PER_SEC << " sec.\n";
   ///system("pause");
   
   cout << "\nThe End\n";
   return 0;
}
