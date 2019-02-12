#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <fstream>
using namespace std;

int main()
{
    ofstream foutc("coord.dat"); //коордтнаты точек
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
    
    int numuc;
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
    
    
    //calculation of neighbors
    double r, max_X_of_centers, max_Y_of_centers, min_Y_of_centers; 
    //double cs = 3/4; //coordination sphere in square
    //cout<<"Input radius of coordination sphere"<<endl;
    //cin>>r;
    
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
    
    vector <vector <double> > hex_array(hex_centers.size()/2);
    vector <vector <double> > bound_array(hex_centers.size()/2);
    
    cout<<"Number of hexagons = " << hex_array.size() <<endl;
    
    int num_spins_in_core=0;
    int num_spins_on_bound=0;
    
    for(unsigned int x=0; x<hex_centers.size(); x+=2)
    {
        for(unsigned int i=0; i<coorm.size(); i+=2)
        {
            //определяем спины, которые входят в гексагон и границу
            r=(hex_centers[x]-coorm[i])*(hex_centers[x]-coorm[i])+
                    (hex_centers[x+1]-coorm[i+1])*(hex_centers[x+1]-coorm[i+1]);
            
            //для гексагона
            if(r<=((double)3/4)+0.0001)
            {
                hex_array[x/2].push_back(i);
                num_spins_in_core++;
                //cout<<"hex["<<x/2<<"] = ";
                //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                if(num_spins_in_core == 6 && num_spins_on_bound == 24)
                {
                    //cout<<"hex["<<x/2<<"] = ";
                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                    break;
                }
            }
            else
            {
                //для границы
                if(r<=((double)27/4)+0.0001)
                {
                    bound_array[x/2].push_back(i);
                    num_spins_on_bound++;
                    if(num_spins_in_core == 6 && num_spins_on_bound == 24)
                    {
                        //cout<<"hex["<<x/2<<"] = ";
                        //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                        break;
                    }
                }
            }
            
            if(i==coorm.size()-2 && (num_spins_in_core<6 || num_spins_on_bound<24))
            {
                for(unsigned int ii=0; ii<coorm.size(); ii+=2)
                {
                    //сдвигаем по X
                    r=(hex_centers[x]-(max_X_of_centers+coorm[ii]))*(hex_centers[x]-(max_X_of_centers+coorm[ii]))+
                            (hex_centers[x+1]-coorm[ii+1])*(hex_centers[x+1]-coorm[ii+1]);
                    
                    if(r<=((double)3/4)+0.0001)
                    {
                        hex_array[x/2].push_back(ii);
                        num_spins_in_core++;
                        if(num_spins_in_core == 6 && num_spins_on_bound == 24)
                        {
                            //cout<<"hex["<<x/2<<"] = ";
                            //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                            break;
                        }
                    }
                    else
                    {
                        if(r<=((double)27/4)+0.0001)
                        {
                            bound_array[x/2].push_back(ii);
                            num_spins_on_bound++;
                            if(num_spins_in_core == 6 && num_spins_on_bound == 24)
                            {
                                //cout<<"hex["<<x/2<<"] = ";
                                //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                break;
                            }
                        }
                    }
                    
                    if(ii==coorm.size()-2 && (num_spins_in_core<6 || num_spins_on_bound<24))
                    {
                        for(unsigned int iii=0; iii<coorm.size(); iii+=2)
                        {
                            //сдвигаем по Y
                            r=(hex_centers[x]-coorm[iii])*(hex_centers[x]-coorm[iii])+
                                    (hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iii+1]))*(hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iii+1]));
                            
                            if(r<=((double)3/4)+0.0001)
                            {
                                hex_array[x/2].push_back(iii);
                                num_spins_in_core++;
                                if(num_spins_in_core == 6 && num_spins_on_bound == 24)
                                {
                                    //cout<<"hex["<<x/2<<"] = ";
                                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                    break;
                                }
                            }
                            else
                            {
                                if(r<=((double)27/4)+0.0001)
                                {
                                    bound_array[x/2].push_back(iii);
                                    num_spins_on_bound++;
                                    if(num_spins_in_core == 6 && num_spins_on_bound == 24)
                                    {
                                        //cout<<"hex["<<x/2<<"] = ";
                                        //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                        break;
                                    }
                                }
                            }
                            
                            if(iii==coorm.size()-2 && (num_spins_in_core<6 || num_spins_on_bound<24))
                            {
                                for(unsigned int iiii=0; iiii<coorm.size(); iiii+=2)
                                {
                                    //сдвигаем по диагонали (X и Y)
                                    r=(hex_centers[x]-(max_X_of_centers+coorm[iiii]))*(hex_centers[x]-(max_X_of_centers+coorm[iiii]))+
                                            (hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iiii+1]))*(hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iiii+1]));
                                    
                                    if(r<=((double)3/4)+0.0001)
                                    {
                                        hex_array[x/2].push_back(iiii);
                                        num_spins_in_core++;
                                        if(num_spins_in_core == 6 && num_spins_on_bound == 24)
                                        {
                                            //cout<<"hex["<<x/2<<"] = ";
                                            //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if(r<=((double)27/4)+0.0001)
                                        {
                                            bound_array[x/2].push_back(iiii);
                                            num_spins_on_bound++;
                                            if(num_spins_in_core == 6 && num_spins_on_bound == 24)
                                            {
                                                //cout<<"hex["<<x/2<<"] = ";
                                                //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                break;
                                            }
                                        }
                                    }
                                    
                                    if(iiii==coorm.size()-2 && num_spins_on_bound<24)
                                    {
                                        for(unsigned int iiiii=0; iiiii<coorm.size(); iiiii+=2)
                                        {
                                            //сдвигаем по -X
                                            r=(hex_centers[x]-(-max_X_of_centers+coorm[iiiii]))*(hex_centers[x]-(-max_X_of_centers+coorm[iiiii]))+
                                                    (hex_centers[x+1]-coorm[iiiii+1])*(hex_centers[x+1]-coorm[iiiii+1]);
                                            
                                            if(r>((double)3/4)+0.0001 && r<=((double)27/4)+0.0001)
                                            {
                                                bound_array[x/2].push_back(iiiii);
                                                num_spins_on_bound++;
                                                if(num_spins_on_bound == 24)
                                                {
                                                    //cout<<"hex["<<x/2<<"] = ";
                                                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                    break;
                                                }
                                            }
                                            
                                            if(iiiii==coorm.size()-2 && num_spins_on_bound<24)
                                            {
                                                for(unsigned int iiiiii=0; iiiiii<coorm.size(); iiiiii+=2)
                                                {
                                                    //сдвигаем по Y
                                                    r=(hex_centers[x]-coorm[iiiiii])*(hex_centers[x]-coorm[iiiiii])+
                                                            (hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiii+1]))*(hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiii+1]));
                                                    
                                                    if(r>((double)3/4)+0.0001 && r<=((double)27/4)+0.0001)
                                                    {
                                                        bound_array[x/2].push_back(iiiiii);
                                                        num_spins_on_bound++;
                                                        if(num_spins_on_bound == 24)
                                                        {
                                                            //cout<<"hex["<<x/2<<"] = ";
                                                            //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                            break;
                                                        }
                                                    }
                                                    
                                                    if(iiiiii==coorm.size()-2 && num_spins_on_bound<24)
                                                    {
                                                        for(unsigned int iiiiiii=0; iiiiiii<coorm.size(); iiiiiii+=2)
                                                        {
                                                            //сдвигаем по диагонали (-X и -Y)
                                                            r=(hex_centers[x]-(-max_X_of_centers+coorm[iiiiiii]))*(hex_centers[x]-(-max_X_of_centers+coorm[iiiiiii]))+
                                                                    (hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiii+1]))*(hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiii+1]));
                                                            
                                                            if(r>((double)3/4)+0.0001 && r<=((double)27/4)+0.0001)
                                                            {
                                                                bound_array[x/2].push_back(iiiiiii);
                                                                num_spins_on_bound++;
                                                                if(num_spins_on_bound == 24)
                                                                {
                                                                    //cout<<"hex["<<x/2<<"] = ";
                                                                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                                    break;
                                                                }
                                                            }
                                                            
                                                            if(iiiiiii==coorm.size()-2 && num_spins_on_bound<24)
                                                            {
                                                                for(unsigned int iiiiiiii=0; iiiiiiii<coorm.size(); iiiiiiii+=2)
                                                                {
                                                                    //сдвигаем по диагонали (-X и Y)
                                                                    r=(hex_centers[x]-(-max_X_of_centers+coorm[iiiiiiii]))*(hex_centers[x]-(-max_X_of_centers+coorm[iiiiiiii]))+
                                                                            (hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iiiiiiii+1]))*(hex_centers[x+1]-(max_Y_of_centers+min_Y_of_centers+coorm[iiiiiiii+1]));
                                                                    
                                                                    if(r>((double)3/4)+0.0001 && r<=((double)27/4)+0.0001)
                                                                    {
                                                                        bound_array[x/2].push_back(iiiiiiii);
                                                                        num_spins_on_bound++;
                                                                        if(num_spins_on_bound == 24)
                                                                        {
                                                                            //cout<<"hex["<<k<<"] = ";
                                                                            //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                                            break;
                                                                        }
                                                                    }
                                                                    
                                                                    if(iiiiiiii==coorm.size()-2 && num_spins_on_bound<24)
                                                                    {
                                                                        for(unsigned int iiiiiiiii=0; iiiiiiiii<coorm.size(); iiiiiiiii+=2)
                                                                        {
                                                                            //сдвигаем по диагонали (X и -Y)
                                                                            r=(hex_centers[x]-(max_X_of_centers+coorm[iiiiiiiii]))*(hex_centers[x]-(max_X_of_centers+coorm[iiiiiiiii]))+
                                                                                    (hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiiiii+1]))*(hex_centers[x+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[iiiiiiiii+1]));
                                                                            
                                                                            if(r>((double)3/4)+0.0001 && r<=((double)27/4)+0.0001)
                                                                            {
                                                                                bound_array[x/2].push_back(iiiiiiiii);
                                                                                num_spins_on_bound++;
                                                                                if(num_spins_on_bound == 24)
                                                                                {
                                                                                    //cout<<"hex["<<k<<"] = ";
                                                                                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                                                                                    break;
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
        num_spins_in_core=0;
        num_spins_on_bound=0;
        
        //cout<<"hex "<<x/2<<" ["<<hex_centers[x]<<","<<hex_centers[x+1]<<"] = ";
        //cout<<hex_array[x/2].size()/2<<endl;
        //cout<<bound_array[x/2].size()/2<<endl;
        
        //проверка на правильное количество спинов в каждом гексагоне
        if(hex_array[x/2].size()!=6)
        {
            cout<<endl<<"Wrong hex size!!!!!!!!!!!!!!!!!!!!!"<<endl;
            cout<<"hex "<<x/2<<" ["<<hex_centers[x]<<", "<<hex_centers[x+1]<<"] = ";
            cout<<hex_array[x/2].size()<<endl<<endl;
        }
        
        //проверка на правильное количество спинов на границе каждого гексагона
        if(bound_array[x/2].size()!=24)
        {
            cout<<endl<<"Wrong border size!!!!!!!!!!!!!!!!!!!!!"<<endl;
            cout<<"hex "<<x/2<<" ["<<hex_centers[x]<<", "<<hex_centers[x+1]<<"] = ";
            cout<<bound_array[x/2].size()<<endl<<endl;
        }
        
        ///вывод на экран всех спинов, входящих в гексагон
        //for(unsigned int c=0; c<hex_array[x/2].size(); c+=2) cout<<hex_array[x/2][c]<<", "<<hex_array[x/2][c+1]<<endl;
    }
    
    
    struct spin_struct
    {
        unsigned int number;
        double x;
        double y;
        
        //конструктор
        spin_struct(unsigned int number, double x, double y){spin_struct::number = number; spin_struct::x = x; spin_struct::y = y; }
    };

    vector <vector <spin_struct> > neighbors(coorm.size()/2); //array of neighbors
    int num_of_neighbors=0;
    
    cout<<"\nneighbors\n";
    //поиск соседей
    for(unsigned int i=0; i<coorm.size(); i+=2)
    {
        for(unsigned int j=0; j<coorm.size(); j+=2)
        {
            r=(coorm[i]-coorm[j])*(coorm[i]-coorm[j])+
                    (coorm[i+1]-coorm[j+1])*(coorm[i+1]-coorm[j+1]);
            
            if(i!=j && r<=4*a+0.0001)
            {
                neighbors[i/2].push_back(spin_struct(j,coorm[j],coorm[j+1]));
                num_of_neighbors++;
                if(num_of_neighbors == 14)
                {
                    //cout<<"neighbors["<<i/2<<"] = ";
                    //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                    break;
                }
            }
            
            if(j==coorm.size()-2 && num_of_neighbors<14)
            {
                for(unsigned int jj=0; jj<coorm.size(); jj+=2)
                {
                    //сдвигаем по X
                    r=(coorm[i]-(max_X_of_centers+coorm[jj]))*(coorm[i]-(max_X_of_centers+coorm[jj]))+
                            (coorm[i+1]-coorm[jj+1])*(coorm[i+1]-coorm[jj+1]);
                    
                    if(i!=jj && r<=4*a+0.0001)
                    {
                        neighbors[i/2].push_back(spin_struct(jj,max_X_of_centers+coorm[jj],coorm[jj+1]));
                        num_of_neighbors++;
                        if(num_of_neighbors == 14)
                        {
                            //cout<<"neighbors["<<i/2<<"] = ";
                            //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                            break;
                        }
                    }
                    
                    if(jj==coorm.size()-2 && num_of_neighbors<14)
                    {
                        for(unsigned int jjj=0; jjj<coorm.size(); jjj+=2)
                        {
                            //сдвигаем по -X
                            r=(coorm[i]-(-max_X_of_centers+coorm[jjj]))*(coorm[i]-(-max_X_of_centers+coorm[jjj]))+
                                    (coorm[i+1]-coorm[jjj+1])*(coorm[i+1]-coorm[jjj+1]);
                            
                            if(i!=jjj && r<=4*a+0.0001)
                            {
                                neighbors[i/2].push_back(spin_struct(jjj, -max_X_of_centers+coorm[jjj], coorm[jjj+1]));
                                num_of_neighbors++;
                                if(num_of_neighbors == 14)
                                {
                                    //cout<<"neighbors["<<i/2<<"] = ";
                                    //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                    break;
                                }
                            }
                            
                            if(jjj==coorm.size()-2 && num_of_neighbors<14)
                            {
                                for(unsigned int jjjj=0; jjjj<coorm.size(); jjjj+=2)
                                {
                                    //сдвигаем по Y
                                    r=(coorm[i]-coorm[jjjj])*(coorm[i]-coorm[jjjj])+
                                            (coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjj+1]))*(coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjj+1]));
                                    
                                    if(i!=jjjj && r<=4*a+0.0001)
                                    {
                                        neighbors[i/2].push_back(spin_struct(jjjj, coorm[jjjj], max_Y_of_centers+min_Y_of_centers+coorm[jjjj+1]));
                                        num_of_neighbors++;
                                        if(num_of_neighbors == 14)
                                        {
                                            //cout<<"neighbors["<<i/2<<"] = ";
                                            //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                            break;
                                        }
                                    }
                                    
                                    
                                    if(jjjj==coorm.size()-2 && num_of_neighbors<14)
                                    {
                                        for(unsigned int jjjjj=0; jjjjj<coorm.size(); jjjjj+=2)
                                        {
                                            //сдвигаем по -Y
                                            r=(coorm[i]-coorm[jjjjj])*(coorm[i]-coorm[jjjjj])+
                                                    (coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjj+1]))*(coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjj+1]));
                                            
                                            if(i!=jjjjj && r<=4*a+0.0001)
                                            {
                                                neighbors[i/2].push_back(spin_struct(jjjjj, coorm[jjjjj], -max_Y_of_centers-min_Y_of_centers+coorm[jjjjj+1]));
                                                num_of_neighbors++;
                                                if(num_of_neighbors == 14)
                                                {
                                                    //cout<<"neighbors["<<i/2<<"] = ";
                                                    //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                                    break;
                                                }
                                            }
                                            
                                            if(jjjjj==coorm.size()-2 && num_of_neighbors<14)
                                            {
                                                for(unsigned int jjjjjj=0; jjjjjj<coorm.size(); jjjjjj+=2)
                                                {
                                                    //сдвигаем по диагонали (X и Y)
                                                    r=(coorm[i]-(max_X_of_centers+coorm[jjjjjj]))*(coorm[i]-(max_X_of_centers+coorm[jjjjjj]))+
                                                            (coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjjjj+1]))*(coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjjjj+1]));
                                                    
                                                    if(i!=jjjjjj && r<=4*a+0.0001)
                                                    {
                                                        neighbors[i/2].push_back(spin_struct(jjjjjj, max_X_of_centers+coorm[jjjjjj], max_Y_of_centers+min_Y_of_centers+coorm[jjjjjj+1]));
                                                        num_of_neighbors++;
                                                        if(num_of_neighbors == 14)
                                                        {
                                                            //cout<<"neighbors["<<i/2<<"] = ";
                                                            //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                                            break;
                                                        }
                                                    }
                                                    
                                                    if(jjjjjj==coorm.size()-2 && num_of_neighbors<14)
                                                    {
                                                        for(unsigned int jjjjjjj=0; jjjjjjj<coorm.size(); jjjjjjj+=2)
                                                        {
                                                            //сдвигаем по диагонали (-X и -Y)
                                                            r=(coorm[i]-(-max_X_of_centers+coorm[jjjjjjj]))*(coorm[i]-(-max_X_of_centers+coorm[jjjjjjj]))+
                                                                    (coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjj+1]))*(coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjj+1]));
                                                            
                                                            if(i!=jjjjjjj && r<=4*a+0.0001)
                                                            {
                                                                neighbors[i/2].push_back(spin_struct(jjjjjjj, -max_X_of_centers+coorm[jjjjjjj], -max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjj+1]));
                                                                num_of_neighbors++;
                                                                if(num_of_neighbors == 14)
                                                                {
                                                                    //cout<<"neighbors["<<i/2<<"] = ";
                                                                    //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                                                    break;
                                                                }
                                                            }
                                                            
                                                            if(jjjjjjj==coorm.size()-2 && num_of_neighbors<14)
                                                            {
                                                                for(unsigned int jjjjjjjj=0; jjjjjjjj<coorm.size(); jjjjjjjj+=2)
                                                                {
                                                                    //сдвигаем по диагонали (X и -Y)
                                                                    r=(coorm[i]-(max_X_of_centers+coorm[jjjjjjjj]))*(coorm[i]-(max_X_of_centers+coorm[jjjjjjjj]))+
                                                                            (coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjjj+1]))*(coorm[i+1]-(-max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjjj+1]));
                                                                    
                                                                    if(i!=jjjjjjjj && r<=4*a+0.0001)
                                                                    {
                                                                        neighbors[i/2].push_back(spin_struct(jjjjjjjj, max_X_of_centers+coorm[jjjjjjjj], -max_Y_of_centers-min_Y_of_centers+coorm[jjjjjjjj+1]));
                                                                        num_of_neighbors++;
                                                                        if(num_of_neighbors == 14)
                                                                        {
                                                                            //cout<<"neighbors["<<i/2<<"] = ";
                                                                            //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                                                            break;
                                                                        }
                                                                    }
                                                                    
                                                                    if(jjjjjjjj==coorm.size()-2 && num_of_neighbors<14)
                                                                    {
                                                                        for(unsigned int jjjjjjjjj=0; jjjjjjjjj<coorm.size(); jjjjjjjjj+=2)
                                                                        {
                                                                            //сдвигаем по диагонали (-X и Y)
                                                                            r=(coorm[i]-(-max_X_of_centers+coorm[jjjjjjjjj]))*(coorm[i]-(-max_X_of_centers+coorm[jjjjjjjjj]))+
                                                                                    (coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjjjjjjj+1]))*(coorm[i+1]-(max_Y_of_centers+min_Y_of_centers+coorm[jjjjjjjjj+1]));
                                                                            
                                                                            if(i!=jjjjjjjjj && r<=4*a+0.0001)
                                                                            {
                                                                                neighbors[i/2].push_back(spin_struct(jjjjjjjjj, -max_X_of_centers+coorm[jjjjjjjjj], max_Y_of_centers+min_Y_of_centers+coorm[jjjjjjjjj+1]));
                                                                                num_of_neighbors++;
                                                                                if(num_of_neighbors == 14)
                                                                                {
                                                                                    //cout<<"neighbors["<<i/2<<"] = ";
                                                                                    //cout<<"num_of_neighbors = "<<num_of_neighbors<<endl;
                                                                                    break;
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
        
        num_of_neighbors=0;
        
        //cout << i/2 << ": " << neighbors[i/2].size()/2 << endl;
        
        //проверка на правильное количество соседей
        if(neighbors[i/2].size()!=14)
        {
            cout<<endl<<"Wrong number of neighbors!!!!!!!!!!!!!!!!!!!!!"<<endl;
            cout<<"spin "<<i/2<<" ["<<coorm[i]<<", "<<coorm[i+1]<<"] = ";
            cout<<neighbors[i/2].size()<<endl<<endl;
        }
    }
    
    //проверяем спин на правильность соседей
    /*
    cout<<"spin "<<47<<" ["<<coorm[47*2]<<", "<<coorm[47*2+1]<<"] = ";
    cout<<neighbors[47].size()<<" neighbors"<<endl<<endl;
    for(unsigned int i=0; i<neighbors[47].size(); ++i)
    {
        //выводим номера соседних спинов
        cout<<i<<": "<< neighbors[47][i].number/2;
        //выводим координаты соседних спинов
        cout<<" ["<<neighbors[47][i].x<<", "<<neighbors[47][i].y<<"]";
        cout<<" (["<<coorm[neighbors[47][i].number]<<", "<<coorm[neighbors[47][i].number+1]<<"])\n";
    }
    //*/
    
    
    
    //расчет энергии
    
    
    
    //system("pause");
    
    
    cout << "\nThe End\n";
    return 0;
}

