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
    
    int period=3;
    
    vector <double> coorm; // array for scaled lattice
    vector <double> neighbors; //array of neighbors
    vector <double> coor(24);//array of coordinates of unit cell число координат
    vector <double> mm(24); //array of magnetic moments of unit cel  число магнитных моментов
    vector <double> mxmy; //array of magnetic moments of scaled lattice
    vector <double> hex_centers(8);
    vector <double> temp_vector;
    
    double sqrt_3 = sqrt(3);
    double a = sqrt_3; //параметр решетки
    
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
    
    cout<<"Number of spins in the unit cell   = "<<coor.size()/2<<endl;
    
    cout<<"Period of lattice = "<<period<<endl;
    
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
    
    
    int numuc;
    cout<<"Input Number of Unit Cels in linear dimension "<<endl;
    cin>>numuc;
    cout<<"Lenear size of system in unit cell's spins ="<<numuc*coor.size()/2<<" spins"<<endl;
    cout<<"Number of spins in system in linear ="<<numuc*coor.size()/2<<endl;
    
    numuc=numuc*numuc;
    cout<<"Total number of unit cells ="<<numuc<<endl;
    
    //if((int)sqrt(numuc*coor.size())!=sqrt(numuc*coor.size()))
    //goto begin1; //case of nonsquare sample 
    
    if(sqrt(numuc)==1)
        coorm=coor;
    vector<double> repit;
    vector<double> repitm;
    
    cout<<"hex_centers.size()_0 = "<<hex_centers.size()<<endl;
    for(int y=0;y<hex_centers.size(); y+=2)
    {
        for(int j=0; j<sqrt(numuc); j++ )
        {
            if(hex_centers[y]!=hex_centers[y+2])
            {
                temp_vector.push_back(hex_centers[y]);
                temp_vector.push_back(hex_centers[y+1]+period*j);
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
                
                for(int i=0; i<repit.size(); i++)
                {
                    temp_vector.push_back(repit[i]);
                }
                
                for(int k=j+1; k<sqrt(numuc); k++)
                {
                    for(int i=0; i<repit.size(); i+=2)
                    { 
                        temp_vector.push_back(repit[i]);  temp_vector.push_back(repit[i+1]+period*k);
                    }
                }
                j=sqrt(numuc);
                repit.erase(repit.begin(),repit.end());
            }
        }
    }
    hex_centers=temp_vector;
    cout<<"hex_centers.size()_1 = "<<hex_centers.size()<<endl;
    
    //scaling over y and ordering of array of coorm and magnetic moments
    for(int y=0;y<coor.size(); y+=2)
    {
        for(int j=0; j<sqrt(numuc); j++ )
        {
            if(coor[y]!=coor[y+2])
            {
                coorm.push_back (coor[y]);
                coorm.push_back (coor[y+1]+period*j);
                
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
                
                for(int i=0; i<repit.size(); i++)
                {
                    coorm.push_back(repit[i]);
                    mxmy.push_back(repitm[i]);
                }
                for(int k=j+1; k<sqrt(numuc); k++)
                    for(int i=0; i<repit.size(); i+=2)
                    { 
                        coorm.push_back(repit[i]);  coorm.push_back(repit[i+1]+period*k);
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
        for(int x=0;x<coor.size(); x+=2)
        {
            coorm.push_back(coor[x]+(period+0.5)*j);
            coorm.push_back(coor[x+1]);
            
            mxmy.push_back(mm[x]);
            mxmy.push_back(mm[x+1]); 
        }
        for(int x=0;x<hex_centers.size(); x+=2)
        {
            temp_vector.push_back(hex_centers[x]+(period+0.5)*j);
            temp_vector.push_back(hex_centers[x+1]);
        }
    }
    hex_centers=temp_vector;
    cout<<"hex_centers.size()_2 = "<<hex_centers.size()<<endl;
    
    cout<<"check of coordinates\n";
    foutc<<"{";
    for(int i=0; i<coorm.size(); i+=2)
    {
        if(i!=coorm.size()-2)
            foutc<<"{"<<coorm[i]<<","<<coorm[i+1]<<"},";
        else
            foutc<<"{"<<coorm[i]<<","<<coorm[i+1]<<"}";
    }
    foutc<<"}\n\n";
    foutc.close();
    
    
    cout<<"last spin="<<coorm[coorm.size()-2]<<endl;
    cout<<"  Number of spins in system ="<<(coorm.size()+1)/2<<endl;
    
    cout<<"check of magnetic moments\n";
    fouts<<"{";
    for(int i=0; i<coorm.size(); i+=2)
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
    double r, max_X, max_Y, min_X, min_Y; 
    //double cs = 3/4; //coordination sphere in square
    //cout<<"Input radius of coordination sphere"<<endl;
    //cin>>r;
    
    vector <vector <int> > hex_array(hex_centers.size()/2);
    cout<<"hex_array.size() = " << hex_array.size() <<endl;
    cout<<"coorm.size() = " << coorm.size() <<endl;
    
    unsigned int k=0;
    int num_spins_in_core=0;
    for(unsigned int x=0; x<hex_centers.size(); x+=2)
    {
        for(unsigned int i=0; i<coorm.size(); i+=2)
        {
            r=(hex_centers[x]-coorm[i])*(hex_centers[x]-coorm[i])+
                    (hex_centers[x+1]-coorm[i+1])*(hex_centers[x+1]-coorm[i+1]);
            if(r<=((double)3/4)+0.0001)
            {
                hex_array[k].push_back(coorm[i]);
                hex_array[k].push_back(coorm[i+1]);
                num_spins_in_core++;
                //cout<<"hex["<<k<<"] = ";
                //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                if(num_spins_in_core == 6)
                {
                    //cout<<"hex["<<k<<"] = ";
                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                    break;
                }
            }
            r=(hex_centers[x]-(hex_centers[x]+coorm[i]))*(hex_centers[x]-(hex_centers[x]+coorm[i]))+
                    (hex_centers[x+1]-coorm[i+1])*(hex_centers[x+1]-coorm[i+1]);
            if(r<=((double)3/4)+0.0001)
            {
                hex_array[k].push_back(coorm[i]);
                hex_array[k].push_back(coorm[i+1]);
                num_spins_in_core++;
                if(num_spins_in_core == 6)
                {
                    //cout<<"hex["<<k<<"] = ";
                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                    break;
                }
            }
            
            r=(hex_centers[x]-coorm[i])*(hex_centers[x]-coorm[i])+
                    (hex_centers[x+1]-(hex_centers[x+1]+coorm[i+1]))*(hex_centers[x+1]-(hex_centers[x+1]+coorm[i+1]));
            if(r<=((double)3/4)+0.0001)
            {
                hex_array[k].push_back(coorm[i]);
                hex_array[k].push_back(coorm[i+1]);
                num_spins_in_core++;
                if(num_spins_in_core == 6)
                {
                    //cout<<"hex["<<k<<"] = ";
                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                    break;
                }
            }
            
            r=(hex_centers[x]-(hex_centers[x]+coorm[i]))*(hex_centers[x]-(hex_centers[x]+coorm[i]))+
                    (hex_centers[x+1]-(hex_centers[x+1]+coorm[i+1]))*(hex_centers[x+1]-(hex_centers[x+1]+coorm[i+1]));
            if(r<=((double)3/4)+0.0001)
            {
                hex_array[k].push_back(coorm[i]);
                hex_array[k].push_back(coorm[i+1]);
                num_spins_in_core++;
                if(num_spins_in_core == 6)
                {
                    //cout<<"hex["<<k<<"] = ";
                    //cout<<"num_spins_in_core = "<<num_spins_in_core<<endl;
                    break;
                }
            }
        }
        num_spins_in_core=0;
        cout<<"hex["<<hex_centers[x]<<","<<hex_centers[x+1]<<"] = ";
        cout<<hex_array[k].size()<<endl;
        k++;
    }
    
    
    system("pause");
    
    
    cout << "\nThe End\n";
    return 0;
}

