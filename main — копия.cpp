#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using namespace std;

int main()
{
    ofstream fout("ASSI_27.csv"); //хранится матрица взаимодействия для Пети
    ofstream foutc("ASSI_27_coord.dat"); //координаты точек
    ofstream fouts("ASSI_27_spins.dat"); //магнитные моменты
    
    int period=3;
    fout << "# Artificial Square Spin Ice Lattice "<<endl; // 
    fout << "# period = "<<period<<endl;
    
    vector <float> coorm; // array for scaled lattice
    vector <float> coor(6);//array of coordinates of unit cell число координат
    vector <float> mm(6); //array of magnetic moments of unit cel  число магнитных моментов
    vector <float> mxmy; //array of magnetic moments of scaled lattice
    
    ///конфигурация с ненулевой энтропией
    coor[0]=0;      coor[1]=1; 
    coor[2]=0.5;    coor[3]=0; 
    coor[4]=1;      coor[5]=1; 
    
    /*
    //SSI
    coor[0]=0.5;      coor[1]=0; 
    coor[2]=1.5;    coor[3]=0; 
    coor[4]=0;      coor[5]=0.5; 
    coor[6]=1;      coor[7]=0.5; 
    coor[8]=2;    coor[9]=0.5; 
    coor[10]=0.5;      coor[11]=1; 
    
    coor[12]=1.5;      coor[13]=1; 
    coor[14]=0;    coor[15]=1.5; 
    coor[16]=1;      coor[17]=1.5; 
    coor[18]=2;      coor[19]=1.5; 
    coor[20]=0.5;    coor[21]=2; 
    coor[22]=1.5;      coor[23]=2; 
    //*/
    
    //unit cell magnetic moments
    ///конфигурация с ненулевой энтропией
    //*
    mm[0]=0.5;  mm[1]=sqrt(1-mm[0]*mm[0]); 
    mm[2]=1;    mm[3]=0; 
    mm[4]=0.5;  mm[5]=-sqrt(1-mm[4]*mm[4]);
    //*/
    
    /*
    //SSI
    mm[0]=1;      mm[1]=0; 
    mm[2]=1;      mm[3]=0; 
    mm[4]=0;      mm[5]=1;
    mm[6]=0;      mm[7]=1; 
    mm[8]=0;      mm[9]=1; 
    mm[10]=1;     mm[11]=0;
    
    mm[12]=1;     mm[13]=0; 
    mm[14]=0;     mm[15]=1; 
    mm[16]=0;     mm[17]=1;
    mm[18]=0;     mm[19]=1; 
    mm[20]=1;     mm[21]=0; 
    mm[22]=1;     mm[23]=0;
    */
    
    cout<<"Number of spins in the unit cell   = "<<coor.size()/2<<endl;
    fout<<"# Number of Spins in the Unit Cell = "<<coor.size()/2<<endl;
    
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
    cout<<"\nInput Number of Unit Cels in linear dimension "<<endl;
    cin>>numuc;
    fout<<"# Number of Unit Cels in linear dimension = "<<numuc<<endl;
    cout<<"Lenear size of system in unit cell's spins ="<<numuc*coor.size()/2<<" spins"<<endl;
    fout<<"# Lenear size of system = "<<(numuc*coor.size()/2)<<endl;
    cout<<"Number of spins in system in linear ="<<numuc*coor.size()/2<<endl;
    
    numuc=numuc*numuc;
    cout<<"Total number of unit cells ="<<numuc<<endl;
    
    //if((int)sqrt(numuc*coor.size())!=sqrt(numuc*coor.size()))
    //goto begin1; //case of nonsquare sample 
    
    if(sqrt(numuc)==1)
        coorm=coor;
    
    vector<float> repit;
    vector<float> repitm;
    //scaling over y and ordering of array of coorm and magnetic moments
    for(unsigned int y=0;y<coor.size(); y+=2)
    {for(int j=0; j<sqrt(numuc); j++ )
        {if(coor[y]!=coor[y+2])
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
                
                for(unsigned int i=0; i<repit.size(); i++)
                {coorm.push_back(repit[i]);
                    mxmy.push_back(repitm[i]);
                }
                for(int k=j+1; k<sqrt(numuc); k++)
                    for(unsigned int i=0; i<repit.size(); i+=2)
                    { coorm.push_back(repit[i]);  coorm.push_back(repit[i+1]+period*k);
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
        for(unsigned int x=0;x<coor.size(); x+=2)
        {
            coorm.push_back (coor[x]+period*j);
            coorm.push_back (coor[x+1]);
            
            mxmy.push_back(mm[x]);
            mxmy.push_back(mm[x+1]); 
        }
    
    cout<<"check of coordinates\n";
    foutc<<"{";
    for(unsigned int i=0; i<coorm.size(); i+=2)
    {if(i!=coorm.size()-2)
            foutc<<"{"<<coorm[i]<<","<<coorm[i+1]<<"},";
        else
            foutc<<"{"<<coorm[i]<<","<<coorm[i+1]<<"}";
    }
    foutc<<"}\n\n";
    
    fout<<"# Number of spins in system ="<<(coorm.size()+1)/2<<endl;
    cout<<"  Number of spins in system ="<<(coorm.size()+1)/2<<endl;
    
    cout<<"check of magnetic moments\n";
    fouts<<"{";
    for(unsigned int i=0; i<coorm.size(); i+=2)
    {if(i!=coorm.size()-2)
            fouts<<"{{"<<coorm[i]-mxmy[i]/2<<","<<coorm[i+1]-mxmy[i+1]/2<<"},{"<<coorm[i]+mxmy[i]/2<<","<<coorm[i+1]+mxmy[i+1]/2<<"}},";
        else
            fouts<<"{{"<<coorm[i]-mxmy[i]/2<<","<<coorm[i+1]-mxmy[i+1]/2<<"},{"<<coorm[i]+mxmy[i]/2<<","<<coorm[i+1]+mxmy[i+1]/2<<"}}";
        
    }
    fouts<<"}\n\n";
    
    
    
    float X, Y; 
    
    //The creation of interaction energy matrix
    //float *e_matrix[(int)(coorm.size()+1)/2][(int)(coorm.size()+1)/2]; // array of pointers on energy
    //cout<<"(int)(coorm.size()+1)/2="<<(int)(coorm.size()+1)/2<<endl;
    //float**  e_m = new float*[(int)(coorm.size()];
    //for(int i = 0; i < coorm.size()/2; ++i)
    //    e_m[i] = new float[coorm.size()/2];
    
    
    cout<<"(coorm.size())="<<(coorm.size())<<endl;
    float  e_m[coorm.size()/2][coorm.size()/2]; // array of pointers on energy
    system("pause");
    
    float en;//zerro=0, one=1.,;
    
    for(unsigned int ix=0; ix<coorm.size()/2; ix++)
        for(unsigned int jy=0; jy<coorm.size()/2; jy++)
        {
            //e_matrix[ix][jy]= &zerro;
            e_m[ix][jy]= 0.;
        }
    
    for(unsigned int j=0; j<coorm.size(); j+=2 )
    {
        for(int i=0; i<coorm.size(); i+=2 )
            if(j!=i)
            {
                //	 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                
                X=(coorm[j]-coorm[i]);
                Y=(coorm[j+1]-coorm[i+1]);
                
                en=(mxmy[j]*mxmy[i]+mxmy[j+1]*mxmy[i+1])/
                        sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                        3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[i]*X+mxmy[i+1]*Y)/
                        sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                
                e_m[j/2][i/2]=en;
                e_m[i/2][j/2]=en;
                
                //cout<< j/2 << " "<<( i+1)/2 <<endl;
                //cout<<"2"<<endl;		
                //	cout<<"en="<<en<<endl;
                //	system("pause");
            }
    }
    
    
    for(unsigned int ix=0; ix<coorm.size()/2; ix++)
    {    for(unsigned int jy=0; jy<coorm.size()/2; jy++)
        {
            if(e_m[ix][jy]!=e_m[jy][ix] || (e_m[ix][jy]!=0 && ix==jy)) 
            {cout<< "matrix is worng!!!"<<endl;
                //system("pause");
            }
            
            if(e_m[ix][jy]!=0)
            {
                //cout<<*e_matrix[ix][jy]<<endl;
                fout<< e_m[ix][jy];
                if(jy<coorm.size()/2-1)
                {
                    fout<<";";
                }
            }
            else
            {//cout<<";";
                if(jy<coorm.size()/2-1)
                {
                    fout<<";";
                }
            }
        }
        fout<<endl;
    }
    
    
    //begin1: //case of nonsquare sample 
    fout.close();
    
    cout << "\nThe End\n";
    return 0;
}
