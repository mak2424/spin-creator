#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using namespace std;

int main()
{
    ofstream fout("matrix_E.csv"); //хранится матрица взаимодействия для Пети
    ofstream foutc("coord.dat"); //координаты точек
    ofstream fouts("spins.dat"); //магнитные моменты
    ofstream fout_centers("centers.dat"); //центры гексагонов
    
    int period=3;
    fout << "# Artificial Honeycobm Spin Ice Lattice "<<endl; // 
    fout << "# period = "<<period<<endl;
    
    vector <double> coorm; // array for scaled lattice
    vector <double> coor(24);//array of coordinates of unit cell число координат
    vector <double> mm(24); //array of magnetic moments of unit cel  число магнитных моментов
    vector <double> mxmy; //array of magnetic moments of scaled lattice
    
    vector <double> hex_centers(8);
    
    double sqrt_3 = sqrt(3);
    
    hex_centers[0]=0;     hex_centers[1]=1.5;
    hex_centers[2]=sqrt_3;     hex_centers[3]=1.5;
    hex_centers[4]=sqrt_3/2;     hex_centers[5]=0;
    hex_centers[6]=(3*sqrt_3)/2;     hex_centers[7]=0;
    
    //hc
    //      X            Y
    //1.73205 = sqrt_3
    //1.29904 = (3*sqrt_3)/4
    //0.433013 = sqrt_3/4
    //0.866025 = sqrt_3/2
    //2.59808 = (3*sqrt_3)/2
    //3.03109 = (7*sqrt_3)/4
    //2.16506 = (5*sqrt_3)/4
    
    coor[0]=sqrt_3;    coor[1]=0; 
    coor[2]=(3*sqrt_3)/4;    coor[3]=0.75; 
    coor[4]=sqrt_3/4;   coor[5]=0.75; 
    coor[6]=0;          coor[7]=0; 
    coor[8]=sqrt_3/4;   coor[9]=-0.75; 
    coor[10]=(3*sqrt_3)/4;   coor[11]=-0.75; 
    
    coor[12]=(7*sqrt_3)/4;   coor[13]=0.75; 
    coor[14]=(5*sqrt_3)/4;   coor[15]=0.75; 
    coor[16]=(5*sqrt_3)/4;   coor[17]=-0.75; 
    coor[18]=(7*sqrt_3)/4;   coor[19]=-0.75; 
    coor[20]=sqrt_3/2;  coor[21]=1.5; 
    coor[22]=(3*sqrt_3)/2;   coor[23]=1.5; 
    
    //unit cell magnetic moments
    
    //hc
    mm[0]=0;            mm[1]=1; 
    mm[2]=-sqrt_3/2;    mm[3]=0.5; 
    mm[4]=-sqrt_3/2;    mm[5]=-0.5;
    mm[6]=0;            mm[7]=-1; 
    mm[8]=sqrt_3/2;     mm[9]=-0.5; 
    mm[10]=sqrt_3/2;    mm[11]=0.5;
    
    mm[12]=-sqrt_3/2;   mm[13]=0.5; 
    mm[14]=-sqrt_3/2;   mm[15]=-0.5; 
    mm[16]=sqrt_3/2;    mm[17]=-0.5;
    mm[18]=sqrt_3/2;    mm[19]=0.5; 
    mm[20]=0;           mm[21]=1; 
    mm[22]=0;           mm[23]=1;
    
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
    
    vector<double> repit;
    vector<double> repitm;
    
    
    unsigned int hex_centers_size = hex_centers.size();
    for(unsigned int y=0;y<hex_centers_size; y+=2)
    {
        for(int j=0; j<sqrt(numuc); j++ )
        {
            if(hex_centers[y]!=hex_centers[y+2])
            {
                hex_centers.push_back(hex_centers[y]);
                hex_centers.push_back(hex_centers[y+1]+period*j);
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
                    hex_centers.push_back(repit[i]);
                }
                for(int k=j+1; k<sqrt(numuc); k++)
                {
                    for(unsigned int i=0; i<repit.size(); i+=2)
                    { 
                        hex_centers.push_back(repit[i]);  hex_centers.push_back(repit[i+1]+period*k);
                    }
                }
                j=sqrt(numuc);
                repit.erase(repit.begin(),repit.end());
            }
        }
    }
    
    //scaling over y and ordering of array of coorm and magnetic moments
    for(unsigned int y=0;y<coor.size(); y+=2)
    {
        for(int j=0; j<sqrt(numuc); j++ )
        {
            if(coor[y]!=coor[y+2])
            {
                coorm.push_back(coor[y]);
                coorm.push_back(coor[y+1]+period*j);
                
                mxmy.push_back(mm[y]);
                mxmy.push_back(mm[y+1]); 
                
            }
            else
            {  
                repit.push_back (coor[y]);
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
                {
                    for(unsigned int i=0; i<repit.size(); i+=2)
                    { 
                        coorm.push_back(repit[i]);  coorm.push_back(repit[i+1]+period*k);
                        mxmy.push_back(repitm[i]);  mxmy.push_back(repitm[i+1]);
                    }
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
    
    
    hex_centers_size = hex_centers.size();
    
    for(int j=1; j<sqrt(numuc); j++ )
    {
        for(unsigned int x=0;x<coor.size(); x+=2)
        {
            coorm.push_back(coor[x]+(period+0.5)*j);
            coorm.push_back(coor[x+1]);
            
            mxmy.push_back(mm[x]);
            mxmy.push_back(mm[x+1]); 
        }
        for(unsigned int x=0;x<hex_centers_size; x+=2)
        {
            hex_centers.push_back(hex_centers[x]+(period+0.5)*j);
            hex_centers.push_back(hex_centers[x+1]);
        }
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
    {
        if(i!=coorm.size()-2)
            fouts<<"{{"<<coorm[i]-mxmy[i]/2<<","<<coorm[i+1]-mxmy[i+1]/2<<"},{"<<coorm[i]+mxmy[i]/2<<","<<coorm[i+1]+mxmy[i+1]/2<<"}},";
        else
            fouts<<"{{"<<coorm[i]-mxmy[i]/2<<","<<coorm[i+1]-mxmy[i+1]/2<<"},{"<<coorm[i]+mxmy[i]/2<<","<<coorm[i+1]+mxmy[i+1]/2<<"}}";
        
    }
    fouts<<"}\n\n";
    
    
    
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
    
    
    double X, Y; 
    
    //The creation of interaction energy matrix
    //double *e_matrix[(int)(coorm.size()+1)/2][(int)(coorm.size()+1)/2]; // array of pointers on energy
    //cout<<"(int)(coorm.size()+1)/2="<<(int)(coorm.size()+1)/2<<endl;
    //double**  e_m = new double*[(int)(coorm.size()];
    //for(int i = 0; i < coorm.size()/2; ++i)
    //    e_m[i] = new double[coorm.size()/2];
    
    
    cout<<"(coorm.size())="<<(coorm.size())<<endl;
    double  e_m[coorm.size()/2][coorm.size()/2]; // array of pointers on energy
    system("pause");
    
    double en;//zerro=0, one=1.,;
    
    for(unsigned int ix=0; ix<coorm.size()/2; ix++)
        for(unsigned int jy=0; jy<coorm.size()/2; jy++)
        {
            //e_matrix[ix][jy]= &zerro;
            e_m[ix][jy]= 0.;
        }
    
    for(unsigned int j=0; j<coorm.size(); j+=2 )
    {
        for(unsigned int i=0; i<coorm.size(); i+=2 )
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
