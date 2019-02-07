#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <fstream>
using namespace std;

int main()
{
    ofstream fout("matrix_E.csv"); //хранится матрица взаимодействия для Пети
    ofstream foutc("coord.dat"); //коордтнаты точек
    ofstream fouts("spins.dat"); //магнитные моменты
    ofstream fout_centers("centers.dat"); //центры гексагонов
    ofstream foutn("neibours.dat"); //соседи
    
    int period=3;
    fout << "# Artificial Honeycobm Spin Ice Lattice "<<endl; // 
    fout << "# period = "<<period<<endl;
    
    vector <double> coorm; // array for scaled lattice
    vector <double> nnn; //array of neighbors
    vector <double> coor(24);//array of coordinates of unit cell число координат
    vector <double> mm(24); //array of magnetic moments of unit cel  число магнитных моментов
    vector <double> mxmy; //array of magnetic moments of scaled lattice
    vector <double> hex_centers(8);
    vector <double> temp_vector;
    
    double sqrt_3 = sqrt(3);
    
    hex_centers[0]=sqrt_3/2;     hex_centers[1]=0.75;
    hex_centers[2]=(3*sqrt_3)/2;     hex_centers[3]=0.75;
    hex_centers[4]=sqrt_3;     hex_centers[5]=2.25;
    hex_centers[6]=2*sqrt_3;     hex_centers[7]=2.25;
    
    
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
    mm[8]=-sqrt_3/2;     mm[9]=0.5; //0.25
    mm[10]=sqrt_3/2;    mm[11]=0.5; //1.25
    
    mm[12]=-sqrt_3/2;   mm[13]=0.5; //1.25
    mm[14]=sqrt_3/2;    mm[15]=0.5; //0.25
    mm[16]=-sqrt_3/2;    mm[17]=0.5; //0.25
    mm[18]=sqrt_3/2;    mm[19]=0.5; //1.25
    mm[20]=0;           mm[21]=1; //1.75
    mm[22]=0;           mm[23]=1; //1.75
    
    cout<<"Number of spins in the unit cell   = "<<coor.size()/2<<endl;
    fout<<"# Number of Spins in the Unit Cell = "<<coor.size()/2<<endl;
    
    cout<<"Period of lattice = "<<period<<endl;
    
    //check of coordinates
    cout<<"{";
    for(int i=0; i<coor.size(); i+=2)
    {
        cout<<"{"<<coor[i]<<","<<coor[i+1]<<"},";
    }
    cout<<"}\n";
    
    //check of magnetic moments
    cout<<"{";
    for(int i=0; i<mm.size(); i+=2)
    {
        cout<<"{"<<mm[i]<<","<<mm[i+1]<<"},";
    }
    cout<<"}\n";
    
    
    int numuc;
    cout<<"Input Number of Unit Cels in linear dimension "<<endl;
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
    cout<<"last spin="<<coorm[coorm.size()-2]<<endl;
    fout<<"# last spin="<<coorm[coorm.size()-2]<<endl;
    fout<<"# Number of spins in system ="<<(coorm.size()+1)/2<<endl;
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
    fouts<<"}\n\n";
    foutc.close();
    fouts.close();
    
    
    fout_centers<<"{";
    for(int i=0; i<hex_centers.size(); i+=2)
    {
        if(i!=hex_centers.size()-2)
            fout_centers<<"{"<<hex_centers[i]<<","<<hex_centers[i+1]<<"},";
        else
            fout_centers<<"{"<<hex_centers[i]<<","<<hex_centers[i+1]<<"}";
        
    }
    fout_centers<<"}\n\n";
    
    fout_centers.close();        
    
    
    //calculation of neighbors
    double r, X, Y; 
    int ii=0, jj=0, it=0, it1=0;
    cout<<"Input radius of coordination sphere"<<endl;
    cin>>r;
    fout<<"# Radius of coordination sphere ="<<r<<endl;
    
    //The creation of interaction energy matrix
    //double *e_matrix[(int)(coorm.size()+1)/2][(int)(coorm.size()+1)/2]; // array of pointers on energy
    //cout<<"(int)(coorm.size()+1)/2="<<(int)(coorm.size()+1)/2<<endl;
    //double**  e_m = new double*[(int)(coorm.size()];
    //for(int i = 0; i < coorm.size()/2; ++i)
    //    e_m[i] = new double[coorm.size()/2];
    
    
    cout<<"(coorm.size())="<<(coorm.size())<<endl;
    double  e_m[coorm.size()/2][coorm.size()/2]; // array of pointers on energy
    system("pause");
    
    double en;
    
    for(int ix=0; ix<coorm.size()/2; ix++)
    {
        for(int jy=0; jy<coorm.size()/2; jy++)
        {
            e_m[ix][jy]= 0.;
        }
    }
    
    cout<<"coorm.size()-1 = " << coorm.size()-1 << endl;
    for(int j=0; j<coorm.size()-1; j+=2)
    {
        cout << "j = " << j<<endl;
        //system("pause");	 
        //cout<<" new point "<<"{"<<coorm[j]<<","<<coorm[j+1]<<"} **************  r="<<r<<endl;
        
        it1=it;
        it=0;
        ii=j;
        
        //cout<<j<<endl;
        while (coorm[j]>=coorm[ii]-r &&  ii<coorm.size())// neighbors from right side over x    
        {
            if((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+
                    (coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1])<=r*r && ii!=j)//circular input condition 
            {
                it++;
                
                //	 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[ii]<<","<<coorm[ii+1]<<"},\n";
                
                en=(mxmy[j]*mxmy[ii]+mxmy[j+1]*mxmy[ii+1])/
                        sqrt((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+(coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1]))/
                        sqrt((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+(coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1]))/
                        sqrt((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+(coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1]))-
                        3*(mxmy[j]*(coorm[j]-coorm[ii])+mxmy[j+1]*(coorm[j+1]-coorm[ii+1]))*(mxmy[ii]*(coorm[j]-coorm[ii])+mxmy[ii+1]*(coorm[j+1]-coorm[ii+1]))/
                        sqrt((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+(coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1]))/
                        sqrt((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+(coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1]))/
                        sqrt((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+(coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1]))/
                        sqrt((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+(coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1]))/
                        sqrt((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+(coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1]));
                
                //		 cout<<"(j+1)/2="<<(j+1)/2<<endl;
                //		 cout<<"(ii+1)/2="<<(ii+1)/2<<endl;
                
                e_m[(j+1)/2][(ii+1)/2]=en;
                e_m[(ii+1)/2][(j+1)/2]=en;
                // cout<<"en="<<en<<endl;
                //cout<<"1"<<endl;
            }
            ii+=2;
        }
        
        //cout<<"01"<<endl;
        //system("pause");	
        ii=j;				
        if(coorm[ii]-r<0)// periodicity over x from left side 
            for(jj=coorm.size()-2; period*sqrt(numuc)-coorm[jj]<=r; jj-=2)		 
                if((coorm[j]+period*sqrt(numuc)-coorm[jj])*(coorm[j]+period*sqrt(numuc)-coorm[jj])+
                        (coorm[j+1]-coorm[jj+1])*(coorm[j+1]-coorm[jj+1])<=r*r && jj!=j)
                {
                    it++;
                    //	 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]+period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"2"<<endl;		
                    //	cout<<"en="<<en<<endl;
                    //	system("pause");
                }
        
        ii=j;						 
        if(coorm[j+1]-r<0)//periodicity over y from bottom
        {
            for(jj=j; coorm[jj]<coorm[j]+r; jj+=2)		 
            {
                if((coorm[j]-coorm[jj])*(coorm[j]-coorm[jj])+
                        (coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))<=r*r && jj!=j)
                {
                    
                    //	cout<<"periodicity over y"<<endl;
                    it++;
                    //	   cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-coorm[jj]);
                    Y=(coorm[j+1]+period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    //cout<<"3"<<endl;		 	
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"en="<<en<<endl;
                    //system("pause");
                }	 
            }
            
            for(jj=j; coorm[jj]>coorm[j]-r && jj>-1; jj-=2)		 
            {
                if((coorm[j]-coorm[jj])*(coorm[j]-coorm[jj])+
                        (coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))<=r*r && jj!=j)
                {
                    //	cout<<"periodicity over y"<<endl;
                    it++;
                    //	   cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-coorm[jj]);
                    Y=(coorm[j+1]+period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    //сout<<"4"<<endl;		 	
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"en="<<en<<endl;
                    //system("pause");		 
                }	 
            }
        }
        
        
        ii=j;
        while (coorm[ii]>=coorm[j]-r && ii>-1)//neighbors from left side over x
        {
            
            if((coorm[j]-coorm[ii])*(coorm[j]-coorm[ii])+
                    (coorm[j+1]-coorm[ii+1])*(coorm[j+1]-coorm[ii+1])<=r*r && ii!=j)
            {
                
                //cout<<"left "<<endl;
                it++;
                //	 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[ii]<<","<<coorm[ii+1]<<"},\n";
                
                X=(coorm[j]-coorm[ii]);
                Y=(coorm[j+1]-coorm[ii+1]);
                
                en=(mxmy[j]*mxmy[ii]+mxmy[j+1]*mxmy[ii+1])/
                        sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                        3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[ii]*X+mxmy[ii+1]*Y)/
                        sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                
                
                e_m[(j+1)/2][(ii+1)/2]=en;
                e_m[(ii+1)/2][(j+1)/2]=en;
                
                //cout<<"5"<<endl;
                //	 cout<<"en="<<en<<endl;
                //	 system("pause");
            }
            ii-=2;
        }
        
        
        if (coorm[j]+r>=period*sqrt(numuc))//cout<<"periodicity over x from right side "<<endl;
        {
            for(jj=0; coorm[j]-period*sqrt(numuc)+coorm[jj]<r; jj+=2)
            {
                if((coorm[j]-period*sqrt(numuc)-coorm[jj])*(coorm[j]-period*sqrt(numuc)-coorm[jj])+
                        (coorm[j+1]-coorm[jj+1])*(coorm[j+1]-coorm[jj+1])<=r*r && jj!=j)
                {
                    
                    //   cout<<"periodicity over x from right side "<<endl;
                    it++;
                    //cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"6"<<endl;
                    //cout<<"en="<<en<<endl;
                    //system("pause");	
                }
            }
        }
        
        
        
        if(coorm[j+1]+r>=period*sqrt(numuc) )//cout<<"periodicity over y from top "<<endl;
        {
            for(jj=j; coorm[j]-r<=coorm[jj] && jj>-1; jj-=2) 
            { //cout<<coorm[jj]<<"  "<<coorm[jj-2]<<endl;
                if((coorm[j]-coorm[jj])*(coorm[j]-coorm[jj])+
                        (coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))<=r*r && jj!=j)
                {
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-coorm[jj]);
                    Y=(coorm[j+1]-period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"6"<<endl;
                    // cout<<"en="<<en<<endl;
                    //system("pause");
                }
            }
            
            for(jj=j; coorm[j]+r>=coorm[jj] && jj<coorm.size(); jj+=2) //
            { //cout<<coorm[j]+r<<"   "<<coorm[jj]<<endl;
                if((coorm[j]-coorm[jj])*(coorm[j]-coorm[jj])+
                        (coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))<=r*r && jj!=j)
                {
                    
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-coorm[jj]);
                    Y=(coorm[j+1]-period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"7"<<endl;
                    //cout<<"en="<<en<<endl;
                    //system("pause");
                }
            }
        }
        
        
        if(coorm[j+1]+r>=period*sqrt(numuc) && coorm[j]-r<0 )//cout<<"periodicity over y from top and x from left "<<endl;
        {
            // cout<<coorm[j]+period*sqrt(numuc)-r<<"  "<<coorm[jj]<<endl;
            for(jj=coorm.size()-2; coorm[j]+period*sqrt(numuc)-r<=coorm[jj]; jj-=2) 
            { 		    
                if((coorm[j]-coorm[jj]+period*sqrt(numuc))*(coorm[j]-coorm[jj]+period*sqrt(numuc))+
                        (coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))<=r*r && jj!=j)
                {
                    
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]+period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]-period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    //cout<<"8"<<endl;
                    // cout<<"en="<<en<<endl;
                    //system("pause");
                }
            } 	   
        }
        
        
        
        if(coorm[j]+r>=sqrt(numuc)*2-1 && coorm[j+1]+r>=sqrt(numuc)*2-1) //cout<<"periodicity over y from top and x from right "<<endl;
        { //cout<<"coorm[j]+r="<<coorm[j]+r<<"  "<<"numuc="<<sqrt(numuc)*2-1<<endl;
            //cout<<"coorm[j+1]+r="<<coorm[j+1]+r<<"  "<<"numuc="<<sqrt(numuc)*2-1<<endl;
            //system("pause");
            for(jj=0; coorm[j]+r-(sqrt(numuc)*2-1)>=coorm[jj]; jj+=2) 
            { 		    
                if((coorm[j]-coorm[jj]-(sqrt(numuc)*2-1))*(coorm[j]-coorm[jj]-(sqrt(numuc)*2-1))+
                        (coorm[j+1]-coorm[jj+1]+(sqrt(numuc)*2-1))*(coorm[j+1]-coorm[jj+1]+(sqrt(numuc)*2-1))<=r*r && jj!=j)
                {
                    
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-(sqrt(numuc)*2-1)-coorm[jj]);
                    Y=(coorm[j+1]+(sqrt(numuc)*2-1)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"9"<<endl;
                    //	cout<<"en="<<en<<endl;
                    //	system("pause");
                }
            }
        }
        
        
        if(coorm[j]-r<=0 && coorm[j+1]-r<=0)//cout<<"periodicity over y from botom and x from left "<<endl;
        {
            //cout<<coorm[j]-period*sqrt(numuc)+r<<"  "<<coorm[jj]<<endl;       
            for(jj=coorm.size()-2; coorm[j]+period*sqrt(numuc)-r<=coorm[jj]; jj-=2) 
            {// cout<<"check point cout ={"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";	
                //system("pause");	    
                if((coorm[j]-coorm[jj]-period*sqrt(numuc))*(coorm[j]-coorm[jj]-period*sqrt(numuc))+
                        (coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))<=r*r && jj!=j)
                {
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]+period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    //cout<<"10"<<endl;		 	
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //		cout<<"en="<<en<<endl;
                    //		system("pause");
                }
                
                if((coorm[j]-coorm[jj]+period*sqrt(numuc))*(coorm[j]-coorm[jj]+period*sqrt(numuc))+
                        (coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))<=r*r && jj!=j)
                {
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]+period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]+period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    //cout<<"11"<<endl;		 	
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"en="<<en<<endl;
                    //system("pause");
                }
            }
            
            for(jj=0; coorm[jj]<=coorm[j]-period*sqrt(numuc)+r; jj+=2) 
            { 		    
                if((coorm[j]-coorm[jj]+period*sqrt(numuc))*(coorm[j]-coorm[jj]+period*sqrt(numuc))+
                        (coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))<=r*r && jj!=j)
                {
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]+period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]-period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    //cout<<"12"<<endl;		 	
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"en="<<en<<endl;
                    //system("pause");
                    
                }
                if((coorm[j]-coorm[jj]+period*sqrt(numuc))*(coorm[j]-coorm[jj]+period*sqrt(numuc))+
                        (coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))<=r*r && jj!=j)
                {
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]+period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]+period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    //cout<<"13"<<endl;		 	
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //cout<<"en="<<en<<endl;
                    //system("pause");
                }
            } 	   
        }
        
        if(coorm[j]+r>=period*sqrt(numuc) && coorm[j+1]+r>=period*sqrt(numuc) )//cout<<"periodicity over y from botom and x from left "<<endl;
        {
            //cout<<coorm[j]-period*sqrt(numuc)+r<<"  "<<coorm[jj]<<endl;
            
            for(jj=0; coorm[jj]<=coorm[j]-period*sqrt(numuc)+r; jj+=2) 
            {
                // cout<<"check point cout ={"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";	
                //system("pause");	    
                if((coorm[j]-coorm[jj]-period*sqrt(numuc))*(coorm[j]-coorm[jj]-period*sqrt(numuc))+
                        (coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]-period*sqrt(numuc))<=r*r && jj!=j)
                {
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]-period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    //cout<<"14"<<endl;		 	
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //	cout<<"en="<<en<<endl;
                    //	system("pause");
                }
            } 	   
        }
        
        
        if(coorm[j]+r>=period*sqrt(numuc) && coorm[j+1]-r<=0)//cout<<"periodicity over y from botom and x from right "<<endl;
        { 
            //cout<<coorm[j]-period*sqrt(numuc)+r<<"  "<<coorm[jj]<<endl;
            
            for(jj=0; coorm[jj]<=coorm[j]-period*sqrt(numuc)+r; jj+=2) 
            {// cout<<"check point cout ={"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";	
                //system("pause");	    
                if((coorm[j]-coorm[jj]-period*sqrt(numuc))*(coorm[j]-coorm[jj]-period*sqrt(numuc))+
                        (coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))*(coorm[j+1]-coorm[jj+1]+period*sqrt(numuc))<=r*r && jj!=j)
                {
                    //cout<<"periodicity over y from top "<<endl;
                    it++;
                    //		 cout<<"{"<<coorm[j]<<","<<coorm[j+1]<<"}," <<"{"<<coorm[jj]<<","<<coorm[jj+1]<<"},\n";
                    
                    X=(coorm[j]-period*sqrt(numuc)-coorm[jj]);
                    Y=(coorm[j+1]+period*sqrt(numuc)-coorm[jj+1]);
                    
                    en=(mxmy[j]*mxmy[jj]+mxmy[j+1]*mxmy[jj+1])/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)-
                            3*(mxmy[j]*X+mxmy[j+1]*Y)*(mxmy[jj]*X+mxmy[jj+1]*Y)/
                            sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y)/sqrt(X*X+Y*Y);
                    
                    //cout<<"15"<<endl;		 	
                    e_m[(j+1)/2][(jj+1)/2]=en;
                    e_m[(jj+1)/2][(j+1)/2]=en;
                    
                    //	cout<<"en="<<en<<endl;
                    //	system("pause");
                }
            } 	   
        }
        
        
        //cout<<"number of neibors="<<it<<endl;
        if(it1!=it)
        {
            foutn<<"# Neibors are ="<<it1<<endl;
            
            //system("pause");
        }
        else
        {
            foutn<<"# Neibors are ="<<it<<endl;
        }
    }
    
    cout<<"HELLO!!!!!!!!!!!!!!!!!!";
    
    for(int ix=0; ix<coorm.size()/2; ix++)
    {    for(int jy=0; jy<coorm.size()/2; jy++)
        {
            if((e_m[ix][jy]!=e_m[jy][ix] || e_m[ix][jy]!=0) && ix==jy) 
            {
                cout<< "matrix is worng!!!"<<endl;
                system("pause");}
            
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
            {
                //cout<<";";
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
    foutn.close();
    
    cout << "\nThe End\n";
    return 0;
}

