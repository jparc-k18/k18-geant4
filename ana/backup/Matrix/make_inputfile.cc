// Author : Hitoshi Sugimura
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <sstream>
#include <fstream>
#include <map>
#include <cmath>

std::vector<double> list10,list11,list12,list13;
void input_mag_curve();
void makeorbitfile(double *mag,char *h);

int main(int argc,char **argv)
{
  if(argc<3)
    {
      std::cout<<"Usage:"
	       <<"./[exefile] [statusfile] [D4field]"
	       <<std::endl;
      exit(-1);
    }

  std::cout<<"magnet status file : "<<argv[1]<<std::endl;
  std::ifstream file(argv[1]);
  
  if(file.fail())
    {
      std::cerr <<argv[1]<<" file open fail ... (>_<)"<<std::endl;
      exit(-1);
    }
  std::map<std::string,double> mag_status;
  while(file.good())
    {
      std::string buf;
      std::getline(file,buf);
      if(buf[0]=='#')continue;
      std::istringstream is(buf);
      std::istream_iterator<std::string> issBegin(is);
      std::istream_iterator<std::string> issEnd;
      std::vector<std::string> param(issBegin,issEnd);

      if(param.empty() || param[0].empty())continue;
      std::cout<<"key = "<<param[0]
	       <<" value = "<<param[1]<<std::endl;
      mag_status[param[0]] = atof(param[1].c_str());
    }
  file.close();
  input_mag_curve();
  int i10 = 0;
  int i11 = 0;
  int i12 = 0;
  int i13 = 0;
  while((fabs(mag_status["Q10"])- i10*100)>0)
    {
      ++i10;
    }
  while((fabs(mag_status["Q11"])- i11*100)>0)
    {
      ++i11;
    }
  while((fabs(mag_status["Q12"])- i12*100)>0)
    {
      ++i12;
    }
  while((fabs(mag_status["Q13"])- i13*100)>0)
    {
      ++i13;
    }

  double diff[4],mag[4],magkg[4];

  diff[0] = list10[21 - i10] - list10[21 - i10 -1];
  mag[0] = fabs(list10[20-i10])+diff[0]*(fabs(mag_status["Q10"])-i10*100)/100;
  magkg[0] = mag[0]/500;

  diff[1] = list11[21 - i11] - list11[21 - i11 -1];
  mag[1] = fabs(list11[20-i11])+diff[1]*(fabs(mag_status["Q11"])-i11*100)/100;
  magkg[1] = mag[1]/500;

  diff[2] = list12[21 - i12] - list12[21 - i12 -1];
  mag[2] = fabs(list12[20-i12])+diff[2]*(fabs(mag_status["Q12"])-i12*100)/100;
  magkg[2] = mag[2]/500;

  diff[3] = list13[21 - i13] - list13[21 - i13 -1];
  mag[3] = fabs(list13[20-i13])+diff[3]*(fabs(mag_status["Q13"])-i13*100)/100;
  magkg[3] = mag[3]/500;

  std::cout<<"///// Magnetic Field /////"<<std::endl;
  std::cout<<"Q10 : "<<magkg[0]<<" [kG]"<<std::endl;
  std::cout<<"Q11 : "<<magkg[1]<<" [kG]"<<std::endl;
  std::cout<<"Q12 : "<<magkg[2]<<" [kG]"<<std::endl;
  std::cout<<"Q13 : "<<magkg[3]<<" [kG]"<<std::endl;
  makeorbitfile(magkg,argv[2]);
  //  makeorbitfile();
  return 0;
}

void input_mag_curve()
{
  std::ifstream file10("mag_curve10");
  std::ifstream file11("mag_curve11");
  std::ifstream file12("mag_curve12");
  std::ifstream file13("mag_curve13");

  if(file10.fail())
    {
      std::cerr <<" file open fail ... (>_<)"<<std::endl;
      exit(-1);
    }
  if(file11.fail())
    {
      std::cerr <<" file open fail ... (>_<)"<<std::endl;
      exit(-1);
    }
  if(file12.fail())
    {
      std::cerr <<" file open fail ... (>_<)"<<std::endl;
      exit(-1);
    }
  if(file13.fail())
    {
      std::cerr <<" file open fail ... (>_<)"<<std::endl;
      exit(-1);
    }

  while(file10.good())
    {
      double param[2];
      std::string buf;
      std::getline(file10,buf);
      if(buf[0]=='#')continue;
      std::istringstream is(buf);
      is>>param[0]>>param[1];
      list10.push_back(param[1]);
    }
  file10.close();

  while(file11.good())
    {
      double param[2];
      std::string buf;
      std::getline(file11,buf);
      if(buf[0]=='#')continue;
      std::istringstream is(buf);
      is>>param[0]>>param[1];
      list11.push_back(param[1]);
    }
  file11.close();

  while(file12.good())
    {
      double param[2];
      std::string buf;
      std::getline(file12,buf);
      if(buf[0]=='#')continue;
      std::istringstream is(buf);
      is>>param[0]>>param[1];
      list12.push_back(param[1]);
    }
  file12.close();

  while(file13.good())
    {
      double param[2];
      std::string buf;
      std::getline(file13,buf);
      if(buf[0]=='#')continue;
      std::istringstream is(buf);
      is>>param[0]>>param[1];
      list13.push_back(param[1]);
    }
  file13.close();
}

void makeorbitfile(double *mag,char *h)
{
  std::ifstream infile("K18MatrixDesignFFpFocus.in");
  std::ofstream outfile("test.out");
  double B = atof(h);
  double Brho = 4*B;
  double fg[4];
  for(int i=0;i<4;++i)
    {
      fg[i] = mag[i]/Brho;
    }
#if 0
  std::cout<<"mag "<<mag[0]<<std::endl;
  std::cout<<"mag "<<mag[1]<<std::endl;
  std::cout<<"mag "<<mag[2]<<std::endl;
  std::cout<<"mag "<<mag[3]<<std::endl;
  std::cout<<h<<std::endl;
  std::cout<<B<<std::endl;
  std::cout<<Brho<<std::endl;
  std::cout<<"fg : "<<fg[0]<<std::endl;
  std::cout<<"fg : "<<fg[1]<<std::endl;
  std::cout<<"fg : "<<fg[2]<<std::endl;
  std::cout<<"fg : "<<fg[3]<<std::endl;
#endif
  std::stringstream ss[4];
  ss[0]<<fg[0];
  ss[1]<<fg[1]*(-1);
  ss[2]<<fg[2];
  ss[3]<<fg[3]*(-1);
  std::vector<std::vector<std::string> > hogehoge;
  hogehoge.resize(32);
  int s=0;
  while(infile.good())
    {
      std::string buf;
      std::getline(infile,buf);
      if(buf[0]=='#')continue;
      std::istringstream is(buf);
      std::istream_iterator<std::string> issBegin(is);
      std::istream_iterator<std::string> issEnd;
      std::vector<std::string> param(issBegin,issEnd);
      for(unsigned int i=0;i<param.size();++i)
	{
	  hogehoge[s].push_back(param[i]);
	}
      ++s;
    }
  infile.close();
  for(int j=0;j<hogehoge.size();++j)
    {
      for(int i=0;i<hogehoge[j].size();++i)
	{
	  if(hogehoge[j][i]=="Q10")
	    hogehoge[j+1][i]=ss[0].str();
	  if(hogehoge[j][i]=="Q11")
	    hogehoge[j+1][i]=ss[1].str();
	  if(hogehoge[j][i]=="Q12")
	    hogehoge[j+1][i]=ss[2].str();
	  if(hogehoge[j][i]=="Q13")
	    hogehoge[j+1][i]=ss[3].str();
	    
	  outfile<<hogehoge[j][i]<<" ";
	}
      outfile<<std::endl;
    }
  outfile.close();
}
