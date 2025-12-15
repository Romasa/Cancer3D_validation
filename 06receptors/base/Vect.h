#ifndef VECT_H
#define VECT_H

#ifndef WWWW
#define WWWW
#define WWWWW(a) {if(((INOUT)?fread(&(a),sizeof(a),1,file):fwrite(&(a),sizeof(a),1,file))==0)return false;}
#endif //WWWW

class Vect3D
{
 public:
  Vect3D(){x=0; y=0; z=0;}
  Vect3D(double x_, double y_, double z_){x=x_; y=y_; z=z_;}
  ~Vect3D(){};

  void operator= (Vect3D v) {x=v.x; y=v.y; z=v.z;}
  void operator= (double t) {x=t; y=t; z=t;}
  Vect3D operator+ (Vect3D v) {return Vect3D(x+v.x,y+v.y,z+v.z);}
  void operator+= (Vect3D v) {x+=v.x; y+=v.y; z+=v.z;}
  Vect3D operator- (Vect3D v) {return Vect3D(x-v.x,y-v.y,z-v.z);}
  void operator-= (Vect3D v) {x-=v.x; y-=v.y; z-=v.z;}
  Vect3D operator* (double lambda) {return Vect3D(lambda*x,lambda*y,lambda*z);}
  Vect3D operator/ (double lambda) {return Vect3D(x/lambda,y/lambda,z/lambda);}

  bool ReadVect(bool INOUT, FILE *file)
  {
    WWWWW(x)
    WWWWW(y)
    WWWWW(z)

}
  
  double x;
  double y;
  double z;
};

#endif //VECT_H