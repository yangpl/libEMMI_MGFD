#include "cstd.h"
  
int main(int argc, char *argv[])
{
  int i1, i2, i3;
  float ox, oy, oz;
  float dx, dy, dz;
  float x, y, z, lenx, leny, lenz, val, depth;
  int nx, ny, nz;
  float ***rho, **bathy;
  FILE *fp;
  static float PI = 3.1415926;


  nx = 100;
  ny = 100;
  nz = 100;
  dx = 200;
  dy = 200;
  dz = 40;
  ox = -10000.;
  oy = -10000.;
  oz = 0.;
  lenx = nx*dx;
  leny = ny*dy;
  lenz = nz*dz;

  rho = alloc3float(nx, ny, nz);
  bathy = alloc2float(nx, ny);

  for(i2=0; i2<ny; i2++){
    y = oy + (i2+0.5)*dy;
    for(i1=0; i1<nx; i1++){
      x = ox + (i1+0.5)*dx;
      depth = 900. + 100.*sin(2*PI*x/(3*lenx)+PI/3.)*sin(3*PI*y/(2*leny)-PI/3.);
      bathy[i2][i1] = depth;
    }
  }
  fp = fopen("fbathy", "wb");
  fwrite(&bathy[0][0], nx*ny*sizeof(float), 1, fp);
  fclose(fp);

  for(i3=0; i3<nz; i3++){
    z = oz + (i3+0.5)*dz;
    for(i2=0; i2<ny; i2++){
      y = oy + (i2+0.5)*dy;
      for(i1=0; i1<nx; i1++){
	x = ox + (i1+0.5)*dx;
	
	depth = 900. + 100.*sin(2*PI*x/(3*lenx)+PI/3.)*sin(3*PI*y/(2*leny)-PI/3.);
	if(z<depth) val = 0.3;
	else val = 1.5 + 5*(z-depth)/lenz;

	//we create two resistor in the depth
	if(z>=2400 &&z<2500 && x>-5500 && x<-2000 && y>-1000 && y<3500) val = 50;
	if(z>=1500 &&z<1600 && sqrt((x-2000)*(x-2000) + y*y)<4000) val = 100;
	
	rho[i3][i2][i1] = val;
      }
    }
  }
  fp = fopen("frho", "wb");
  fwrite(&rho[0][0][0], nx*ny*nz*sizeof(float), 1, fp);
  fclose(fp);

  for(i3=0; i3<nz; i3++){
    z = oz + (i3+0.5)*dz;
    for(i2=0; i2<ny; i2++){
      y = oy + (i2+0.5)*dy;
      for(i1=0; i1<nx; i1++){
	x = ox + (i1+0.5)*dx;
	
	depth = 900. + 100.*sin(2*PI*x/(3*lenx)+PI/3.)*sin(3*PI*y/(2*leny)-PI/3.);
	if(z<depth) val = 0.3;
	else val = 1.5;
	
	rho[i3][i2][i1] = val;
      }
    }
  }
  fp = fopen("frho_init", "wb");
  fwrite(&rho[0][0][0], nx*ny*nz*sizeof(float), 1, fp);
  fclose(fp);
  
  free3float(rho);
  free2float(bathy);
}
