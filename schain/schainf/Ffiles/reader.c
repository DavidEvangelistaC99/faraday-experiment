#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include "complex.h"

#define NFREQ 512
#define LNES 12
#define LTES 16
#define NANG 19

#define MAX(x,y) ((x)>(y) ?  (x) : (y))
#define MIN(x,y) ((x)<(y) ?  (x) : (y))

fcomplex ****exlib;
int first=1;

void swab_test(void);
//float swab(signed char *a);
void read_exlib(char *s,int nfreq, int lnes, int lTes, int nang,
                fcomplex ****exlib);

void initialize(void);
//void collision_(float *np, float *tp, float *fp, float *ap, float *yr, float *yi)// step 1
void collision_(float *n, float *t, float *f, float *a, fcomplex *y);

/*
void main(void){
  float n=4.9e11,t=1100.0,f=300.0,a=1.0;
  int i;
  fcomplex j;
  for(i=0;i<10;i++){
    a+=0.05;
    collision_(&n,&t,&f,&a,&j);
    printf("%f %f %f \n",a,j.r,j.i);
  }
}
*/

//void collision_(float *np, float *tp, float *fp, float *ap, float *yr, float *yi){ // step 2
void collision_(float *np, float *tp, float *fp, float *ap, fcomplex *y){
  float ne[]={5.0e10, 1.0e11, 2.0e11, 3.5e11, 5.0e11, 7.5e11,
	      1.0e12, 1.5e12, 2.0e12, 2.5e12, 3.0e12, 3.5e12};
  float te[]={600.0, 700.0, 800.0, 900.0, 1000.0, 1150.0, 1300.0, 1500.0,
	      1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3500.0, 4000.0};
  float alpha[]={.125, .1875, .25, .3125, .375, .5, .625, .75, 1.0, 1.25,
		 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0};
  float freq[2],fmax=6237.79;
  float n,t,f,a;
  float v,w1,w2,w3,w4,weight;
  fcomplex g;
  int i,j,k,l,m,i1,j1,k1,l1;

  if(first==1){
	//printf("initi");
    initialize();
    first=0;
  }
//printf(" inputs in collision_ \n y.r %f    y.i %f   \n *************\n",(*y).r,(*y).i);

//printf("After:print\n");
a=*ap; t=*tp; f=*fp; n=*np;

  a=MAX(MIN(a,6.0),0.125);
  t=MAX(MIN(t,4000.0),600.0);
  f=MIN(f,fmax-1.0);
  n=MAX(MIN(n,3.5e12),5.0e10);

  i1=((float)(NFREQ-1))*f/fmax;
  i1=MIN(i1,NFREQ-1);
  freq[0]=i1*fmax/(float)(NFREQ-1);
  freq[1]=(i1+1)*fmax/(float)(NFREQ-1);

  j1=0;
  while(n>ne[j1+1]) j1++;
  j1=MIN(j1,LNES-1);

  k1=0;
  while(t>te[k1+1]) k1++;
  k1=MIN(k1,LTES);

  l1=0;
  while(a>alpha[l1+1]) l1++;
  l1=MIN(l1,NANG);

  // printf("i,j,k,l: %i %i %i %i \n",i1,j1,k1,l1);

  v=(fmax/(float)(NFREQ-1))*
    (ne[j1+1]-ne[j1])*(te[k1+1]-te[k1])*(alpha[l1+1]-alpha[l1]);
  g=Complex(0.0,0.0);

  for(i=0;i<2;i++){
    if(i==0)
      w1=freq[1]-f;
    else
      w1=f-freq[0];
    for(j=0;j<2;j++){
      if(j==0)
	w2=ne[j1+1]-n;
      else
	w2=n-ne[j1];
      for(k=0;k<2;k++){
	if(k==0)
	  w3=te[k1+1]-t;
	else
	  w3=t-te[k1];
	for(l=0;l<2;l++){
	  if(l==0)
	    w4=alpha[l1+1]-a;
	  else
	    w4=a-alpha[l1];
	  weight=w1*w2*w3*w4;
	  g.r+=exlib[i1+i][j1+j][k1+k][l1+l].r*weight;
	  g.i+=exlib[i1+i][j1+j][k1+k][l1+l].i*weight;
 //        printf(" in collision    i1+i %d j1+j %d k1+k %d l1+l %d  exlib.r %f  exlib.i %f \n",i1+i,j1+j,k1+k,l1+l,exlib[i1+i][j1+j][k1+k][l1+l].r,exlib[i1+i][j1+j][k1+k][l1+l].i);
 //        printf(" in collision    g.r %f  g.i %f \n",g.r,g.i);
	}
      }
    }
  }
  g.r/=v; g.i/=v;

  *y=g;
//*yr=g.r
//*yi=g.i
//printf(" outputs in collision_ \n g.r %f    g.i %f   \n *************\n",g.r,g.i);
//exit(-2);
//getchar();
return;
}

void initialize(void){
//  char s[]="/usr/local/lib/faraday/jlib26feb2002";
//  char s[]="/usr/local/lib/faraday/jlib26feb2002";
//  FILE *s=fopen("/usr/local/lib/faraday/jlib26feb2002","r");
  int nfreq=NFREQ, lnes=LNES, lTes=LTES, nang=NANG;
  int i,j,k;
  int crr,crr2;
  int the_number=1025;

  char the_path[1025];

  get_path_reader_(the_path,&the_number);
  //printf("C the_number: %i",the_number);
  //printf("After\n");
  //printf("The Path: %s", the_path);
  //printf("END\n");
  //printf("%s%n",the_path,&crr);
  //printf("=%d",crr);
  //printf("END2\n");
  char the_path_true[the_number];
  memcpy(the_path_true,the_path,the_number);
  the_path_true[the_number] = '\0';
  //printf("%s%n",the_path_true,&crr2);
  //printf("=%d",crr2);
  //printf("END3\n");
  //printf("PATH: %s",s);
  //fprintf(s);
  exlib=malloc(nfreq*8);
  for(i=0;i<nfreq;i++){
    exlib[i]=malloc(lnes*8);
    for(j=0;j<lnes;j++){
      exlib[i][j]=malloc(lTes*8);
      for(k=0;k<lTes;k++){
	exlib[i][j][k]=(fcomplex *)malloc(nang*sizeof(fcomplex));
      }
    }
  }

  read_exlib(the_path_true,nfreq,lnes,lTes,nang,exlib);
}

float swab(signed char *a)
// table is in big endian format
{
  float mant;
  unsigned char *amant;
  printf(" in swab 0 ***  a[0] %c   a[1] %c   a[2] %c    a[3] %c\n",a[0],a[1],a[2],a[3]);
  amant=(unsigned char *) &mant;
  amant[3]=a[0]; amant[2]=a[1]; amant[1]=a[2]; amant[0]=a[3];
  //printf(" in swab 1 ***\n");
  return(mant);
}

void swab_test(void)
// table is in big endian format
{
  printf("*********in swab stest ********\n");
  getchar();
}


void read_exlib(char *the_path_true,int nfreq, int lnes, int lTes, int nang,
                fcomplex ****exlib)
{
      FILE *fp;
      int  i,j,k,l;
      float a;
      float t1, t2;
      signed char *t1a,*t2b;
	      float mant;
              unsigned char *amant;
      printf("Start reading exlib file \n");
      if( (fp=fopen(the_path_true,"r")) == NULL){
        printf("Je Library file opening error\n");
        //printf("%s\n", strerror(errno));
        exit(-2);
      }
      for(l=0;l<nang;l++)
	for(k=0;k<lTes;k++)
	  for(j=0;j<lnes;j++)
            for(i=0;i<nfreq;i++){
	      //printf(" in read_exlib 0***\n");
	      fread(&t1,sizeof(float),1,fp);
	      fread(&t2,sizeof(float),1,fp);
	      //printf(" in read_exlib 1***\n");
              //printf(" i %d j %d k %d l %d  t1 %f  t2 %f \n",i,j,k,l,t1,t2);
	      t1a=(signed char *)&t1;
              t2b=(signed char *)&t2;
	      //printf(" in read_exlib 1.5  t1a %d  t2b %d ***\n",t1a,t2b);
	      ////t1=swab(t1a);
              //printf(" t1 -- in readexlib 1.8  0 ***  t1a[0] %c   t1a[1] %c   t1a[2] %c    t1a[3] %c\n",t1a[0],t1a[1],t1a[2],t1a[3]);
              amant=(unsigned char *) &mant;
              amant[3]=t1a[0]; amant[2]=t1a[1]; amant[1]=t1a[2]; amant[0]=t1a[3];
              t1=mant;

              ////t2=swab(t2b);
              //printf(" t2 -- in readexlib 1.8  0 ***  t2b[0] %c   t2b[1] %c   t2b[2] %c    t2b[3] %c\n",t2b[0],t2b[1],t2b[2],t2b[3]);
              amant=(unsigned char *) &mant;
              amant[3]=t2b[0]; amant[2]=t2b[1]; amant[1]=t2b[2]; amant[0]=t2b[3];
              t2=mant;

	      ////t1=swab((signed char *)&t1);
              ////t2=swab((signed char *)&t2);
	      //printf(" in read_exlib 2   i %d j %d k %d l %d  t1 %f  t2 %f \n",i,j,k,l,t1,t2);
	      exlib[i][j][k][l]=Complex(t1,t2);
              //printf(" in read_exlib 3   i %d j %d k %d l %d  t1 %f  t2 %f \n",i,j,k,l,exlib[i][j][k][l].r,exlib[i][j][k][l].i);
	      //getchar();
	    }
          // getchar();

      fclose(fp);
      printf("Done reading exlib file \n");
}
/*
void read_exlib(char *s,int nfreq, int lnes, int lTes, int nang,
                fcomplex ****exlib)
{
      FILE *fp;
      int  i,j,k,l;
      float a;
      float t1, t2;

      printf("Start reading exlib fileee \n");

      if( (fp=fopen(s,"r")) == NULL){
        printf("Je Library file opening error\n");
        exit(-2);
      }

      for(l=0;l<nang;l++)
	for(k=0;k<lTes;k++)
	  for(j=0;j<lnes;j++)
            for(i=0;i<nfreq;i++){
	      printf(" in read_exlib 0***\n");
	      fread(&t1,sizeof(float),1,fp);
	      fread(&t2,sizeof(float),1,fp);
	      printf(" in read_exlib 1***\n");
              printf(" i %d j %d k %d l %d  t1 %f  t2 %f \n",i,j,k,l,t1,t2);
	      t1=swab((char *)&t1); t2=swab((char *)&t2);
	      printf(" in read_exlib 2   i %d j %d k %d l %d  t1 %f  t2 %f \n",i,j,k,l,t1,t2);
	      exlib[i][j][k][l]=Complex(t1,t2);
              printf(" in read_exlib 3   i %d j %d k %d l %d  t1 %f  t2 %f \n",i,j,k,l,exlib[i][j][k][l].r,exlib[i][j][k][l].i);
	      getchar();
	    }

      fclose(fp);
      printf("Done reading exlib file \n");
}
*/
