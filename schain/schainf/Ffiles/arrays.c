float **array2(int m, int n){

  int i,j;
  float **temp;

  temp=(float **)malloc(sizeof(double)*m);

  temp[0]=(float *)malloc(sizeof(double)*m*n);
  for(j=1;j<m;j++){
    temp[j]=temp[j-1]+n;
  }

  return temp;
}


float ***array3(int m, int n, int l){

  int i,j;
  float ***temp;

  temp=(float ***)malloc(sizeof(double)*m);
  temp[0]=(float **)malloc(sizeof(double)*m*n);
  for(j=1;j<m;j++){
    temp[j]=temp[j-1]+n;
  }
  temp[0][0]=(float *)malloc(sizeof(double)*l*n*m);;
  for(j=0;j<m;j++)
    for(i=0;i<n;i++){      
      temp[j][i]=temp[0][0]+ l*(i+n*j);
    }
  return temp;
}

fcomplex **carray2(int m, int n){

  int i,j;
  fcomplex **temp;

  temp=(fcomplex **)malloc(sizeof(fcomplex)*m);

  temp[0]=(fcomplex *)malloc(sizeof(fcomplex)*m*n);
  for(j=1;j<m;j++){
    temp[j]=temp[j-1]+n;
  }

  return temp;
}


fcomplex ***carray3(int m, int n, int l){

  int i,j;
  fcomplex ***temp;

  temp=(fcomplex ***)malloc(sizeof(fcomplex)*m);
  temp[0]=(fcomplex **)malloc(sizeof(fcomplex)*m*n);
  for(j=1;j<m;j++){
    temp[j]=temp[j-1]+n;
  }
  temp[0][0]=(fcomplex *)malloc(sizeof(fcomplex)*l*n*m);;
  for(j=0;j<m;j++)
    for(i=0;i<n;i++){      
      temp[j][i]=temp[0][0]+ l*(i+n*j);
    }
  return temp;
}
