/* file: dccag.c	  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FILE *fentra1,*fentra2,*fsaida1, *fsaida2, *fsaida3, *fsaida4;

/* Function prototypes. */
long input(void);
void error(char error_text[]);



/* Global variables. */
char *pname;	/* this program's name (for use in error messages) */
char c, string1[40], string2[40], string3[40], string4[40], string5[40], 
        string6[40]; /*imput file*/
double *seq;	/* input data buffer; allocated and filled by input() */
int iflag = 1;	/* integrate the input data if non-zero */

int main(int argc, char **argv){
    
unsigned long int sw=0, ir=0, temp, i , j, n, valor1, valor2, count;
double nmin, nmax, Nn, npts, xmean1, xmean2;
double *x1, *x2, *y1, *y2, *sigma ;
double SUMx1, SUMy1, SUMxy1, SUMxx1, SUMx2, SUMy2, SUMxy2, SUMxx2;
double slope1, slope2, y_intercept1,  y_intercept2, y_estimate1, y_estimate2 ;
double p_f2dcca, p_F2dcca, F2dcca, sqf2dcca;
double p_f2dfa1, p_f2dfa2, p_F2dfa1, p_F2dfa2, sqf2dfa1, 
       sqf2dfa2, F2dfa, varsigma, esigma, sigmamedio;
    
    
  puts("Enter the name for the First Time Series to be read:");
  scanf("%s",string1);
  fentra1 = fopen(string1,"r"); /* Arquivo ASCII, para leitura */
  if(!fentra1){
    printf( "Error opening in the first file\n");
    system("pause");
    exit(0);
  }
  fclose(fentra1);
 
  puts("Enter the name for the Second Time Series to be read:");
  puts("(....with the same number of points N):");
  scanf("%s",string2);
  fentra2 = fopen(string2,"r"); /* Arquivo ASCII, para leitura */
  if(!fentra2){
    printf( "Error opening in the second file\n");
    system("pause");
    exit(0);
  }
  fclose(fentra2);

/* Allocate and fill the input data array seq[]. */
  
  npts = input(); /*Acha o número de pontos da Serie 1*/
  
  x1 = (double *) malloc ((npts+1)*sizeof(double));
  x2 = (double *) malloc ((npts+1)*sizeof(double));
  y1 = (double *) malloc ((npts+1)*sizeof(double));
  y2 = (double *) malloc ((npts+1)*sizeof(double));
  sigma = (double *) malloc ((npts+1)*sizeof(double));

/* Le os dois arquivos, para calcular a soma dos pontos*/
  
  fentra1 = fopen(string1,"r");
  fentra2 = fopen(string2,"r");

  for (i=1; i<= npts; i++) {
    fscanf (fentra1, "%lf", &x1[i]);
    fscanf (fentra2, "%lf", &x2[i]);
  }

  fclose(fentra1);
  fclose(fentra2);

  xmean1=xmean2=0;
  
  for (i=1; i<= npts; i++) {
  	xmean1+=x1[i];
    xmean2+=x2[i];
  }
  xmean1=xmean1/npts; /* Calcula a media global */
  xmean2=xmean2/npts;

  y1[0]=y2[0]=0;

  for (i=1; i<=npts; i++) {
    y1[i]=y1[i-1]+(x1[i]-xmean1); /* Calcula a soma integrada */
    y2[i]=y2[i-1]+(x2[i]-xmean2);
  }

  nmin=4; /* tamanho minimo de box */
  nmax=npts/4; /* tamanho maximo de box */
  n=nmin;

/* Abre o arquivo de saida para o coeficiente rho_DCCA */  
  puts("Put the name for Detrended Cross-Correlation coefficient?");
  scanf("%s",string6);
  fsaida4= fopen(string6,"w"); /* Arquivo ASCII, para escrita */


/* Variacao de do box ate seu valor maximo */
  while (n<=nmax){

    count=0;
    valor1=1;
    valor2=n+1;

    p_F2dcca=0; 
    p_F2dfa1=0;
    p_F2dfa2=0;

    /* em cada janela para um valor de n*/
    while (++count < (npts-n)){


      /* Interpolação de polinomio*/
      SUMx1 = SUMy1 = SUMxy1 = SUMxx1 = 0;
      SUMx2 = SUMy2 = SUMxy2 = SUMxx2 = 0;

      for (i=valor1; i<=valor2; i++) {
        SUMx1 +=  (double)(i);
        SUMy1 += y1[i];
        SUMxy1 += (double)(i)*y1[i];
        SUMxx1 += (double)(i)*(double)(i);

        SUMx2 += (double)(i);
        SUMy2 += y2[i];
        SUMxy2 += (double)(i)*y2[i];
        SUMxx2 += (double)(i)*(double)(i);
      }

/* coeficiente angular e linear dentro do box (i,n) s�rie 1*/
      slope1 = ( SUMx1*SUMy1 - (double)(n+1)*SUMxy1 ) / ( SUMx1*SUMx1 - (double)(n+1)*SUMxx1 );
      y_intercept1 = ( SUMy1 - slope1*SUMx1 ) / (double)(n+1);

/* coeficiente angular e linear dentro do box (i,n) s�rie 2*/
      slope2 = ( SUMx2*SUMy2 - (double)(n+1)*SUMxy2 ) / ( SUMx2*SUMx2 - (double)(n+1)*SUMxx2 );
      y_intercept2 = ( SUMy2 - slope2*SUMx2 ) / (double)(n+1);

      p_f2dcca=0;
      p_f2dfa1=0;
      p_f2dfa2=0;
      
/* Calculo da ordenada do ajuste */  
      for (i=valor1; i<=valor2; i++) {
        y_estimate1 = slope1*(double)(i) + y_intercept1;
        y_estimate2 = slope2*(double)(i) + y_intercept2;

        p_f2dcca+=(y1[i]-y_estimate1)*(y2[i]-y_estimate2);
        
        p_f2dfa1+=(y1[i]-y_estimate1)*(y1[i]-y_estimate1);
        p_f2dfa2+=(y2[i]-y_estimate2)*(y2[i]-y_estimate2);
      }
/* covariancia entre as s�ries 1 e 2 */
      
      p_F2dcca+=p_f2dcca/(double)(n);
      
      p_F2dfa1+=p_f2dfa1/(double)(n);
      p_F2dfa2+=p_f2dfa2/(double)(n);
      
      sigma[count]=p_F2dcca/(sqrt(p_F2dfa1)*sqrt(p_F2dfa2));

      
      valor1+=1; /* janelas sobrepostas por uma unidade*/
      valor2+=1;
    }

    Nn=count;
/* Detrended covariance para o tamanho n
   F2DCCA(n)= 1/(N-n)*sum(i=1;i=N-n)(f2DCCA)*/
    F2dcca=p_F2dcca/Nn;
    
    sqf2dfa1=sqrt(p_F2dfa1/Nn);
    sqf2dfa2=sqrt(p_F2dfa2/Nn);
    F2dfa=sqf2dfa1*sqf2dfa2;
     
    sigmamedio= F2dcca/F2dfa;

    varsigma=0;
    
     for (j=1; j<Nn; j++) {
       varsigma+=(sigma[j]-sigmamedio)*(sigma[j]-sigmamedio);
    }
 
    esigma=sqrt(varsigma/Nn);
 
    fprintf(fsaida4,"%lu %lf\n", n, sigmamedio); 
/*   fprintf(fsaida4,"%u %lf %lf\n", n, sigmamedio, esigma);  */

    ir++; 
    n = (nmin+ir)*pow(pow(2,1.0/8.0), ir); 
  
  } /* loop do n */

  fclose(fsaida4); /*fecha arquivo de saida */

  exit(0);

} /* end main */


long input(){
   
    long maxdat = 0L, npts = 0L;
    double y,yp = 0.0;

    fentra1= fopen(string1,"r"); /* Arquivo ASCII, para leitura */
    
     while(fscanf(fentra1,"%lf",&y)==1){
         if (++npts >= maxdat) {
	       double *s;
           maxdat += 50000;	/* allow the input buffer to grow (the
				   increment is arbitrary) */
	       if ((s = realloc(seq, maxdat * sizeof(double))) == NULL) {
		     fprintf(stderr,
		      "%s: insufficient memory, truncating input at row %ld\n",
		      pname, npts);
	         break;
          }
	      seq = s;
	    }
	    seq[npts] = iflag ? (yp += y) : y;
     }
    
     fclose(fentra1);
 
     if (npts < 1) error("no data read");
    return (npts);
}




void error(char error_text[])
{
    fprintf(stderr, "%s: %s\n", pname, error_text);
    exit(1);
}
