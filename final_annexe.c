/* Projet Visual Sound Group 10 Polytech Nice Sophia
   LE Trung Hieu
   LIBERATO Enzo
   PLISSON Louis
   GUERAND Arthur
   Ait-Ali Youssef
   Dernière Modification:25 Mai 2018
   
   Pour bien utilisé ce programme, il faut utiliser émulateur de Terminal picocom
   sudo apt-get install picocom
   
   Pour initialisé l'environnement
   sudo picocom -b 115200 /dev/ttyACM0
   
   Cross Compilation Command : $CC -std=c99 final_annexe.c -lm -lsndfile -lncurses -o final_color

   Version testé sur l'ordinateur et sur la carte SAMA Xplained en Couleur
   @2018
*/


#include <sndfile.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <ncurses.h>
#define M_PI		3.14159265358979323846
#define SWAP(a,b) ctmp=(a); (a)=(b); (b)=ctmp
#define ARRAY_LEN(x)	((int) (sizeof (x) / sizeof (x [0])))
#define log2(x) log(x)/log(2)
#define MAX(x,y)	((x)>(y) ? (x) : (y))
#define MIN(x,y)	((x)<(y) ? (x) : (y))
#define ABS(x)		((x)>=0 ? (x) : (-x))
#define complex _Complex //Compilateur C99

#define BUFFER_LEN 4096 
//Nous avons choisi la taille de FFT 4096
#define FREQUENCE_MAX 2048
// Nous prenons que des fréquences normalisé.

static int MAXLINE=15;
// Définition de la hauteur
static int MAXCOLONNE=40;
// Définition de la largeur
static int color_check=1;

static int division = 3;


void init_graphic()
{
  // initialisation de l'interface graphic l'environnement curses
  initscr();
  cbreak();
  nonl;
  noecho();
  curs_set(0);
  if((has_colors() == FALSE)&&color_check)
	{
		printf("Your terminal does not support color dommage\n");
    color_check = 0;
	}
  clear();
}

double complex TW[BUFFER_LEN];

void ASK_NAME(char *s) 
{
	printf("Nom du fichier(.wav <100 caractere):\n");
	scanf("%s",s);
}

int ASK_TYPE()
{
	int i = 0, r;
	while (i == 0)
	{
		printf("\n TYPE DE FFT ?\n 1- Iterative\n 2- Recursive\n 3- Quitter le programme\n\n(La reponse attendue est le nombre devant la réponse de votre choix)\n");
		scanf("%1d",&r);
		switch(r) 
			{
				case 1 : i=1; return 0; break;
				case 2 : i=1; return 1; break;
				case 3 : i=1; return 2; break;
				default : i=0; break;
			}
	}
	return EXIT_SUCCESS;
}

sf_count_t sfx_mix_mono_read_double (SNDFILE * file, double * data, sf_count_t datalen)
{
  SF_INFO info;
  static double multi_data [2048] ;
  int k, ch, frames_read ;
  sf_count_t dataout = 0 ;

  sf_command (file, SFC_GET_CURRENT_SF_INFO, &info, sizeof (info)) ;
  if (info.channels == 1)
    return sf_read_double (file, data, datalen) ;

  while (dataout < datalen)
    {	int this_read ;

      this_read = MIN (ARRAY_LEN (multi_data) / info.channels, datalen) ;
      frames_read = sf_readf_double (file, multi_data, this_read) ;
      if (frames_read == 0) break ;

      for (k = 0 ; k < frames_read ; k++)
		{	double mix = 0.0 ;

	  	for (ch = 0 ; ch < info.channels ; ch++)
	    	mix += multi_data [k * info.channels + ch] ;
	  		data [dataout + k] = mix / info.channels ;
		} 
      dataout += this_read ;
    } 
  return dataout ;
}

void fftrec(double complex *data, double complex *result, unsigned int size, int log2n)
{
  double complex ypair [size],yimpair [size],Fimpair [size], Fpair [size] ;
  int n,k,N2;

  if (size > 1)
    {
      N2=size/2;
      for (n=0;n<N2;n++)
	{
	  ypair[n]=data[n]+data[n+N2];
	  yimpair[n]=(data[n]-data[n+N2])*TW[n*BUFFER_LEN/size];
 	}
      fftrec(ypair,Fpair,N2,log2n);
      fftrec(yimpair,Fimpair,N2,log2n);
      for (n=0;n<N2;n++)
	{
	  result[2*n]=Fpair[n];
	  result[2*n+1]=Fimpair[n];
	}
    }
  else 
    {    
     result[0]=data[0];
     return ;
    } 
}

//Pour éviter une erreur de compilation on crée la fonction bitrev
int bitrev(int inp, int numbits){
  int i, rev = 0;
  for(i = 0; i<numbits; i++){
    rev = (rev<<1)|(inp&1);
    inp >>=1;
  }
  return rev;
}

void fftiterTW(double complex *data, unsigned int size, int log2n)
{
  int j,N2, Bpair, Bimpair,Bp=1,N=size;
  double complex impair, pair,ctmp;

 for (int k=0;k<log2n;k++)
   {
    N2=N/2;
    Bpair=0;
    for (int b=0;b<Bp;b++)
      {
        Bimpair=Bpair+N2;
	for (int n=0;n<N2;n++) 
	  { 
	    impair=data[Bpair+n]+data[Bimpair+n];
            pair=(data[Bpair+n]-data[Bimpair+n])*TW[n*size/N];
            data[Bpair+n]=pair;
            data[Bimpair+n]=impair;
	  }
        Bpair=Bpair+N;
      }
    Bp=Bp*2;
    N=N/2;
   }
 for (int i=0;i<size;i++)
    { 
      j=bitrev(i,log2n);
      if (j>i)
	{ 
            SWAP(data[j],data[i]);
        }
     }
 for (int i=size-1;i>0;i--)
   data[i]=data[i-1];
 data[0]=ctmp;
 return ;
}

void twiddle(double complex *TW, unsigned int size)
{
  double complex phi=cexp(-2*I*M_PI/size);

  TW[0]=1;
  for (int i=1;i<size;i++)
      TW[i]=TW[i-1]*phi;
}

void module(double complex *data,double *rslt, int longueur)
{ //On enregistre le module de data dans rslt
  for (int i=0; i<longueur; i++){
    rslt[i]=creal(data[i])*creal(data[i])+cimag(data[i])*cimag(data[i]);//Le module est égal a Re^2+Im^2 
  }
}

void gain(double *data, double *rslt, int longueur)
{//On enregistre le gain en dB de data dans rslt
  for (int i=0; i<longueur; i++){
    if (data[i] > 1) rslt[i]= 10*log(data[i]);
    else rslt[i]=0;
  }
}


void dessiner(double GAIN[],int frequence_ech,int count,float pasdB)
{
  int xi,yi; 
  // coordonnée pour dessiner
  int temps;
  int nb_division; 
  // le nombre de division
  int pos = 0; 
  // la position de la fréquence où on va ramasser pour faire la moyenne
  double Gain_moyenne; 
  // le Gain moyenne
  double data_affiche[MAXCOLONNE]; 
  // le tableau de valeur moyenne
  nb_division = FREQUENCE_MAX/MAXCOLONNE; 
  // le nombre de fréquence nous allons faire la moyenne
  
  // Nous faisons l'algorithme pour regrouper des fréquences.
  for (int i = 0;i<MAXCOLONNE;i++)
    {
      Gain_moyenne = 0;
      for (int j=0;j<nb_division;j++)
	{
	  Gain_moyenne += GAIN[pos];
	  pos = pos +1;

	}
      Gain_moyenne = Gain_moyenne /nb_division;
      data_affiche[i]=Gain_moyenne;
    }
  // Affichage de temps de la piste en lecture
  mvprintw(0,0,"Temps Actuel:%0.1f second(s)",(double)(count*BUFFER_LEN/frequence_ech));
  mvprintw(1,0,"%0.1f(dB) par div horizontale",pasdB);
  // Le temps est calculé en fonction de fichier en lecture, pas à l'horloge
  // Affichage du spectre
  for (xi = 0; xi<MAXCOLONNE;xi++)  
    for (yi = 0;yi<MIN((int) (data_affiche[xi])/3+1,MAXLINE-2);yi++)
      // Si amplitudes dépasse la taille de la fenetre, on bloque par la limit MAXLINE
      {
	mvprintw(MAXLINE-yi,xi,"%c",'-');
      }
  refresh();
  // rafraichir l'écran, afficher toutes les points
}

void set_color(void)
{
  // Définition des couleurs
  int color_back;
  int color_front;
  int choix;
  if(color_check)
    {
      mvprintw(0,0,"Vous parametrer des couleurs ? 1:Oui 0:Par defaut choix:");
      scanw("%d",&choix);
      if(choix)
      {
        mvprintw(1,0,"\nCOLOR_BLACK 0\nCOLOR_RED 1\nCOLOR_GREEN 2\nCOLOR_YELLOW 3\nCOLOR_BLUE 4\nCOLOR_MAGENTA 5\nCOLOR_CYAN 6\nCOLOR_WHITE  7\n");
        mvprintw(10,0,"Couleur Background:");
        scanw("%d",&color_back);
        mvprintw(12,0,"Couleur Frontground:");
        scanw("%d",&color_front);
        if (color_back >= 0&&color_back<8||color_front>=0&&color_front<8)
        {
          start_color();
          init_pair(1,color_front,color_back);
          attron(COLOR_PAIR(1));
        }

      }
    }
}

int main(int argc, char const *argv[])
{
  int ter_lines,ter_cols;
  init_graphic(); 
  ter_cols = COLS-1;
  ter_lines = LINES -2;
  clear();
  endwin();
  clock_t ti,tf;
  long int temps;
	SF_INFO      sfinfo;
	SNDFILE      *infile, *outfile ;
	char s[100]; 
  //Les titres des fichiers ne peuvent pas dépasser 100 caractères
	int type;
	int HAUTEUR = 40;
	int BANDE = 128;
	sf_count_t readcount;
	double maximum = LONG_MAX;
	float pasdB;
	float pasHz;
  int choix;
  printf("Choisi mode de configuration:\n 0.Par défaut \n 1.Manuelle \n");
  scanf("%d",&choix);
  if (choix)
  {
      printf("Votre fenetre maximal dectecte par le système est: %d x %d",ter_lines-2,ter_cols-1);
      printf("\nEntrer la taille de la fenetre d'affichage (supérieur à 0 et de type ENTIER)\n");
      printf("HAUTEUR:");
      scanf("%d",&MAXLINE);
      printf("LARGEUR:");
      scanf("%d",&MAXCOLONNE);
      if (MAXLINE>ter_lines||MAXCOLONNE>ter_cols||MAXLINE<1||MAXCOLONNE<1)
      {
        printf("ERREUR DIMENSION: PARAMETRE PAR DEFAUT ACTIVE !\n");
        MAXLINE=15;
        MAXCOLONNE = 40;
      }
      printf("Pas de division pour %d",MAXLINE);
      scanf("%d",&choix);
      if (choix>1)
      {
        division=choix;
      }
      else
      {
        printf("ERR:Paramère par défaut sélectioné\n");
      }
  }
	ASK_NAME(s);
	type = ASK_TYPE();
	pasdB = MAXLINE/division;
  // Calculer les échelles
	infile = sf_open(s, SFM_READ, &sfinfo);
  	if ( infile == NULL)
    	{   /* Open failed so print an error message. */
      	printf ("Not able to open input file %s.\n", s) ;
      	/* Print the error message fron libsndfile. */
      	sf_perror (NULL) ;
      	return  1 ;
    	} 
	int count = 0;
	int samplerate=sfinfo.samplerate;
 	int channels=sfinfo.channels;
 	double data[2*BUFFER_LEN]; //On y enregistre la data
 	double complex data_compl[2*BUFFER_LEN]; //On y enregistre la data en complexe
	int log2n;
	double complex fft[BUFFER_LEN]; //On y enregistrera la fft
	double fft_mod[BUFFER_LEN]; // On y enregistrera le module de la fft
	double GAIN[BUFFER_LEN]; //On y enregistre le gain en dB
    init_graphic();
    set_color();
  	while ((readcount = sfx_mix_mono_read_double (infile, data, BUFFER_LEN))>0)
    	{
	      ti = clock();
			  log2n = log2(BUFFER_LEN);
			  twiddle(TW,readcount);
			  count ++;

  			for (int i =0;i<2*BUFFER_LEN;i++) data_compl[i] = (double complex)data[i];
        switch(type)
        {
          case 2 : return EXIT_SUCCESS ; break;
          case 1 : fftrec(data_compl,fft,readcount,log2(readcount));break;
          case 0 : 	fftiterTW(data_compl,readcount,log2(readcount));for (int j=0; j<readcount; j++) fft[j] = data_compl[j];break;
          defaut: fprintf(stderr, "Mauvaise Selection");return EXIT_FAILURE ;
        } 

        module(fft,fft_mod,readcount);
        gain(fft_mod,GAIN,readcount);
        dessiner(GAIN,samplerate,count,pasdB);
        refresh();
        /*Calculer le temps de calculer du microprocesseur */
        tf=clock();
        temps = 93091 - (long int)((tf-ti)*1000000/CLOCKS_PER_SEC);
        usleep(temps);
        clear();
    	}
	endwin();
  	/* Close input and output files. */
  sf_close (infile) ;
}
