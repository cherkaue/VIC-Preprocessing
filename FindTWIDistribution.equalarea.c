/******************************************************************************
   SUMMARY:  
   This program is designed to create a lake parameter file for use with the VIC
   version 4.1.1 and higher.  It reads in both an ascii dem of a VIC grid cell 
   (in arc/info ascii format with standard 6 line header) and a classified land 
   use map (in arc/info ascii format with standard 6 line header) to calculate 
   either the elevation distribution versus area (VIC 4.1.1 lake parameter 
   format) or the topographic wetness index distribution versus area 
   (experimental VIC 4.1.2_SEA format) for wetland areas only.

   The classification code in the land cover file for wetland and open water 
   areas must be specified here (the codes associated with all other 
   vegetation types do not matter): */

#define WETLANDCLASS 2  /* modified 12 to 141 to fit ACRE*/
#define WATERCLASS 1    /* no water class in ACRE  */
   
/*****
 This code assumes that the raster file has been projected into an equal 
  area projection before processing.
*****/

/******************************************************************************
   NOTES: DEM and vegetation file should correspond to one VIC model grid cell.  
   They must hve the exact same domain and resolution as specified in the header.

   The program uses the method presented by Pelletier (2008) to fill sinks and 
   eliminate errors in flat areas by calculating flow accumulation using 
   a multiple flow direction algorithm. 

   REFERENCES: Jon Pelletier (2008) Quantitative Modeling of Earth Surface Processes. 
  
   USAGE: CreatLakeParam <DEM file> <Grid no> <vegetation file> <SEA flag> ;
     DEM file: Name of DEM (elevation) floating point grid with arcinfo header
     Gridno: integer - the number of the VIC grid cell for the parameter file
     veg file: Name of land cover integer grid with arcinfo header
     SEA flag: "SEA" for output in SEA code file format; "LAKE" for original lake model format

   AUTHOR:       Chun-Mei Chiu / Laura Bowling
   DESCRIPTION:                  
   Usage: 
   Compile with: gcc CreateLakeParam.c -lm -o CreateLakeParam
                 
   COMMENTS:
   Modified: 4/22/2011


*******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<malloc.h>

/* This part is used to check time spent in various functions. */ 
#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>
struct tms tt,uu;

#define MAXSTRING 500
#define NNEIGHBORS  8
#define VERTRES 2.3       /* assumed vertical resolution of the dem (m) */
#define OUTSIDEBASIN -99

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 10000000

#define FREE_ARG char*
#define NR_END 1
#define fillincrement 0.01
#define oneoversqrt2 0.707106781187
double **topo,**flow,**flow1,**flow2,**flow3,**flow4,**flow5,**flow6,**flow7,**flow8;
int *iup,*idown,*jup,*jdown;

typedef struct 
{
  double Rank;
  double x;
  double y;
}ITEM;


/*--- Function Declaration---*/ 
void quick(ITEM *item, int count);
void qs(ITEM *item, int left, int right);
void Topindex(double **dem, double **sink, double **flowacc,int columns, int rows, double xorig, double yorig, double delta, double nodata, char gridno[], char option[], FILE *fo);
void VICcalculation(double** iniarray, int n, int m, char gridno[], float wetlandVeg, float waterVeg, int totalVeg, double delta, char option[]);
double correlation(double *AREASUM, double *DEMSUM, int counter7);
void PrintResult(FILE *file, int columns, int rows, double xorig, double yorig, double delta, double nodata);
double **Memoryalloc(int columns, int rows);

/* for contributing area */
void fillin(double **dem, int columns, int rows, double delta, double nodata);
int *ivector(long nl, long nh);
double *vector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_vector(double *v, long nl, long nh);
double **matrix(int nrl, int nrh, int ncl,int nch);
void indexx(int n,double arr[], int indx[]);
void setupgridneighbors(int lattice_size_x,int lattice_size_y);
void fillinpitsandflats(int i,int j, int lattice_size_x,int lattice_size_y, double nodata);
void mfdflowroute(int i,int j, double nodata);


int main(int argc ,char *argv[])  
{
  FILE   *fdem, *fo, *fv;
  char tempstr[MAXSTRING];         
  char   demfile[1000], option[1000], gridno[1000], outfile[1000];
  int    columns, rows, lattice_size_x, lattice_size_y;   
  int    i, j, cnt;  /* counters */
  double xorig, yorig, delta, nodata;  
  double **dem, **flowacc, **sink; 
  int **veg;
  double min_wetland_elev, max_wetland_elev, min_elev;

  /* test time in each step */
  long time_begin,time_end;
  float elapsed_time; 
  long time_begin2,time_end2;
  float elapsed_time2; 


  /*-------------print the usage ------------------*/  
  if (argc < 3 || argc > 4 )
    {
      printf("Usage: CreatLakeParam <DEM file> <output file> [<min elevation>]\n");
      printf("\t\t DEM file : DEM(elevation) floating point grid with arcinfo header;\n");
      printf("\t\t output file : XYZ style file with coordinates and elevation, TWI and sink depth;\n");
      printf("\t\t min elevation : Minimum elevation to process, helps remove empty pixels created by projection (default = 0.1);\n");
      exit(0);
    }
  strcpy(demfile, argv[1]);
  strcpy(outfile, argv[2]);
  if ( argc == 4 ) min_elev = atof( argv[3] );
  else min_elev = 0.1;
  
  /*-----------------------------------------------*/
  /*	 OPEN FILES*/
  /*-----------------------------------------------*/  
  if((fdem=fopen(demfile,"r"))==NULL)
    { 
      fprintf(stderr, "cannot open/read dem file,%s\n",demfile);
      exit(1);
    }
  else {
    fprintf(stderr, "Opening %s\n",demfile);
  }

  if((fo=fopen(outfile,"w"))==NULL)
    { 
      fprintf(stderr, "cannot open/read dem file,%s\n",demfile);
      exit(1);
    }
  else {
    fprintf(stderr, "Opening %s\n",demfile);
  }

  /* check data file has data inside */ 
  if ( getc(fdem) == EOF) {
    fprintf(stderr, "DEM is empty\n");
    exit(0);
  }
  
 

  /*----------------------------------------------*/
  /*Scan and read in DEM's and VEG's header*/
  /*----------------------------------------------*/  
  fscanf(fdem,"%s %d",tempstr,&columns);
  fscanf(fdem,"%s %d",tempstr,&rows);
  fscanf(fdem,"%s %lf",tempstr,&xorig);
  fscanf(fdem,"%s %lf",tempstr,&yorig);
  fscanf(fdem,"%s %lf",tempstr,&delta);
  fscanf(fdem,"%s %lf",tempstr,&nodata);

  fprintf( stdout, "Header: %i %i %lf %lf %lf %lf\n", columns, rows, xorig, yorig, delta, nodata );

 
  lattice_size_x = columns;
  lattice_size_y = rows;
  fprintf(stderr, "DEM is %d by %d\n", rows, columns);

  /*----------------------------------------------*/
  /* Allocate memory to arrays for handling huge data */
  /*----------------------------------------------*/  
 
  dem = Memoryalloc(columns, rows);
  sink = Memoryalloc(columns, rows);
  flowacc = Memoryalloc(columns, rows);
  fprintf(stderr, "Memory allocated.\n");
  fflush(stderr);

  /*-----------------------------------------------*/
  /* READ IN DEM's Mask FILES                      */
  /*---------------------------------------------- */
  for(i=0; i<rows;i++) 
    {
      if(i==1000 || i==5000 || i==10000 || i==15000 || i==20000)
	fprintf(stderr, "i=%d\n",i);
      for(j=0; j<columns; j++)
	{
	  fscanf(fdem,"%lf",&dem[i][j]);
	  if(dem[i][j] < min_elev)
	    {
	      dem[i][j] = nodata;  //check the dem file
	    }
	} 
    }
  
  fprintf(stderr, "DEM read.\n");
  
  /***********************************/
  /*  fill and calculate multi flow accumulation from dem.    */
  /*  Creates filled dem (topo) and accumulation grid (flow). */
  /***********************************/
  
  fillin(dem, columns, rows, delta, nodata);
  fprintf(stderr, "DEM filled\n");
 
  /*************************************/
  /* wetness index calculation         */
  /*************************************/

  /* Replace dem with filled dem. Plletier code indexes arrays starting at 1, so offset is needed. */
  for (i = 0; i < rows; i++) {
    for (j = 0; j < columns; j++){      
  
      sink[i][j] = topo[j+1][i+1] - dem[i][j];
      dem[i][j] = topo[j+1][i+1];
      flowacc[i][j] = flow[j+1][i+1];
      //  fprintf(stderr,"%.2lf ",flowacc[i][j]);
      //      fprintf(stderr,"%.0lf ",dem[i][j]);
      if (dem[i][j] == nodata )
	cnt++;
    }
    //         fprintf(stderr,"\n");
  }
  
  /* Check to make sure dem contains some data. */
  if(cnt < rows*columns)
    { 
      /* This will generate the 526x526 grid lake paramater */
      //      time_begin2 = times(&tt);
      Topindex(dem, sink, flowacc, columns, rows, xorig, yorig, delta, nodata, gridno, option, fo);
      //  time_end2 = times(&uu);
      //  elapsed_time2= (float)(time_end2-time_begin2)/HZ ;
     }
  else {
    printf("No valid value in this grid %s\n", gridno);
  }
  cnt = 0;


  /*  free memory */
  for (i = 0; i < rows; i++)
    {
      free(dem[i]);
    }
    free(dem);

    return (0);
} /*END OF MAIN FUNCTION*/

/*****************************************************************************/
/*   Topindex Function                                                       */
/*****************************************************************************/
void Topindex(double **dem, double **sink, double **flowacc, int columns, 
	      int rows, double xorig, double yorig, double delta, double nodata, 
	      char gridno[], char option[], FILE *fo)
{ 
  int xneighbor[NNEIGHBORS] = { -1, 0, 1, 1, 1, 0, -1, -1 }; /*8 neighbor*/
  int yneighbor[NNEIGHBORS] = { 1, 1, 1, 0, -1, -1, -1, 0 }; /*8 neighbor*/
  int    i, j, k, x, y, n, lower, count;  /* counters */
  double  dx, dy;  
  double  celev;  /*celev =center elevation, Delev = the difference of elevation */
  double  neighbor_elev[NNEIGHBORS], neighbor_twi[NNEIGHBORS], temp_slope[NNEIGHBORS];
  double  length_diagonal;
  double  **tanbeta, **tanbeta_pixel;
  double  **contour_length, **Delev, **AveDelev;
  double  **wetnessindex, **mask;  
  ITEM    *OrderedCellsDEM;
  ITEM    *OrderedCellsTWI;
  int     Norow, VICrow;

  long time_begin,time_end;
  float elapsed_time; 
  int t95, t90, t85, t80, t75, t70;
  int VICcolumn = 5;  

  double **VIC;

  /*-------------- allocate memory------------*/
  /* mask */
  mask = Memoryalloc(columns, rows);

  /* wetnessindex */
  wetnessindex = Memoryalloc(columns, rows);

  /* tanbeta */
  tanbeta = Memoryalloc(columns, rows);

  /* tanbeta_pixel */
  tanbeta_pixel = Memoryalloc(columns, rows);

  /*contour_length */
  contour_length  = Memoryalloc(columns, rows);

  /*Delev */
  Delev = Memoryalloc(columns, rows);

  /* AveDelev */
  AveDelev = Memoryalloc(columns, rows);

  /* VIC output storage */
  VIC = Memoryalloc(rows*columns, VICcolumn);


  /******  exclude the nodata  *********/
   /* This was already done in main. */
   Norow = 0;
   for(i=0; i<rows;i++) {
     for(j=0; j<columns; j++) {
       if(dem[i][j] == nodata)     {
	 Norow++;
       }}}
   
   VICrow = rows*columns - Norow; 
   fprintf(stderr, "Active cells = %d\n", VICrow);

   /*----------------------------------------------- */
   /*          Start calculation                     */
   /* Set 1d array to put elevation value for ranking*/
   /* This is not really needed since flow accumulation */
   /*   was already calculated, but it doesn't hurt. */
   /*------------------------------------------------*/

   count = (int)(rows*columns - Norow);   /*the size of 1d array */

   if(count==0)    /* There is no data in this cell*/
     { 
       fprintf(stderr, "ERROR: No data in current cell: %s. \n", gridno );
       exit(0);
     }

   /** Allocate memory **/
   if(!(OrderedCellsDEM=(ITEM*) calloc(count,sizeof(ITEM)))) 
     { 
       printf("Cannot allocate memory to first record: OrderedCellsDEM\n");
       exit(1); 
     } 
   if(!(OrderedCellsTWI=(ITEM*) calloc(count,sizeof(ITEM)))) 
     { 
       printf("Cannot allocate memory to first record: OrderedCellsTWI\n");
       exit(1); 
     } 

   /* Go through each row and column,and assign the elevation = dem[row][column] */
   count =0;
   for(i=0; i<rows;i++) {
     for(j=0; j<columns; j++)
       {
	 if(dem[i][j] != nodata) {
	   OrderedCellsDEM[count].Rank = dem[i][j];
	   OrderedCellsDEM[count].y = i;
	   OrderedCellsDEM[count].x = j;
	   count++;  
	 }
       }}

   /* Sort OrderedCellsfine/dems into ascending order (from low to high) */
   quick(OrderedCellsDEM, count);
   

   fprintf(stderr, "count=%d\n",count);
   /* -----------------------------------------------*/
   /* Loop through all cells in descending order(from high to low) of elevation */
   /* -----------------------------------------------*/ 
   for (k = count-1; k >-1; k--) 
     { 
       y = (int)OrderedCellsDEM[k].y;/* y is rows*/
       x = (int)OrderedCellsDEM[k].x;/* x is columns*/	    
       dx = delta;  /* meter to km */
       dy = delta;  /* meter to km */   
       length_diagonal = sqrt((pow(dx, 2)) + (pow(dy, 2)));
      
       /* fill neighbor array*/
       for (n = 0; n < NNEIGHBORS; n++) 
	{
	  int xn = x + xneighbor[n]; /* calculate the x-axis of the neighbor cell */ 
	  int yn = y + yneighbor[n]; /* calculate the y-axis of the neighbor cell */
	  
	  /*Initialize neighbor_elev */
	  neighbor_elev[n] = (double) OUTSIDEBASIN;
	  
	  /*check to see if xn and yn are with in dem boundries*/
	  if(xn>=0 && yn>=0 && xn<columns && yn<rows)
	    {
	      neighbor_elev[n] = ((dem[yn][xn]!=nodata) ?   dem[yn][xn] :(double) OUTSIDEBASIN);
	    }
	}
      
       celev = dem[y][x]; /* the elevation of the center cell */
       lower = 0;         /* determine landscape position     */
       for (n = 0; n < NNEIGHBORS; n++)
	 {
	   if(neighbor_elev[n] == OUTSIDEBASIN) 
	     { /* If the neighbor_elev[n] still equal OUTSIDEBASIN,that means 
		  that it doesn't initialize  */ 
	       neighbor_elev[n] = nodata; 
	     }
	   
	   /* Calculating tanbeta as tanbeta * length of cell boundary between 
	      the cell of interest and downsloping neighbor. */
	   if(neighbor_elev[n] < celev && neighbor_elev[n] != nodata )
	     {
	       if(n==0 || n==2 || n==4 || n==6)
		 {
		   temp_slope[n] = (celev - neighbor_elev[n])/length_diagonal; /*slope               */
		   contour_length[y][x] += 0.2*dx+0.2*dy;                      /*contour length      */
		   tanbeta[y][x] += temp_slope[n]*(0.2*dx+0.2*dy);             /*tan beta            */
		 }
	       else if(n==1||n==5)
		 {
		   temp_slope[n] = (celev - neighbor_elev[n])/dy; /*slope               */
		   contour_length[y][x] += 0.6*dx ;
		   tanbeta[y][x] += temp_slope[n]*0.6*dx;
		 }
	       else if (n==3||n==7)
		 {
		   temp_slope[n] = (celev - neighbor_elev[n])/dx; /*slope               */
		   contour_length[y][x] += 0.6*dy;
		   tanbeta[y][x] += temp_slope[n]*0.6*dy;
		 }
	       /* Count how many neighbors are lower than current pixel. */
	       lower++;
	     }
	 }/* end for (n = 0; n < NNEIGHBORS; n++)*/
     
       /* if this is a flat area(slope doesn't change), then tanbeta = sum of 
	  (0.5 * vertical delta of elevation data)/ horizontal distance between 
	  centers of neighboring grid cells  --------------------------------------*/
       if (lower == 0) 
	 { 
	   tanbeta[y][x] = (4.*((0.5 * VERTRES)/length_diagonal) + 
			    (2.0*((0.5 * VERTRES)/dx)) + (2.0*((0.5 * VERTRES)/dy)))/NNEIGHBORS ;
	   tanbeta_pixel[y][x]= (4.*((0.5 * VERTRES)/length_diagonal) + 
			    (2.0*((0.5 * VERTRES)/dx)) + (2.0*((0.5 * VERTRES)/dy)))/NNEIGHBORS ;
	   contour_length[y][x] = 2.*dx + 2.*dy;
	}            
       else {
	 /* Calculate weighted average tanbeta at end of loop through pixel neighbors. */
	 tanbeta_pixel[y][x] =  tanbeta[y][x]/contour_length[y][x];
	 contour_length[y][x] /= (double)lower;
       }
     
       	 /* Add in an extra safety check. */ 
       if(tanbeta_pixel[y][x] <  (4.*((0.5 * VERTRES)/length_diagonal) + 
			    (2.0*((0.5 * VERTRES)/dx)) + 
				  (2.0*((0.5 * VERTRES)/dy)))/NNEIGHBORS ) {

       tanbeta_pixel[y][x] = (4.*((0.5 * VERTRES)/length_diagonal) + (2.0*((0.5 * VERTRES)/dx)) + (2.0*((0.5 * VERTRES)/dy)))/NNEIGHBORS; 

       }

	 /* Calculate general topographic index TI = A/(C.L.* tan B) */
       wetnessindex[y][x] = (double)(flowacc[y][x])/(contour_length[y][x]*tanbeta_pixel[y][x]);

     } /* end  for (k = 0; k < count-1; k++) { */
  

     
  /* ----------------------------------------------- */
  /* Rank the wetness index order. */
  /* ----------------------------------------------- */
  
   count = (int)(rows*columns - Norow);   /*the size of 1d array */
   for(i=0; i<count;i++) 
      {
	      OrderedCellsTWI[count].Rank = 0.0;
	      OrderedCellsTWI[count].y = 0;
	      OrderedCellsTWI[count].x = 0;
      }

    count =0;
    for(i=0; i<rows;i++) 
      {
	for(j=0; j<columns; j++)
	  {
	      OrderedCellsTWI[count].Rank = wetnessindex[i][j];
	      OrderedCellsTWI[count].y = i;
	      OrderedCellsTWI[count].x = j;

	      OrderedCellsDEM[count].Rank = dem[i][j];
	      OrderedCellsDEM[count].y = i;
	      OrderedCellsDEM[count].x = j;
	      count++;
	  }
      }
 
    /* Sort OrderedCellsfine/wetnessindex into ascending order 
       (from low to high) */
    quick(OrderedCellsTWI, count);
    quick(OrderedCellsDEM, count);

    for (k =0; k<count; k++)
      { 

	/***** What is going on here?  It does seem to be trying to output a sorted list, but x and y are rows and columns, not coordinates.  So what is actually in the output file?  VIC[1][k] is not the same values as wetnessindex[y][k]!!!!!! *****/
	/* Assign values to VIC array in descending TWI order (from high to low.) */
	y = OrderedCellsTWI[count-1-k].y;
	x = OrderedCellsTWI[count-1-k].x;
 
	VIC[0][k]= flowacc[y][x];
	VIC[1][k]= wetnessindex[y][x];
	VIC[2][k]= tanbeta_pixel[y][x];
	VIC[3][k]= AveDelev[y][x];

	y = OrderedCellsDEM[k].y;
	x = OrderedCellsDEM[k].x;
	VIC[4][k]= dem[y][x];

	//fprintf(fo, "%lf %lf %lf\n", xorig+x*delta, yorig+y*delta, VIC[1][k]);
      }

    for ( y = 0; y < rows; y++ ) {
      for ( x = 0; x < columns; x++ ) {
	  fprintf(fo, "%lf %lf %lf %lf %lf\n", xorig+x*delta, yorig+y*delta, dem[y][x], wetnessindex[y][x], sink[y][x] );
      }
    }

    t95 = (int) 0.05*count;
    t90 = (int) 0.1*count;
    t85 = (int) 0.15*count;
    t80 = (int) 0.2*count;
    t75 = (int) 0.25*count;
    t70 = (int) 0.3*count;
    fprintf(stdout, "Thresholds: %lf, %lf, %lf, %lf, %lf, %lf\n", VIC[1][t95],VIC[1][t90],VIC[1][t85],VIC[1][t80],VIC[1][t75],VIC[1][t70]);

  //  time_begin = times(&tt);

  //  time_end = times(&uu);
  // elapsed_time= (float)(time_end-time_begin)/HZ ;
  

  /*This is used to free the memory that have been allocated*/ 
  for(i=0; i< VICcolumn-1; i++)
    {
      free(VIC[i]);     
   }
  free(VIC); 
  for(i=0; i<rows; i++)  
    { 
      free(mask[i]);
      free(wetnessindex[i]);
      free(tanbeta[i]);
      free(tanbeta_pixel[i]);
      free(contour_length[i]);
      free(Delev[i]);
      free(AveDelev[i]);
    }  
  free(mask);
  free(wetnessindex); 
  free(tanbeta);
  free(tanbeta_pixel);
  free(contour_length);
  free(Delev);
  free(AveDelev);
  //fprintf(stdout, " here here here2 elapsed_time =%f  %d %d %d %f\n", elapsed_time, i, VICcolumn, count, VIC[3][100]); 	
  //  free(OrderedCellsfine);


  return;
}/* END wetness FUNCTION*/


/* ----------------------  
  Allocate the memory 
 ------------------------*/
double **Memoryalloc(int columns, int rows)
{
  int i;
  double **arr2;
  /*-------------- allocate memory------------*/
  if(!(arr2 = (double**) calloc(rows,sizeof(double*))))
    { printf("Cannot allocate memory to first record: arr2\n");
      exit(8); 
    }
  for(i=0; i<rows;i++)
      if(!(arr2[i] = (double*) calloc(columns,sizeof(double)))) 
	{ printf("Cannot allocate memory to first record: arr2\n");
	  exit(8); 
	}

  return arr2;
}


/*-----------------------------------------------------------------
  Quick Sort Function: this subroutine starts the quick sort
 ------------------------------------------------------------------*/
void quick(ITEM *item, int count)
{
  qs(item,0,count-1);
  return;
}

/*----------------------------------------------------------------
 this is the quick sort subroutine - it returns the values in
 an array from high to low.
 -----------------------------------------------------------------*/
void qs(ITEM *item, int left,  int right)  
{
  register int i,j;
  ITEM x,y;
  
  i=left;
  j=right;
  x=item[(left+right)/2];
  
  do {
    while(item[i].Rank < x.Rank && i<right) i++;
    while(x.Rank < item[j].Rank && j>left) j--;
    
    if (i<=j) 
      {
	y=item[i];
	item[i]=item[j];
	item[j]=y;
	i++;
	j--;
      }
  } while (i<=j);
  
  if(left<j) qs(item,left,j);
  if(i<right) qs(item,i,right); 

  return;
}







/* ------------------------------------------------------------------------
 * Function: least-squares.c and correlation
 * This program computes a linear model for a set of given data.
 *
 * PROBLEM DESCRIPTION:
 *  The method of least squares is a standard technique used to find
 *  the equation of a straight line from a set of data. Equation for a
 *  straight line is given by 
 *	 y = mx + b
 *  where m is the slope of the line and b is the y-intercept.
 *
 *  Given a set of n points {(x1,y1), x2,y2),...,xn,yn)}, let
 *      SUMx = x1 + x2 + ... + xn
 *      SUMy = y1 + y2 + ... + yn
 *      SUMxy = x1*y1 + x2*y2 + ... + xn*yn
 *      SUMxx = x1*x1 + x2*x2 + ... + xn*xn
 *
 *  The slope and y-intercept for the least-squares line can be 
 *  calculated using the following equations:
 *        slope (m) = ( n*SUMxy -SUMx*SUMy ) / ( n*SUMxx - SUMx*SUMx ) 
 *  y-intercept (b) = ( SUMy - slope*SUMx ) / n
 *  R  = ( n*SUMxy - SUMx*SUMy ) / sqrt(( n*SUMxx - SUMx*SUMx)*(n*SUMyy - SUMy*SUMy));
 *
 * AUTHOR: Dora Abdullah (Fortran version, 11/96)
 * REVISED: RYL (converted to C, 12/11/96)
 * ADAPTED: CMC (add correlation calculation)
 * ---------------------------------------------------------------------- */
double correlation (double *x, double *y, int n)
{
  double SUMx, SUMy, SUMxy, SUMxx, SUMyy, SUMres, res, slope, y_intercept, y_estimate, cor;
  int i;
   
  SUMx = 0; 
  SUMy = 0; 
  SUMxy = 0; 
  SUMxx = 0;
  SUMyy = 0;

  for (i=0; i<n; i++) 
    {
      SUMx = SUMx + x[i];
      SUMy = SUMy + y[i];
      SUMxy = SUMxy + x[i]*y[i];
      SUMxx = SUMxx + x[i]*x[i];
      SUMyy = SUMyy + y[i]*y[i];
    }
  slope = ( n*SUMxy -  SUMx*SUMy ) / ( n*SUMxx - SUMx*SUMx );
  y_intercept = ( SUMy - slope*SUMx ) / n;
  cor =   ( n*SUMxy - SUMx*SUMy ) / sqrt(( n*SUMxx - SUMx*SUMx)*(n*SUMyy - SUMy*SUMy));
  //printf ("\n The linear equation that best fits the given data:\n");
  //printf (" y = %6.2lfx + %6.2lf\n", slope, y_intercept);
  //printf ("--------------------------------------------------\n");
  //printf ("   Original (x,y)     Estimated y     Residual\n");
      
  SUMres = 0;
  for (i=0; i<n ; i++) 
    {
      y_estimate = slope*x[i] + y_intercept;
      res = y[i] - y_estimate;
      SUMres = SUMres + res*res;
      //printf ("(%6.2lf %6.2lf) %6.2lf %6.2lf\n", x[i], y[i], y_estimate, res);
    }
  //printf("--------------------------------------------------\n");
  //printf("Residual sum = %6.2lf  correlation = %6.2lf \n", SUMres, cor);

  return (cor);
}









/***************************************************************************/
/*                     Fill increment                                     */
/* Creates the global filled dem matrix (topo) and flow accumulation (flow) */
/* size_x = columns, size_y = rows [i][j] rows:columns
/**************************************************************************/
void fillin(double **dem, int lattice_size_x, int lattice_size_y, double delta, double nodata)
{
  int i,j,t,*topovecind;
  double *topovec;

  setupgridneighbors(lattice_size_x, lattice_size_y); /* The neighbor setting */

  topo=matrix(1,lattice_size_x,1,lattice_size_y);
  topovec=vector(1,lattice_size_x*lattice_size_y);
  topovecind=ivector(1,lattice_size_x*lattice_size_y);
  flow=matrix(1,lattice_size_x,1,lattice_size_y);
  flow1=matrix(1,lattice_size_x,1,lattice_size_y);
  flow2=matrix(1,lattice_size_x,1,lattice_size_y);
  flow3=matrix(1,lattice_size_x,1,lattice_size_y);
  flow4=matrix(1,lattice_size_x,1,lattice_size_y);
  flow5=matrix(1,lattice_size_x,1,lattice_size_y);
  flow6=matrix(1,lattice_size_x,1,lattice_size_y);
  flow7=matrix(1,lattice_size_x,1,lattice_size_y);
  flow8=matrix(1,lattice_size_x,1,lattice_size_y);

  /* Indexing conventions are different than we typically use. */
  /* topo begins indexing at 1 not zero and topo[col][row] vs. dem[row][col] */

  for (j=1;j<=lattice_size_y;j++) {
    for (i=1;i<=lattice_size_x;i++)
      {
	topo[i][j] = dem[j-1][i-1];
        flow[i][j]= delta*delta;
      } }

  for (j=1;j<=lattice_size_y;j++) {
    for (i=1;i<=lattice_size_x;i++)
      {

	fillinpitsandflats(i,j,lattice_size_x, lattice_size_y, nodata);
      } }

  fprintf(stderr, "Done with fill...\n");
    
  for (j=1; j<=lattice_size_y; j++){
    for (i=1; i<=lattice_size_x; i++){
      topovec[(j-1)*lattice_size_x+i]=topo[i][j];
    }}
  
  indexx(lattice_size_x*lattice_size_y,topovec,topovecind);
  t=lattice_size_x*lattice_size_y+1;

  while (t>1)
    {t--;
      i=(topovecind[t])%lattice_size_x;
      if (i==0) i=lattice_size_x;
      j=(topovecind[t])/lattice_size_x+1;
      if (i==lattice_size_x) j--;
      mfdflowroute(i,j, nodata);
    }

} /* End of fillin() */


int *ivector(long nl,long nh)
{ /* allocate an int vector with subscript range v[nl..nh] */
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

double *vector(long nl, long nh)
{ /* allocate a double vector with subscript range v[nl..nh] */
        double *v;

        v=(double *)malloc((unsigned) ((nh-nl+1+NR_END)*sizeof(double)));
        return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
{ /* free an int vector allocated with ivector() */
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(double *v, long nl, long nh)
{ /* free an int vector allocated with ivector() */
        free((FREE_ARG) (v+nl-NR_END));
}

double **matrix(int nrl,int nrh,int ncl,int nch)
{  /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
  int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

    /*allocate pointers to rows */
    m=(double **) malloc((unsigned) (nrow+1)*sizeof(double*));
    m+=1;
    m -= nrl;
    
    m[nrl]=(double *) malloc((unsigned)((nrow*ncol+1)*sizeof(double)));
    m[nrl] += 1;
    m[nrl] -= ncl;

   /*allocate rows and set pointers to them */
    for(i=nrl+1;i<=nrh;i++) {
      m[i]= m[i]=m[i-1]+ncol;
    }
    /* return pointer to array of pointers to rows */
    return m;
}



void indexx(int n,double arr[], int indx[])
{
        unsigned long i,indxt,ir=n,itemp,j,k,l=1;
        int jstack=0,*istack;
        double a;

        istack=ivector(1,NSTACK);
        for (j=1;j<=n;j++) indx[j]=j;
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                indxt=indx[j];
                                a=arr[indxt];
                                for (i=j-1;i>=1;i--) {
                                        if (arr[indx[i]] <= a) break;
                                        indx[i+1]=indx[i];
                                }
                                indx[i+1]=indxt;
                        }
                        if (jstack == 0) break;
                        ir=istack[jstack--];
                        l=istack[jstack--];
                } else {
                        k=(l+ir) >> 1;
                        SWAP(indx[k],indx[l+1]);
                        if (arr[indx[l+1]] > arr[indx[ir]]) {
                                SWAP(indx[l+1],indx[ir])
                        }
                        if (arr[indx[l]] > arr[indx[ir]]) {
                                SWAP(indx[l],indx[ir])
                        }
                        if (arr[indx[l+1]] > arr[indx[l]]) {
                                SWAP(indx[l+1],indx[l])
                        }
                        i=l+1;
                        j=ir;
                        indxt=indx[l];
                        a=arr[indxt];
                        for (;;) {
                                do i++; while (arr[indx[i]] < a);
                                do j--; while (arr[indx[j]] > a);
                                if (j < i) break;
                                SWAP(indx[i],indx[j])
                        }
                        indx[l]=indx[j];
                        indx[j]=indxt;
                        jstack += 2;
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                        } else {
                                istack[jstack]=j-1;
                                istack[jstack-1]=l;
                                l=i;
                        }
                }
        }
        free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP


void setupgridneighbors(int lattice_size_x, int lattice_size_y)
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     
     for (i=1;i<=lattice_size_x;i++)
      {
	idown[i]=i-1;
	iup[i]=i+1;
      }
     idown[1]=1;

     iup[lattice_size_x]=lattice_size_x;

     for (j=1;j<=lattice_size_y;j++)
      {
	jdown[j]=j-1;
	jup[j]=j+1;
      }
     jdown[1]=1;
     jup[lattice_size_y]=lattice_size_y;
}

void fillinpitsandflats(int i,int j, int lattice_size_x, int lattice_size_y, double nodata)
{    double min;
     
  // I don't think anything should happen if topo cell is equal to nodata, so added 
  // brackets to capture the statements after the nodata check.    KAC
  if (topo[i][j] != nodata) {
    min=topo[i][j]; 
    if (topo[iup[i]][j] < min && topo[iup[i]][j] != nodata ) min=topo[iup[i]][j];
    if (topo[idown[i]][j]<min && topo[idown[i]][j] != nodata) min=topo[idown[i]][j];
    if (topo[i][jup[j]]<min && topo[i][jup[j]]!= nodata) min=topo[i][jup[j]];
    if (topo[i][jdown[j]]<min && topo[i][jdown[j]] != nodata) min=topo[i][jdown[j]];
    if (topo[iup[i]][jup[j]]<min && topo[iup[i]][jup[j]] != nodata) min=topo[iup[i]][jup[j]];
    if (topo[idown[i]][jup[j]]<min && topo[idown[i]][jup[j]] != nodata) min=topo[idown[i]][jup[j]];
    if (topo[idown[i]][jdown[j]]<min && topo[idown[i]][jdown[j]] != nodata) min=topo[idown[i]][jdown[j]];
    if (topo[iup[i]][jdown[j]]<min && topo[iup[i]][jdown[j]] != nodata) min=topo[iup[i]][jdown[j]];
    
    if ((topo[i][j] <= min)&& (topo[i][j]!=nodata)&&(i>1)&&(j>1)&&(i<lattice_size_x)&&(j<lattice_size_y))
      {
	topo[i][j]=min+fillincrement;
	fillinpitsandflats(i,j, lattice_size_x, lattice_size_y, nodata);
	fillinpitsandflats(iup[i],j, lattice_size_x, lattice_size_y, nodata);
	fillinpitsandflats(idown[i],j, lattice_size_x, lattice_size_y, nodata);
	fillinpitsandflats(i,jup[j], lattice_size_x, lattice_size_y, nodata);
	fillinpitsandflats(i,jdown[j], lattice_size_x, lattice_size_y, nodata);
	fillinpitsandflats(iup[i],jup[j], lattice_size_x, lattice_size_y, nodata);
	fillinpitsandflats(idown[i],jup[j], lattice_size_x, lattice_size_y, nodata);
	fillinpitsandflats(idown[i],jdown[j], lattice_size_x, lattice_size_y, nodata);
	fillinpitsandflats(iup[i],jdown[j], lattice_size_x, lattice_size_y, nodata);
      }
  }
}

void mfdflowroute(int i,int j, double nodata)
{ 
  double tot;
 
  if(topo[i][j] == nodata)
    flow1[i][j]=flow2[i][j]=flow3[i][j]=flow4[i][j]=flow5[i][j]=flow6[i][j]=flow7[i][j]=flow8[i][j]=0.0;
  else {
     tot=0.;
     if (topo[i][j]>topo[iup[i]][j] && topo[iup[i]][j]!= nodata) 
      tot+=pow(topo[i][j]-topo[iup[i]][j],1.1);
     if (topo[i][j]>topo[idown[i]][j] && topo[idown[i]][j]!= nodata) 
      tot+=pow(topo[i][j]-topo[idown[i]][j],1.1);
     if (topo[i][j]>topo[i][jup[j]] && topo[i][jup[j]]!= nodata) 
      tot+=pow(topo[i][j]-topo[i][jup[j]],1.1);
     if (topo[i][j]>topo[i][jdown[j]] && topo[i][jdown[j]]!= nodata) 
      tot+=pow(topo[i][j]-topo[i][jdown[j]],1.1);
     if (topo[i][j]>topo[iup[i]][jup[j]] && topo[iup[i]][jup[j]]!= nodata) 
      tot+=pow((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[iup[i]][jdown[j]] && topo[iup[i]][jdown[j]]!= nodata) 
      tot+=pow((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[idown[i]][jup[j]] && topo[idown[i]][jup[j]]!= nodata) 
      tot+=pow((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[idown[i]][jdown[j]] && topo[idown[i]][jdown[j]]!= nodata) 
      tot+=pow((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2,1.1);
    
     if (topo[i][j]>topo[iup[i]][j] && topo[iup[i]][j]!= nodata) 
       flow1[i][j]=pow(topo[i][j]-topo[iup[i]][j],1.1)/tot; 
     else flow1[i][j]=0;

     if (topo[i][j]>topo[idown[i]][j] && topo[idown[i]][j]!= nodata) 
      flow2[i][j]=pow(topo[i][j]-topo[idown[i]][j],1.1)/tot; 
     else flow2[i][j]=0;
     
     if (topo[i][j]>topo[i][jup[j]] && topo[i][jup[j]]!= nodata) 
      flow3[i][j]=pow(topo[i][j]-topo[i][jup[j]],1.1)/tot; 
     else flow3[i][j]=0;
     
     if (topo[i][j]>topo[i][jdown[j]] && topo[i][jdown[j]]!= nodata) 
       flow4[i][j]=pow(topo[i][j]-topo[i][jdown[j]],1.1)/tot; 
     else flow4[i][j]=0;
     
     if (topo[i][j]>topo[iup[i]][jup[j]] && topo[iup[i]][jup[j]]!= nodata) 
       flow5[i][j]=pow((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2,1.1)/tot;
     else flow5[i][j]=0;
     
     if (topo[i][j]>topo[iup[i]][jdown[j]] && topo[iup[i]][jdown[j]]!= nodata) 
       flow6[i][j]=pow((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2,1.1)/tot;
     else flow6[i][j]=0;
     
     if (topo[i][j]>topo[idown[i]][jup[j]] && topo[idown[i]][jup[j]]!= nodata) 
       flow7[i][j]=pow((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2,1.1)/tot;
     else flow7[i][j]=0;
     
     if (topo[i][j]>topo[idown[i]][jdown[j]] && topo[idown[i]][jdown[j]]!= nodata) 
       flow8[i][j]=pow((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2,1.1)/tot;
     else flow8[i][j]=0;
  }
     
     flow[iup[i]][j]+=flow[i][j]*flow1[i][j];
     flow[idown[i]][j]+=flow[i][j]*flow2[i][j];
     flow[i][jup[j]]+=flow[i][j]*flow3[i][j];
     flow[i][jdown[j]]+=flow[i][j]*flow4[i][j];
     flow[iup[i]][jup[j]]+=flow[i][j]*flow5[i][j];
     flow[iup[i]][jdown[j]]+=flow[i][j]*flow6[i][j];
     flow[idown[i]][jup[j]]+=flow[i][j]*flow7[i][j];
     flow[idown[i]][jdown[j]]+=flow[i][j]*flow8[i][j];
}
