#include"stdio.h"
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include"math.h"
#include <ctype.h>
//#include"openacc.h"
//float *px,*py;

#pragma acc routine vector
int distD(int i,int j,float *x,float*y)
{
	float dx=x[i]-x[j];
	float dy=y[i]-y[j]; 
	return(sqrtf( (dx*dx) + (dy*dy) ));
}

void setCoord(int *r,float *posx,float *posy,float *px,float *py,int cities)
{
	int i;
	for(i=0;i<cities;i++)
	{
		px[i]=posx[r[i]];
		py[i]=posy[r[i]];
	}
}

int *twoOpt(int x,int y,int *route,int city)
{
	int *tmp;
	tmp=(int *)malloc(sizeof(int )*city);
	int i,j;

	for (i = 0; i <=x; ++i)
	{
		tmp[i] = route[i];
	}

	for (i = x+1, j = y; i <= y; ++i, --j)
	{
		tmp[i] = route[j];
	}


	for (i = y+1; i < city; ++i)
	{
		tmp[i] = route[i];
	}


	return tmp;

}


int distH(float *px,float *py,int cit)
{
	float dx,dy,cost=0;
	int i;
	for(i=0;i<(cit-1);i++)
	{
		dx=px[i]-px[i+1];
		dy=py[i]-py[i+1]; 
		cost+=sqrtf( (dx*dx) + (dy*dy) );

	}
	dx=px[i]-px[0];
	dy=py[i]-py[0]; 
	cost+=sqrtf( (dx*dx) + (dy*dy) );

	return cost;

}

/*int* randRoute(int dim)
{
	int i,j;
	int *seq_city;
	int *a;
	seq_city = (int*) malloc(sizeof(int)*dim);
	a = (int*) malloc(sizeof(int)*dim);
	
	//srand ( time(NULL) );
	for (i = 0; i < dim; ++i)
	{
	a[i]=0;
	}
	i=0;
	while(i<dim)
	{
		j=rand() % dim;
		if(a[j])
		{
			continue;
		}
		else
		{
			seq_city[i]=j;
			a[j]=1;
			i++;	
		}

	}
	return seq_city;

}*/


int main(int argc, char *argv[])
{
	int ch, cnt, in1, cities;
	float in2, in3;
	FILE *f;
	float *posx, *posy;
	char str[256];  
	int *r;//*route;
	int dst,d,tid=0;
	int sol,i,j;// x,y;
	float *px,*py;
	clock_t start,end;

	f = fopen(argv[1], "r");
	if (f == NULL) {fprintf(stderr, "could not open file \n");  exit(-1);}

	ch = getc(f);  while ((ch != EOF) && (ch != '\n')) ch = getc(f);
	ch = getc(f);  while ((ch != EOF) && (ch != '\n')) ch = getc(f);
	ch = getc(f);  while ((ch != EOF) && (ch != '\n')) ch = getc(f);

	ch = getc(f);  while ((ch != EOF) && (ch != ':')) ch = getc(f);
	fscanf(f, "%s\n", str);
	cities = atoi(str);
	if (cities <= 2) {fprintf(stderr, "only %d cities\n", cities);  exit(-1);}

	sol=cities*(cities-1)/2;
	posx = (float *)malloc(sizeof(float) * cities);  if (posx == NULL) {fprintf(stderr, "cannot allocate posx\n");  exit(-1);}
	posy = (float *)malloc(sizeof(float) * cities);  if (posy == NULL) {fprintf(stderr, "cannot allocate posy\n");  exit(-1);}
	px = (float *)malloc(sizeof(float) * cities);  if (posy == NULL) {fprintf(stderr, "cannot allocate posy\n");  exit(-1);}
	py = (float *)malloc(sizeof(float) * cities);  if (posy == NULL) {fprintf(stderr, "cannot allocate posy\n");  exit(-1);}
	
	
	r = (int *)malloc(sizeof(int) * cities);
	/*{r[i]=i;}*/
	//r=randRoute(cities);	
		
	ch = getc(f);  while ((ch != EOF) && (ch != '\n')) ch = getc(f);
	fscanf(f, "%s\n", str);
	if (strcmp(str, "NODE_COORD_SECTION") != 0) {fprintf(stderr, "wrong file format\n");  exit(-1);}

	cnt = 0;

	while (fscanf(f, "%d %f %f\n", &in1, &in2, &in3)) 
	{
		posx[cnt] = in2;
		posy[cnt] = in3;
		cnt++;
		if (cnt > cities) {fprintf(stderr, "input too long\n");  exit(-1);}
		if (cnt != in1) {fprintf(stderr, "input line mismatch: expected %d instead of %d\n", cnt, in1);  exit(-1);}
	}

	if (cnt != cities) {fprintf(stderr, "read %d instead of %d cities\n", cnt, cities);  exit(-1);}
	fscanf(f, "%s", str);
	if (strcmp(str, "EOF") != 0) {fprintf(stderr, "didn't see 'EOF' at end of file\n");  exit(-1);}
	
//----------------------------------------------------
	r[0]=0;
	int k=1;i=0;float min;int minj,mini,count=1,flag=0;dst=0;
	int *v=(int*)calloc(cities,sizeof(int));
	v[0]=1;
	while(count!=cities)
	{
		flag=0;
		for(j=1;j<cities;j++)
		{
			if(i!=j && !v[j])
			{
				min=distD(i,j,posx,posy);
				minj=j;
				break;	
			}
		}

		for(j=minj+1;j<cities;j++)
		{
			
				 if( !v[j])
				{
					if(min>distD(i,j,posx,posy))
					{
						min=distD(i,j,posx,posy);
						mini=j;
						//printf("\nmin =%f minI=%d\n",min,mini);
						flag=1;				
					}
				}
		}
		if(flag==0)
			i=minj;
		else
			i=mini;
		dst+=min;
		r[k++]=i;v[i]=1;
		count++;
	}
	free(v);
	dst+=distD(r[0],r[cities-1],posx,posy);

//----------------------------------------------------
	printf("\ninitial cost : %d\n",dst);
	start = clock();
	int cost=0,dist=dst;
	int x=0,y=0;
	register int change=0;
	setCoord(r,posx,posy,px,py,cities);
	count=0;	
	do{
		cost=0;
		dist=dst;
		for(i=0;i<(cities-1);i++)
		{	
		
			for(j=i+1; j<(cities);j++ )
			{
					cost=dist;			
					change=distD(i,j,px,py)+distD(i+1,(j+1)%cities,px,py)-distD(i,(i+1)%cities,px,py)-distD(j,(j+1)%cities,px,py);
					cost+=change;	
					if(cost<dst)
					{
						x=i;
						y=j;
						dst=cost;
					}
			}

		}

	
		if(dst<dist)
		{
			float *tmp_x,*tmp_y;
			tmp_x=(float*)malloc(sizeof(float)*(y-x));	
			tmp_y=(float*)malloc(sizeof(float)*(y-x));	
			for(j=0,i=y;i>x;i--,j++)
			{
				tmp_x[j]=px[i];
				tmp_y[j]=py[i];
			}
			for(j=0,i=x+1;i<=y;i++,j++)
			{
				px[i]=tmp_x[j];
				py[i]=tmp_y[j];
			}
			free(tmp_x);
			free(tmp_y);
		}
	count++;
	}while(dst<dist);
	printf("\nMinimal distance found %d\n",dst);
	printf("\nnumber of time hill climbed %d\n",count);
	end = clock();
	printf("\ntime : %f\n",((double) (end - start)) / CLOCKS_PER_SEC);

	free(posx);
	free(posy);
	free(px);
	free(py);
	free(r);
	return 0;
}

