#include"stdio.h"
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include"math.h"
#include <ctype.h>

int distD(int i,int j,float **dist_mat)
{
	return(dist_mat[i][j]);
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

int* randRoute(int dim)
{
	int i,j;
	int *seq_city;
	int *a;
	seq_city = (int*) malloc(sizeof(int)*dim);
	a = (int*) malloc(sizeof(int)*dim);
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

}

int DistH(int*route,int N,float **dist_mat)
{
	int i,cost=0;
	for(i=0;i<N-1;i++)
	{
		cost+=dist_mat[route[i]][route[i+1]];
	}
	cost+=dist_mat[route[i]][route[0]];
	return cost;
}

int main(int argc, char *argv[])
{
	int ch, cnt, in1, cities;
	float in2, in3;
	FILE *f;
	float *posx, *posy;
	char str[256];  
	int *r;
	int dst,d,tid=0;
	int sol,i,j;// x,y;
	
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
	
	r = (int *)malloc(sizeof(int) * cities);
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

	float **dist_mat;
	dist_mat = (float **)malloc(sizeof(float*) * cities);
	for(i=0;i<cities;i++)
	{dist_mat[i] = (float *)malloc(sizeof(float) * cities);}

	for (i = 0; i < cities; ++i)
	{
		for (j = 0; j < cities; ++j)
		{
		dist_mat[i][j] = sqrtf(pow(posx[i] - posx[j], 2)
		             +powf(posy[i] - posy[j], 2));
		
		}
	}
//---------------------route construction-------------------------
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
				min=dist_mat[i][j];
				minj=j;
				break;	
			}
		}

		for(j=minj+1;j<cities;j++)
		{
			
				 if( !v[j])
				{
					if(min>dist_mat[i][j])
					{
						min=dist_mat[i][j];
						mini=j;
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
	dst+=dist_mat[r[0]][r[cities-1]];
//-------------------------------------------------------------------
	printf("\ninitial cost : %d\n",dst);
	start = clock();
	float cost=0,dist=dst;
	float x=0,y=0;
	register int change=0;
	count=0;	
	do{
		cost=0;
		dist=dst;
			for(i=0;i<(cities-1);i++)
			{	
		
				for(j=i+1; j<(cities);j++ )
				{
					cost=dist;			
					change=distD(r[i],r[j],dist_mat)+distD(r[i+1],r[(j+1)%cities],dist_mat)-distD(r[i],r[(i+1)%cities],dist_mat)-distD(r[j],r[(j+1)%cities],dist_mat);
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
				float *tmp_x;
				tmp_x=(float*)malloc(sizeof(float)*(y-x));	
				for(j=0,i=y;i>x;i--,j++)
				{
					tmp_x[j]=r[i];
				}
				for(j=0,i=x+1;i<=y;i++,j++)
				{
					r[i]=tmp_x[j];
				}
				free(tmp_x);
			}
		count++;
	}while(dst<dist);

	printf("\nMinimal distance found %d\n",dst);
	printf("\nnumber of time hill climbed %d\n",count);
	end = clock();
	printf("\ntime : %f\n",((double) (end - start)) / CLOCKS_PER_SEC);

	free(posx);
	free(posy);
	return 0;
}

