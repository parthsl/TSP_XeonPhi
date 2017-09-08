#include"stdio.h"
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include"math.h"
#include <ctype.h>
//int *tmp;
__host__  __device__ int distD(int x,int y,int N,float*dt)
{
	int id;
	id=x*N+y;
	return(dt[id]);
}

/*kernel function to generate the feasible solution and determine the best improvement among them */
__global__ void tsp(int *rt,int cost,unsigned long long *dst_tid,int cit,float *dt)
    {

    	int i,j;
    	register int change=0;
    	int sol=(cit)*(cit-1)/2;
    	int id=threadIdx.x+blockIdx.x*blockDim.x;
	
    	if(id<sol)
    	{
    		i=cit-2-floorf(((int)__dsqrt_rn(8*(sol-id-1)+1)-1)/2);
    		j=id-i*(cit-1)+(i*(i+1)/2)+1;
		if(i<j && i<cit && j<cit && i!=j-1)
		{
		
			change=distD(rt[i],rt[j],cit,dt)+distD(rt[(i+1)%cit],rt[(j+1)%cit],cit,dt)-distD(rt[i],rt[(i+1)%cit],cit,dt)-distD(rt[j],rt[(j+1)%cit],cit,dt);
			cost+=change;
			if (change < 0)
				atomicMin(dst_tid, ((unsigned long long)cost << 32) | id);
    		}
	}
    	
    }
/*calculate cost of specified route*/
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

 int* twoOpt(int x,int y,int *route,int city)
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



void setCoord(int *r,float *posx,float *posy,float *px,float *py,int cities)
{

	for(int i=0;i<cities;i++)
	{
		px[i]=posx[r[i]];

		py[i]=posy[r[i]];
	
	}
}
 
int main(int argc, char *argv[])
{
	int ch, cnt, in1, cities;
	float in2, in3;
	FILE *f;
	float *posx, *posy;

	char str[256];  
	int dst,d,tid;
        unsigned long long *d_dst_tid;
	int x,y;
	int blk,thrd;
	clock_t start,end;
	long sol;
	int *r,i,j;
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
	r = (int *)malloc(sizeof(int) * cities);  if (posy == NULL) {fprintf(stderr, "cannot allocate posy\n");  exit(-1);}
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
    	fflush(f);
	fclose(f);
	/*thread and block setting up*/
	if(sol<=50000)
	{
		blk=(sol-1)/512+1;
		thrd=512;
	}
	else
	{
		blk=(sol-1)/1024+1;
		thrd=1024;
	}
	/*generate distance matrix*/
	float *dist_mat;
	dist_mat = (float *)malloc(sizeof(float) * (cities*cities));
	for (int i = 0; i < cities; ++i)
	{
		for (int j = 0; j < cities; ++j)
		{
		dist_mat[i*cities+j] = sqrtf(pow(posx[i] - posx[j], 2)
		             +powf(posy[i] - posy[j], 2));
		//k++;		
		}
	}
//----------------------------------------------------
/*initial route generation and its cost calculation*/
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
				int id;
				if(i>j)
				{id=j*cities+i;}
				else{id=i*cities+j;}	
				min=dist_mat[id];
				minj=j;
				break;	
			}
		}

		for(j=minj+1;j<cities;j++)
		{
			
				 if( !v[j])
				{
					int id;
				if(i>j)
				{id=j*cities+i;}
				else{id=i*cities+j;}	
					if(min>dist_mat[id])
					{
						min=dist_mat[id];
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
	dst+=dist_mat[r[cities-1]];
//-----------------------------------------------------------
	start = clock();
 	unsigned long long dst_tid = (((long)dst+1) << 32) -1;
        unsigned long long dtid;
	int *d_r;
    	float *d_mt;
	printf("\ninitial cost : %d\n",dst);
	/*Allocating memory on GPU */
	if(cudaSuccess!=cudaMalloc((void**)&d_dst_tid,sizeof(unsigned long long)))printf("\nAllocating memory for dst_tid on GPU");
	if(cudaSuccess!=cudaMalloc((void**)&d_mt,sizeof(float)*(cities*cities)))printf("\nAllocating memory for thread id on GPU");
    	if(cudaSuccess!=cudaMalloc((void**)&d_r,sizeof(int)*cities))printf("\nAllocating memory for thread id on GPU");
	/*Data transfer on GPU */
    	if(cudaSuccess!=cudaMemcpy(d_dst_tid,&dst_tid,sizeof(unsigned long long),cudaMemcpyHostToDevice))printf("\ntransfer on GPU");
	if(cudaSuccess!=cudaMemcpy(d_mt,dist_mat,sizeof(float)*(cities*cities),cudaMemcpyHostToDevice))printf("\ntransfer on GPU 1");
    	if(cudaSuccess!=cudaMemcpy(d_r,r,sizeof(int)*cities,cudaMemcpyHostToDevice))printf("\ntransfer on GPU 1");

	tsp<<<blk,thrd>>>(d_r,dst,d_dst_tid,cities,d_mt);

    	if(cudaSuccess!=cudaMemcpy(&dtid,d_dst_tid,sizeof(unsigned long long),cudaMemcpyDeviceToHost))
	printf("\nCan't transfer minimal cost back to CPU");
        d = dtid >> 32;
    	count =0;
	/*This loop contineous until no improvement possible*/
    	while( d < dst )
    	{
    		dst=d;
	        tid = dtid & ((1ull<<32)-1); 
	    	x=cities-2-floor((sqrt(8*(sol-tid-1)+1)-1)/2);
	    	y=tid-x*(cities-1)+(x*(x+1)/2)+1;
	    	r=twoOpt(x,y,r,cities);

                unsigned long long dst_tid = (((long)dst+1) << 32) -1;
    		if(cudaSuccess!=cudaMemcpy(d_r,r,sizeof(int)*cities,cudaMemcpyHostToDevice))printf("\ntransfer on GPU 1");
    	        if(cudaSuccess!=cudaMemcpy(d_dst_tid,&dst_tid,sizeof(unsigned long long),cudaMemcpyHostToDevice))
		printf("\ntransfer on GPU");

    		tsp<<<blk,thrd>>>(d_r,dst,d_dst_tid,cities,d_mt);

    	        if(cudaSuccess!=cudaMemcpy(&dtid,d_dst_tid,sizeof(unsigned long long),cudaMemcpyDeviceToHost))
                printf("\nCan't transfer minimal cost back to CPU");
    		d = dtid >> 32;
                count++;
    	}
printf("\nMinimal Distance : %d\n",d);

printf("\nnumber of time climbed %d\n",count);
end = clock();
double t=((double) (end - start)) / CLOCKS_PER_SEC;
printf("\ntime : %f\n",t);

cudaFree(d_dst_tid);
free(posx);
free(posy);
free(r);
return 0;
}
