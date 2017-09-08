#include"stdio.h"
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include"math.h"
#include <ctype.h>
//int *tmp;
texture<float, 1, cudaReadModeElementType> t_posx;
texture<float, 1, cudaReadModeElementType> t_posy;

int distNN(int i,int j,float *x,float*y)
{
	unsigned int dx=x[i]-x[j];
	unsigned int dy=y[i]-y[j]; 
	return(sqrtf( (dx*dx) + (dy*dy) ));
}

__device__ int distD(int i,int j)
{
int dx=tex1D(t_posx,i)-tex1D(t_posx,j);
int dy=tex1D(t_posy,i)-tex1D(t_posy,j); 
return(sqrtf( (dx*dx) + (dy*dy) ));
}


__global__ void tsp(int cost,unsigned long long *dst_tid,int cit)
{

	int i,j;
	register int change=0;
	long sol=(cit)*(cit-1)/2;
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	if(id<sol)
	{
		
		i=cit-2-floorf(((int)__dsqrt_rn(8*(sol-id-1)+1)-1)/2);
		j=id-i*(cit-1)+(i*(i+1)/2)+1;
		if(i<j && i<cit && j<cit && i!=j-1)
		{
		change=distD(i,j)+distD((i+1)%cit,(j+1)%cit)-distD(i,(i+1)%cit)-distD(j,(j+1)%cit);
		cost+=change;	
		if(change < 0)
			 atomicMin(dst_tid, ((unsigned long long)cost << 32) | id);
		}
		if(i<0 || j <0 || i>=cit ||j>=cit)
			printf("\nD id=%d i=%d j=%d cost=%d change:%d\n",id,i,j,cost,change);
	}
	
}

int distH(float *px,float *py,int cit)
{
	int dx,dy,cost=0;
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

void twoOpt(int x,int y,float *pox,float *poy)
{
		float *tmp_x,*tmp_y;
		int i,j;
		tmp_x=(float*)malloc(sizeof(float)*(y-x));	
		tmp_y=(float*)malloc(sizeof(float)*(y-x));	
		for(j=0,i=y;i>x;i--,j++)
		{
			tmp_x[j]=pox[i];
			tmp_y[j]=poy[i];
		}
		for(j=0,i=x+1;i<=y;i++,j++)
		{
			pox[i]=tmp_x[j];
			poy[i]=tmp_y[j];
		}
		free(tmp_x);
		free(tmp_y);

}


void setCoord(int *r,float *posx,float *posy,float *px,float *py,int cities)
{

	for(int i=0;i<cities;i++)
	{
	px[i]=posx[r[i]];

	py[i]=posy[r[i]];
	
	}
}

static void CudaTest(char *msg)
{
	cudaError_t e;
	cudaThreadSynchronize();
	if (cudaSuccess != (e = cudaGetLastError()))
	{
		fprintf(stderr, "%s: %d\n", msg, e);
		fprintf(stderr, "%s\n", cudaGetErrorString(e));
		exit(-1);
	}
}
 
int main(int argc, char *argv[])
{
	int ch, cnt, in1, cities;
	float in2, in3;
	FILE *f;
	float *posx, *posy;
	float *px, *py;
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
	px = (float *)malloc(sizeof(float) * cities);  if (posx == NULL) {fprintf(stderr, "cannot allocate posx\n");  exit(-1);}
	py = (float *)malloc(sizeof(float) * cities);  if (posy == NULL) {fprintf(stderr, "cannot allocate posy\n");  exit(-1);}
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
				min=distNN(i,j,posx,posy);
				minj=j;
				break;	
			}
		}

		for(j=minj+1;j<cities;j++)
		{
			
				 if( !v[j])
				{
					if(min>distNN(i,j,posx,posy))
					{
						min=distNN(i,j,posx,posy);
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
	dst+=distNN(0,r[cities-1],posx,posy);

//----------------------------------------------------
	setCoord(r,posx,posy,px,py,cities);
	count=0;
	start = clock();
	cudaEvent_t strt, stp;
	cudaEventCreate(&strt);
	cudaEventCreate(&stp);
 	unsigned long long dst_tid = (((long)dst+1) << 32) -1;
        unsigned long long dtid;

	cudaArray* d_posx = NULL; 
	cudaMallocArray(&d_posx, &t_posx.channelDesc, cities, 1); 
	cudaMemcpyToArray(d_posx, 0, 0, px, cities * sizeof(float), cudaMemcpyHostToDevice); 
	cudaBindTextureToArray(t_posx, d_posx); 
	t_posx.normalized = false; 
	t_posx.addressMode[0] = cudaAddressModeClamp;

	cudaArray* d_posy = NULL; 
	cudaMallocArray(&d_posy, &t_posy.channelDesc, cities, 1); 
	cudaMemcpyToArray(d_posy, 0, 0, py, cities * sizeof(float), cudaMemcpyHostToDevice); 
	cudaBindTextureToArray(t_posy, d_posy); 
	t_posy.normalized = false; 
	t_posy.addressMode[0] = cudaAddressModeClamp;

	printf("\ninitial cost : %d\n",dst);
	if(cudaSuccess!=cudaMalloc((void**)&d_dst_tid,sizeof(unsigned long long)))printf("\nAllocating memory for dst_tid on GPU");
    	if(cudaSuccess!=cudaMemcpy(d_dst_tid,&dst_tid,sizeof(unsigned long long),cudaMemcpyHostToDevice))printf("\ntransfer on GPU");
	cudaEventRecord(strt,0);	
	tsp<<<blk,thrd>>>(dst,d_dst_tid,cities);
	CudaTest("kernel launch failed");
	cudaEventRecord(stp,0);		cudaEventSynchronize(stp);
	if(cudaSuccess!=cudaMemcpy(&dtid,d_dst_tid,sizeof(unsigned long long),cudaMemcpyDeviceToHost))
	printf("\nCan't transfer minimal cost back to CPU");

	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, strt, stp);
  	d = dtid >> 32;
	printf("\nfirst cost found %d ",d);	
	while( d < dst )
	{
		dst=d;
		tid = dtid & ((1ull<<32)-1); 
		x=cities-2-floor((sqrt(8*(sol-tid-1)+1)-1)/2);
		y=tid-x*(cities-1)+(x*(x+1)/2)+1;
		twoOpt(x,y,px,py);

		unsigned long long dst_tid = (((long)dst+1) << 32) -1;
		if(cudaSuccess!=cudaMemcpy(d_dst_tid,&dst_tid,sizeof(unsigned long long),cudaMemcpyHostToDevice))
		printf("\ntransfer on GPU");
		
		cudaMemcpyToArray(d_posx, 0, 0, px, cities * sizeof(float), cudaMemcpyHostToDevice); 
		cudaBindTextureToArray(t_posx, d_posx); 

		cudaMemcpyToArray(d_posy, 0, 0, py, cities * sizeof(float), cudaMemcpyHostToDevice); 
		cudaBindTextureToArray(t_posy, d_posy); 

		tsp<<<blk,thrd>>>(dst,d_dst_tid,cities);
		CudaTest("kernel launch failed");
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
	long long climb=1LL*count * sol;
	printf("\ngmoves/sec :%f\n",climb * 0.000000001 / t);
	cudaFree(d_posy);
	cudaFree(d_posx);
	cudaFree(d_dst_tid);
	free(posx);
	free(posy);
	free(px);
	free(py);
	free(r);
	return 0;
}
