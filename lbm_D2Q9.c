#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mt19937ar.h"

#define Lx 400//lattice sites size along x
#define Ly 100//lattice sites size along y
#define N (Lx*Ly)
#define Q 9//number of velocities
#define D 2//lattice dimension
#define tau 0.6//characteristic timescale
#define rho0 100//initial density at each site
#define rad 15.0//radius of the obstacle
#define GAP 100
#define T_EQ 100000
#define boundary_flag 0//0 for periodic along both x and y, 1 for periodic along x and reflecting along y
#define obstacle_flag 1// 0 for circular, 1 for rctangular obstacle
#define gnu_flag 2//0 for density plot, 1 for velocity plot, 2 for vorticity plot

int nbr[N][Q-1],lat[N];
double fi[N*Q],temp_fi[N*Q],vel[Q][D],vorticity[N];

FILE *pipe;

void take_input()
{
	/*reads in seed for random number generator  and initialises random number generator*/
	long seedval;
	//printf("\nEnter value of seed : ");
	//scanf("%ld",&seedval);
	seedval=time(NULL);
	//seedval=9656546;	
	init_genrand(seedval);
}
//given x, y it gives site index
int find_site(int x,int y)
{
    x=(x+Lx)%Lx;
    y=(y+Ly)%Ly;
    return x+y*Lx;
}
//neighbour index
void neighbour_index()
{
	int x,y,s;	
	for(x=0;x<Lx;x++)
   	{
   		for(y=0;y<Ly;y++)
      		{
      			s=find_site(x,y);
         		nbr[s][0]=find_site(x+1,y);
         		nbr[s][1]=find_site(x,y+1);
         		nbr[s][2]=find_site(x-1,y);
         		nbr[s][3]=find_site(x,y-1);
			nbr[s][4]=find_site(x+1,y+1);
         		nbr[s][5]=find_site(x-1,y+1);
			nbr[s][6]=find_site(x-1,y-1);
			nbr[s][7]=find_site(x+1,y-1);
   		}
	}
}//initialize lattice with obstacle
void initialize_lattice_c(int site, double r)
{
	int i,x,y,x1,y1;
	x=site%Lx;
	y=site/Lx;
	for(i=0;i<N;i++)
	{
		x1=i%Lx;
		y1=i/Lx;
		//circle
		if(sqrt((x-x1)*(x-x1)*1.0+(y-y1)*(y-y1)*1.0)<=r)
		{
			lat[i]=1;
		}
		else
		{
			lat[i]=0;
		}
	}
}
void initialize_lattice_r(int site, int len)
{
	int i,x,y,j,s;
	x=site%Lx;
	y=site/Lx;
	for(i=0;i<N;i++)
	{
		lat[i]=0;
	}
	//rectangle
	for(j=0;j<len;j++)
	{
		s=find_site(x,y-len/2+j);
		lat[s]=1;
		lat[nbr[s][0]]=1;
	}
}
//make initial random density distribution
void make_initial_density_distribution()
{
	int i,j;
	double sm;
	for(i=0;i<N*Q;i++)
	{
		fi[i]=1.0+genrand_real3();
	}
	for(i=0;i<N;i++)
	{
		sm=0.0;
		for(j=0;j<Q;j++)
		{
			sm+=fi[i*Q+j];
		}
		for(j=0;j<Q;j++)
		{
			fi[i*Q+j]=fi[i*Q+j]*rho0/sm;
		}
	}
}
void make_initial_density_distribution_drift_x()
{
	int i,x,j;
	double sm;
	for(i=0;i<N*Q;i++)
	{
		fi[i]=1.0+(0.02*genrand_real3()-0.01);
	}
	for(i=0;i<N;i++)
	{
		x=i%Lx;
		fi[i*Q+1]+=2.0*(1.0+0.2*cos(8.0*M_PI*x/(1.0*Lx)));
		//printf("i=%d\tfi=%.8e\n",i,fi[i]);
		sm=0.0;
		for(j=0;j<Q;j++)
		{
			sm+=fi[i*Q+j];
		}
		if(lat[i]==0)
		{
			for(j=0;j<Q;j++)
			{
				fi[i*Q+j]=fi[i*Q+j]*(1.0*rho0/sm);
			}
		}
		else
		{
			for(j=0;j<Q;j++)
			{
				fi[i*Q+j]=0.0;
			}
		}
	}
}
//initialization
void initialize_velocity_parameters()
{
	vel[0][0]=0.0;
	vel[0][1]=0.0;
	vel[1][0]=1.0;
	vel[1][1]=0.0;
	vel[2][0]=0.0;
	vel[2][1]=1.0;
	vel[3][0]=-1.0;
	vel[3][1]=0.0;
	vel[4][0]=0.0;
	vel[4][1]=-1.0;
	vel[5][0]=1.0;
	vel[5][1]=1.0;
	vel[6][0]=-1.0;
	vel[6][1]=1.0;
	vel[7][0]=-1.0;
	vel[7][1]=-1.0;
	vel[8][0]=1.0;
	vel[8][1]=-1.0;
}
//get local equilibrium density
void get_local_equilibrium_density(int site, double av_density, double av_velocity[D], double feq[Q])
{
	double usq,cdotu;
	usq=av_velocity[0]*av_velocity[0]+av_velocity[1]*av_velocity[1];
	feq[0]=av_density*(4.0/9.0)*(1.0-1.5*usq);
	cdotu=vel[1][0]*av_velocity[0];
	feq[1]=av_density*(1.0/9.0)*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*usq);
	cdotu=vel[2][1]*av_velocity[1];
	feq[2]=av_density*(1.0/9.0)*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*usq);
	cdotu=vel[3][0]*av_velocity[0];	
	feq[3]=av_density*(1.0/9.0)*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*usq);
	cdotu=vel[4][1]*av_velocity[1];	
	feq[4]=av_density*(1.0/9.0)*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*usq);
	cdotu=vel[5][0]*av_velocity[0]+vel[5][1]*av_velocity[1];	
	feq[5]=av_density*(1.0/36.0)*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*usq);
	cdotu=vel[6][0]*av_velocity[0]+vel[6][1]*av_velocity[1];	
	feq[6]=av_density*(1.0/36.0)*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*usq);
	cdotu=vel[7][0]*av_velocity[0]+vel[7][1]*av_velocity[1];	
	feq[7]=av_density*(1.0/36.0)*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*usq);
	cdotu=vel[8][0]*av_velocity[0]+vel[8][1]*av_velocity[1];	
	feq[8]=av_density*(1.0/36.0)*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*usq);
}
//calculate delta fi
void get_delta_fi(int site, double feq[Q], double delta_fi[Q])
{
	int i,s;
	for(i=0;i<Q;i++)
	{
		s=site*Q+i;
		delta_fi[i]=(feq[i]-fi[s])/(1.0*tau);
	}
}
//streaming 
void streaming_step_periodic_xy()
{
	int i,s;
	for(i=0;i<N;i++)
	{
		if(lat[i]==0)
		{
			fi[i*Q]=temp_fi[i*Q];
		}
		s=nbr[i][1];
		if(lat[s]==1)
		{
			fi[i*Q+4]=temp_fi[i*Q+2];
		}
		else
		{
			fi[i*Q+4]=temp_fi[s*Q+4];
		}
		s=nbr[i][4];
		if(lat[s]==1)
		{
			fi[i*Q+7]=temp_fi[i*Q+5];
		}
		else
		{
			fi[i*Q+7]=temp_fi[s*Q+7];
		}
		s=nbr[i][5];
		if(lat[s]==1)
		{
			fi[i*Q+8]=temp_fi[i*Q+6];
		}
		else
		{
			fi[i*Q+8]=temp_fi[s*Q+8];
		}
		s=nbr[i][3];
		if(lat[s]==1)
		{
			fi[i*Q+2]=temp_fi[i*Q+4];
		}
		else
		{
			fi[i*Q+2]=temp_fi[s*Q+2];
		}
		s=nbr[i][6];
		if(lat[s]==1)
		{
			fi[i*Q+5]=temp_fi[i*Q+7];
		}
		else
		{
			fi[i*Q+5]=temp_fi[s*Q+5];
		}
		s=nbr[i][7];
		if(lat[s]==1)
		{
			fi[i*Q+6]=temp_fi[i*Q+8];
		}
		else
		{
			fi[i*Q+6]=temp_fi[s*Q+6];
		}
		s=nbr[i][2];
		if(lat[s]==1)
		{
			fi[i*Q+1]=temp_fi[i*Q+3];
		}
		else
		{
			fi[i*Q+1]=temp_fi[s*Q+1];
		}
		s=nbr[i][0];
		if(lat[s]==1)
		{
			fi[i*Q+3]=temp_fi[i*Q+1];
		}
		else
		{
			fi[i*Q+3]=temp_fi[s*Q+3];
		}
	}
}
void streaming_step_reflecting_y()
{
	int i,s,x,y;
	for(i=0;i<N;i++)
	{
		x=i%Lx;
		y=i/Lx;
		if(lat[i]==0)
		{
			fi[i*Q]=temp_fi[i*Q];
					
			if(y==0)
			{
				fi[i*Q+2]=temp_fi[i*Q+4];
				fi[i*Q+5]=temp_fi[i*Q+7];
				fi[i*Q+6]=temp_fi[i*Q+8];
				s=nbr[i][1];
				fi[i*Q+4]=temp_fi[s*Q+4];
				s=nbr[i][4];
				fi[i*Q+7]=temp_fi[s*Q+7];
				s=nbr[i][5];
				fi[i*Q+8]=temp_fi[s*Q+8];

				s=nbr[i][2];
				fi[i*Q+1]=temp_fi[s*Q+1];
				s=nbr[i][0];
				fi[i*Q+3]=temp_fi[s*Q+3];
			}
			else if(y==Ly-1)
			{
				fi[i*Q+4]=temp_fi[i*Q+2];
				fi[i*Q+7]=temp_fi[i*Q+5];
				fi[i*Q+8]=temp_fi[i*Q+6];
				s=nbr[i][3];
				fi[i*Q+2]=temp_fi[s*Q+2];
				s=nbr[i][6];
				fi[i*Q+5]=temp_fi[s*Q+5];
				s=nbr[i][7];
				fi[i*Q+6]=temp_fi[s*Q+6];
		
				s=nbr[i][2];
				fi[i*Q+1]=temp_fi[s*Q+1];
				s=nbr[i][0];
				fi[i*Q+3]=temp_fi[s*Q+3];
			}
			else
			{
				s=nbr[i][1];
				if(lat[s]==1)
				{
					fi[i*Q+4]=temp_fi[i*Q+2];
				}
				else
				{
					fi[i*Q+4]=temp_fi[s*Q+4];
				}
				s=nbr[i][4];
				if(lat[s]==1)
				{
					fi[i*Q+7]=temp_fi[i*Q+5];
				}
				else
				{
					fi[i*Q+7]=temp_fi[s*Q+7];
				}
				s=nbr[i][5];
				if(lat[s]==1)
				{
					fi[i*Q+8]=temp_fi[i*Q+6];
				}
				else
				{
					fi[i*Q+8]=temp_fi[s*Q+8];
				}
				s=nbr[i][3];
				if(lat[s]==1)
				{
					fi[i*Q+2]=temp_fi[i*Q+4];
				}
				else
				{
					fi[i*Q+2]=temp_fi[s*Q+2];
				}
				s=nbr[i][6];
				if(lat[s]==1)
				{
					fi[i*Q+5]=temp_fi[i*Q+7];
				}
				else
				{
					fi[i*Q+5]=temp_fi[s*Q+5];
				}
				s=nbr[i][7];
				if(lat[s]==1)
				{
					fi[i*Q+6]=temp_fi[i*Q+8];
				}
				else
				{
					fi[i*Q+6]=temp_fi[s*Q+6];
				}
				s=nbr[i][2];
				if(lat[s]==1)
				{
					fi[i*Q+1]=temp_fi[i*Q+3];
				}
				else
				{
					fi[i*Q+1]=temp_fi[s*Q+1];
				}
				s=nbr[i][0];
				if(lat[s]==1)
				{
					fi[i*Q+3]=temp_fi[i*Q+1];
				}
				else
				{
					fi[i*Q+3]=temp_fi[s*Q+3];
				}
			}
		}				
	}
}
//evolve
void evolve()
{
	int i,j,x;
	double den,av_velocity[D],feq[Q],delta_fi[Q];
	double sm;
	for(i=0;i<N;i++)
	{
		x=i%Lx;
		if(lat[i]==0)
		{
			av_velocity[0]=0.0;
			av_velocity[1]=0.0;
			den=0.0;
			
			for(j=0;j<Q;j++)
			{
				av_velocity[0]+=fi[i*Q+j]*vel[j][0];
				av_velocity[1]+=fi[i*Q+j]*vel[j][1];
				den+=fi[i*Q+j];
			}
			av_velocity[0]=av_velocity[0]/(1.0*den);
			av_velocity[1]=av_velocity[1]/(1.0*den);
			get_local_equilibrium_density(i,den,av_velocity,feq);
			sm=0.0;
			for(j=0;j<Q;j++)
			{
				sm+=feq[j];
			}
			get_delta_fi(i,feq,delta_fi);
			
			for(j=0;j<Q;j++)
			{
				temp_fi[i*Q+j]=fi[i*Q+j]+delta_fi[j];
			}
			
		}
	}
	//streaming to neighbours
	if(boundary_flag==0)
	{
		streaming_step_periodic_xy();
	}
	else if(boundary_flag==1)
	{
		streaming_step_reflecting_y();
	}
}
//data file
void create_datafile(int time)
{
	int i,j;
	double sm,av_velocity[D],vel_x,vel_y,vel_amp;
	double av_vel[N][D],site_den[N];
	char outfile[100];

	sprintf(outfile,"datafileLx%dLy%dt%d",Lx,Ly,time);
	FILE *fp;
	fp=fopen(outfile,"w");
   	fprintf(fp,"#x\ty\tdensity\n");
   	fclose(fp);
	
	fp=fopen(outfile,"a");
	for(i=0;i<N;i++)
	{
		sm=0.0;
		av_velocity[0]=0.0;
		av_velocity[1]=0.0;
		for(j=0;j<Q;j++)
		{
			sm+=fi[i*Q+j];
			av_velocity[0]+=fi[i*Q+j]*vel[j][0];
			av_velocity[1]+=fi[i*Q+j]*vel[j][1];
		}
		vel_x=av_velocity[0]/(1.0*sm);
		vel_y=av_velocity[1]/(1.0*sm);
		
		av_vel[i][0]=vel_x;
		av_vel[i][1]=vel_y;
		site_den[i]=sm;
	}
	for(i=0;i<N;i++)
	{
		if(lat[i]==0)
		{
			vorticity[i]=(av_vel[nbr[i][0]][1]-av_vel[i][1])-(av_vel[nbr[i][1]][0]-av_vel[i][0]);
		}
		else
		{
			vorticity[i]=0.0;
		}
		vel_amp=sqrt(1.0*av_vel[i][0]*av_vel[i][0]+1.0*av_vel[i][1]*av_vel[i][1]);
		fprintf(fp,"%d\t%d\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",(i%Lx),(i/Lx),fi[i*Q+0],fi[i*Q+1],fi[i*Q+2],fi[i*Q+3],fi[i*Q+4],fi[i*Q+5],fi[i*Q+6],fi[i*Q+7],fi[i*Q+8],site_den[i],av_vel[i][0]/(1.0*vel_amp),av_vel[i][1]/(1.0*vel_amp),vorticity[i]);
	}
	fclose(fp);
}
void call_gnuplot_vorticity(int t)
{
	fprintf(pipe,"set xrange [%d:%d]\n",0,Lx-1);
	fprintf(pipe,"set yrange [%d:%d]\n",0,Ly-1);
	fprintf(pipe,"set cbrange [-0.01:0.01]\n");
	fprintf(pipe,"set size ratio %3.4lf\n",(1.0*Ly/Lx));
	fprintf(pipe,"set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')\n");
	fprintf(pipe, "set title 'time=%d'\n",t);

	fprintf(pipe,"p 'datafileLx%dLy%dt%d' u 1:2:15 w image tit ''\n",Lx,Ly,t);
	fflush(pipe);
}
void call_gnuplot_density(int t)
{
	fprintf(pipe,"set xrange [%d:%d]\n",0,Lx-1);
	fprintf(pipe,"set yrange [%d:%d]\n",0,Ly-1);
	fprintf(pipe,"set size ratio %3.4lf\n",(1.0*Ly/Lx));
	fprintf(pipe,"set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')\n");
	fprintf(pipe, "set title 'time=%d'\n",t);

	fprintf(pipe,"p 'datafileLx%dLy%dt%d' u 1:2:4 w image tit ''\n",Lx,Ly,t);
	fflush(pipe);
}
void call_gnuplot_velocity(int t)
{
	fprintf(pipe,"set xrange [%d:%d]\n",0,Lx-1);
	fprintf(pipe,"set yrange [%d:%d]\n",0,Ly-1);
	fprintf(pipe,"set size ratio %3.4lf\n",(1.0*Ly/Lx));
	fprintf(pipe, "set title 'time=%d'\n",t);

	fprintf(pipe,"p 'datafileLx%dLy%dt%d' u 1:2:13:14 w vectors lc rgb 'red' tit ''\n",Lx,Ly,t);
	fflush(pipe);
}
int main()
{
	take_input();
	neighbour_index();

	if(obstacle_flag==0)
	{
		initialize_lattice_c(40+Lx*Ly/2,rad);
	}
	else if(obstacle_flag==1)
	{
		initialize_lattice_r(20+Lx*Ly/2,20);
	}
	
	//make_initial_density_distribution();
	make_initial_density_distribution_drift_x();
	initialize_velocity_parameters();
	
	pipe=popen("\\gnuplot","w");
	fprintf(pipe, "unset xtics\n");
	fprintf(pipe, "unset ytics\n");
	
	int t;
	for(t=1;t<=T_EQ;t++)
	{
		if(t%GAP==0)
		{
			create_datafile(t);
			if(gnu_flag==0)
			{
				call_gnuplot_density(t);
			}
			else if(gnu_flag==1)
			{
			 	call_gnuplot_velocity(t);
			}
			else if(gnu_flag==2)
			{
			 	call_gnuplot_vorticity(t);
			}
		}
		evolve();
	}
	fclose(pipe);
}
