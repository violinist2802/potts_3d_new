// file: init.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
int* alloc_attach(int NRc)
{
	int i;
	int * attached;

	attached = calloc(NRc+1, sizeof(int));

	for(i=0;i<NRc+1;i++)
		attached[i]=0;
	
	return attached;
}

////////////////////////////////////////////////////////////////////////////////
VOX* init_voxels(int NVX, int NVY, int NVZ)
{
	VOX* pv;
	int v, vx, vy, vz;
	int i;

   	pv = calloc(NV, sizeof(VOX));

	// set voxel information
   	for(vy=0; vy<NVY; vy++)
   	for(vx=0; vx<NVX; vx++)
	for(vz=0; vz<NVZ; vz++)
   	{
   		v = vx + vy*NVX+vz*NVX*NVY;

		//pv[v].x = vx * VOXSIZE; pv[v].y = vy * VOXSIZE;
		pv[v].ctag = 0;
		pv[v].type = 0;
		pv[v].contact = 0;
		pv[v].bond = 0;
	}
	return pv;
}

////////////////////////////////////////////////////////////////////////////////
int init_cells(VOX* pv, int * types, BOX* pb, int NCX, int NCY, int NCZ, double PART, int shifts, double TARGETVOLUME_FB, double VOXSIZE, int NVX, int NVY, int NVZ)
{
	int v, vx, vy, vz, i, j, k, ix, iy, iz;
	int NRc;
	double r01;
	double d;
	double dx, dy, dz, dvx,dvy, dvz; // distance to center
	int r;

	NRc = 0;

	dx = (double) (NVX - 2 * MARGINX) / (NCX);
	dy = (double) (NVY - 2 * MARGINY) / (NCY);
	dz = (double) (NVZ - 2 * MARGINY) / (NCZ);
	r  = (int) (cbrt(STARTVOLUME));
	if(dx<2*r || dy<2*r || dz<2*r)
		printf("Too dense!");
	for (iy = 0; iy < NCY; iy++){
		for (ix = 0; ix < NCX; ix++){
			for (iz = 0; iz < NCZ; iz++){
				dvx = (mt_random()%((int) dx-2*r+1)) -(dx/2 - r);
				dvy = (mt_random()%((int) dy-2*r+1)) -(dy/2 - r);
				dvz = (mt_random()%((int) dz-2*r+1)) -(dz/2 - r);
				vx = MARGINX + (int) (((double) ix + 0.5) * dx + shifts*dvx);
				vy = MARGINY + (int) (((double) iy + 0.5) * dy + shifts*dvy);
				//vz = MARGINZ + (int) (((double) iz + 0.5) * dz + shifts*dvz);
				//vx = (int)(NVX/2);
				//vy = (int)(NVY/2);
				vz = (int)(NVZ/2);

				NRc++;
				types[NRc] = (PART<(rand()/(double)RAND_MAX) ? 1 : 2);
				pb[NRc].x1 = vx-r;
				pb[NRc].x2 = vx+r;
				pb[NRc].y1 = vy-r;
				pb[NRc].y2 = vy+r;
				pb[NRc].z1 = vz-r;
				pb[NRc].z2 = vz+r;
				for(i = -r; i<=r; i++){
					for (j = -r; j<=r; j++){
						for (k = -r; k<=r; k++){
							v = (vx + i) + (vy + j)*NVX + (vz + k)*NVX*NVY;
							if (v<NV){
								pv[v].ctag = NRc;
								pv[v].type = types[NRc]; 
							}
							else
								printf("Cell out of area: (%d,%d,%d)\n",vx+i-NVX,vy+j-NVY,vz+k-NVZ);
						}
					}
				}
			}
		}
		}
			
	

	return NRc;
}

////////////////////////////////////////////////////////////////////////////////
FIBERS* set_fibers(double distanceF, double VOXSIZE, int NVX, int NVY, int NVZ) //idk with fibers angle
{
	FIBERS* pf;
	
	int v, vx, vy,vz, vd, fd, kc;
	int i;
	double dx,dy,dz;
	double k, k0;

	dx = F_ANGLE!=0    ? F_DISTANCE / sin(F_ANGLE) : 0.0;
	dy = F_ANGLE!=0 ? F_DISTANCE / sin(F_ANGLE) : 0.0;
	dz = F_ANGLE!=PI/2 ? F_DISTANCE / cos(F_ANGLE) : 0.0;

   	pf = calloc(NV, sizeof(FIBERS));

	// set voxel information
    for(v=0; v<NV; v++)
    	pf[v].Q = 0;

   /*	if(F_ANGLE!=PI/2){
   		fd = round(fiberD / VOXSIZE / cos(F_ANGLE));
	   	fd = fd<1 ? 1 : fd;
	   	for(vx=0; vx<NVX; vx++){
	   		k0 = fmod(vx*tan(F_ANGLE),dy);
	   		k = k0;
	   		for(vy=0; vy<=(NVY-k0)/dy; vy++){
	   			for(vd=0; vd<fd; vd++){
	   				kc = round(k)+vd;
	   				v = vx + kc*NVX;
	   				if(kc<NVY)
	   					pf[v].Q = 1;
	   			}
	   			k += dy;
	   		}
	   	}
	}

	if(F_ANGLE!=0){
		fd = round(fiberD / VOXSIZE / sin(F_ANGLE));
	   	fd = fd<1 ? 1: fd;
	   	printf("fiberdX: %d\n", fd);
	   	for(vy=0; vy<NVY; vy++){
	   		k0 = fmod(vy/tan(F_ANGLE),dx);
	   		k = k0;
	   		for(vx=0; vx<=(NVX-k0)/dx; vx++){
	   			for(vd=0; vd<fd; vd++){
	   				kc = round(k)+vd;
		   			v = kc + vy*NVX;
		   			if(kc<NVX)
		   				pf[v].Q = 1;
	   			}
	   			k += dx;
	   		}
	   	}
	}*/
	for(vx=0; vx<NVX; vx+=F_DISTANCE){
		for(vy=0; vy<=NVY; vy+=F_DISTANCE){
			for(vz=0; vz<NVZ; vz++){
				v= vx+vy*NVX+vz*NVX*NVY;
				pf[v].Q = 1;
			}
		}
		}
	return pf;
}
