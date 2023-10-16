// file: CM.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
CM* allocCM(int NRc)
// alloc CM
{
	CM* CMs;

	CMs = calloc(NRc+1, sizeof(CM));
	
	return CMs;
}

BOX* allocBOX(int NRc)
{
	BOX* box;

	box = calloc(NRc+1,sizeof(BOX));

	return box;
}

////////////////////////////////////////////////////////////////////////////////
void findCM(VOX* pv, CM* CMs, int NRc, int NVX, int NVY, int NVZ)
// find center of mass
{
	int vx,vy,vz,v,i;
	long int x[NRc+1], y[NRc+1], z[NRc+1], n[NRc+1];

	for(i=0;i<=NRc;i++){
		x[i]=0;
		y[i]=0;
		z[i]=0;
		n[i]=0;
		
	}
	
	for(vx=0; vx<NVX; vx++) {
        for (vy=0; vy<NVY; vy++) {
			for (vz=0; vz<NVZ; vz++) {
            	v = vx + vy*NVX+vz*NVX*NVY;
				
        		if (pv[v].ctag != 0 ){
					x[pv[v].ctag] += vx;
					y[pv[v].ctag] += vy;
					z[pv[v].ctag] += vz;
					n[pv[v].ctag] ++;
				}
			}
    	}
	}
    for(i=0;i<=NRc;i++){
		
    	CMs[i].x = (double) x[i] / (double) n[i];
    	CMs[i].y = (double) y[i] / (double) n[i];
		CMs[i].z = (double) z[i] / (double) n[i];
		
    }
}

////////////////////////////////////////////////////////////////////////////////
