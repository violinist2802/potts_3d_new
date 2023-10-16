// file: cellmoves.c
#include "functions.h"

#define MAX_FOCALS_T(a) (a==1 ? MAX_FOCALS_CM : MAX_FOCALS_FB)
#define LMAX(a) (a==1 ? LMAX_CM : LMAX_FB)

////////////////////////////////////////////////////////////////////////////////
double CPM_moves(VOX* pv, short * CCAlabels, BOX* pb, FIBERS* pf, CM* CMs, 
int* attached, int* csize, double MAX_FOCALS_CM, double MAX_FOCALS_FB, double TARGETVOLUME_CM, double TARGETVOLUME_FB, double INELASTICITY_CM, double INELASTICITY_FB, double LMAX_CM, double LMAX_FB, double GN_CM, double GN_FB, double UNLEASH_CM, double UNLEASH_FB, double DETACH_CM, double DETACH_FB, double VOXSIZE, int NVX, int NVY, int NVZ, double JCMCM, double JCMMD, double JFBFB, double JFBMD,  double JFBCM, char CONT, char CONT_INHIB)
// cellular potts model: one Monte Carlo step
{
	int i,j,NRsteps = NV;
	int xs, xt; // source and target pixel
	int xtx,xty,xtz,xsx,xsy,xsz; // x and y position of target pixel
	int ttag, stag, type; // target and source label
	int nbs[26],pick, reject, accept; // neighbors of target pixel
	BOOL go_on;
	double dH, prob, phi;

	accept=0; reject=0;

	for(i=0;i<NRsteps;i++)
	{
		
		xt = mt_random()%NV; // pick random element
		
		xty = (xt%(NVX*NVY))/NVX; xtx = (xt%(NVX*NVY))%NVX; xtz = xt/(NVX*NVY);

		if((xtx>0)&&(xtx<NVX-1)&&(xty>0)&&(xty<NVY-1)&&(xtz>0)&&(xtz<NVZ-1)) // exclude outer rim
		{
			nbs[8]=xt-1+NVX+NVX*NVY; nbs[9]=xt+NVX+NVX*NVY; nbs[10]=xt+1+NVX+NVX*NVY;
			nbs[11]=xt-1+NVX*NVY;	 nbs[12]=xt+NVX*NVY;    nbs[13]=xt+1+NVX*NVY;
			nbs[14]=xt-1-NVX+NVX*NVY;nbs[15]=xt-NVX+NVX*NVY;nbs[16]=xt+1-NVX+NVX*NVY;

			nbs[0]=xt-1+NVX; nbs[1]=xt+NVX; nbs[2]=xt+1+NVX;
			nbs[7]=xt-1;                    nbs[3]=xt+1;
			nbs[6]=xt-1-NVX; nbs[5]=xt-NVX; nbs[4]=xt+1-NVX;

			nbs[17]=xt-1+NVX-NVX*NVY; nbs[18]=xt+NVX-NVX*NVY; nbs[19]=xt+1+NVX-NVX*NVY;
			nbs[20]=xt-1-NVX*NVY;	 nbs[21]=xt-NVX*NVY;    nbs[22]=xt+1-NVX*NVY;
			nbs[23]=xt-1-NVX-NVX*NVY;nbs[23]=xt-NVX-NVX*NVY;nbs[25]=xt+1-NVX-NVX*NVY;
			pick = mt_random()%26;
			xs = nbs[pick]; // pick random neighbor
			xsy = (xs%(NVX*NVY))/NVX; xsx = (xs%(NVX*NVY))%NVX; xsz = xs/(NVX*NVY);

			ttag = pv[xt].ctag;
			stag = pv[xs].ctag;
			
			
			go_on = 0;
			
			
			if(ttag!=stag || (CONT && pv[xs].contact==1) )//don't bother if no difference
			{
				
        		go_on = 1;
        		//if(!pv[xt].contact && !pv[xt].contact && mt_random()%10!=0)			//if contact is not involved -- pick up again
        		//	go_on = 0; 
        		//if(go_on && ttag) // if a cell in xt (retracting)
				
        		if(ttag) // if a cell in xt (retracting)
				{	
					
            		if (splitcheckCCR(pv,CCAlabels,pb,csize,xt,ttag, NVX, NVY, NVZ)){
						
                	go_on = 0;
					}
					
            		if(csize[ttag-1]==1) // cell cannot disappear (constraint may be removed)
                		go_on = 0;
					
				}
				// if(go_on && shifts==1 && distanceF<0.1 && (xtx<MARGINX || xtx>NVX-MARGINX))
				// 	go_on = 0;
				
				if(E_bond && ttag!=stag && pv[xt].type==1 && pv[xs].type==1 && pv[xt].contact==1 && pv[xs].contact==1 && pv[xt].bond==0 && pv[xs].bond==0){
					printf("\nn %d", (int) 2);
					pv[xt].bond=xs;
					pv[xs].bond=xt;
					
				}

			}
			
			if(go_on)
			{
				dH = calcdH(pv,pf,CMs,csize,xt,xs,pick,ttag,stag,TARGETVOLUME_CM, TARGETVOLUME_FB, INELASTICITY_CM, INELASTICITY_FB, LMAX_CM, LMAX_FB, GN_CM, GN_FB, UNLEASH_CM, UNLEASH_FB, DETACH_CM, DETACH_FB, VOXSIZE, NVX, NVY, JCMCM, JCMMD, JFBFB, JFBMD, JFBCM);
        		prob = exp(-IMMOTILITY*dH);
        		if (prob>(rand()/(double)RAND_MAX))
        		{
            		pv[xt].ctag = stag; // a move is made
            		pv[xt].type = pv[xs].type;

					if(stag){									//box update
						if(xtx<pb[stag].x1) pb[stag].x1=xtx;
						if(xtx>pb[stag].x2) pb[stag].x2=xtx;
						if(xty<pb[stag].y1) pb[stag].y1=xty;
						if(xty>pb[stag].y2) pb[stag].y2=xty;
						if(xtz<pb[stag].z1) pb[stag].z1=xtz;
						if(xtz>pb[stag].z2) pb[stag].z2=xtz;
					}

            		if(pv[xs].contact){		//contact moves
						if(pv[xt].contact)
							attached[ttag]--;
            			pv[xt].contact=1;
            			pv[xs].contact=0;
            		}else
	            		if(pv[xt].contact){
	            			pv[xt].contact = 0;
	            			attached[ttag]--;
	            		}

            		if(pv[xt].bond!=0){							//bonds break
            			/*printf("\n%d-%d breaks",xt,pv[xt].bond);*/
            			pv[pv[xt].bond].bond = 0;
            			pv[xt].bond = 0;
            		}

            		if(pv[xs].bond!=0){
            			/*printf("\n%d-%d breaks",xs,pv[xs].bond);*/
            			pv[pv[xs].bond].bond = 0;
            			pv[xs].bond = 0;
            		}

            		if((!CONT_INHIB || ttag == 0) && stag && pv[xt].contact==0 && attached[stag]<MAX_FOCALS_T(pv[xs].type)){    //if there was no contact that moved and there are still less than MAX contacts -> create a new one
            			//ttag == 0 adds sort of contact inhibition, because cell can not make new adhesion sites with the substrate if it is facing another cell there
            			// so if CONT_INHIB == 1, we apply this extra contact inhibition (ttag == 0). If CONT_INHIB == 0 we skip it.
						pv[xt].contact = 1;
						attached[stag]++;
					}
            		if(ttag) {csize[ttag-1]--;}
            		if(stag) {csize[stag-1]++;}
            		accept++;
				}else{
					reject++;
				}
				

			}
			/*else
				i--;*/
		}
		
	}
	

	return ((double) accept / (double) (reject + accept));
}

////////////////////////////////////////////////////////////////////////////////
double CH_moves(VOX* pv, CM* CMs, double k, double VOXSIZE, int NVX, int NVY, int NVZ)
// cellular potts model: one Monte Carlo step
{
	int i,j,NRsteps = NV;
	int xs, xt; // source and target pixel
	int xtx,xty,xtz,xsx,xsy,xsz; // x and y position of target pixel
	int nbs[26],pick, reject, accept; // neighbors of target pixel
	double dH, prob, phi;

	accept=0; reject=0;
	
	for(i=0;i<NRsteps;i++)
	{
		//xt = (rand()*NV/RAND_MAX); // pick random element
		xt = mt_random()%NV; // pick random element
		xty = (xt%(NVX*NVY))/NVX; xtx = (xt%(NVX*NVY))%NVX; xtz=xt/(NVX*NVY);

		if((xtx>0)&&(xtx<NVX-1)&&(xty>0)&&(xty<NVY-1)&&(xtz>0)&&(xtz<NVZ-1)) // exclude outer rim
		{
			nbs[8]=xt-1+NVX+NVX*NVY; nbs[9]=xt+NVX+NVX*NVY; nbs[10]=xt+1+NVX+NVX*NVY;
			nbs[11]=xt-1+NVX*NVY;	 nbs[12]=xt+NVX*NVY;    nbs[13]=xt+1+NVX*NVY;
			nbs[14]=xt-1-NVX+NVX*NVY;nbs[15]=xt-NVX+NVX*NVY;nbs[16]=xt+1-NVX+NVX*NVY;

			nbs[0]=xt-1+NVX; nbs[1]=xt+NVX; nbs[2]=xt+1+NVX;
			nbs[7]=xt-1;                    nbs[3]=xt+1;
			nbs[6]=xt-1-NVX; nbs[5]=xt-NVX; nbs[4]=xt+1-NVX;

			nbs[17]=xt-1+NVX-NVX*NVY; nbs[18]=xt+NVX-NVX*NVY; nbs[19]=xt+1+NVX-NVX*NVY;
			nbs[20]=xt-1-NVX*NVY;	 nbs[21]=xt-NVX*NVY;    nbs[22]=xt+1-NVX*NVY;
			nbs[23]=xt-1-NVX-NVX*NVY;nbs[24]=xt-NVX-NVX*NVY;nbs[25]=xt+1-NVX-NVX*NVY;
			pick = mt_random()%26;
			xs = nbs[pick]; // pick random neighbor
			xsy = (xs%(NVX*NVY))/NVX; xsx = (xs%(NVX*NVY))%NVX; xsz = xs/(NVX*NVY);
			
			
			if(pv[xt].contact!=pv[xs].contact && pv[xs].contact!=0)
			{
				
        		dH = calcdH_CH(pv, CMs, xt, xs, NVX, NVY);
        		prob = exp(-k*IMMOTILITY_CH*dH);
        		if (prob>(rand()/(double)RAND_MAX))
				{
            		pv[xt].contact = pv[xs].contact; // a move is made
            		accept++;
				}else{
					reject++;
				}
			}
			
		}
	}
	
	return ((double) accept / (double) (reject + accept));
}

////////////////////////////////////////////////////////////////////////////////
BOOL splitcheckCCR(VOX* pv, short * CCAlabels, BOX * pb, int* csize, int xt, int ttag, int NVX, int NVY, int NVZ)
{
	
	BOOL split;
	int nbs[26],n,nb,prev,curr,in;
	int v, nrblue, nrgrey, startnb;
	int greys[csize[ttag-1]];
	int i, nrgrey0, g, nbsg[26];
	int vx,vy,vz;
	
	nbs[8]=xt-1+NVX+NVX*NVY; nbs[9]=xt+NVX+NVX*NVY; nbs[10]=xt+1+NVX+NVX*NVY;
	nbs[11]=xt-1+NVX*NVY;	 nbs[12]=xt+NVX*NVY;    nbs[13]=xt+1+NVX*NVY;
	nbs[14]=xt-1-NVX+NVX*NVY;nbs[15]=xt-NVX+NVX*NVY;nbs[16]=xt+1-NVX+NVX*NVY;

	nbs[0]=xt-1+NVX; nbs[1]=xt+NVX; nbs[2]=xt+1+NVX;
	nbs[7]=xt-1;                    nbs[3]=xt+1;
	nbs[6]=xt-1-NVX; nbs[5]=xt-NVX; nbs[4]=xt+1-NVX;

	nbs[17]=xt-1+NVX-NVX*NVY; nbs[18]=xt+NVX-NVX*NVY; nbs[19]=xt+1+NVX-NVX*NVY;
	nbs[20]=xt-1-NVX*NVY;	 nbs[21]=xt-NVX*NVY;    nbs[22]=xt+1-NVX*NVY;
	nbs[23]=xt-1-NVX-NVX*NVY;nbs[24]=xt-NVX-NVX*NVY;nbs[25]=xt+1-NVX-NVX*NVY;
	
	prev = pv[nbs[25]].ctag; in = 0;
	for(n=0;n<26;n++)
	{
		curr = pv[nbs[n]].ctag;
		if((prev!=ttag)&&(curr==ttag))
			in++;
		prev = curr;
	}
	
	split = FALSE;
	
	if(in>1)
	{
		// CONNECTED COMPONENT ALGORITHM Rene-style (CCR)
    	// connected checking "label":
    	// 0: blue;    neighbors of retracted element
    	// 1: white;   undiscovered
    	// 2: grey;    discovered but not finished processing
    	// 3: black;   finished processing
		
		for(vx=pb[ttag].x1-1; vx<=pb[ttag].x2+1; vx++) {
			
      		for (vy=pb[ttag].y1-1; vy<=pb[ttag].y2+1; vy++) {
				
				for (vz=pb[ttag].z1-1; vz<=pb[ttag].z2+1; vz++) {
					
					v = vx + vy * NVX + vz * NVX * NVY;
					CCAlabels[v] = 1;
					
				}
      		}
      	}
		
		CCAlabels[xt] = 3;
		
		nrblue = -1;
		for(n=0;n<26;n++)
		{
			nb = nbs[n];
			if(pv[nb].ctag==ttag)
			{
				CCAlabels[nb]=0; nrblue++;
				startnb = nb;
			}
		}
		CCAlabels[startnb]=2; nrgrey=1; greys[0]=startnb;

		while(nrgrey&&nrblue)
		{
			nrgrey0 = nrgrey;
			// make neighbors of discovered greys grey
			for(i=0;i<nrgrey0;i++)
			{
				g = greys[i];
				nbsg[8]=g-1+NVX+NVX*NVY; nbsg[9]=g+NVX+NVX*NVY; nbsg[10]=g+1+NVX+NVX*NVY;
				nbsg[11]=g-1+NVX*NVY;	 nbsg[12]=g+NVX*NVY;    nbsg[13]=g+1+NVX*NVY;
				nbsg[14]=g-1-NVX+NVX*NVY;nbsg[15]=g-NVX+NVX*NVY;nbsg[16]=g+1-NVX+NVX*NVY;

				nbsg[0]=g-1+NVX; nbsg[1]=g+NVX; nbsg[2]=g+1+NVX;
				nbsg[7]=g-1;                    nbsg[3]=g+1;
				nbsg[6]=g-1-NVX; nbsg[5]=g-NVX; nbsg[4]=g+1-NVX;

				nbsg[17]=g-1+NVX-NVX*NVY; nbsg[18]=g+NVX-NVX*NVY; nbsg[19]=g+1+NVX-NVX*NVY;
				nbsg[20]=g-1-NVX*NVY;	 nbsg[21]=g-NVX*NVY;    nbsg[22]=g+1-NVX*NVY;
				nbsg[23]=g-1-NVX-NVX*NVY;nbsg[24]=g-NVX-NVX*NVY;nbsg[25]=g+1-NVX-NVX*NVY;
				for(n=0;n<26;n++)
				{
					nb = nbsg[n];
					if((pv[nb].ctag==ttag)&&(CCAlabels[nb]<2))
					{
						if(CCAlabels[nb]==0) {nrblue--;}
						CCAlabels[nb]=2; nrgrey++; greys[nrgrey-1]=nb;
					}
				}
			}

			// make processed greys black
			for(i=0;i<nrgrey0;i++)
			{
				g = greys[i];
				CCAlabels[g]=3;
				greys[i]=greys[nrgrey-1]; nrgrey--;
			}

		}
		if(nrblue) {split = TRUE;}

	}

	return split;
}
