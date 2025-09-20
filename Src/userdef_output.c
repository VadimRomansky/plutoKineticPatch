#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;  
  
  double ***T;
  double ***Shock;
  double ***Bturb;
  double ***Jmc;
  double ***Fkin;
  double ***Pkin;
  double ***rho;
  double ***prs;
  double ***Vjump;
  //double ***Vdx;
  //double ***Vux;
  //double ***Vdy;
  //double ***Vuy;
  double mu;
  double us[256];
  int nv;

  double UNIT_MASS = UNIT_DENSITY*pow(UNIT_LENGTH, 3.0);
  double UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY;
  double UNIT_ENERGY = UNIT_MASS*UNIT_VELOCITY*UNIT_VELOCITY;
  double UNIT_MFIELD = (UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY));

  T = GetUserVar("T");
  Shock = GetUserVar("Shock");
#if PARTICLES == PARTICLES_KIN
  Fkin = GetUserVar("Fkin");
  Pkin = GetUserVar("Pkin");
#endif
  Vjump = GetUserVar("Vjump");
  //Vdx = GetUserVar("Vdx");
  //Vux = GetUserVar("Vux");
  //Vdy = GetUserVar("Vdy");
  //Vuy = GetUserVar("Vuy");

#if TURBULENT_FIELD == YES
  Bturb = GetUserVar("Bturb");
  Jmc = GetUserVar("Jmc");
#endif

  rho = d->Vc[RHO];
  prs = d->Vc[PRS];


  DOM_LOOP(k,j,i){
      NVAR_LOOP(nv) us[nv]=d->Vc[nv][k][j][i];
      mu = MeanMolecularWeight(us)*CONST_mp;
      T[k][j][i] = (prs[k][j][i]*mu/rho[k][j][i])*UNIT_VELOCITY*UNIT_VELOCITY/(1.6E-12);

      Shock[k][j][i] = 0;
      if(!(d->flag[k][j][i] & FLAG_ENTROPY)){
          Shock[k][j][i] = 1.0;
      }
      Vjump[k][j][i] = d->velocityJump[k][j][i];
      //Vdx[k][j][i] = d->downstreamV1[k][j][i];
      //Vux[k][j][i] = d->upstreamV1[k][j][i];
      //Vdy[k][j][i] = d->downstreamV2[k][j][i];
      //Vuy[k][j][i] = d->upstreamV2[k][j][i];
#if TURBULENT_FIELD == YES
      double B2 = 0;
      double Jx = 0;
      double Jy = 0;
      double Jz = 0;
      for(int l = 0; l < NTURB; ++l){
          double dk;
          if(l == 0){
              dk = d->k_turb[1] - d->k_turb[0];
          } else {
              dk = d->k_turb[l] - d->k_turb[l-1];
          }
          B2 += 8*CONST_PI*d->Wt[k][j][i][l]*dk;
      }
      Bturb[k][j][i] = sqrt(B2);

      for(int l = 0; l < NMOMENTUM; ++l){
          double dp;
          if(l == 0){
              dp = d->p_grid[1] - d->p_grid[0];
          } else {
              dp = d->p_grid[l] - d->p_grid[l-1];
          }
          Jx += d->Jkin1[k][j][i][l]*dp;
          Jy += d->Jkin2[k][j][i][l]*dp;
          Jz += d->Jkin3[k][j][i][l]*dp;
      }
      Jmc[k][j][i] = sqrt(Jx*Jx + Jy*Jy + Jz*Jz);
#endif
#if PARTICLES == PARTICLES_KIN
      Pkin[k][j][i] = d->Pkin[k][j][i];
      Fkin[k][j][i] = 0;
      for(int l = 0; l < NMOMENTUM; ++l){
          double dp;
          if(l == 0){
              dp = d->p_grid[1] - d->p_grid[0];
          } else {
              dp = d->p_grid[l] - d->p_grid[l-1];
          }
          Fkin[k][j][i] += d->Fkin[k][j][i][l]*dp/(d->p_grid[l]);
      }
#endif
}
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
    SetOutputVar ("T", VTK_OUTPUT, YES);
    SetOutputVar ("Shock", VTK_OUTPUT, YES);
    SetOutputVar ("Bturb", VTK_OUTPUT, YES);
    SetOutputVar ("Fkin", VTK_OUTPUT, YES);
    SetOutputVar ("Jmc", VTK_OUTPUT, YES);
    SetOutputVar ("Pkin", VTK_OUTPUT, YES);
  Image *image;

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





