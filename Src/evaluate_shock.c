#include "pluto.h"

void updateShockFront(Data* d, Grid* grid){
    int i,j,k;
    FlagShock(d, grid);
    TOT_LOOP(k,j,i){
        d->shockWidth[k][j][i] = grid->dx[0][i];
        d->velocityJump[k][j][i] = 0.0;
        d->upstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
        d->downstreamDensity[k][j][i] = d->Vc[RHO][k][j][i];
    }

    DOM_LOOP(k,j,i){
        if(!(d->flag[k][j][i] & FLAG_ENTROPY)){

            int upstreami = i;
            int upstreamj = j;
            int upstreamk = k;

            double x = grid->x[0][i];
            double y = grid->x[1][j];
            double z = grid->x[2][k];

            while(!(d->flag[upstreamk][upstreamj][upstreami] & FLAG_ENTROPY)){
                double pgradx = 0;
                double pgrady = 0;
                double pgradz = 0;

#if INCLUDE_IDIR
                pgradx = grid->dx_dl[IDIR][j][i]*(d->Vc[PRS][k][j][i+1] - d->Vc[PRS][k][j][i-1])/(grid->x[0][i+1] - grid->x[0][i-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][j][i]*(d->Vc[PRS][k][j+1][i] - d->Vc[PRS][k][j-1][i])/(grid->x[1][j+1] - grid->x[1][j-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][j][i]*(d->Vc[PRS][k+1][j][i] - d->Vc[PRS][k-1][j][i])/(grid->x[2][k+1] - grid->x[2][k-1]);
#endif
                traceNextCell(grid, &x, &y, &z, -pgradx, -pgrady, -pgradz, &upstreami, &upstreamj, &upstreamk);
            }

            int downstreami = i;
            int downstreamj = j;
            int downstreamk = k;

            x = grid->x[0][i];
            y = grid->x[1][j];
            z = grid->x[2][k];

            while(!(d->flag[downstreamk][downstreamj][downstreami] & FLAG_ENTROPY)){
                double pgradx = 0;
                double pgrady = 0;
                double pgradz = 0;

#if INCLUDE_IDIR
                pgradx = grid->dx_dl[IDIR][j][i]*(d->Vc[PRS][k][j][i+1] - d->Vc[PRS][k][j][i-1])/(grid->x[0][i+1] - grid->x[0][i-1]);
#endif
#if INCLUDE_JDIR
                pgrady = grid->dx_dl[JDIR][j][i]*(d->Vc[PRS][k][j+1][i] - d->Vc[PRS][k][j-1][i])/(grid->x[1][j+1] - grid->x[1][j-1]);
#endif
#if INCLUDE_KDIR
                pgradz = grid->dx_dl[KDIR][j][i]*(d->Vc[PRS][k+1][j][i] - d->Vc[PRS][k-1][j][i])/(grid->x[2][k+1] - grid->x[2][k-1]);
#endif
                traceNextCell(grid, &x, &y, &z, pgradx, pgrady, pgradz, &downstreami, &downstreamj, &downstreamk);
            }

            double Vd1 = d->Vc[VX1][downstreamk][downstreamj][downstreami];
            double Vd2 = d->Vc[VX2][downstreamk][downstreamj][downstreami];
            double Vd3 = d->Vc[VX3][downstreamk][downstreamj][downstreami];

            double Vu1 = d->Vc[VX1][upstreamk][upstreamj][upstreami];
            double Vu2 = d->Vc[VX2][upstreamk][upstreamj][upstreami];
            double Vu3 = d->Vc[VX3][upstreamk][upstreamj][upstreami];

            double rhod = d->Vc[RHO][downstreamk][downstreamj][downstreami];
            double rhou = d->Vc[RHO][upstreamk][upstreamj][upstreami];

            double xd, yd, zd;
            double xu, yu, zu;
            double Vdx, Vdy, Vdz;
            double Vux, Vuy, Vuz;

            double width = 1E100;

#if GEOMETRY == CARTESIAN
            xd = grid->x[0][downstreami];
            yd = grid->x[1][downstreamj];
            zd = grid->x[2][downstreamk];

            xu = grid->x[0][upstreami];
            yu = grid->x[1][upstreamj];
            zu = grid->x[2][upstreamk];

            Vdx = Vd1;
            Vdy = Vd2;
            Vdz = Vd3;

            Vux = Vu1;
            Vuy = Vu2;
            Vuz = Vu3;
#elif GEOMETRY == CYLINDRICAL
            xd = grid->x[0][downstreami]*cos(grid->x[2][downstreamk]);
            yd = grid->x[0][downstreami]*sin(grid->x[2][downstreamk]);
            zd = grid->x[1][downstreamj];

            xu = grid->x[0][upstreami]*cos(grid->x[2][upstreamk]);
            yu = grid->x[0][upstreami]*sin(grid->x[2][upstreamk]);
            zu = grid->x[1][upstreamj];
#elif GEOMETRY == POLAR
            xd = grid->x[0][downstreami]*cos(grid->x[1][downstreamj]);
            yd = grid->x[0][downstreami]*sin(grid->x[1][downstreamj]);
            zd = grid->x[2][downstreamk];

            xu = grid->x[0][upstreami]*cos(grid->x[1][upstreamj]);
            yu = grid->x[0][upstreami]*sin(grid->x[1][upstreamj]);
            zu = grid->x[2][upstreamk];
#elif GEOMETRY == SPHERICAL
            xd = grid->x[0][downstreami]*sin(grid->x[1][downstreamj])*cos(grid->x[2][downstreamk]);
            yd = grid->x[0][downstreami]*sin(grid->x[1][downstreamj])*sin(grid->x[2][downstreamk]);
            zd = grid->x[0][downstreami]*cos(grid->x[1][downstreamj]);

            xu = grid->x[0][upstreami]*sin(grid->x[1][upstreamj])*cos(grid->x[2][upstreamk]);
            yu = grid->x[0][upstreami]*sin(grid->x[1][upstreamj])*sin(grid->x[2][upstreamk]);
            zu = grid->x[0][upstreami]*cos(grid->x[1][upstreamj]);
#else
#endif
            width = sqrt((xd-xu)*(xd-xu) + (yd-yu)*(yd-yu) + (zd-zu)*(zd-zu));
            double V = sqrt((Vdx - Vux)*(Vdx - Vux) + (Vdy - Vuy)*(Vdy - Vuy) + (Vdz - Vuz)*(Vdz - Vuz));

            d->shockWidth[k][j][i] = width;
            d->velocityJump[k][j][i] = V;
            d->upstreamDensity[k][j][i] = rhou;
            d->downstreamDensity[k][j][i] = rhod;
        }
    }
}

void traceNextCell(Grid* grid, double* x1, double* x2, double* x3, double vx, double vy, double vz, int* i, int* j, int* k){
//todo proper stright lines for other geometries
#if INCLUDE_KDIR
    double lx = grid->xl[0][*i];
    double rx = grid->xr[0][*i];
    double ly = grid->xl[1][*j];
    double ry = grid->xr[1][*j];
    double lz = grid->xl[2][*k];
    double rz = grid->xr[2][*k];

    double dx;
    double dy;
    double dz;
    if(vx > 0){
        dx = (rx - *x1)/grid->dx_dl[IDIR][*j][*i];
    } else {
        dx = (*x1 - lx)/grid->dx_dl[IDIR][*j][*i];
    }
    if(vy > 0){
        dy = (ry - *x2)/grid->dx_dl[JDIR][*j][*i];
    } else {
        dy = (*x2 - ly)/grid->dx_dl[JDIR][*j][*i];
    }
    if(vz > 0){
        dz = (rz - *x3)/grid->dx_dl[KDIR][*j][*i];
    } else {
        dz = (*x3 - lz)/grid->dx_dl[KDIR][*j][*i];
    }

    if(fabs(dx*vy) > fabs(dy*vx)){
        if(fabs(dz*vy) > fabs(dy*vz)){
            double dt = fabs(dy/vy);
            *x1 = *x1 + dt*vx;
            *x3 = *x3 + dt*vz;
            if(vy > 0){
                *x2 = ry;
                *j = (*j)+1;
            } else {
                *x2 = ly;
                *j = (*j)-1;
            }
        } else {
            double dt = fabs(dz/vz);
            *x1 = *x1 + dt*vx;
            *x2 = *x2 + dt*vy;
            if(vz > 0){
                *x3 = rz;
                *k = (*k)+1;
            } else {
                *x3 = lz;
                *k = (*k)-1;
            }
        }
    } else {
        if(fabs(dz*vx) > fabs(dx*vz)){
            double dt = fabs(dx/vx);
            *x2 = *x2 + dt*vy;
            *x3 = *x3 + dt*vz;
            if(vx > 0){
                *x1 = rx;
                *i = (*i)+1;
            } else {
                *x1 = lx;
                *i = (*i)-1;
            }
        } else {
            double dt = fabs(dz/vz);
            *x1 = *x1 + dt*vx;
            *x2 = *x2 + dt*vy;
            if(vz > 0){
                *x3 = rz;
                *k = (*k)+1;
            } else {
                *x3 = lz;
                *k = (*k)-1;
            }
        }
    }
#elif INCLUDE_JDIR
    int a = *i;
    double lx = grid->xl[0][*i];
    double rx = grid->xr[0][*i];
    double ly = grid->xl[1][*j];
    double ry = grid->xr[1][*j];

    double dx;
    double dy;
    if(vx > 0){
        dx = (rx - *x1)/grid->dx_dl[IDIR][*j][*i];
    } else {
        dx = (*x1 - lx)/grid->dx_dl[IDIR][*j][*i];
    }
    if(vy > 0){
        dy = (ry - *x2)/grid->dx_dl[JDIR][*j][*i];
    } else {
        dy = (*x2 - ly)/grid->dx_dl[JDIR][*j][*i];
    }

    if(fabs(dx*vy) > fabs(dy*vx)){
        double dt = fabs(dy/vy);
        *x1 = *x1 + dt*vx;
        if(vy > 0){
            *x2 = ry;
            *j = (*j)+1;
        } else {
            *x2 = ly;
            *j = (*j)-1;
        }
    } else {
        double dt = fabs(dx/vx);
        *x2 = *x2 + dt*vy;
        if(vx > 0){
            *x1 = rx;
            *i = (*i)+1;
        } else {
            *x1 = lx;
            *i = (*i)-1;
        }
    }
#else
    if(v1 > 0){
        *x1 = grid->xr[0][*i];
        *i = (*i) + 1;
    } else if (v1 < 0){
        *x1 = grid->xl[0][*i];
        *i = (*i) - 1;
    } else {
        printLog("vx = 0 in trace cell\n");
        exit(0);
    }
#endif
}
