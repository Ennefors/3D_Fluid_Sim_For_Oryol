//
// Created by bill on 18/01/18.
//

#ifndef ORYOL_MACGRID_H
#define ORYOL_MACGRID_H

#pragma once

#include "glm/mat4x4.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "random"
#include "Support.h"
#include <cuda.h>
#include <time.h>
#include <map>
#include <math.h>

#define PARTICLE_MASS 0.02
#define PARTICLE_RADIUS M_PI/17
#define CELL_WIDTH 0.5
#define MAGIC 1942735123




    const unsigned char cell0 = 0x01;
    const unsigned char cell1 = 0x02;
    const unsigned char cell2 = 0x04;
    const unsigned char cell3 = 0x08;
    const unsigned char cell4 = 0x10;
    const unsigned char cell5 = 0x20;
    const unsigned char cell6 = 0x40;


enum types
{
    fluid = 0,
    air = 1,
    solid = 2,
};

struct gridCell
{
    glm::vec3 velocity = glm::vec3(0.f);
    glm::vec3 convection = glm::vec3(0.f);
    glm::vec3 position = glm::vec3(0.f);
    glm::vec3 viscosity = glm::vec3(0.f);
    gridCell* cells[7] = {nullptr};


    int cNeighbours[6] = {MAGIC};
    gridCell* cNeigh[6] = {nullptr};

    unsigned char neighbours;

    double pressure;
    double width;
    int layer;

    int spot = 0;

    int type;
    gridCell()
    {

    }
    gridCell(glm::vec3 loc, int t, int lay, int id)
    {
        pressure = 0.f;
        velocity = glm::vec3(0.f);
        convection = glm::vec3(0.f);
        viscosity = glm::vec3(0.f);

        position = loc;
        layer = lay;
        type = t;
    }


    inline gridCell& operator=(const gridCell c)
    {
        this->velocity = c.velocity;
        this->position = c.position;
        for (int i = 0; i < 6; ++i)
        {
            this->cNeighbours[i] = c.cNeighbours[i];
        }

        this->neighbours = c.neighbours;

        this->pressure = c.pressure;
        this->width = c.width;
        this->layer = c.layer;

        this->spot = c.spot;

        this->type = c.type;

        return *this;
    }
};

struct particle
{
    glm::vec3 pos;
    particle(int i) {
        pos = glm::vec3(0);
    }
};

class MacGrid
{

public:
    MacGrid(const glm::vec3 ori, const glm::vec3 dim);

    ~MacGrid(){};

    float density;
    float cellWidth;
    float viscousity = 0.3;

    float maxvel = 15.f;
    float uMaxx = 15.f;
    float uMaxy = 15.f;
    float uMaxz = 15.f;
    float specConstMin = 2000;
    float specConstMax = 3;

    time_t  timev;

    bool firstTime = true;

    //bounds = max
    //-bounds = min, min.z = 0;
    glm::vec3 bounds;
    glm::vec3 origin;


    inline glm::vec3 getMinB() {glm::vec3 min = glm::vec3(-CELL_WIDTH, -CELL_WIDTH, -CELL_WIDTH); return min;};

    std::vector<particle> particles;
    std::map<float, gridCell> cells;
    std::map<float, glm::vec3> emptyCells;

    float timeStep();

    inline glm::vec3 cellOwner(glm::vec3 pPos)
    {

        float cx = floorf(pPos.x+CELL_WIDTH);
        float cy = floorf(pPos.y+CELL_WIDTH);
        float cz = floorf(pPos.z+CELL_WIDTH);
        return glm::vec3(cx, cy, cz);
    };
    inline glm::vec3 getUx(glm::vec3 v)
    {

        glm::vec3 Ux(v.x-1, v.y, v.z);
        return Ux;
    };
    inline glm::vec3 getUy(glm::vec3 v)
    {

        glm::vec3 Uy(v.x, v.y-1, v.z);
        return Uy;
    };
    inline glm::vec3 getUz(glm::vec3 v)
    {

        glm::vec3 Uz(v.x, v.y, v.z-1);
        return Uz;
    };

    ////Randomly spawns a particle within the given bounds and adds it to the grid.
    inline void addParticle(glm::vec3 max, glm::vec3 min)
    {
        double x = (rand()%(int)max.x) + min.x;
        double y = (rand()%(int)max.y) + min.y;
        double z = (rand()%(int)max.z) + min.z;
        glm::vec3 pPos(x, y, z);
        particle newP(particles.size());
        newP.pos = pPos;


        float rx = (float)rand()/(float)RAND_MAX;
        float ry = (float)rand()/(float)RAND_MAX;
        float rz = (float)rand()/(float)RAND_MAX;

        newP.pos.x += rx;
        newP.pos.y += ry;
        newP.pos.z += rz;

        particles.push_back(newP);
        return;
    };

    //The hashing function for cells
    inline float hashFunc(glm::vec3 v){return (v.x*541 + 79*v.y + 31*v.z);};

    int update(float dt);
    int fluidUpdate(float dt);

    //Particle backwards tracing algorithms
    glm::vec3 traceParticle(double x, double y, double z, float dt);
    glm::vec3 getVelocity(double x, double y, double z);
    glm::vec3 interpolateVel(double x, double y, double z);
    float interPolVal(double x, double y, double z, int ind);
    float interPolVal2(double x, double y, double z, int ind);


    inline float calc_laplace2(glm::vec3 v, int ind, gridCell c)
    {

        glm::vec3 vf1;
        glm::vec3 vf2;
        glm::vec3 vf3;

        glm::vec3 vb1;
        glm::vec3 vb2;
        glm::vec3 vb3;

        auto cfr1 = cells.find(hashFunc(glm::vec3(v.x+1, v.y, v.z)));
        auto cfr2 = cells.find(hashFunc(glm::vec3(v.x, v.y+1, v.z)));
        auto cfr3 = cells.find(hashFunc(glm::vec3(v.x, v.y, v.z+1)));

        auto cbr1 = cells.find(hashFunc(glm::vec3(v.x-1, v.y, v.z)));
        auto cbr2 = cells.find(hashFunc(glm::vec3(v.x, v.y-1, v.z)));
        auto cbr3 = cells.find(hashFunc(glm::vec3(v.x, v.y, v.z-1)));

        auto cellsEnd = cells.end();



        gridCell cf1 = cfr1->second;
        gridCell cf2 = cfr2->second;
        gridCell cf3 = cfr3->second;

        gridCell cb1 = cbr1->second;
        gridCell cb2 = cbr2->second;
        gridCell cb3 = cbr3->second;

        int cf1reltype = -1;
        int cf2reltype = -1;
        int cf3reltype = -1;

        int cb1reltype = -1;
        int cb2reltype = -1;
        int cb3reltype = -1;

        int creltype = -1;

        if (ind == 0)
        {
            if(cfr1 != cellsEnd)
            {
                if(cf1.cNeighbours[1] != MAGIC)
                    cf1reltype = cf1.cNeigh[1]->type;
            }
            if(cbr1 != cellsEnd)
            {
                if(cb1.cNeighbours[1] != MAGIC)
                    cb1reltype = cb1.cNeigh[1]->type;
            }

            if(cfr2 != cellsEnd)
            {
                if(cf2.cNeighbours[1] != MAGIC)
                    cf2reltype = cf2.cNeigh[1]->type;
            }
            if(cbr2 != cellsEnd)
            {
                if(cb2.cNeighbours[1] != MAGIC)
                    cb2reltype = cb2.cNeigh[1]->type;
            }

            if(cfr3 != cellsEnd)
            {
                if(cf3.cNeighbours[1] != MAGIC)
                    cf3reltype = cf3.cNeigh[1]->type;
            }
            if(cbr3 != cellsEnd)
            {
                if(cb3.cNeighbours[1] != MAGIC)
                    cb3reltype = cb3.cNeigh[1]->type;
            }

            if(c.cNeighbours[1] != MAGIC)
                creltype = c.cNeigh[1]->type;
        }
        else if(ind == 1)
        {
            if(cfr1 != cellsEnd)
            {
                if(cf1.cNeighbours[3] != MAGIC)
                    cf1reltype = cf1.cNeigh[3]->type;
            }
            if(cbr1 != cellsEnd)
            {
                if(cb1.cNeighbours[3] != MAGIC)
                    cb1reltype = cb1.cNeigh[3]->type;
            }

            if(cfr2 != cellsEnd)
            {
                if(cf2.cNeighbours[3] != MAGIC)
                    cf2reltype = cf2.cNeigh[3]->type;
            }
            if(cbr2 != cellsEnd)
            {
                if(cb2.cNeighbours[3] != MAGIC)
                    cb2reltype = cb2.cNeigh[3]->type;
            }

            if(cfr3 != cellsEnd)
            {
                if(cf3.cNeighbours[3] != MAGIC)
                    cf3reltype = cf3.cNeigh[3]->type;
            }
            if(cbr3 != cellsEnd)
            {
                if(cb3.cNeighbours[3] != MAGIC)
                    cb3reltype = cb3.cNeigh[3]->type;
            }

            if(c.cNeighbours[3] != MAGIC)
                creltype = c.cNeigh[3]->type;
        }
        else if(ind == 2)
        {
            if(cfr1 != cellsEnd)
            {
                if(cf1.cNeighbours[5] != MAGIC)
                    cf1reltype = cf1.cNeigh[5]->type;
            }
            if(cbr1 != cellsEnd)
            {
                if(cb1.cNeighbours[5] != MAGIC)
                    cb1reltype = cb1.cNeigh[5]->type;
            }

            if(cfr2 != cellsEnd)
            {
                if(cf2.cNeighbours[5] != MAGIC)
                    cf2reltype = cf2.cNeigh[5]->type;
            }
            if(cbr2 != cellsEnd)
            {
                if(cb2.cNeighbours[5] != MAGIC)
                    cb2reltype = cb2.cNeigh[5]->type;
            }

            if(cfr3 != cellsEnd)
            {
                if(cf3.cNeighbours[5] != MAGIC)
                    cf3reltype = cf3.cNeigh[5]->type;
            }
            if(cbr3 != cellsEnd)
            {
                if(cb3.cNeighbours[5] != MAGIC)
                    cb3reltype = cb3.cNeigh[5]->type;
            }

            if(c.cNeighbours[5] != MAGIC)
                creltype = c.cNeigh[5]->type;

        }

        float  ret = 0;
        float mult_nr = 0;
        if((cf1.type == fluid && cfr1 != cellsEnd) || cf1reltype == fluid)
        {
            ret += cf1.velocity[ind];
//            ret += getVelocity(cf1.position.x, cf1.position.y, cf1.position.z)[ind];
            mult_nr++;
        }
        if((cb1.type == fluid && cbr1 != cellsEnd) || cb1reltype == fluid)
        {
            ret += cb1.velocity[ind];
//            ret += getVelocity(cb1.position.x, cb1.position.y, cb1.position.z)[ind];
            mult_nr++;
        }
        if((cf2.type == fluid && cfr2 != cellsEnd) || cf2reltype == fluid)
        {
            ret += cf2.velocity[ind];
//            ret += getVelocity(cf2.position.x, cf2.position.y, cf2.position.z)[ind];
            mult_nr++;
        }
        if((cb2.type == fluid && cbr2 != cellsEnd) || cb2reltype == fluid)
        {
            ret += cb2.velocity[ind];
//            ret += getVelocity(cb2.position.x, cb2.position.y, cb2.position.z)[ind];
            mult_nr++;
        }
        if((cf3.type == fluid && cfr3 != cellsEnd) || cf3reltype == fluid)
        {
            ret += cf3.velocity[ind];
//            ret += getVelocity(cf3.position.x, cf3.position.y, cf3.position.z)[ind];
            mult_nr++;
        }
        if((cb3.type == fluid && cbr3 != cellsEnd) || cb3reltype == fluid)
        {
            ret += cb3.velocity[ind];
//            ret += getVelocity(cb3.position.x, cb3.position.y, cb3.position.z)[ind];
            mult_nr++;
        }

        ret -= c.velocity[ind]*mult_nr;
//            ret -= getVelocity(c.position.x, c.position.y, c.position.z)[ind]*6.f;

        if (ret>1000)
        {
            printf("Wat");
        }
        return ret;
    }

    inline float calc_diver(glm::vec3 v, gridCell c)
    {

        glm::vec3 vf1(0.f);
        glm::vec3 vf2(0.f);
        glm::vec3 vf3(0.f);

        glm::vec3 v1;

        gridCell cx = cells.find(hashFunc(glm::vec3(v.x+1, v.y, v.z)))->second;
        gridCell cy = cells.find(hashFunc(glm::vec3(v.x, v.y+1, v.z)))->second;
        gridCell cz = cells.find(hashFunc(glm::vec3(v.x, v.y, v.z+1)))->second;

//        vf1 = getVelocity(v.x+1, v.y, v.z);
//        vf2 = getVelocity(v.x, v.y+1, v.z);
//        vf3 = getVelocity(v.x, v.y, v.z+1);

        vf1 = cx.velocity;
        vf2 = cy.velocity;
        vf3 = cz.velocity;

        v1 = c.velocity;
//        v1 = getVelocity(v.x, v.y, v.z);

        float vx = 0;
        float vy = 0;
        float vz = 0;

        if(cx.type != solid && c.type != solid)
        {
            vx = vf1[0] - v1[0];
        }else if(cx.type == solid && c.type == fluid){
//            vx = 0;
            vx = -v1[0];
        }else if(cx.type == fluid && c.type == solid){
//            vx = 0;
            vx = vf1[0];
        }
        if(cy.type != solid && c.type != solid)
        {
            vy = vf2[1] - v1[1];
        }else if(cy.type == solid && c.type == fluid){
//            vy = 0;
            vy = -v1[1];
        }else if(cy.type == fluid && c.type == solid){
//            vy = 0;
            vy = vf2[1];
        }
        if(cz.type != solid && c.type != solid)
        {
            vz = vf3[2] - v1[2];
        }else if(cz.type == solid && c.type == fluid){
//            vz = 0;
            vz = -v1[2];
        }else if(cz.type == fluid && c.type == solid){
//            vz = 0;
            vz = vf3[2];
        }

        return (vx) + (vy) + (vz);
    }

    inline glm::vec3 calc_grad(glm::vec3 v, gridCell c)
    {

        gridCell vb1;
        gridCell vb2;
        gridCell vb3;

        vb1 = cells.find(hashFunc(glm::vec3(v.x-1, v.y, v.z)))->second;
        vb2 = cells.find(hashFunc(glm::vec3(v.x, v.y-1, v.z)))->second;
        vb3 = cells.find(hashFunc(glm::vec3(v.x, v.y, v.z-1)))->second;


        float vx = 0;
        float vy = 0;
        float vz = 0;

        if(vb1.type != 9)
        {
            vx = c.pressure - vb1.pressure;
        } else {
            vx = -c.velocity.x;
        }
        if(vb2.type != 9)
        {
            vy = c.pressure - vb2.pressure;
        } else {
            vy = -c.velocity.y;
        }
        if(vb3.type != 9)
        {
            vz = c.pressure - vb3.pressure;
        } else {
            vz = -c.velocity.z;
        }

        return glm::vec3((vx), (vy), (vz));
    }

    void NavStokes(float dt);
    ////*Functions used by NavierStokes*////

    glm::vec3 applyConvection(gridCell c, float dt);
    glm::vec3 applyViscousity(gridCell c, float dt);
    void applyPressure(gridCell** c, int nrf,float dt);

    void calcDiffuse();
    void calcProj();
    void calcAdv();

    void ChangeGravity(int input);
    glm::vec3 getGravity();
};





#endif //ORYOL_MACGRID_H
