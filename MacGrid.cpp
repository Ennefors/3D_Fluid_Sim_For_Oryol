//
// Created by bill on 18/01/18.
//

#include "MacGrid.h"

glm::vec3 GRAVITY = glm::vec3(0.f, -9.82f,0.f );
float VISC = 0.614f;

glm::vec3 MacGrid::getGravity()
{
    return GRAVITY;
}
void MacGrid::ChangeGravity(int input)
{
    switch(input)
    {
        case 0:
            GRAVITY.x += 0.2f;
            break;
        case 1:
            GRAVITY.x -= 0.2f;
            break;
        case 2:
            GRAVITY.y += 0.2f;
            break;
        case 3:
            GRAVITY.y -= 0.2f;
            break;
        case 4:
            GRAVITY.z += 0.2f;
            break;
        case 5:
            GRAVITY.z -= 0.2f;
            break;
    }

    if(GRAVITY.x <= -9.82f)
    {
        GRAVITY.x = -9.82f;
    }
    if(GRAVITY.x > 9.82f)
    {
        GRAVITY.x = 9.82f;
    }
    if(GRAVITY.y <= -9.82f)
    {
        GRAVITY.y = -9.82f;
    }
    if(GRAVITY.y > 9.82f)
    {
        GRAVITY.y = 9.82f;
    }
    if(GRAVITY.z <= -9.82f)
    {
        GRAVITY.z = -9.82f;
    }
    if(GRAVITY.z > 9.82f)
    {
        GRAVITY.z = 9.82f;
    }
};

double rateRK4(double x, double y)
{
    return x * sqrt(y);
}

MacGrid::MacGrid(const glm::vec3 ori, const glm::vec3 dim)
{
    srand(0);
    this->bounds = dim;
    this->bounds.x += CELL_WIDTH;
    this->bounds.y += CELL_WIDTH;
    this->bounds.z += CELL_WIDTH;
    this->origin = ori;
};



int MacGrid::update(float dt)
{
    int cellNR = 0;

    for (auto j = cells.begin(); j != cells.end() ; ++j)
    {
        j->second.layer = -1;
    }

    for (unsigned int i = 0; i < particles.size(); ++i)
    {
        float gi = hashFunc(cellOwner(particles[i].pos));
        auto it = cells.find(gi);
        if(it != cells.end())
        {
            if(it->second.type != solid)
            {
                it->second.layer = 0;
                it->second.type = fluid;
            }
        } else
        {
            if((it->second.position.x <= bounds.x || it->second.position.x >= getMinB().x) ||
               (it->second.position.y <= bounds.y || it->second.position.y >= getMinB().y) ||
               (it->second.position.z <= bounds.z || it->second.position.z >= getMinB().z)  ) {
                gridCell c = gridCell(cellOwner(particles[i].pos), fluid, 0, gi);
                std::pair<int, gridCell> p;
                p.first = gi;
                p.second = c;
                cells.insert(p);
                cellNR++;
            }
        }
    }

    for (int i = 1; i <= 2 ; ++i)
    {

        for (auto j = cells.begin(); j != cells.end(); ++j)
        {
            if ( (j->second.type == fluid && j->second.layer == i - 1) || (j->second.layer == i - 1 && j->second.type == air))
            {

                glm::ivec3 c = j->second.position;
                //get all cell neighbours
                glm::vec3 neigh[6] =
                        {
                                glm::vec3(c.x+1, c.y, c.z),
                                glm::vec3(c.x-1, c.y, c.z),
                                glm::vec3(c.x, c.y+1, c.z),
                                glm::vec3(c.x, c.y-1, c.z),
                                glm::vec3(c.x, c.y, c.z+1),
                                glm::vec3(c.x, c.y, c.z-1),
                        };
                //see if neighbour exsist, if not add them as a air cell;
                for (int k = 0; k < 6; ++k)
                {
                    float gi = hashFunc(neigh[k]);
                    auto it = cells.find(gi);
                    if (it != cells.end())
                    {
                        if (it->second.layer == -1 && it->second.type != solid)
                        {
                            it->second.layer = i;
                            it->second.type = air;
                            it->second.pressure = 1.0f;
                        } else if (it->second.type == solid)
                        {
                            it->second.layer = i;
                        }

                    } else {
                        gridCell nc(gridCell(neigh[k], air, i, gi));
                        std::pair<int, gridCell> p;
                        if((nc.position.x >= bounds.x || nc.position.x <= getMinB().x) ||
                           (nc.position.y >= bounds.y || nc.position.y <= getMinB().y) ||
                           (nc.position.z >= bounds.z || nc.position.z <= getMinB().z)  )
                        {
                            nc.type = solid;
                        } else {
                            nc.type = air;
                            nc.pressure = 0.986;
                        }
                        nc.layer = i;
                        p.first = gi;
                        p.second = nc;
                        cells.insert(p);

                        cellNR++;
                    }
                }
            }
        }
    }
    //remove unused cells
    for (auto j = cells.begin(); j != cells.end();) {


        if (j->second.layer == -1) {
            j = cells.erase(j);
        } else {
            ++j;
        }
    }

    for (auto j = cells.begin(); j != cells.end(); j++) {
        glm::ivec3 c = j->second.position;
        //get all cell neighbours
        glm::vec3 neigh[6] =
                {
                        glm::vec3(c.x + 1, c.y, c.z),
                        glm::vec3(c.x - 1, c.y, c.z),
                        glm::vec3(c.x, c.y + 1, c.z),
                        glm::vec3(c.x, c.y - 1, c.z),
                        glm::vec3(c.x, c.y, c.z + 1),
                        glm::vec3(c.x, c.y, c.z - 1),
                };
        for (int k = 0; k < 6; ++k)
        {
            float gi = hashFunc(neigh[k]);
            auto it = cells.find(gi);

            if (it != cells.end()) {
                j->second.cNeighbours[k] = gi;
                j->second.cNeigh[k] = &it->second;
            } else
            {
                j->second.cNeighbours[k] = MAGIC;
                j->second.cNeigh[k] = nullptr;
            }
        }

        glm::vec3 bNeigh[7] =
            {
                    glm::vec3(c.x+1, c.y, c.z),
                    glm::vec3(c.x, c.y+1, c.z),
                    glm::vec3(c.x+1, c.y+1, c.z),
                    glm::vec3(c.x, c.y, c.z+1),
                    glm::vec3(c.x+1, c.y, c.z+1),
                    glm::vec3(c.x, c.y+1, c.z+1),
                    glm::vec3(c.x+1, c.y+1, c.z+1)

            };

        for (int k = 0; k < 7; ++k)
        {
            float gi = hashFunc(bNeigh[k]);
            auto it = cells.find(gi);

            if (it != cells.end())
            {

                switch (k)
                {
                    case 0:
                        j->second.neighbours |= cell0;
                        j->second.cells[0] = &it->second;
                        break;
                    case 1:
                        j->second.neighbours |= cell1;
                        j->second.cells[1] = &it->second;

                        break;
                    case 2:
                        j->second.neighbours |= cell2;
                        j->second.cells[2] = &it->second;

                        break;
                    case 3:
                        j->second.neighbours |= cell3;
                        j->second.cells[3] = &it->second;

                        break;
                    case 4:
                        j->second.neighbours |= cell4;
                        j->second.cells[4] = &it->second;

                        break;
                    case 5:
                        j->second.neighbours |= cell5;
                        j->second.cells[5] = &it->second;

                        break;
                    case 6:
                        j->second.neighbours |= cell6;
                        j->second.cells[6] = &it->second;

                        break;
                }
            }else
            {
                switch (k)
                {
                    case 0:
                        j->second.neighbours ^= cell0;
                        j->second.cells[0] = nullptr;
                        break;
                    case 1:
                        j->second.neighbours ^= cell1;
                        j->second.cells[1] = nullptr;

                        break;
                    case 2:
                        j->second.neighbours ^= cell2;
                        j->second.cells[2] = nullptr;

                        break;
                    case 3:
                        j->second.neighbours ^= cell3;
                        j->second.cells[3] = nullptr;

                        break;
                    case 4:
                        j->second.neighbours ^= cell4;
                        j->second.cells[4] = nullptr;

                        break;
                    case 5:
                        j->second.neighbours ^= cell5;
                        j->second.cells[5] = nullptr;

                        break;
                    case 6:
                        j->second.neighbours ^= cell6;
                        j->second.cells[6] = nullptr;

                        break;
                }
            }
        }
    }

    return cellNR;
};

int MacGrid::fluidUpdate(float dt)
{
    float dTime = timeStep();

    int ret = update(dt);
    float adv = 0;
    for (adv; adv < dt; adv += dTime)
    {
        NavStokes(dTime);

        for(auto it = particles.begin(); it != particles.end(); ++it)
        {
    //        glm::vec3 pos(0.f);

    //        glm::vec3 vel = cells.find(hashFunc(cellOwner(it->pos)))->second.velocity;
            glm::vec3 vel;


                vel = interpolateVel(it->pos.x, it->pos.y, it->pos.z);
                vel = dTime*interpolateVel(it->pos.x + 0.5*dTime*vel.x, it->pos.y + 0.5*dTime*vel.y, it->pos.z + 0.5*dTime*vel.z);

    //            auto nc = cells.find(hashFunc(cellOwner(it->pos+0.5f*dTime*vel)));
    //            vel = dTime*nc->second.velocity;
                it->pos.x += vel.x;
                it->pos.y += vel.y;
                it->pos.z += vel.z;



        }
    }

    return ret;
};

float MacGrid::timeStep()
{
 return 0.1f*1.f/20;
// return 0.25f*(2.f/uMaxx)*(2/uMaxy)*(2/uMaxz)/0.05;
};

void MacGrid::NavStokes(float dt)
{
    int nFluidCells = 0;
    gridCell* fluidCells[cells.size()];

    for(auto it = cells.begin(); it != cells.end(); ++it) {

        if (it->second.layer <= 1)
        {
            it->second.convection = applyConvection(it->second, -dt);

        }



        if(it->second.type == fluid)
        {
            if(it->second.convection != it->second.convection)
                printf("waow");
            it->second.spot = nFluidCells;
            fluidCells[nFluidCells] = &it->second;
            nFluidCells++;
        }

    }

    for(auto it = cells.begin(); it != cells.end(); ++it)
    {
        if(it->second.type == fluid)
        {
            it->second.velocity = it->second.convection;
            it->second.velocity += dt*GRAVITY;
        }
        else
        {
            if(it->second.cNeigh[1] != nullptr)
            {
                if(it->second.cNeigh[1]->type == fluid)
                {
                    it->second.velocity.x = it->second.convection.x;
                    it->second.velocity.x += dt*GRAVITY.x;
                }
            }
            if(it->second.cNeigh[3] != nullptr)
            {
                if(it->second.cNeigh[3]->type == fluid)
                {
                    it->second.velocity.y = it->second.convection.y;
                    it->second.velocity.y += dt*GRAVITY.y;
                }
            }
            if(it->second.cNeigh[5] != nullptr)
            {
                if(it->second.cNeigh[5]->type == fluid)
                {
                    it->second.velocity.z = it->second.convection.z;
                    it->second.velocity.z += dt*GRAVITY.z;
                }
            }

            if(it->second.velocity != it->second.velocity)
                printf("Convection Blew");
        }
    }

    for(auto it = cells.begin(); it != cells.end(); ++it)
    {
        if (it->second.layer <= 2)
        {
            it->second.viscosity = applyViscousity(it->second, dt);
            int obongo = 12;
        }

    }

    for(auto it = cells.begin(); it != cells.end(); ++it)
    {
        if(it->second.type == fluid)
        {
            it->second.velocity += VISC*it->second.viscosity;
        }
        else
        {
            if(it->second.cNeigh[1] != nullptr)
            {
                if(it->second.cNeigh[1]->type == fluid)
                {
                    it->second.velocity.x += VISC*it->second.viscosity.x;
                }
            }
            if(it->second.cNeigh[3] != nullptr)
            {
                if(it->second.cNeigh[3]->type == fluid)
                {
                     it->second.velocity.y += VISC*it->second.viscosity.y;
                }
            }
            if(it->second.cNeigh[5] != nullptr)
            {
                if(it->second.cNeigh[5]->type == fluid)
                {
                    it->second.velocity.z += VISC*it->second.viscosity.z;
                }
            }
        }
            if(it->second.velocity != it->second.velocity)
                printf("Viscousity blew");

    }

    applyPressure(fluidCells, nFluidCells, dt);

    for(auto it = cells.begin(); it != cells.end(); ++it)
    {


        if (it->second.type == fluid)
        {
            glm::vec3 pressureGrad = calc_grad(it->second.position, it->second);
//            it->second.velocity -= dt*pressureGrad;
            if(it->second.cNeigh[1] != nullptr)
            {
                if(it->second.cNeigh[1]->type != solid)
                {
                    it->second.velocity.x -= dt*pressureGrad.x;
                }
            }
            if(it->second.cNeigh[3] != nullptr)
            {
                if(it->second.cNeigh[3]->type != solid)
                {
                    it->second.velocity.y -= dt*pressureGrad.y;
                }
            }
            if(it->second.cNeigh[5] != nullptr)
            {
                if(it->second.cNeigh[5]->type != solid)
                {
                    it->second.velocity.z -= dt*pressureGrad.z;
                }
            }
            it->second.layer = 0;

        } else if(it->second.type == air)
        {
            glm::vec3 pressureGrad = calc_grad(it->second.position, it->second);
            if(it->second.cNeigh[1] != nullptr)
            {
                if(it->second.cNeigh[1]->type == fluid)
                {
                    it->second.velocity.x -= dt*pressureGrad.x;
                }
            }
            if(it->second.cNeigh[3] != nullptr)
            {
                if(it->second.cNeigh[3]->type == fluid)
                {
                    it->second.velocity.y -= dt*pressureGrad.y;
                }
            }
            if(it->second.cNeigh[5] != nullptr)
            {
                if(it->second.cNeigh[5]->type == fluid)
                {
                    it->second.velocity.z -= dt*pressureGrad.z;
                }
            }
            it->second.layer = -1;
        }
        else
        {
            it->second.layer = -1;
        }
        if(it->second.velocity != it->second.velocity)
            printf("presure grad blew");
    }

    for (int i = 1; i <= 2; ++i)
    {
        for (auto it = cells.begin(); it != cells.end(); ++it)
        {
            if(it->second.layer == -1)
            {
                glm::vec3 avg = it->second.velocity;
//                glm::vec3 avg = glm::vec3(0.f);
                int divider = 1;
                int dividerx = 1;
                int dividery = 1;
                int dividerz = 1;

                bool hasNX = false;
                bool hasNY = false;
                bool hasNZ = false;

                bool hasN = false;

                for (int j = 0; j < 6; ++j)
                {
                    if(it->second.cNeighbours[j] != MAGIC)
                    {
                        gridCell fc = *it->second.cNeigh[j];

                        if(fc.layer == i-1)
                        {
                            hasN = true;
                        }
                        if (j == 1)
                        {
                            if(fc.type == fluid)
                                hasNX = true;

                        }
                        else if(j == 3)
                        {
                            if(fc.type == fluid)
                                hasNY = true;
                        }
                        else if(j == 5)
                        {
                            if(fc.type == fluid)
                                hasNZ = true;
                        }

                    }
                }
                for (int j = 0; j < 6; ++j)
                {
                    if(it->second.cNeighbours[j] != MAGIC)
                    {
                        gridCell fc = *it->second.cNeigh[j];
                        if(fc.layer == i-1)
                        {
                            if(!hasNX)
                            {
                                avg.x += fc.velocity.x;
                                dividerx++;
                            }
                            if(!hasNY)
                            {
                                avg.y += fc.velocity.y;
                                dividery++;
                            }
                            if(!hasNZ)
                            {
                                avg.z += fc.velocity.z;
                                dividerz++;
                            }

                            divider++;
                        }
                    }
                }

                if(hasN)
                {

                    it->second.velocity.x = avg.x/(float)dividerx;
                    it->second.velocity.y = avg.y/(float)dividery;
                    it->second.velocity.z = avg.z/(float)dividerz;
//                    it->second.velocity = avg/(float)divider;
                    it->second.layer = i;
                } else
                {
//                    if(it->second.velocity.x >= 100000 || it->second.velocity.y >= 100000 || it->second.velocity.z >= 100000)
//                    {
//                        printf("extrapolation blew\n");
//                    };
                }
            }
        }
    }
    float max_vel = 0;
    float max_velx = 0;
    float max_vely = 0;
    float max_velz = 0;
    for (auto it = cells.begin(); it != cells.end(); ++it) {
        if (it->second.type == air || it->second.type == fluid)
        {
            glm::vec3 direction = glm::normalize(it->second.velocity);

            float xDir = 0;
            float yDir = 0;
            float zDir = 0;
            if (direction.x <= 0.f) {
                xDir = floorf(direction.x);
            } else if (direction.x > 0.f) {
                xDir = ceilf(direction.x);
            } else {
                it->second.velocity.x = 0;
            }

            if (direction.y <= 0.0f) {
                yDir = floorf(direction.y);
            } else if (direction.y > 0.f) {
                yDir = ceilf(direction.y);
            } else {
                it->second.velocity.y = 0;
            }

            if (direction.z <= 0.f) {
                zDir = floorf(direction.z);
            } else if (direction.z > 0.f) {
                zDir = ceilf(direction.z);
            } else {
                it->second.velocity.z = 0;
            }


            auto cex = cells.find(hashFunc(it->second.position + glm::vec3(xDir, 0.f, 0.f)));
            auto cey = cells.find(hashFunc(it->second.position + glm::vec3(0.f, yDir, 0.f)));
            auto cez = cells.find(hashFunc(it->second.position + glm::vec3(0.f, 0.f, zDir)));

            if (cex != cells.end())
            {
                if (cex->second.type == solid)
                {
                    it->second.velocity.x = 0.f;
                    cex->second.velocity.x = -it->second.velocity.x;
                }
            }
            if (cey != cells.end())
            {
                if (cey->second.type == solid)
                {
                    it->second.velocity.y = 0.f;
                    cey->second.velocity.y = -it->second.velocity.y;
                }
            }
            if (cez != cells.end())
            {
                if (cez->second.type == solid)
                {
                    it->second.velocity.z = 0.f;
                    cez->second.velocity.z = -it->second.velocity.z;
                }
            }

        }

//            float vecL = sqrt(pow(it->second.velocity.x, 2)+ pow(it->second.velocity.y, 2) + pow(it->second.velocity.z, 2));
        float vecL = fmax(fmax(fabs(it->second.velocity.x), fabs(it->second.velocity.y)), it->second.velocity.z);
        if(vecL > max_vel)
        {
            max_vel = vecL;
        }
        float ux = fabs(it->second.velocity.x);
        float uy = fabs(it->second.velocity.y);
        float uz = fabs(it->second.velocity.z);
        if(ux > max_velx)
        {
            max_velx = vecL;
        }
        if(uy > max_vely)
        {
            max_vely = vecL;
        }
        if(uz > max_velz)
        {
            max_velz = vecL;
        }

//        if (it->second.type != fluid)
//        {
//            it->second.velocity = glm::vec3(0.f);
//        }
       it->second.pressure = 1.f;
    }


    uMaxx = max_velx;
    uMaxy = max_vely;
    uMaxz = max_velz;

    maxvel = fmax(fmax(max_velx, max_vely), max_velz);
};

glm::vec3 MacGrid::applyConvection(gridCell c, float dt)
{
    glm::vec3 tracedPosition = traceParticle(c.position.x, c.position.y, c.position.z, dt);
//    return getVelocity(tracedPosition.x, tracedPosition.y, tracedPosition.z);
    return cells.find(hashFunc(cellOwner(tracedPosition)))->second.velocity;

};

glm::vec3 MacGrid::applyViscousity(gridCell c, float dt)
{
    return dt*glm::vec3(calc_laplace2(c.position, 0, c), calc_laplace2(c.position, 1, c), calc_laplace2(c.position, 2, c));
};

void MacGrid::applyPressure(gridCell** c, int nrf, float dt)
{
    Eigen::SparseMatrix<double> mat(nrf, nrf);
    mat.setZero();
    Eigen::VectorXd vsc(nrf);
    std::vector<Eigen::Triplet<double>> coeff;
    for (int i = 0; i < nrf; ++i)
    {
        float ac = 0.f;
        float nsc = 0.f;
        for (int j = 0; j < 6; ++j)
        {
            if(c[i]->cNeighbours[j] != MAGIC)
            {
                auto NC = c[i]->cNeigh[j];
                if(NC->type == fluid)
                {
                    coeff.push_back(Eigen::Triplet<double>(c[i]->spot, NC->spot, 1.0));
                    nsc += 1.f;
                }else if(NC->type == air)
                {
                    nsc += 1.f;
                    ac += 1.f;
                }
            }
        }
        coeff.push_back(Eigen::Triplet<double>(c[i]->spot, c[i]->spot, -nsc));
        float bVal = (1.f/dt)*calc_diver(c[i]->position, *c[i]) - 1.f*ac;
        if(bVal != bVal)
            printf("pressure1 Blew");
//        float bVal2 = 1.f/dt*calc_diver(c[i]->position, *c[i]) - 0.986*ac;
        vsc[c[i]->spot] = bVal;
    }

    mat.setFromTriplets(coeff.begin(), coeff.end());
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;

//    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(mat);
    Eigen::VectorXd B = solver.compute(mat).solve(vsc);

    for (int k = 0; k < nrf; ++k)
    {
        c[k]->pressure = B[k];
        if(c[k]->pressure != c[k]->pressure)
            printf("Pressure Blew");
    }

};

//////

glm::vec3 MacGrid::traceParticle(double x, double y, double z, float dt)
{
    glm::vec3 V;

    V = getVelocity(x, y, z);
    V = getVelocity(x + 0.5f*dt*V.x, y  + 0.5f*dt*V.y, z  + 0.5f*dt*V.z);
    return (glm::vec3(x, y, z) + (dt*V));

};

glm::vec3 MacGrid::getVelocity(double x, double y, double z)
{
    glm::vec3 V;
    V.x = interPolVal(x, y-0.5, z-0.5, 0);
    V.y = interPolVal(x-0.5, y, z-0.5, 1);
    V.z = interPolVal(x-0.5, y-0.5, z, 2);
    return V;
};

float MacGrid::interPolVal(double x, double y, double z, int ind)
{
    float i;
    float j;
    float k;




    i = floorf(x+0.00001);
    j = floorf(y+0.00001);
    k = floorf(z+0.00001);

    float divWeights =  (i + 1.f - x) * (j + 1.f - y) * (k + 1.f - z) +
                        (x - i) * (j + 1.f - y) * (k + 1.f - z) +
                        (i + 1.f - x) * (y - j) * (k + 1.f - z) +
                        (x - i) * (y - j) * (k + 1.f - z) +
                        (i + 1.f - x) * (j + 1.f - y) * (z - k) +
                        (x - i) * (j + 1.f - y) * (z - k) +
                        (i + 1.f - x) * (y - j) * (z - k) +
                        (x - i) * (y - j) * (z - k);

    float gi = hashFunc(glm::vec3(i, j, k));
    auto che = cells.find(gi);
    gridCell ce = che->second;
    float ret = 0;

//    bool lost0 = false;
//    bool lost1 = false;
//    bool lost2 = false;
//    bool lost3 = false;
//    bool lost4 = false;
//    bool lost5 = false;
//    bool lost6 = false;

    if (che != cells.end()) {
        ret = (i + 1.f - x) * (j + 1.f - y) * (k + 1.f - z) * ce.velocity[ind];

        if (ce.cells[0] != nullptr) {
            ret +=        (x - i) * (j + 1.f - y) * (k + 1.f - z) * ce.cells[0]->velocity[ind];
        } else{
            divWeights -= (x - i) * (j + 1.f - y) * (k + 1.f - z);
        }

        if (ce.cells[1] != nullptr) {
            ret +=        (i + 1.f - x) * (y - j) * (k + 1.f - z) * ce.cells[1]->velocity[ind];
        } else{
            divWeights -= (i + 1.f - x) * (y - j) * (k + 1.f - z);
        }

        if (ce.cells[2] != nullptr) {
            ret +=        (x - i) * (y - j) * (k + 1.f - z) * ce.cells[2]->velocity[ind];
        } else{

            divWeights -= (x - i) * (y - j) * (k + 1.f - z);
        }

        if (ce.cells[3] != nullptr) {
            ret +=        (i + 1.f - x) * (j + 1.f - y) * (z - k) * ce.cells[3]->velocity[ind];
        } else{

            divWeights -= (i + 1.f - x) * (j + 1.f - y) * (z - k);
        }

        if (ce.cells[4] != nullptr) {
            ret +=        (x - i) * (j + 1.f - y) * (z - k) * ce.cells[4]->velocity[ind];
        } else{

            divWeights -= (x - i) * (j + 1.f - y) * (z - k);
        }

        if (ce.cells[5] != nullptr) {
            ret +=        (i + 1 - x) * (y - j) * (z - k) * ce.cells[5]->velocity[ind];
        } else{

            divWeights -= (i + 1 - x) * (y - j) * (z - k);
        }

        if (ce.cells[6] != nullptr) {
            ret +=        (x - i) * (y - j) * (z - k) * ce.cells[6]->velocity[ind];
        } else{

            divWeights -= (x - i) * (y - j) * (z - k);
        }
        if (divWeights == 0)
            divWeights = 1;

        return ret / divWeights;
    }
    return 0;
};

glm::vec3 MacGrid::interpolateVel(double x, double y, double z)
{
    glm::vec3 V;
    V.x = interPolVal(x, y, z, 0);
    V.y = interPolVal(x, y, z, 1);
    V.z = interPolVal(x, y, z, 2);
    return V;
};


float MacGrid::interPolVal2(double x, double y, double z, int ind)
{
    float i;
    float j;
    float k;




    i = roundf(x);
    j = roundf(y);
    k = roundf(z);

    float divWeights =  (i + 1.f - x) * (j + 1.f - y) * (k + 1.f - z) +
                        (x - i) * (j + 1.f - y) * (k + 1.f - z) +
                        (i + 1.f - x) * (y - j) * (k + 1.f - z) +
                        (x - i) * (y - j) * (k + 1.f - z) +
                        (i + 1.f - x) * (j + 1.f - y) * (z - k) +
                        (x - i) * (j + 1.f - y) * (z - k) +
                        (i + 1.f - x) * (y - j) * (z - k) +
                        (x - i) * (y - j) * (z - k);

    float gi = hashFunc(glm::vec3(i, j, k));
    auto che = cells.find(gi);
    gridCell ce = che->second;
    float ret = 0;

//    bool lost0 = false;
//    bool lost1 = false;
//    bool lost2 = false;
//    bool lost3 = false;
//    bool lost4 = false;
//    bool lost5 = false;
//    bool lost6 = false;

    if (che != cells.end()) {
        ret = (i + 1.f - x) * (j + 1.f - y) * (k + 1.f - z) * ce.velocity[ind];

        if (ce.cells[0] != nullptr) {
            ret +=        (x - i) * (j + 1.f - y) * (k + 1.f - z) * ce.cells[0]->velocity[ind];
        } else{
            divWeights -= (x - i) * (j + 1.f - y) * (k + 1.f - z);
        }

        if (ce.cells[1] != nullptr) {
            ret +=        (i + 1.f - x) * (y - j) * (k + 1.f - z) * ce.cells[1]->velocity[ind];
        } else{
            divWeights -= (i + 1.f - x) * (y - j) * (k + 1.f - z);
        }

        if (ce.cells[2] != nullptr) {
            ret +=        (x - i) * (y - j) * (k + 1.f - z) * ce.cells[2]->velocity[ind];
        } else{

            divWeights -= (x - i) * (y - j) * (k + 1.f - z);
        }

        if (ce.cells[3] != nullptr) {
            ret +=        (i + 1.f - x) * (j + 1.f - y) * (z - k) * ce.cells[3]->velocity[ind];
        } else{

            divWeights -= (i + 1.f - x) * (j + 1.f - y) * (z - k);
        }

        if (ce.cells[4] != nullptr) {
            ret +=        (x - i) * (j + 1.f - y) * (z - k) * ce.cells[4]->velocity[ind];
        } else{

            divWeights -= (x - i) * (j + 1.f - y) * (z - k);
        }

        if (ce.cells[5] != nullptr) {
            ret +=        (i + 1 - x) * (y - j) * (z - k) * ce.cells[5]->velocity[ind];
        } else{

            divWeights -= (i + 1 - x) * (y - j) * (z - k);
        }

        if (ce.cells[6] != nullptr) {
            ret +=        (x - i) * (y - j) * (z - k) * ce.cells[6]->velocity[ind];
        } else{

            divWeights -= (x - i) * (y - j) * (z - k);
        }
        if (divWeights == 0)
            divWeights = 1;

        return ret / divWeights;
    }
    return 0;
};
