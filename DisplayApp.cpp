//
// Created by bill on 16/01/18.
//




#include <Modules/Dbg/Dbg.h>
#include "DisplayApp.h"


const int nrOfParticles = 100;

class DisplayApp : public App
{
public:
    bool hasChanged = false;

    GfxSetup win = GfxSetup::WindowMSAA4(1024, 1024, "fluidSim");
    Id shd;

    DbgSetup dbg;

    double dt = 0.016;

    double resetTimer = 0;

    TimePoint fps;

    MeshBuilder red;
    MeshBuilder green;
    MeshBuilder black;

    DrawState drawState;
    DrawState drawLines;
    DrawState drawVelVecs;
    DrawState drawRCells;
    DrawState drawGCells;
    DrawState drawBCells;

    Shader::params params;
    glm::mat4 view;
    glm::mat4 proj;

    MacGrid grid = MacGrid(glm::vec3(5.f, 5.f, 2.5f), glm::vec3(10, 10, 5));

    float angleX = 0.0f;
    float angleY = 0.0f;

    int frameCount = 0;

    int prevNRCells = 0;


    glm::vec2 startMousePos;
    glm::vec3 pointOfInterest;

    bool pointerLock = false;

    bool dCells = false;
    bool updateFluid = false;
    bool dVelField = false;

    bool resetField = false;

    bool firstRun = true;

    bool drop2 = false;

    glm::vec3 optVec1;
    glm::vec3 optVec2;


    void testMouseButton(MouseButton::Code btn, const char* name);
    void printMouseState();
    void handleMouseInput();

    void drawGridBounds(GfxSetup win, Id shd);
    void GridCell(GfxSetup win, Id shd, MeshBuilder* mesh, int type);
    void drawUI();

    AppState::Code OnInit();

    AppState::Code OnRunning();

    AppState::Code OnCleanup();
};


OryolMain(DisplayApp);

Id
createIndexMesh(int numIndices, const void* data, int dataSize) {
    auto setup = MeshSetup::FromData(Usage::InvalidUsage, Usage::Immutable);
    setup.NumVertices = 0;
    setup.NumIndices = numIndices;
    setup.IndicesType = IndexType::Index16;
    setup.VertexDataOffset = InvalidIndex;
    setup.IndexDataOffset = 0;
    return Gfx::CreateResource(setup, data, dataSize);
}

Id
createPipeline(PrimitiveType::Code primType, int layoutSlot, const VertexLayout& layout, Id shd, int sampleCount) {
    auto pipSetup = PipelineSetup::FromShader(shd);
    pipSetup.DepthStencilState.DepthWriteEnabled = true;
    pipSetup.DepthStencilState.DepthCmpFunc = CompareFunc::LessEqual;
    pipSetup.RasterizerState.SampleCount = sampleCount;
    pipSetup.Layouts[layoutSlot] = layout;
    pipSetup.PrimType = primType;
    return Gfx::CreateResource(pipSetup);
}

AppState::Code
DisplayApp::OnInit()
{

    win.DefaultPassAction = PassAction::Clear(glm::vec4(1.f, 1.f, 1.f, 0.96f));
    Gfx::Setup(win);
    Dbg::Setup(dbg);
    Input::Setup();

    Input::SetPointerLockHandler([this](const InputEvent& event) -> PointerLockMode::Code {
        if (event.Button == MouseButton::Left) {
            if (event.Type == InputEvent::MouseButtonDown) {
                this->pointerLock = true;
                return PointerLockMode::Enable;
            }
            else if (event.Type == InputEvent::MouseButtonUp) {
                this->pointerLock = false;
                return PointerLockMode::Disable;
            }
        }
        return PointerLockMode::DontCare;
    });


    /*the particle shape builder setup*/
    ShapeBuilder shapeBuilder;
    shapeBuilder.Color(glm::vec4(0.23f, 0.65f, 0.97f, 1.0f));
    shapeBuilder.Layout = {
            { VertexAttr::Position, VertexFormat::Float3 },
            { VertexAttr::Color0, VertexFormat::UByte4N }
    };
    shapeBuilder.Sphere(PARTICLE_RADIUS, 5, 3);


    this->drawState.Mesh[0] = Gfx::CreateResource(shapeBuilder.Build());


    MeshBuilder vecLines = MeshBuilder();
    vecLines.NumVertices = 2;
    vecLines.IndicesType = IndexType::None;
    vecLines.Layout = {
            {VertexAttr::Position, VertexFormat::Float3},
            {VertexAttr::Color0,   VertexFormat::UByte4N}
    };

    vecLines.Begin();


    float may = 0;
    float miy = CELL_WIDTH;

    //line1
    vecLines.Vertex(0, VertexAttr::Position, may, 0, 0);
    vecLines.Vertex(1, VertexAttr::Position, miy, 0, 0);
    vecLines.Vertex(1, VertexAttr::Color0, 0.97f, 0.0f, 0.1f, 1.0f);
    vecLines.Vertex(0, VertexAttr::Color0, 0.0f, 0.97f, 0.1f, 1.0f);

    StaticArray<uint16_t, 2> indices;
    indices[0] = 0;
    indices[1] = 1;

    Id ad = Gfx::CreateResource(vecLines.Build());



    //grid.bounds = glm::vec3(10+CELL_WIDTH,10+CELL_WIDTH,10+CELL_WIDTH);

    optVec1 = glm::vec3(4 , 4, 2);
    optVec2 = glm::vec3(2, 2, 1);
    for (int i = 0; i < nrOfParticles*nrOfParticles; ++i)
    {
        grid.addParticle(optVec1, optVec2);
    }
    ////////////////////////////////////

    /*Shader and pipeline setup*/
    shd = Gfx::CreateResource(Shader::Setup());
    auto ps = PipelineSetup::FromLayoutAndShader(shapeBuilder.Layout, shd);
    ps.DepthStencilState.DepthWriteEnabled = true;
    ps.DepthStencilState.DepthCmpFunc = CompareFunc::LessEqual;
    ps.RasterizerState.SampleCount = win.SampleCount;
    this->drawState.Pipeline = Gfx::CreateResource(ps);

    ////////////////////////////////////

    GridCell(win, shd, &red, air);
    GridCell(win, shd, &green, fluid);
    GridCell(win, shd, &black, solid);

    const float fbWidth = (const float) Gfx::DisplayAttrs().FramebufferWidth;
    const float fbHeight = (const float) Gfx::DisplayAttrs().FramebufferHeight;
    this->proj = glm::perspectiveFov(glm::radians(45.0f), fbWidth, fbHeight, 0.01f, 100.0f);
    this->view = glm::lookAt(grid.origin,glm::vec3(0.0f, 0.0f, -1.0f),glm::vec3(0.0f, 1.0f, 0.0f));

    drawVelVecs.Pipeline = createPipeline(PrimitiveType::Lines, 1, vecLines.Layout, shd, win.SampleCount);
    drawVelVecs.Mesh[0] = createIndexMesh(2, &indices[0], 2*2);
    drawVelVecs.Mesh[1] = ad;

    drawGridBounds(win, shd);

    return App::OnInit();
};

AppState::Code
DisplayApp::OnRunning()
{

    if(Input::MouseAttached())
    {
        const glm::vec2& mouseScroll = Input::MouseScroll();

        if(mouseScroll[1] > 0)
        {

        }
        else if(mouseScroll[1] < 0)
        {

        }

        const glm::vec2& mouseMove = Input::MouseMovement();
        if(pointerLock)
        {
            this->angleX += mouseMove[0] / 50.f;
            this->angleY += mouseMove[1] / 50.f;
        }
    }

    if (Input::KeyboardAttached()) {

        if (Input::KeyDown(Key::Space)) {
            if(updateFluid)
            {
                this->updateFluid = false;
                firstRun = false;
            } else
            {
                this->updateFluid = true;

            }
        }

        if (Input::KeyDown(Key::R)) {
            if(dCells)
            {
                this->dCells = false;
            } else
            {
                this->dCells = true;
            }
        }

        if (Input::KeyDown(Key::F)) {
            if(dVelField)
            {
                this->dVelField = false;
            } else
            {
                this->dVelField = true;
            }
        }

        if (Input::KeyDown(Key::V)) {
            resetField = true;
            updateFluid = false;
            firstRun = true;
            resetTimer = 0;
        }

        if (Input::KeyDown(Key::N1)) {
            resetField = true;
            optVec1 = glm::vec3(4 , 4, 2);
            optVec2 = glm::vec3(2, 2, 1);
            drop2 = false;
            resetTimer = 0;
        }
        if (Input::KeyDown(Key::N2)) {
            resetField = true;
            optVec1 = glm::vec3(4 , 4, 2);
            optVec2 = glm::vec3(2, 5, 1);
            drop2 = false;
            resetTimer = 0;
        }
        if (Input::KeyDown(Key::N3)) {
            resetField = true;
            optVec1 = glm::vec3(1 , 4, 2);
            optVec2 = glm::vec3(5, 5, 2);
            drop2 = false;
            resetTimer = 0;
        }
        if (Input::KeyDown(Key::N4)) {
            resetField = true;
            optVec1 = glm::vec3(1 , 4, 2);
            optVec2 = glm::vec3(1, 1, 2);
            drop2 = true;
            resetTimer = 0;
        }
        if (Input::KeyDown(Key::N5)) {
            resetField = true;
            optVec1 = glm::vec3(4 , 4, 2);
            optVec2 = glm::vec3(1, 1, 2);
            drop2 = true;
            resetTimer = 0;
        }


        if (Input::KeyPressed(Key::A)) {
            if(!updateFluid && firstRun)
                grid.ChangeGravity(0);
        }
        if (Input::KeyPressed(Key::Z)) {
            if(!updateFluid && firstRun)
                grid.ChangeGravity(1);
        }
        if (Input::KeyPressed(Key::S)) {
            if(!updateFluid && firstRun)
                grid.ChangeGravity(2);
        }
        if (Input::KeyPressed(Key::X)) {
            if(!updateFluid && firstRun)
                grid.ChangeGravity(3);
        }
        if (Input::KeyPressed(Key::D)) {
            if(!updateFluid && firstRun)
                grid.ChangeGravity(4);
        }
        if (Input::KeyPressed(Key::C)) {
            if(!updateFluid && firstRun)
                grid.ChangeGravity(5);
        }

    }

    this->frameCount++;

    if (updateFluid)
    {
        grid.fluidUpdate(this->dt);
    }

    int curNRCells = grid.cells.size();

    if (curNRCells != prevNRCells)
    {
        prevNRCells = curNRCells;
        hasChanged = true;
    } else
    {
        hasChanged = false;
    }

    if(resetField)
    {
        grid.particles.clear();
        grid.cells.clear();
        srand(0);
        if (!drop2)
        {
            for (int i = 0; i < nrOfParticles*nrOfParticles; ++i)
            {
                grid.addParticle(optVec1, optVec2);
            }
        } else
        {
            for (int i = 0; i < nrOfParticles*nrOfParticles/2; ++i)
            {
                grid.addParticle(optVec1, optVec2);
            }
            for (int i = 0; i < nrOfParticles*nrOfParticles/2; ++i)
            {
                glm::vec3 vec = glm::vec3(grid.origin.x, 0, 0);
                grid.addParticle(optVec1, vec+optVec2);
            }
        }
        resetField = false;
    }

    Gfx::BeginPass();
    Gfx::ApplyDrawState(this->drawState);

    float angle = this->frameCount * 0.003f;
    glm::vec3 pos(sin(angle)*45.f, 10.f, cos(angle)*45.f);

    glm::vec3 camPos = glm::vec3(sin(angleX)*45.f, 10, cos(angleX)*45.f);

    this->view = glm::lookAt(camPos, glm::vec3(5,5, 2.5), glm::vec3(0.0f, 1.0f, 0.0f));
    this->params.view = this->view;
    this->params.projection = this->proj;
    //int primGroupIndex = 0;
    for (unsigned int i = 0; i < grid.particles.size(); ++i) {
        glm::mat4 modelTform;
        modelTform = glm::translate(glm::mat4(), grid.particles[i].pos);
        this->params.model = modelTform;

        Gfx::ApplyUniformBlock(this->params);
        Gfx::Draw(0);
        Gfx::Draw(1);
        Gfx::Draw(2);

    }

    Gfx::ApplyDrawState(this->drawLines);
    this->params.model = glm::mat4();
    Gfx::ApplyUniformBlock(this->params);
    Gfx::Draw(PrimitiveGroup(0, 24));

    if (dCells)
    {
        for (auto it = grid.cells.begin(); it != grid.cells.end(); ++it) {

            this->params.model = glm::translate(glm::mat4(), it->second.position);
            Gfx::ApplyUniformBlock(this->params);
            if (it->second.type == fluid)
            {
                Gfx::ApplyDrawState(this->drawBCells);
                Gfx::Draw(PrimitiveGroup(0, 24));
            }
            else if (it->second.type == air)
            {
                Gfx::ApplyDrawState(this->drawGCells);
                Gfx::Draw(PrimitiveGroup(0, 24));
            } else
            {
                Gfx::ApplyDrawState(this->drawRCells);
                Gfx::Draw(PrimitiveGroup(0, 24));
            }
        }
    }

    if(dVelField)
    {
        Gfx::ApplyDrawState(this->drawVelVecs);
        for (auto it = grid.cells.begin(); it != grid.cells.end(); ++it)
        {
            glm::vec3 new_y;
            glm::vec3 new_z;
            glm::vec3 new_x;
            glm::vec3 vec;
            vec = it->second.position;
            if(it->second.velocity[0] != 0 || it->second.velocity[1] != 0 || it->second.velocity[2] != 0)
            {
                if (it->second.velocity.x == 0 && it->second.velocity.z == 0)
                {
                    if (it->second.velocity.y < 0) // rotate 180 degrees
                        vec = glm::vec3(-vec.x, -vec.y, vec.z);

                    // else if direction.y >= 0, leave `vec` as it is.
                }
                else
                {
                    new_y = glm::normalize(it->second.velocity);
                    new_z = glm::normalize(glm::cross(new_y, glm::vec3(0, 1, 0)));

                    // code below will go here.
                }
                new_y = glm::normalize(it->second.velocity);
                new_z = glm::normalize(glm::cross(new_y,glm::vec3(0,1,0)));
                new_x = glm::normalize(glm::cross(new_y,new_z));

                this->params.model = glm::mat4(glm::vec4(new_x, 0), glm::vec4(new_z,0), glm::vec4(new_y,0), glm::vec4(it->second.position, 1));

                Gfx::ApplyUniformBlock(this->params);
                Gfx::Draw(PrimitiveGroup(1, 2));
            }
        }
    }


    Dbg::DrawTextBuffer();
    Dbg::TextColor(0,0,0,1);

    Gfx::EndPass();
    Gfx::CommitFrame();

    Duration frameTime = Clock::LapTime(this->fps);
    this->dt = 1/frameTime.AsMilliSeconds();
    if(updateFluid)
    {
        resetTimer += this->dt;
    }
    if(resetTimer > 25)
    {
        resetField = true;
        firstRun = true;
        resetTimer = 0;
    }
    Dbg::PrintF("\n %d Cells\n\r %d Particles\n\r fps=%f\n\r Gravity = (%f , %f , %f)\n\r",
                curNRCells,
                nrOfParticles*nrOfParticles,
                1/frameTime.AsSeconds(),
                grid.getGravity().x,
                grid.getGravity().y,
                grid.getGravity().z);

    return Gfx::QuitRequested() ? AppState::Cleanup : AppState::Running;
};

void DisplayApp::drawGridBounds(GfxSetup win, Id shd)
{
    MeshBuilder meshBuilder;
    meshBuilder.NumVertices = 24;
    meshBuilder.IndicesType = IndexType::None;
    meshBuilder.Layout = {
            { VertexAttr::Position, VertexFormat::Float3 },
            { VertexAttr::Color0, VertexFormat::UByte4N }
    };

    meshBuilder.Begin();
    //line1
    meshBuilder.Vertex(0, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.bounds.y-0.5, this->grid.bounds.z-0.5);
    meshBuilder.Vertex(1, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.bounds.y-0.5, this->grid.bounds.z-0.5);
    //line2
    meshBuilder.Vertex(2, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.bounds.y-0.5, this->grid.bounds.z-0.5);
    meshBuilder.Vertex(3, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.getMinB().y+0.5f, this->grid.bounds.z-0.5);
    //line3
    meshBuilder.Vertex(4, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.getMinB().y+0.5f, this->grid.bounds.z-0.5);
    meshBuilder.Vertex(5, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.getMinB().y+0.5f, this->grid.getMinB().z+0.5f);
    //line4
    meshBuilder.Vertex(6, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.getMinB().y+0.5f, this->grid.bounds.z-0.5);
    meshBuilder.Vertex(7, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.getMinB().y+0.5f, this->grid.bounds.z-0.5);
    //line5
    meshBuilder.Vertex(8, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.bounds.y-0.5, this->grid.bounds.z-0.5);
    meshBuilder.Vertex(9, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.bounds.y-0.5, this->grid.getMinB().z+0.5f);


    //line4
    meshBuilder.Vertex(10, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.getMinB().y+0.5f, this->grid.getMinB().z+0.5f);
    meshBuilder.Vertex(11, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.getMinB().y+0.5f, this->grid.getMinB().z+0.5f);
    //line6
    meshBuilder.Vertex(12, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.getMinB().y+0.5f, this->grid.getMinB().z+0.5f);
    meshBuilder.Vertex(13, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.bounds.y-0.5, this->grid.getMinB().z+0.5f);
    //line6
    meshBuilder.Vertex(14, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.bounds.y-0.5, this->grid.getMinB().z+0.5f);
    meshBuilder.Vertex(15, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.bounds.y-0.5, this->grid.bounds.z-0.5);
    //line6
    meshBuilder.Vertex(16, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.bounds.y-0.5, this->grid.getMinB().z+0.5f);
    meshBuilder.Vertex(17, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.bounds.y-0.5, this->grid.getMinB().z+0.5f);
    //line7
    meshBuilder.Vertex(18, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.getMinB().y+0.5f, this->grid.getMinB().z+0.5f);
    meshBuilder.Vertex(19, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.getMinB().y+0.5f, this->grid.bounds.z-0.5);

    meshBuilder.Vertex(20, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.getMinB().y+0.5f, this->grid.bounds.z-0.5);
    meshBuilder.Vertex(21, VertexAttr::Position, this->grid.getMinB().x+0.5f, this->grid.bounds.y-0.5, this->grid.bounds.z-0.5);

    meshBuilder.Vertex(22, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.getMinB().y+0.5f, this->grid.getMinB().z+0.5f);
    meshBuilder.Vertex(23, VertexAttr::Position, this->grid.bounds.x-0.5, this->grid.bounds.y-0.5, this->grid.getMinB().z+0.5f);

    for (int j = 0; j < 24; ++j) {
        meshBuilder.Vertex(j, VertexAttr::Color0, 0.0f, 0.47f, 0.74f, 1.0f);
    }

    const int numIndices = 24;
    StaticArray<uint16_t, numIndices> indices;

    for (int i = 0;  i < 24; i+=2) {
            indices[i] = i;
            indices[i+1] = i + 1;

    }

    Id vertexMesh = Gfx::CreateResource(meshBuilder.Build());
    auto& ds = this->drawLines;
    ds.Pipeline = createPipeline(PrimitiveType::Lines, 1, meshBuilder.Layout, shd, win.SampleCount);
    ds.Mesh[0] = createIndexMesh(numIndices, &indices[0], numIndices*2);
    ds.Mesh[1] = vertexMesh;
};

void DisplayApp::GridCell(GfxSetup win, Id shd, MeshBuilder* cellMesh, int type) {
    cellMesh->NumVertices = 24;
    cellMesh->IndicesType = IndexType::None;
    cellMesh->Layout = {
            {VertexAttr::Position, VertexFormat::Float3},
            {VertexAttr::Color0,   VertexFormat::UByte4N}
    };

    cellMesh->Begin();

    float max = CELL_WIDTH;
    float may = CELL_WIDTH;
    float maz = CELL_WIDTH;

    float mix = 0 - CELL_WIDTH;
    float miy = 0 - CELL_WIDTH;
    float miz = 0 - CELL_WIDTH;



    //line1
    cellMesh->Vertex(0, VertexAttr::Position, max, may, maz);
    cellMesh->Vertex(1, VertexAttr::Position, mix, may, maz);
    //line2
    cellMesh->Vertex(2, VertexAttr::Position, max, may, maz);
    cellMesh->Vertex(3, VertexAttr::Position, max, miy, maz);
    //line3
    cellMesh->Vertex(4, VertexAttr::Position, max, miy, maz);
    cellMesh->Vertex(5, VertexAttr::Position, max, miy, miz);
    //line4
    cellMesh->Vertex(6, VertexAttr::Position, max, miy, maz);
    cellMesh->Vertex(7, VertexAttr::Position, mix, miy, maz);
    //line5
    cellMesh->Vertex(8, VertexAttr::Position, max, may, maz);
    cellMesh->Vertex(9, VertexAttr::Position, max, may, miz);


    //line4
    cellMesh->Vertex(10, VertexAttr::Position, mix, miy, miz);
    cellMesh->Vertex(11, VertexAttr::Position, max, miy, miz);
    //line6
    cellMesh->Vertex(12, VertexAttr::Position, mix, miy, miz);
    cellMesh->Vertex(13, VertexAttr::Position, mix, may, miz);
    //line6
    cellMesh->Vertex(14, VertexAttr::Position, mix, may, miz);
    cellMesh->Vertex(15, VertexAttr::Position, mix, may, maz);
    //line6
    cellMesh->Vertex(16, VertexAttr::Position, mix, may, miz);
    cellMesh->Vertex(17, VertexAttr::Position, max, may, miz);
    //line7
    cellMesh->Vertex(18, VertexAttr::Position, mix, miy, miz);
    cellMesh->Vertex(19, VertexAttr::Position, mix, miy, maz);

    cellMesh->Vertex(20, VertexAttr::Position, mix, miy, maz);
    cellMesh->Vertex(21, VertexAttr::Position, mix, may, maz);

    cellMesh->Vertex(22, VertexAttr::Position, max, miy, miz);
    cellMesh->Vertex(23, VertexAttr::Position, max, may, miz);

    const int numIndices = 24;
    StaticArray<uint16_t, numIndices> indices;

    if (type == 0) {
        for (int i = 0; i < 24; ++i) {
            cellMesh->Vertex(i, VertexAttr::Color0, 0.23f, 0.23f, 0.23f, 0.16f);
        }

        for (int i = 0; i < 24; i += 2) {
            indices[i] = i;
            indices[i + 1] = i + 1;
        }

        Id vertexMesh = Gfx::CreateResource(cellMesh->Build());
        auto &ds = this->drawGCells;
        ds.Pipeline = createPipeline(PrimitiveType::Lines, 1, cellMesh->Layout, shd, win.SampleCount);
        ds.Mesh[0] = createIndexMesh(numIndices, &indices[0], numIndices*2);
        ds.Mesh[1] = vertexMesh;
    }

    if (type == 1) {
        for (int i = 0; i < 24; ++i) {
            cellMesh->Vertex(i, VertexAttr::Color0, 0.77f, 0.24f, 0.24f, 1.0f);
        };

        for (int i = 0; i < 24; i += 2) {
            indices[i] = i;
            indices[i + 1] = i + 1;
        }

        Id vertexMesh = Gfx::CreateResource(cellMesh->Build());
        auto &ds = this->drawRCells;
        ds.Pipeline = createPipeline(PrimitiveType::Lines, 1, cellMesh->Layout, shd, win.SampleCount);
        ds.Mesh[0] = createIndexMesh(numIndices, &indices[0], numIndices*2);
        ds.Mesh[1] = vertexMesh;
    }

    if (type == 2) {
        for (int i = 0; i < 24; ++i) {
            cellMesh->Vertex(i, VertexAttr::Color0, 0.1f, 0.32f, 0.87f, 1.0f);
        }


        for (int i = 0; i < 24; i += 2) {
            indices[i] = i;
            indices[i + 1] = i + 1;
        }

        Id vertexMesh = Gfx::CreateResource(cellMesh->Build());

        auto &ds = this->drawBCells;
        ds.Pipeline = createPipeline(PrimitiveType::Lines, 1, cellMesh->Layout, shd, win.SampleCount);
        ds.Mesh[0] = createIndexMesh(numIndices, &indices[0], numIndices*2);
        ds.Mesh[1] = vertexMesh;
    }
};

AppState::Code
DisplayApp::OnCleanup()
{
    Gfx::Discard();
    Input::Discard();
    return App::OnCleanup();
};
