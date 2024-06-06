#include <time.h>
#include <unistd.h>

#define NUM_PARTICLES 1

typedef struct{
  float x;
  float y;
} Vector2 ;

typedef struct{
  Vector2 position;
  Vector2 velocity;
  float mass;
} Particle ;


Particle particles[NUM_PARTICLES] ;

void PrintParticles(){
  for (int i = 0; i < NUM_PARTICLES ; i++) {
    Particle *particle = &particles[i];
    printf("particle[%i] (%.2f, %.2f)\n", i, particle->position.x, particle->position.y);
  }
}


// init all Paritcles with random positions

void InitializeParticles() {
    // Seed the random number generator
    srand(time(NULL));

    for (int i = 0; i < NUM_PARTICLES; ++i) {
        particles[i].position = (Vector2){rand() % 50, rand() % 50}; // Random positions between 0 and 49
        particles[i].velocity = (Vector2){0, 0};
        particles[i].mass = 1;
    }
}

// Just applies Earth's gravity force (mass times gravity acceleration 9.81 m/s^2) to each particle.

Vector2 ComputeForce(Particle *particle) {
    return (Vector2){0, particle->mass * -9.81};  // -ve becauce we are assuming force against earth's gravity as +ve 
}



void RunSimulation() {
    float totalSimulationTime = 10; // The simulation will run for 10 seconds.
    float currentTime = 0; // This accumulates the time that has passed.
    float dt = 0; // Each step's time delta.

    InitializeParticles();
    PrintParticles();

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    while (currentTime < totalSimulationTime) {
        clock_gettime(CLOCK_MONOTONIC, &end);

        // Calculate dt in seconds
        dt = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9f;

        for (int i = 0; i < NUM_PARTICLES; ++i) {
            Particle *particle = &particles[i];
            Vector2 force = ComputeForce(particle);
            Vector2 acceleration = (Vector2){force.x / particle->mass, force.y / particle->mass};
            particle->velocity.x += acceleration.x * dt;
            particle->velocity.y += acceleration.y * dt;
            particle->position.x += particle->velocity.x * dt;
            particle->position.y += particle->velocity.y * dt;
        }

        PrintParticles();
        currentTime += dt;

        // Update start time for the next iteration
        clock_gettime(CLOCK_MONOTONIC, &start);

        // Sleep for a short duration to avoid busy-waiting (optional)
        usleep(10000); // 10 milliseconds
    }
}



// ------------------------------------------------------------
// -------------------- Rigid Body Dynamics -------------------
// ------------------------------------------------------------


#define NUM_RIGID_BODIES 1

// 2D box shape. Physics engines usually have a couple different classes of shapes
// such as circles, spheres (3D), cylinders, capsules, polygons, polyhedrons (3D)...
typedef struct {
    float width;
    float height;
    float mass;
    float momentOfInertia;
} BoxShape;

// Calculates the inertia of a box shape and stores it in the momentOfInertia variable.
void CalculateBoxInertia(BoxShape *boxShape) {
    float m = boxShape->mass;
    float w = boxShape->width;
    float h = boxShape->height;
    boxShape->momentOfInertia = m * (w * w + h * h) / 12;
}

// Two dimensional rigid body
typedef struct {
    Vector2 position;
    Vector2 linearVelocity;
    float angle;
    float angularVelocity;
    Vector2 force;
    float torque;
    BoxShape shape;
} RigidBody;

// Global array of rigid bodies.
RigidBody rigidBodies[NUM_RIGID_BODIES];

// Prints the position and angle of each body on the output.
// We could instead draw them on screen.
void PrintRigidBodies() {
    for (int i = 0; i < NUM_RIGID_BODIES; ++i) {
        RigidBody *rigidBody = &rigidBodies[i];
        printf("body[%i] p = (%.2f, %.2f), a = %.2f\n", i, rigidBody->position.x, rigidBody->position.y, rigidBody->angle);
    }
}

// Initializes rigid bodies with random positions and angles and zero linear and angular velocities.
// They're all initialized with a box shape of random dimensions.
void InitializeRigidBodies() {
    for (int i = 0; i < NUM_RIGID_BODIES; ++i) {
        RigidBody *rigidBody = &rigidBodies[i];
        rigidBody->position = (Vector2){arc4random_uniform(50), arc4random_uniform(50)};
        rigidBody->angle = arc4random_uniform(360)/360.f * M_PI * 2;
        rigidBody->linearVelocity = (Vector2){0, 0};
        rigidBody->angularVelocity = 0;
        
        BoxShape shape;
        shape.mass = 10;
        shape.width = 1 + arc4random_uniform(2);
        shape.height = 1 + arc4random_uniform(2);
        CalculateBoxInertia(&shape);
        rigidBody->shape = shape;
    }
}

// Applies a force at a point in the body, inducing some torque.
void ComputeForceAndTorque(RigidBody *rigidBody) {
    Vector2 f = (Vector2){0, 100};
    rigidBody->force = f;
    // r is the 'arm vector' that goes from the center of mass to the point of force application
    Vector2 r = (Vector2){rigidBody->shape.width / 2, rigidBody->shape.height / 2};
    rigidBody->torque = r.x * f.y - r.y * f.x;
}

void RunRigidBodySimulation() {
    float totalSimulationTime = 10; // The simulation will run for 10 seconds.
    float currentTime = 0; // This accumulates the time that has passed.
    float dt = 1; // Each step will take one second.
    
    InitializeRigidBodies();
    PrintRigidBodies();
    
    while (currentTime < totalSimulationTime) {
        sleep(dt);
        
        for (int i = 0; i < NUM_RIGID_BODIES; ++i) {
            RigidBody *rigidBody = &rigidBodies[i];
            ComputeForceAndTorque(rigidBody);
            Vector2 linearAcceleration = (Vector2){rigidBody->force.x / rigidBody->shape.mass, rigidBody->force.y / rigidBody->shape.mass};
            rigidBody->linearVelocity.x += linearAcceleration.x * dt;
            rigidBody->linearVelocity.y += linearAcceleration.y * dt;
            rigidBody->position.x += rigidBody->linearVelocity.x * dt;
            rigidBody->position.y += rigidBody->linearVelocity.y * dt;
            float angularAcceleration = rigidBody->torque / rigidBody->shape.momentOfInertia;
            rigidBody->angularVelocity += angularAcceleration * dt;
            rigidBody->angle += rigidBody->angularVelocity * dt;
        }
        
        PrintRigidBodies();
        currentTime += dt;
    }
}
