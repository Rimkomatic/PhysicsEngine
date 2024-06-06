# Rigid Body Dynamics Simulation

This project demonstrates basic concepts of rigid body dynamics and particle simulation in a 2D space. The simulation includes gravitational force acting on particles and rigid body physics including force, torque, and inertia.

## Project Structure

- `main.c`: The main file to run the simulation.
- `physics.h`: Header file containing the simulation logic and data structures for particles and rigid bodies.

## Simulation Overview

### Particle Simulation

The particle simulation initializes particles with random positions and simulates the effect of gravity on them.

#### Key Functions

- `InitializeParticles()`: Initializes particles with random positions and zero velocity.
- `ComputeForce(Particle *particle)`: Computes the force acting on a particle due to gravity.
- `RunSimulation()`: Runs the particle simulation for a given duration, updating particle positions based on computed forces.

### Rigid Body Simulation

The rigid body simulation initializes rigid bodies with random positions and angles and applies forces and torques to them, simulating their motion.

#### Key Functions

- `InitializeRigidBodies()`: Initializes rigid bodies with random positions, angles, and box shapes.
- `CalculateBoxInertia(BoxShape *boxShape)`: Calculates the moment of inertia for a box shape.
- `ComputeForceAndTorque(RigidBody *rigidBody)`: Applies a force to a rigid body and computes the resulting torque.
- `RunRigidBodySimulation()`: Runs the rigid body simulation for a given duration, updating rigid body positions and angles based on computed forces and torques.

## How to Compile and Run

### Prerequisites

- GCC compiler

### Compilation

To compile the project, run the following command:

``` sh
gcc -o simulation main.c -lm
```

### Execution
To run the particle simulation, execute the compiled binary:

``` sh
./simulation
```

# Physics Behind the Simulation

This project simulates basic 2D particle dynamics and rigid body physics. The simulation is divided into two main parts: particle simulation and rigid body simulation. Each part uses fundamental principles of mechanics to model motion under various forces.

## Particle Simulation

### 1. Initialization
Particles are initialized with random positions within a specified range. Each particle has properties such as position, velocity, and mass.

#### Code:
```c
void InitializeParticles() {
    srand(time(NULL)); // Seed the random number generator
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        particles[i].position = (Vector2){rand() % 50, rand() % 50}; // Random positions
        particles[i].velocity = (Vector2){0, 0}; // Initial velocity is zero
        particles[i].mass = 1; // Uniform mass for simplicity
    }
}
```
### 2. Force Calculation

The force acting on each particle is computed based on gravity. In this case, we simulate Earth's gravitational force.

#### Code:

```c
Vector2 ComputeForce(Particle *particle) {
    return (Vector2){0, particle->mass * -9.81}; // Gravity force (negative as it's downward)
}
```

### 3. Simulation Loop

The simulation runs for a specified total time, updating particle positions and velocities at each time step based on the forces acting on them.

#### Code:

```c
void RunSimulation() {
    float totalSimulationTime = 10; // Simulation runs for 10 seconds
    float currentTime = 0;
    float dt = 0;
    InitializeParticles();
    PrintParticles();
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    while (currentTime < totalSimulationTime) {
        clock_gettime(CLOCK_MONOTONIC, &end);
        dt = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9f; // Calculate time step

        for (int i = 0; i < NUM_PARTICLES; ++i) {
            Particle *particle = &particles[i];
            Vector2 force = ComputeForce(particle);
            Vector2 acceleration = (Vector2){force.x / particle->mass, force.y / particle->mass}; // F = ma
            particle->velocity.x += acceleration.x * dt;
            particle->velocity.y += acceleration.y * dt;
            particle->position.x += particle->velocity.x * dt;
            particle->position.y += particle->velocity.y * dt;
        }

        PrintParticles();
        currentTime += dt;
        clock_gettime(CLOCK_MONOTONIC, &start);
        usleep(10000); // Sleep for 10 milliseconds
    }
}
```

## Rigid Body Simulation

### 1. Initialization

Rigid bodies are initialized with random positions, angles, and box shapes. Each rigid body has properties such as position, linear velocity, angle, angular velocity, force, and torque.

#### Code:

```c
void InitializeRigidBodies() {
    srand(time(NULL)); // Seed the random number generator
    for (int i = 0; i < NUM_RIGID_BODIES; ++i) {
        RigidBody *rigidBody = &rigidBodies[i];
        rigidBody->position = (Vector2){rand() % 50, rand() % 50}; // Random positions
        rigidBody->angle = rand() % 360 / 360.f * M_PI * 2; // Random angles
        rigidBody->linearVelocity = (Vector2){0, 0}; // Initial velocity is zero
        rigidBody->angularVelocity = 0; // Initial angular velocity is zero

        BoxShape shape;
        shape.mass = 10; // Uniform mass for simplicity
        shape.width = 1 + rand() % 2; // Random width
        shape.height = 1 + rand() % 2; // Random height
        CalculateBoxInertia(&shape); // Calculate moment of inertia
        rigidBody->shape = shape;
    }
}
```

### 2. Moment of Inertia Calculation

The moment of inertia for a box shape is calculated based on its dimensions and mass.

#### Code:

```c
void CalculateBoxInertia(BoxShape *boxShape) {
    float m = boxShape->mass;
    float w = boxShape->width;
    float h = boxShape->height;
    boxShape->momentOfInertia = m * (w * w + h * h) / 12; // Formula for a rectangle's moment of inertia
}
``` 

### 3. Force and Torque Calculation

Forces and torques acting on the rigid bodies are computed. The force is applied at a specific point, creating a torque.

#### Code:

```c
void ComputeForceAndTorque(RigidBody *rigidBody) {
    Vector2 f = (Vector2){0, 100}; // Applied force
    rigidBody->force = f;
    Vector2 r = (Vector2){rigidBody->shape.width / 2, rigidBody->shape.height / 2}; // Arm vector
    rigidBody->torque = r.x * f.y - r.y * f.x; // Torque calculation
}
```

### 4. Simulation Loop

The simulation runs for a specified total time, updating rigid body positions, angles, velocities, and angular velocities at each time step.

#### Code:

```c
void RunRigidBodySimulation() {
    float totalSimulationTime = 10; // Simulation runs for 10 seconds
    float currentTime = 0;
    float dt = 1; // Each time step is 1 second
    InitializeRigidBodies();
    PrintRigidBodies();

    while (currentTime < totalSimulationTime) {
        sleep(dt);

        for (int i = 0; i < NUM_RIGID_BODIES; ++i) {
            RigidBody *rigidBody = &rigidBodies[i];
            ComputeForceAndTorque(rigidBody);
            Vector2 linearAcceleration = (Vector2){rigidBody->force.x / rigidBody->shape.mass, rigidBody->force.y / rigidBody->shape.mass}; // F = ma
            rigidBody->linearVelocity.x += linearAcceleration.x * dt;
            rigidBody->linearVelocity.y += linearAcceleration.y * dt;
            rigidBody->position.x += rigidBody->linearVelocity.x * dt;
            rigidBody->position.y += rigidBody->linearVelocity.y * dt;
            float angularAcceleration = rigidBody->torque / rigidBody->shape.momentOfInertia; // τ = Iα
            rigidBody->angularVelocity += angularAcceleration * dt;
            rigidBody->angle += rigidBody->angularVelocity * dt;
        }

        PrintRigidBodies();
        currentTime += dt;
    }
}
```

## Summary

This simulation provides a basic understanding of particle dynamics and rigid body physics using fundamental principles such as Newton's laws of motion, force, torque, and moment of inertia. By simulating these physical properties, we can model and visualize how particles and rigid bodies behave under different forces in a 2D space.
